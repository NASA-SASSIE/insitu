#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 18:15:18 2023

@author: suzanne
"""

import sys,os
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import datetime as dt
# from scipy import interpolate, stats, signal
# from collections import deque
# import itertools
# import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# from matplotlib.widgets import Cursor

# from matplotlib.widgets import Button
# from matplotlib.text import Annotation
# from mpl_point_clicker import clicker
from matplotlib.dates import (MONTHLY,DateFormatter,rrulewrapper,RRuleLocator,drange)
# import CalculateDepths as CD
import UpTempO_BuoyMaster as BM
# import UpTempO_PlotsLevel2 as uplotsL2

sys.path.append('/Users/suzanne/git_repos/polarstereo-lonlat-convert-py/polar_convert')
from scipy.signal import medfilt
from scipy.optimize import curve_fit
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import timezone
from waterIce import getL1
import uuid

cList=['k','purple','blue','deepskyblue','cyan','limegreen','darkorange','red'] #'lime','yellow','darkorange','orangered','red','saddlebrown','darkgreen','olive','goldenrod','tan','slategrey']
colorList=['k','purple','blue','deepskyblue','cyan','limegreen','lime','yellow','darkorange','orangered','red','saddlebrown','darkgreen','olive','goldenrod','tan','slategrey']

# assign paths. L2 data brought in for comparison
L1path = '/Users/suzanne/Google Drive/UpTempO/UPTEMPO/WebData'
SASSIEpath = '/Users/suzanne/SASSIE/dataConvert'
figspath = '/Users/suzanne/SASSIE/dataConvert'

bids = ['300534062898720', # 2022 01
        '300534062897730', # 2022 02
        '300534063704980', # 2022 03
        '300534063807110', # 2022 04
        '300534063803100', # 2022 05
        '300534062892700', # 2022 06
        '300534062894700', # 2022 07
        '300534062894740', # 2022 08
        '300534062896730', # 2022 09
        '300534062894730', # 2022 10
        '300534062893700'  # 2022 11
       ]
bb=2
bid=bids[bb]

binf = BM.BuoyMaster(bid)
# make unique buoy figures path
figspath = os.path.join(figspath,f'{binf["name"][0]}_{int(binf["name"][1]):02d}')
if not os.path.exists(figspath):
    os.mkdir(figspath)

level1File = f'{L1path}/UpTempO_{binf["name"][0]}_{int(binf["name"][1]):02d}_{binf["vessel"]}-Last.dat'
# print(level1File)
df1,pdepths,tdepths,sdepths,ddepths,tiltdepths = getL1(level1File,bid,figspath)
print(df1.columns)
Pcols = [col for col in df1.columns if col.startswith('P')]
Tcols = [col for col in df1.columns if col.startswith('T')]
Scols = [col for col in df1.columns if col.startswith('S') and 'SUB' not in col]
print(Pcols)
print(Tcols)
print(Scols)

df1.loc[(df1['Lon']>180),'Lon'] -= 360

df1.fillna(-999.)

# excel file with instrument attributes
fileInst = f'{SASSIEpath}/Instruments_Details.xlsx'
dfi = pd.read_excel(fileInst,sheet_name=f'Buoy{bb+1:02d}')
dfi.set_index('sensor',inplace=True)

# create and open netcdf file
fn = f'{SASSIEpath}/SASSIE_Fall_2022_UpTempO_{int(binf["name"][1]):02d}_L1.nc'
# print(fn)
ds = nc.Dataset(fn,'w',format='NETCDF4')

# get established global attributes
global_attrsFile = pd.read_excel(f'{SASSIEpath}/SASSIE_attributes.xlsx')
dfGA = pd.DataFrame(global_attrsFile)
GA = dfGA[dfGA['title'].str.contains('Drifter Hydrography')==True].to_dict(orient='dict')
# print(GA.keys())
# 
# write global attrs to fn
for key,value in GA.items():
    if key == 'title':
        value[2] = value[2].replace('<, XXXX>',f', Buoy {int(binf["name"][1]):02d}')
    if key == 'summary':
        value[2] = value[2] # Mike will add
    if key == 'references':
        value[2] = value[2] # Mike will add
    if key == 'acknowledgement':
        value[2] = value[2].split('[')[0]
    if key == 'date_created':
        value[2] = dt.datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    if key == 'geospatial_lat_min':
        value[2] = df1['Lat'].min()
    if key == 'geospatial_lat_max':
        value[2] = df1['Lat'].max()
    if key == 'geospatial_lon_min':
        value[2] = df1['Lon'].min()
    if key == 'geospatial_lon_max':
        value[2] = df1['Lon'].max()
    if key == 'time_coverage_start':
        value[2] = df1['Dates'].iloc[0].strftime('%Y-%m-%dT%H:%M:%SZ')
    if key == 'time_coverage_end':
        value[2] = df1['Dates'].iloc[-1].strftime('%Y-%m-%dT%H:%M:%SZ')
    if key == 'time_coverage_duration':
        duration = df1['Dates'].iloc[-1] - df1['Dates'].iloc[0]
        hr, rem = divmod(duration.seconds,3600)
        mn,sec = divmod(rem,60)
        value[2] = f'P{duration.days}DT{hr:02d}H{mn:02d}M{np.round(sec):02d}S'
    if key == 'uuid':
        value[2] = str(uuid.uuid1())
    # print(key,value[2])
    if not key.startswith('instrument'):
        setattr(ds,key,value[2])   # dont know why there is a '2' in the value dict, or why value is a dict
setattr(ds,'BuoyID',bid)

# create time dimension 
time = ds.createDimension('time',None)

# dictionary for netcdf file variable names  {UpTempO: SASSIE}
varNames = {'Lat':'latitude',
            'Lon':'longitude'}
for col in Pcols:
    varNames[col] = f'pressure{col[1:]}'
for col in Tcols:
    varNames[col] = f'temperature{col[1:]}'
for col in Scols:
    varNames[col] = f'salinity{col[1:]}'
varNames['SUB'] = 'submergence_percentage'
varNames['BATT'] = 'battery_voltage'

# create variables
time = ds.createVariable('time','f4',('time',))
time.long_name = 'time of drifting buoy measurement'
time.axis = 'T'
time.standard_name = 'time'
time.units = 'days since 1950-01-01'
time.coverage_content_type = 'coordinate'
time.axis = 'T'
# fill. timedelta objects only have .days and .seconds (called durations)
time[:] = df1['Dates'].apply(lambda x: (x-dt.datetime(1950,1,1)).days + (x-dt.datetime(1950,1,1)).seconds/86400)

for col in df1.columns:
    if col not in ['Year','Month','Day','Hour','Dates']:
        print(varNames[col])
        print(col)
        varNames[col] = ds.createVariable(varNames[col],'f4',('time',),fill_value=-999.)  # no longer a dict elem, but a ds variable! Hmmm
        varNames[col].coordinate = 'time'
        if col == 'Lat':
            varNames[col].long_name = 'latitude of drifting buoy measurement'
            varNames[col].standard_name = 'latitude'
            varNames[col].units = 'degrees_north'
        if col == 'Lon':
            varNames[col].long_name = 'longitude of drifting buoy measurement'
            varNames[col].standard_name = 'longitude'
            varNames[col].units = 'degrees_east'
        if col.startswith('P'):
            long_name = f'pressure measurement from drifting buoy at nominal depth of {pdepths[int(col[1:])-1]}m'
            varNames[col].long_name = long_name
            varNames[col].standard_name = f'ocean pressure{int(col[1:])}'
            varNames[col].units = 'dbar'  
            # varNames[col].nominal_depth = f'{pdepths[int(col[1:])-1]}m'
            varNames[col].valid_min = 0
            varNames[col].valid_max = 100
        if col.startswith('T'):
            long_name = f'temperature measurement from drifting buoy at nominal depth of {tdepths[int(col[1:])]}m'
            varNames[col].long_name = long_name
            varNames[col].standard_name = f'sea water temperature{int(col[1:])}'
            varNames[col].units = 'degrees C'
            # varNames[col].nominal_depth = f'{tdepths[int(col[1:])]}m'
            varNames[col].valid_min = -40
            varNames[col].valid_max = 30
        if col.startswith('S') and col != 'SUB':
            long_name = f'salinity measurement from drifting buoy at nominal depth of {sdepths[int(col[1:])]}m'
            varNames[col].long_name = long_name
            varNames[col].standard_name = f'sea water salinity{int(col[1:])}'
            varNames[col].units = '1'
            # varNames[col].nominal_depth = f'{sdepths[int(col[1:])]}m'
            varNames[col].valid_min = 15
            varNames[col].valid_max = 40
        if col =='BATT':
            varNames[col].long_name = 'battery voltage measurement on buoy'
            varNames[col].standard_name = 'battery_voltage'
            varNames[col].units = 'volts'
        if col =='SUB':
            varNames[col].long_name = 'submergence of buoy'
            varNames[col].standard_name = 'submergence_percentage'
            varNames[col].units = 'percent'
            varNames[col].valid_min = 0
            varNames[col].valid_max = 100
            
        varNames[col][:] = df1[col]
        # varNames[col].valid_min = df1.loc[(df1[col]>-999.),col].min()
        # varNames[col].valid_max = df1.loc[(df1[col]>-999.),col].max()
        varNames[col].coverage_content_type = 'physicalMeasurement'
        if col in dfi.index:
            for attri in dfi.columns:
                if attri != 'depth' and attri != 'abbrev':
                    setattr(varNames[col],f'instrument_{attri}',dfi.loc[col][attri])
            
            
ds.close()

# plot to ensure
fn = f'{SASSIEpath}/SASSIE_Fall_2022_UpTempO_{int(binf["name"][1]):02d}_L1.nc'
ds = xr.open_dataset(fn)
df = ds.to_dataframe()
df.reset_index(inplace=True)  # time was default index, but can't plot like that
# change SASSIE variable names back to UpTempO names
#  must re-establish dict (it got converted to netcdf variables)
varNames = {'Lat':'latitude',
            'Lon':'longitude'}
for col in Pcols:
    varNames[col] = f'pressure{col[1:]}'
for col in Tcols:
    varNames[col] = f'temperature{col[1:]}'
for col in Scols:
    varNames[col] = f'salinity{col[1:]}'
varNames['SUB'] = 'submergence_percentage'
varNames['BATT'] = 'battery_voltage'

df.rename(columns = {v:k for k,v in varNames.items()},inplace=True)  # invert key/value
print(df.columns)  # time is the row index column

buoylab=f"SASSIE {binf['name'][0]} #{binf['name'][1]}"
abbv=binf['imeiabbv']

rule=rrulewrapper(MONTHLY,interval=1)
loc=RRuleLocator(rule)
formatter=DateFormatter('%m/%d/%y')

plt.rcParams['figure.figsize']=(20,8)
plt.rcParams['font.size']=18

for quan in ['Temp','Pressure','Salinity']:
    fig,ax=plt.subplots()
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_tick_params(rotation=85,labelsize=14)
    
    if quan == 'Temp':
        tcols = [col for col in df.columns if col.startswith('T') and not col.startswith('Tilt')]
        yaxlab='Temperature (C)'
        ax.set_ylim(-2.0,10.0)
    
        for ii,tcol in enumerate(tcols):
            ax.plot(df['time'],df[tcol],'-o',color=colorList[ii],ms=1)
    
    if quan == 'Pressure':
        pcols = [col for col in df.columns if col.startswith('P')]
        yaxlab='Pressure (dB)'
        ax.set_ylim(max(pdepths)+10.,0)
    
        for ii,pcol in enumerate(pcols):
            ax.plot(df['time'],df[pcol],'-o',color=colorList[ii],ms=1)
    
    if quan == 'Salinity':
        scols = [col for col in df.columns if col.startswith('S') and not col.startswith('SUB')]
        yaxlab='Salinity'
        ax.set_ylim(20.,40.)
    
        for ii,scol in enumerate(scols):
            ax.plot(df['time'],df[scol],'-o',color=colorList[ii],ms=1)
    
    # print(df['time'].dt.month.iloc[0])
    # print(df['time'].dt.day.iloc[0])
    # print(df['time'].dt.year.iloc[0])
    # print(df['time'].dt.month.iloc[-1])
    # print(df['time'].dt.day.iloc[-1])
    # print(df['time'].dt.year.iloc[-1])
    # exit(-1)
    datelab=f"{int(df['time'].dt.month.iloc[0]):02}/{int(df['time'].dt.day.iloc[0]):02}/{int(df['time'].dt.year.iloc[0])} to {int(df['time'].dt.month.iloc[-1]):02}/{int(df['time'].dt.day.iloc[-1]):02}/{int(df['time'].dt.year.iloc[-1])}"
    plt.title(buoylab+' ('+bid+') '+yaxlab+' Time Series: '+datelab,fontsize=20)
    ax.set_ylabel(yaxlab)
    plt.subplots_adjust(bottom=0.15)
    plt.grid(True)
    plt.savefig(f'{figspath}/{quan}Series{abbv}L2.png')
    plt.show()

