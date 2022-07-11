#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2/7/2022

@author: suzanne
"""
import os, glob, re, sys
import matplotlib.pyplot as plt  # ver 3.5
import matplotlib as mpl  # ver 3.5
import numpy as np
import datetime as dt
import cartopy.crs as ccrs
import ProcessSatFields as pfields
from satelliteBase import BeaufortSatelliteMap
import matlab.engine
import UpTempO_BuoyMaster as BM

# sys.path.append('/Users/suzanne/git_repos/polarstereo-lonlat-convert-py/polar_convert')
# import polar_convert


# # strdate = '20211108'  # or None for 'day before' strdat format mmddyyyy
strdate = '20180929'  

surface = 'SST'
# # get domain in x,y
# xsw,ysw = polar_convert.polar_lonlat_to_xy(-155,68,70,6378.273,0.081816153,'North')
ax0,fig0,figstr0 = BeaufortSatelliteMap(strdate,surface,extent=[-160,-138,69,81])
strdateSSS = '20210901'
# ax1,fig1,figstr1 = smapBeaufort(strdateSSS,level='L3',extent=[-160,-138,69,81])
ax1,fig1,figstr1 = BeaufortSatelliteMap(strdateSSS,surface='SSS',extent=[-160,-138,69,81],src='JPL-L3')  # options: 'JPL-L3', 'RSS-L3', 'JPL-L2B-NRT','JPL-L2B-del'


sstcolors=['k','darkslateblue','blue','deepskyblue','cyan','limegreen','lime','yellow','darkorange','orangered','red','saddlebrown']
sstbounds=[-100,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,2.0,3.0,4.0,5.0]   # temperature
sstcmap = mpl.colors.ListedColormap(sstcolors)
sstnorm = mpl.colors.BoundaryNorm(sstbounds, len(sstcolors))

ssscolors=['k','darkslateblue','blue','deepskyblue','cyan','limegreen','lime','yellow','darkorange','orangered','red','saddlebrown']
sssbounds=list(np.arange(22,34))   # salinity
sssbounds.insert(0,0)
# print(sssbounds)
# print(len(sssbounds),len(ssscolors))
ssscmap = mpl.colors.ListedColormap(ssscolors)
sssnorm = mpl.colors.BoundaryNorm(sssbounds, len(ssscolors))

newline = '\n'
degree = '\u00B0'

# get buoy data
bids={'300534060649670':'2021-01',
      '300534062158480':'2021-04'}
buoyTpts=[]
buoySpts=[]
for ii,bid in enumerate(bids.keys()):
    df = pfields.getPGbuoy(bid,'pscapluw','20210929')
    print(bid,df.tail())
    binf = BM.BuoyMaster(bid)
    if 'Ts' in df and not np.isnan(df['Ts'].iloc[-1]):
        buoyTlabel = f"{bids[bid]} SST={df['Ts'].iloc[-1]:.1f}{degree}C, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N" 
        buoyTpts.append(ax0.scatter(df['Lon'].iloc[-1],df['Lat'].iloc[-1], s=50, c=df['Ts'].iloc[-1],
                    cmap=sstcmap, norm=sstnorm, transform=ccrs.PlateCarree(),
                    edgecolor='k', label=buoyTlabel))
    if 'T1' in df and not np.isnan(df['T1'].iloc[-5]):
        buoyTlabel = f"{bids[bid]} T1={df['T1'].iloc[-5]:.1f}{degree}C at {binf['tdepths'][1]}m, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N" 
        buoyTpts.append(ax0.scatter(df['Lon'].iloc[-5], df['Lat'].iloc[-5], s=50, c=df['T1'].iloc[-5],
                    cmap=sstcmap, norm=sstnorm, transform=ccrs.PlateCarree(),
                    edgecolor='k', label=buoyTlabel))
    if 'CTD-S2' in df and 'D2' in df and not np.isnan(df['CTD-S2'].iloc[-1]):
        buoySlabel = f"{bids[bid]} SSS={df['CTD-S2'].iloc[-1]:.1f} at {df['D2'].iloc[-1]:.2f}m, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N"
        buoySpts.append(ax1.scatter(df['Lon'].iloc[-1], df['Lat'].iloc[-1], s=50, c=df['CTD-S2'].iloc[-1],
                    cmap=ssscmap, norm=sssnorm, transform=ccrs.PlateCarree(),
                    edgecolor='k', label=buoySlabel))
    if 'S1' in df and not np.isnan(df['S1'].iloc[-5]):
        buoySlabel = f"{bids[bid]} S1={df['S1'].iloc[-5]:.2f} at {binf['tdepths'][1]:.2f}m, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N" 
        buoySpts.append(ax1.scatter(df['Lon'].iloc[-5], df['Lat'].iloc[-5], s=50, c=df['S1'].iloc[-5],
                    cmap=ssscmap, norm=sssnorm, transform=ccrs.PlateCarree(),
                    edgecolor='k', label=buoySlabel))

legend10 = ax0.legend(handles=buoyTpts,bbox_to_anchor=(1.1,1),loc=2,borderaxespad=0.,fontsize=9,title='HydroBuoy Data')
frame10 = legend10.get_frame()
frame10.set_facecolor('lightgray')
frame10.set_edgecolor('black')
leg = ax0.get_legend()
for ii in range(len(bids)):
    leg.legendHandles[ii].set_color('k')
    
legend11 = ax1.legend(handles=buoySpts,bbox_to_anchor=(1.1,1),loc=2,borderaxespad=0.,fontsize=9,title='HydroBuoy Data')
frame11 = legend11.get_frame()
frame11.set_facecolor('lightgray')
frame11.set_edgecolor('black')
leg = ax1.get_legend()
for ii in range(len(bids)):
    leg.legendHandles[ii].set_color('k')
    
# get swift data
eng = matlab.engine.start_matlab()
eng.addpath('../swift_telemetry')
# IDs = ['12','16','17']  # '09', ,'13','15'  # matlab format with the semi colons
IDs = ['09']  # '09', ,'13','15'  # matlab format with the semi colons
starttime = '2018-09-01T22:00:00'
endtime = '2018-09-30T22:00:00' # leave empty to retrieve data up to present time
# can do eng.cd(r'../swift_telemetry') if you need to.

swiftTpts=[]
swiftSpts=[]
for ID in IDs:
    dfSwift = pfields.getSWIFT(ID, starttime, endtime, eng)
    # fig5,ax5 = plt.subplots(1,1)
    # ax5.scatter(dfSwift['Lon'],dfSwift['Lat'],20,dfSwift['WaterTemp-0'],cmap=sstcmap)
    # ax5.scatter(dfSwift['Lon'].iloc[-1],dfSwift['Lat'].iloc[-1],150,dfSwift['WaterTemp-0'].iloc[-1],cmap=sstcmap)
    # print(dfSwift.tail(30))
    # plt.show()
    # exit(-1)
    # reqCols=['Lon','Lat]
    # if dfSwfit[]  IF STATEMENT LOOKING FOR NON-NANS IN lon/lat/data
    idx = -1
    plotTemp = False
    plotSali = False
    # plot latest lon/lat/temp/sali. If lon/lat are null, check back a maximum of 5 data times. 
    while idx >= -5:
        if np.isnan(dfSwift['Lon'].iloc[idx]) or np.isnan(dfSwift['Lon'].iloc[idx]):
            idx -= 1
        else:
            Tcols = [col for col in dfSwift.columns if col.startswith('WaterTemp')]
            Scols = [col for col in dfSwift.columns if col.startswith('Salinity')]
            Dcols = [col for col in dfSwift.columns if col.startswith('CTdepth')]
            for tcol,dcol in zip(Tcols,Dcols):
                if not np.isnan(dfSwift[tcol].iloc[idx]):            
                    swiftTlabel=f"SWIFT {ID} SST={dfSwift[tcol].iloc[idx]:.1f}{degree}C at {dfSwift[dcol].iloc[idx]:.2f}m, {dfSwift['Lon'].iloc[idx]:.2f}W, {dfSwift['Lat'].iloc[idx]:.2f}N"
                    swiftTpts.append(ax0.scatter(dfSwift['Lon'].iloc[idx], dfSwift['Lat'].iloc[idx], 50,dfSwift[tcol].iloc[idx],
                               cmap=sstcmap, norm=sstnorm, marker='s', edgecolor='k',transform=ccrs.PlateCarree(), label=swiftTlabel))
                    plotTemp=True
                    break            
            if not plotTemp:
                    swiftTlabel=f"SWIFT {ID} No SST data at {dfSwift['Lon'].iloc[idx]:.2f}W, {dfSwift['Lat'].iloc[idx]:.2f}N"
                
            for scol,dcol in zip(Scols,Dcols):
                if not np.isnan(dfSwift[scol].iloc[idx]):                    
                    swiftSlabel=f"SWIFT {ID} SSS={dfSwift[scol].iloc[idx]:.2f} at {dfSwift[dcol].iloc[idx]:.2f}m, {dfSwift['Lon'].iloc[idx]:.2f}W, {dfSwift['Lat'].iloc[idx]:.2f}N"
                    swiftSpts.append(ax1.scatter(dfSwift['Lon'].iloc[idx], dfSwift['Lat'].iloc[idx], 50,dfSwift[scol].iloc[idx],
                               cmap=ssscmap, norm=sssnorm, marker='s', edgecolor='k',transform=ccrs.PlateCarree(), label=swiftSlabel))
                    plotSali=True
                    break            
            if not plotSali:
                    swiftSlabel=f"SWIFT {ID} No SSS data at {dfSwift['Lon'].iloc[idx]:.2f}W, {dfSwift['Lat'].iloc[idx]:.2f}N"
    if idx >= -5:
        swiftTlabel=f"SWIFT {ID} No SST data."
        swiftSlabel=f"SWIFT {ID} No SSS data."
        
legend20 = ax0.legend(handles=swiftTpts,bbox_to_anchor=(1.1, 0.8), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
frame20 = legend20.get_frame()
frame20.set_facecolor('lightblue')
frame20.set_edgecolor('black')
leg = ax0.get_legend()
for ii in range(len(IDs)):
    leg.legendHandles[ii].set_color('k')

legend21 = ax1.legend(handles=swiftSpts,bbox_to_anchor=(1.1, 0.8), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
frame21 = legend21.get_frame()
frame21.set_facecolor('lightblue')
frame21.set_edgecolor('black')
leg = ax1.get_legend()
for ii in range(len(IDs)):
    leg.legendHandles[ii].set_color('k')

ax0.add_artist(legend10) # july sea ice content, need to add legends back, as by default only last legend is shown
ax0.add_artist(legend20) # drift track month
fig0.savefig(figstr0)

ax1.add_artist(legend11) # july sea ice content, need to add legends back, as by default only last legend is shown
ax1.add_artist(legend21) # drift track month
fig1.savefig(figstr1)

plt.show()
exit(-1)





