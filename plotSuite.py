#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2/7/2022

@author: suzanne
"""
import os, glob, re, sys
import matplotlib.pyplot as plt  # ver 3.5
import matplotlib as mpl  # ver 3.5
import matplotlib.colors as colors
import numpy as np
import datetime as dt
import pandas as pd
import cartopy.crs as ccrs
import ProcessFields as pfields
from satelliteBase import BeaufortSatelliteMap
import matlab.engine
import UpTempO_BuoyMaster as BM
import openpyxl

# sys.path.append('/Users/suzanne/git_repos/polarstereo-lonlat-convert-py/polar_convert')
# import polar_convert

def plotSuite(args):
    
    today = dt.datetime.now()
    # BECAUSE we don't have data from all sources for 'now', set today = strdate in .args file.
    today = dt.datetime.strptime(args.strdate,'%Y%m%d')
    starttimeObj = today - dt.timedelta(hours = args.hourstoPlot)
    print(starttimeObj,today)
    # exit(-1)
    # # strdate = '20211108'  # or None for 'day before' strdat format mmddyyyy
    # # get domain in x,y
    # xsw,ysw = polar_convert.polar_lonlat_to_xy(-155,68,70,6378.273,0.081816153,'North')
    extent=[int(item.strip()) for item in args.mapDomain.split(',')]

    ax0,fig0,figstr0 = BeaufortSatelliteMap(args,surface='SST')
    ax1,fig1,figstr1 = BeaufortSatelliteMap(args,surface='SSS') 
    
    ax10,fig10, figstr10 = BeaufortSatelliteMap(args,surface='SST',zoom=True)
    ax11,fig11, figstr11 = BeaufortSatelliteMap(args,surface='SSS',zoom=True)
    # plt.show()
    # exit(-1)
    
    # sstcolors=['k','darkslateblue','blue','deepskyblue','cyan','limegreen','lime','yellow','darkorange','orangered','red','saddlebrown']
    # sstbounds=[-100,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,2.0,3.0,4.0,5.0]   # temperature
    # sstcmap = mpl.colors.ListedColormap(sstcolors)
    # sstnorm = mpl.colors.BoundaryNorm(sstbounds, len(sstcolors))
    
    # ssscolors=['k','darkslateblue','blue','deepskyblue','cyan','limegreen','lime','yellow','darkorange','orangered','red','saddlebrown']
    # sssbounds=list(np.arange(22,34))   # salinity
    # sssbounds.insert(0,0)
    # print(sssbounds)
    # print(len(sssbounds),len(ssscolors))
    # ssscmap = mpl.colors.ListedColormap(ssscolors)
    # sssnorm = mpl.colors.BoundaryNorm(sssbounds, len(ssscolors))
    
    cmap = plt.cm.turbo
    normsst = colors.BoundaryNorm(np.arange(-2,6,0.5),cmap.N)
    normsss = colors.BoundaryNorm(np.arange(22,31,0.25),cmap.N)

    newline = '\n'
    degree = '\u00B0'

    filexls = f'InsituData_{today.year}{today.month:02}{today.day:02}.xlsx'
    with pd.ExcelWriter(filexls) as writer:
        outHeaders = {'Date':'Date',
           'DateTimeStr':'Date',
           'CTD-S2':'Salinity',
           'S1':'Salinity',
           'Salinity-0':'Salinity',
           'Ts':'Temperature',
           'T1':'Temperature',
           'WaterTemp-0':'Temperature'}        
        
        if args.UTOids:
            # get buoy data
            bids = [item.strip() for item in args.UTOids.split(',')]
            print()
            print(bids)
            print()
            # exit(-1)
            # bids={'300534060649670':'2021-01',
            #       '300534062158480':'2021-04'}
            buoyTpts=[]
            buoyTptsZ=[]
            buoySpts=[]
            buoySptsZ=[]
            for ii,bid in enumerate(bids):
                print(ii,bid)
                print()
                if bid=='300534060649670':  # THIS PARA SHOULD BE REMOVED DURING EXPERIMENT, WHEN STRDATE = TODAY
                    df = pfields.getPGbuoy(args,bid,'pscapluw','20210929')
                elif bid=='300534062158480':
                    df = pfields.getPGbuoy(args,bid,'pscapluw','20211013')
                    
                df.reset_index(inplace=True)
                binf = BM.BuoyMaster(bid)
                columnsWrite = ['Date','Lat','Lon']
                
                # if 'Ts' in df and not np.isnan(df['Ts']).all():
                if 'Ts' in df and not df['Ts'].isnull().all():                # buoyTlabel = f"{bids[bid]} SST={df['Ts'].iloc[-1]:.1f}{degree}C, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N" 
                    buoyTlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: {df['Ts'].iloc[-1]:.1f}{degree}C, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N" 
                    buoyTpts.append(ax0.scatter(df['Lon'],df['Lat'], 2*df['index'], c=df['Ts'],
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoyTlabel))
                    if args.smallDomain is not None:
                        buoyTptsZ.append(ax10.scatter(df['Lon'],df['Lat'], 2*df['index'], c=df['Ts'],
                                    cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),  #cmap=sstcmap, norm=sstnorm, 
                                    edgecolor='face', label=buoyTlabel))
                    columnsWrite.append('Ts')
                    
                elif 'T1' in df and not df['T1'].isnull().all():
                    buoyTlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: {df['T1'].iloc[-5]:.1f}{degree}C at {binf['tdepths'][1]}m, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N" 
                    buoyTpts.append(ax0.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df['T1'],
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoyTlabel))
                    if args.smallDomain is not None:
                        buoyTptsZ.append(ax10.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df['T1'],
                                    cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                    edgecolor='face', label=buoyTlabel))
                    columnsWrite.append('T1')
                            
                if 'CTD-S2' in df and 'D2' in df and not df['CTD-S2'].isnull().all():
                    buoySlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: {df['CTD-S2'].iloc[-1]:.1f} at {df['D2'].iloc[-1]:.2f}m, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N"
                    buoySpts.append(ax1.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df['CTD-S2'],
                                cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoySlabel))
                    if args.smallDomain is not None:
                        buoySptsZ.append(ax11.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df['CTD-S2'],
                                    cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                    edgecolor='face', label=buoySlabel))
                    columnsWrite.append('CTD-S2')
    
                if 'S1' in df and not df['S1'].isnull().all():
                    buoySlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: {df['S1'].iloc[-1]:.2f} at {binf['tdepths'][1]:.2f}m, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N" 
                    buoySpts.append(ax1.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df['S1'],
                                cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoySlabel))
                    if args.smallDomain is not None:
                        buoySptsZ.append(ax11.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df['S1'],
                                    cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                    edgecolor='face', label=buoySlabel))
                    columnsWrite.append('S1')
           
                # reorder to descending dates for writing to ouput
                df.sort_values(by='DateObj',ascending=False,inplace=True)
                df = df.reset_index(drop=True)
                df.drop(columns=(['index','DateObj']),axis=1, inplace=True)
                dfWrite = df[columnsWrite]
                dfWrite.rename(columns=outHeaders,inplace=True)
                print('dfWrite')
                print(dfWrite.head())
                sheetname = f'UTO_{binf["name"][0]}-{int(binf["name"][1]):02d}'
                dfWrite.to_excel(writer, sheet_name=sheetname)
                utoFile = f'UTO_{binf["name"][0]}-{int(binf["name"][1]):02d}.xlsx'
                dfWrite.to_excel(utoFile)
                
            legend10 = ax0.legend(handles=buoyTpts,bbox_to_anchor=(1.1,1),loc=2,borderaxespad=0.,fontsize=9,title='HydroBuoy Data')
            frame10 = legend10.get_frame()
            frame10.set_facecolor('lightgray')
            frame10.set_edgecolor('black')
            leg = ax0.get_legend()
            for ii in range(len(bids)):
                leg.legendHandles[ii].set_color('k')
    
            legend110 = ax10.legend(handles=buoyTptsZ,bbox_to_anchor=(1.1,1),loc=2,borderaxespad=0.,fontsize=9,title='HydroBuoy Data')
            frame10 = legend110.get_frame()
            frame10.set_facecolor('lightgray')
            frame10.set_edgecolor('black')
            leg = ax10.get_legend()
            for ii in range(len(bids)):
                leg.legendHandles[ii].set_color('k')
                
            legend11 = ax1.legend(handles=buoySpts,bbox_to_anchor=(1.1,1),loc=2,borderaxespad=0.,fontsize=9,title='HydroBuoy Data')
            frame11 = legend11.get_frame()
            frame11.set_facecolor('lightgray')
            frame11.set_edgecolor('black')
            leg = ax1.get_legend()
            for ii in range(len(bids)):
                leg.legendHandles[ii].set_color('k')
                        
            legend111 = ax11.legend(handles=buoySptsZ,bbox_to_anchor=(1.1,1),loc=2,borderaxespad=0.,fontsize=9,title='HydroBuoy Data')
            frame11 = legend111.get_frame()
            frame11.set_facecolor('lightgray')
            frame11.set_edgecolor('black')
            leg = ax11.get_legend()
            for ii in range(len(bids)):
                leg.legendHandles[ii].set_color('k')
    
            # plt.show()
            # exit(-1)
            
        if args.swiftIDS:
            # get swift data
            eng = matlab.engine.start_matlab()
            eng.addpath(f'{args.base_dir}/swift_telemetry')
            # IDs = ['12','16','17']  # '09', ,'13','15'  # matlab format with the semi colons
            IDs = [item for item in args.swiftIDS.split()]  # '09', ,'13','15'  # matlab format with the semi colons
            print('swift ids',IDs)
            # exit(-1)
            # WHEN data are current, use the next two lines
            # starttime = starttimeObj.strftime('%Y-%m-%dT%H:%M:%S')
            # endtime = today.strftime('%Y-%m-%dT%H:%M:%S')
    
            starttime = '2018-09-12T00:00:00'
            endtime = '2018-09-30T00:00:00' 
            # can do eng.cd(r'../swift_telemetry') if you need to.
            
            swiftTpts=[]
            swiftTptsZ=[]
            swiftSpts=[]
            swiftSptsZ=[]
            for ID in IDs:
                dfSwift = pfields.getSWIFT(args, ID, starttime, endtime, eng) 
                print(dfSwift.head())
                dfSwift.reset_index(inplace=True)  # used for plotting
                print(dfSwift.tail())
                columnsWrite = ['DateTimeStr','Lat','Lon','WaterTemp-0','Salinity-0']
                # fig5,ax5 = plt.subplots(1,1)
                # ax5.scatter(dfSwift['Lon'],dfSwift['Lat'],20,dfSwift['WaterTemp-0'],cmap=sstcmap)
                # ax5.scatter(dfSwift['Lon'].iloc[-1],dfSwift['Lat'].iloc[-1],150,dfSwift['WaterTemp-0'].iloc[-1],cmap=sstcmap)
                # print(dfSwift.tail(30))
                # plt.show()
                # exit(-1)
                # reqCols=['Lon','Lat]
                # if dfSwfit[]  IF STATEMENT LOOKING FOR NON-NANS IN lon/lat/data
                if not dfSwift['Lon'].isnull().all():
                    Tcols = [col for col in dfSwift.columns if col.startswith('WaterTemp')]
                    Scols = [col for col in dfSwift.columns if col.startswith('Salinity')]
                    Dcols = [col for col in dfSwift.columns if col.startswith('CTdepth')]
                    for tcol,dcol in zip(Tcols,Dcols):
                        if not dfSwift[tcol].isnull().all():            
                            swiftTlabel=f"SWIFT {ID}: {dfSwift[tcol].iloc[-1]:.1f}{degree}C at {dfSwift[dcol].iloc[-1]:.2f}m, {dfSwift['Lon'].iloc[-1]:.2f}W, {dfSwift['Lat'].iloc[-1]:.2f}N"
                            swiftTpts.append(ax0.scatter(dfSwift['Lon'], dfSwift['Lat'], 2*dfSwift['index'], dfSwift[tcol],
                                       cmap=cmap, norm=normsst, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftTlabel))
                            if args.smallDomain is not None:
                                swiftTptsZ.append(ax10.scatter(dfSwift['Lon'], dfSwift['Lat'], 2*dfSwift['index'], dfSwift[tcol],
                                           cmap=cmap, norm=normsst, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftTlabel))
                            
                        else:
                            swiftTlabel=f"SWIFT {ID} No SST data."
                                                                   
                    for scol,dcol in zip(Scols,Dcols):
                        if not dfSwift[scol].isnull().all():                    
                            swiftSlabel=f"SWIFT {ID}: {dfSwift[scol].iloc[-1]:.2f} at {dfSwift[dcol].iloc[-1]:.2f}m, {dfSwift['Lon'].iloc[-1]:.2f}W, {dfSwift['Lat'].iloc[-1]:.2f}N"
                            swiftSpts.append(ax1.scatter(dfSwift['Lon'], dfSwift['Lat'],2*dfSwift['index'],dfSwift[scol],
                                       cmap=cmap, norm=normsss, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftSlabel))
                            if args.smallDomain is not None:
                                swiftSptsZ.append(ax11.scatter(dfSwift['Lon'], dfSwift['Lat'],2*dfSwift['index'],dfSwift[scol],
                                           cmap=cmap, norm=normsss, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftSlabel))
                        else:
                            swiftSlabel=f"SWIFT {ID} No SSS data."
                else:
                    swiftTlabel=f"SWIFT {ID} No lon/lat."
                    swiftSlabel=f"SWIFT {ID} No lon/lat."
            
                # reorder to descending dates for writing to ouput
                dfSwift.sort_values(by='DateObj',ascending=False,inplace=True)
                dfSwift = dfSwift.reset_index(drop=True)
                dfSwift.drop(columns=(['index','DateObj']),axis=1, inplace=True)
                dfSwiftWrite = dfSwift[columnsWrite]
                dfSwiftWrite.rename(columns=outHeaders,inplace=True)
                print('dfSwiftWrite')
                print(dfSwiftWrite.head())
                sheetname = f'Swift{ID}'
                dfSwiftWrite.to_excel(writer, sheet_name=sheetname)
                SwiftFile = f'Swift{ID}.xlsx'
                dfSwiftWrite.to_excel(SwiftFile)

            legend20 = ax0.legend(handles=swiftTpts,bbox_to_anchor=(1.1, 0.8), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
            frame20 = legend20.get_frame()
            frame20.set_facecolor('lightblue')
            frame20.set_edgecolor('black')
            leg = ax0.get_legend()
            for ii in range(len(IDs)):
                leg.legendHandles[ii].set_color('k')
                
            legend120 = ax10.legend(handles=swiftTpts,bbox_to_anchor=(1.1, 0.8), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
            frame20 = legend120.get_frame()
            frame20.set_facecolor('lightblue')
            frame20.set_edgecolor('black')
            leg = ax10.get_legend()
            for ii in range(len(IDs)):
                leg.legendHandles[ii].set_color('k')
            
            legend21 = ax1.legend(handles=swiftSpts,bbox_to_anchor=(1.1, 0.8), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
            frame21 = legend21.get_frame()
            frame21.set_facecolor('lightblue')
            frame21.set_edgecolor('black')
            leg = ax1.get_legend()
            for ii in range(len(IDs)):
                leg.legendHandles[ii].set_color('k')
        
            legend121 = ax11.legend(handles=swiftSpts,bbox_to_anchor=(1.1, 0.8), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
            frame21 = legend121.get_frame()
            frame21.set_facecolor('lightblue')
            frame21.set_edgecolor('black')
            leg = ax11.get_legend()
            for ii in range(len(IDs)):
                leg.legendHandles[ii].set_color('k')
        
        if args.UTOids:
            ax0.add_artist(legend10) # UTO legend: need to add legends back, as by default only last legend is shown
            ax10.add_artist(legend110)
        if args.swiftIDS:
            ax0.add_artist(legend20) # swift legend
            ax10.add_artist(legend120) 
        fig0.savefig(figstr0)
        fig10.savefig(figstr10)
        
        
        if args.UTOids:
            ax1.add_artist(legend11) # UTO legend 
            ax11.add_artist(legend111) 
        if args.swiftIDS:
            ax1.add_artist(legend21) # swift legend
            ax11.add_artist(legend121)
        fig1.savefig(figstr1)
        fig11.savefig(figstr11)
        print(figstr11)
        
        plt.show(block=False)
        plt.pause(0.001)
        input('Press enter to close figures.')
        plt.close('all')
        exit(-1)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,  fromfile_prefix_chars='@')
    parser.add_argument('--base_dir', type=str, default=os.getcwd(), help='root directory')
    # parser.add_argument('--minLon', type=int, help='minimum Longitude for plot domain')  # cant seem to pass list of ints...
    # parser.add_argument('--maxLon', type=int, help='maximum Longitude for plot domain')
    # parser.add_argument('--minLat', type=int, help='minimum Latitude for plot domain')
    # parser.add_argument('--maxLat', type=int, help='maximum Latitude for plot domain')
    parser.add_argument('--strdate', type=str, default='None', help='date for which to get data, in yyyymmdd format')
    parser.add_argument('--mapDomain', type=str, help='list of lons/lats: W, E, S, N')
    parser.add_argument('--shipLocation', type=str, help='lon, lat of Ship')
    parser.add_argument('--smallDomain', type=int, help='+/- m from shipLocation')
    
    # parser.add_argument('--strdateSSS',type=str,default='None',help='date for which to get , in yyyymmdd format')
    parser.add_argument('--satelliteICE',type=str,default='g02202', help='Ice Concentration source: \n'
                                                          '\t g02202 \n'
                                                          '\t nsidc-0081')
    parser.add_argument('--satelliteSSS',type=str,default='JPL-L3', help='Sea Surface Salinity source: \n'
                                                          '\t JPL-L3 \n'
                                                          '\t RSS-L3 \n'
                                                          '\t JPL-L2B-NRT \n'
                                                          '\t JPL-L2B-del \n'
                                                          '\t these have not all been tested')
    parser.add_argument('--UTOids', type=str, help='list of 15-digit UpTempO buoy IDs. Options: \n'
                                                   '\t 300534060649670 (2021-01) \n'
                                                   '\t 300534062158480 (2021-04)')
    parser.add_argument('--swiftIDS', type=str, help='list of 2-digit Swift buoy IDs. Options: \n'
                                                   '\t 09 \n'
                                                   '\t 13 \n'
                                                   '\t 15')
    parser.add_argument('--hourstoPlot', type=int, help='plot data this number of hours leading up to latest')


    # parser.add_argument('-c','--cycles', type=str, help="4-digit number specifying first/last cycles for output filename")
    # parser.add_argument('-R','--Release', type=str, help="3-digit release number for output filename")
    # parser.add_argument('-v','--version', type=str, help="2-digit version number for output filename")
    # parser.add_argument('--Hemisphere','-H', type=int, default=1, help='1 for Northern, -1 for Southern')
    # parser.add_argument('--mosaic', '-m', type=str)
    # parser.add_argument('--out_path', '-o', type=str, help='default is ATL15_file path')
    # parser.add_argument('--pdf', action='store_true', default=False, help='write images to .pdf file')
    # parser.add_argument('--nolog', action='store_true', default=False, help='no writing errors to .log file')
    args, unknown = parser.parse_known_args()
    print(args)    
    # exit(-1)
    plotSuite(args) 





