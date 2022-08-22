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

def plotSuite(args):

    today = dt.datetime.now()
    args.strdate = "%d%.2d%.2d" % (today.year,today.month,today.day)

    extent=[int(item.strip()) for item in args.mapDomain.split(',')]

    # get the latest ship location for the plots (will show up as a red star)
    if os.path.exists(f'{args.base_dir}/Ship_track.csv'):
        os.remove(f'{args.base_dir}/Ship_track.csv')  # I get an error if I try to overwrite. :(
    os.system(f'/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/; get -e Ship_track.csv -o {args.base_dir}/Ship_track.csv; bye"')
    dfship = pd.read_csv(f'{args.base_dir}/Ship_track.csv',header=None)
    args.shipLon = dfship.iloc[-1][6]
    args.shipLat = dfship.iloc[-1][5]
    print(f'ship location: {args.shipLat:.2f}N, {-1*args.shipLon:.2f}W')

    # make the base for the plots, one for Temperature, one for Salinity
    ax0,fig0,figstr0 = BeaufortSatelliteMap(args,today,surface='SST')
    ax1,fig1,figstr1 = BeaufortSatelliteMap(args,today,surface='SSS')

    # same as above but for zoomed in view
    ax10,fig10, figstr10 = BeaufortSatelliteMap(args,today,surface='SST',zoom=True)
    ax11,fig11, figstr11 = BeaufortSatelliteMap(args,today,surface='SSS',zoom=True)

    cmap = plt.cm.turbo
    normsst = colors.BoundaryNorm(np.arange(-2,6,0.5),cmap.N)
    normsss = colors.BoundaryNorm(np.arange(22,31,0.25),cmap.N)

    newline = '\n'
    degree = '\u00B0'

    # filexlsx = f'{args.base_dir}/InsituData_{today.year}{today.month:02}{today.day:02}T{today.hour:02}:{today.minute:02}:{today.second:02}.xlsx'
    # print('filexlsx',filexlsx)
    # with pd.ExcelWriter(filexlsx) as writer:

    # catch all dictionary for column names
    outHeaders = {'Date':'Date','Lat':'Lat','Lon':'Lon',
       'DateTimeStr':'Date',
       'DateTime':'Date',
       'Depth':'Depth',
       'CTDepth-0':'Depth',
       'CTD-S2':'Salinity',
       'S1':'Salinity',
       'Salinity-0':'Salinity',
       'Salinity':'Salinity',
       'Ts':'Temperature',
       'T1':'Temperature',
       'WaterTemp':'Temperature'}

    ##################################### DRIFTING BUOY DATA ###################
    if bool(args.buoyIDs):
        # get buoy data
        bids = [item.strip() for item in args.buoyIDs.split(',')]
        buoyTpts=[]
        buoyTptsZ=[]
        buoySpts=[]
        buoySptsZ=[]
        for ii,bid in enumerate(bids):
            print(f'Buoy ID: {bid}')
            # make pandas dataframe
            df = pfields.getPGbuoy(args,bid,'pscapluw')
            binf = BM.BuoyMaster(bid)
            columnsWrite = ['Date','Lat','Lon']

            if len(df)>0:  # if there are data
                # temperature data on plots
                try:
                    (tcol,) = [col for col in df.columns if col.startswith('T') or col.startswith('CTD-T')]
                    buoyTlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: {df[tcol].iloc[-1]:.1f}{degree}C, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N"
                    buoyTpts.append(ax0.scatter(df['Lon'],df['Lat'], 2*df['index'], c=df[tcol],
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoyTlabel))
                    if args.smallDomain is not None:
                        buoyTptsZ.append(ax10.scatter(df['Lon'],df['Lat'], 2*df['index'], c=df[tcol],
                                    cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),  #cmap=sstcmap, norm=sstnorm,
                                    edgecolor='face', label=buoyTlabel))
                    columnsWrite.append(tcol)
                except:
                    df['Ts'] = np.nan
                    buoyTlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: no temperature data"
                    buoyTpts.append(ax1.scatter(df['Lon'], df['Lat'], df['Ts'], c=df['Ts'],  # <- dummy vars
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoyTlabel))
                    buoyTptsZ.append(ax11.scatter(df['Lon'], df['Lat'], df['Ts'], c=df['Ts'],
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoyTlabel))
                    columnsWrite.append('Ts')

                # salinity data on plots
                try:
                    (scol,) = [col for col in df.columns if col.startswith('S') or col.startswith('CTD-S')]
                    buoySlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: {df[scol].iloc[-1]:.1f} psu, {df['Lon'].iloc[-1]:.2f}W, {df['Lat'].iloc[-1]:.2f}N"
                    buoySpts.append(ax1.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df[scol],
                                cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoySlabel))
                    if args.smallDomain is not None:
                        buoySptsZ.append(ax11.scatter(df['Lon'], df['Lat'], 2*df['index'], c=df[scol],
                                    cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                    edgecolor='face', label=buoySlabel))
                    columnsWrite.append(scol)
                except:
                    df['S1'] = np.nan
                    buoySlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: no salinity data"
                    buoySpts.append(ax1.scatter(df['Lon'], df['Lat'], df['S1'], c=df['S1'],  # <- dummy vars
                                cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoySlabel))
                    buoySptsZ.append(ax11.scatter(df['Lon'], df['Lat'], df['S1'], c=df['S1'],
                                cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=buoySlabel))
                    columnsWrite.append('S1')

            else:   # if datafiles from website are empty.

                df['Ts'] = np.nan
                buoyTlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: no temperature data"
                buoyTpts.append(ax1.scatter(df['Lon'], df['Lat'], df['Ts'], c=df['Ts'],  # <- dummy vars
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=buoyTlabel))
                buoyTptsZ.append(ax11.scatter(df['Lon'], df['Lat'], df['Ts'], c=df['Ts'],
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=buoyTlabel))
                columnsWrite.append('Ts')

                df['S1'] = np.nan
                buoySlabel = f"{binf['name'][0]}-{int(binf['name'][1]):02d}: no salinity data"
                buoySpts.append(ax1.scatter(df['Lon'], df['Lat'], df['S1'], c=df['S1'],  # <- dummy vars
                            cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=buoySlabel))
                buoySptsZ.append(ax11.scatter(df['Lon'], df['Lat'], df['S1'], c=df['S1'],
                            cmap=cmap, norm=normsss, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=buoySlabel))
                columnsWrite.append('S1')

            print('line 134 in plotSuite',bid,columnsWrite)

            # # reorder to descending dates for writing to ouput
            # if len(df)>0:
            #     df.sort_values(by='Date',ascending=False,inplace=True)
            #     df = df.reset_index(drop=True)
            #     print(df.head())

            dfNew = df[columnsWrite]
            dfNew.rename(columns=outHeaders,inplace=True)
            # sheetname = f'UTO_{binf["name"][0]}-{int(binf["name"][1]):02d}'
            # print('sheetname',sheetname)
            # dfWrite.to_excel(writer, sheet_name=sheetname)

            # writes/append  data to csv file
            utoFile = f'{args.base_dir}/csv/UTO_{binf["name"][0]}-{int(binf["name"][1]):02d}.csv'
            print('utoFile',utoFile)
            if os.path.exists(utoFile):
                dfPrev = pd.read_csv(utoFile)
                if dfPrev.empty is True:
                    print('line 150 dfPrev is empty I guess',dfPrev.head())
                    dfNew.to_csv(utoFile,index=False)
                else:
                    dfPrev = pd.concat([dfPrev,dfNew],axis=0,ignore_index=True)
                    dfPrev.drop_duplicates(inplace=True)
                    dfPrev.to_csv(utoFile,index=False)
            else:
                dfNew.to_csv(utoFile,float_format='%.3f',index=False)

        # legends for temperatures
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

        # legends for salinities
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

    ##################################### SWIFT FLOAT DATA #####################
    if bool(args.swiftIDs):
        # get swift data, this requires matlab
        eng = matlab.engine.start_matlab()
        # eng.addpath(f'{args.base_dir}/swift_telemetry')
        eng.addpath(f'{args.base_dir}/pyfiles')

        IDs = [item for item in args.swiftIDs.split(',')]  # '09', ,'13','15'  # matlab format with the semi colons

        swiftTpts=[]
        swiftTptsZ=[]
        swiftSpts=[]
        swiftSptsZ=[]
        for ID in IDs:
            dfSwift = pfields.getSWIFT(args, ID, eng)
            dfSwift.reset_index(inplace=True)  # used for plotting
            print(dfSwift.tail())
            print(dfSwift.columns)
            columnsWrite = ['DateTime','Lat','Lon','Temperature','Salinity','Depth']

            if not dfSwift['Lon'].isnull().all():
                if not dfSwift['Temperature'].isnull().all():
                    swiftTlabel=f"{ID}: {dfSwift['Temperature'].iloc[-1]:.1f}{degree}C, {dfSwift['Lon'].iloc[-1]:.2f}W, {dfSwift['Lat'].iloc[-1]:.2f}N"
                    swiftTpts.append(ax0.scatter(dfSwift['Lon'], dfSwift['Lat'], 2*dfSwift['index'], dfSwift['Temperature'],
                               cmap=cmap, norm=normsst, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftTlabel))
                    if args.smallDomain is not None:
                        swiftTptsZ.append(ax10.scatter(dfSwift['Lon'], dfSwift['Lat'], 2*dfSwift['index'], dfSwift['Temperature'],
                                   cmap=cmap, norm=normsst, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftTlabel))

                else:
                    swiftTlabel=f"{ID} No SST data."

                if not dfSwift['Salinity'].isnull().all():
                    swiftSlabel=f"{ID}: {dfSwift['Salinity'].iloc[-1]:.2f}, {dfSwift['Lon'].iloc[-1]:.2f}W, {dfSwift['Lat'].iloc[-1]:.2f}N"
                    swiftSpts.append(ax1.scatter(dfSwift['Lon'], dfSwift['Lat'],2*dfSwift['index'],dfSwift['Salinity'],
                               cmap=cmap, norm=normsss, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftSlabel))
                    if args.smallDomain is not None:
                        swiftSptsZ.append(ax11.scatter(dfSwift['Lon'], dfSwift['Lat'],2*dfSwift['index'],dfSwift['Salinity'],
                                   cmap=cmap, norm=normsss, marker='s', edgecolor='face',transform=ccrs.PlateCarree(), label=swiftSlabel))
                else:
                    swiftSlabel=f"{ID} No SSS data."
            else:
                dfSwift['Temperature'] = np.nan
                swiftTlabel = f"{ID}: no temperature data"
                swiftTpts.append(ax1.scatter(dfSwift['Lon'], dfSwift['Lat'], dfSwift['Temperature'], c=dfSwift['Temperature'],  # <- dummy vars
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=swiftTlabel))
                swiftTptsZ.append(ax11.scatter(dfSwift['Lon'], dfSwift['Lat'], dfSwift['Temperature'], c=dfSwift['Temperature'],
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=swiftTlabel))

                dfSwift['Salinity'] = np.nan
                swiftSlabel = f"{ID}: no salinity data"
                swiftSpts.append(ax1.scatter(dfSwift['Lon'], dfSwift['Lat'], dfSwift['Salinity'], c=dfSwift['Salinity'],  # <- dummy vars
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=swiftSlabel))
                swiftSptsZ.append(ax11.scatter(dfSwift['Lon'], dfSwift['Lat'], dfSwift['Salinity'], c=dfSwift['Salinity'],
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=swiftSlabel))

            # reorder to descending dates for writing to ouput
            # dfSwift.sort_values(by='DateObj',ascending=False,inplace=True)
            # dfSwift = dfSwift.reset_index(drop=True)
            # dfSwift.drop(columns=(['index','DateObj']),axis=1, inplace=True)
            print('line 285 polotSuite')
            print(dfSwift.tail())
            dfSwift = dfSwift[columnsWrite]
            print(dfSwift.tail())
            swiftFile = f'{args.base_dir}/csv/Swift{ID}.csv'
            dfSwift.to_csv(swiftFile,float_format='%.3f',index=False)

        # work the legends on the plots
        legend20 = ax0.legend(handles=swiftTpts,bbox_to_anchor=(1.1, 0.6), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
        frame20 = legend20.get_frame()
        frame20.set_facecolor('lightgray')
        frame20.set_edgecolor('black')
        leg = ax0.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')

        legend120 = ax10.legend(handles=swiftTpts,bbox_to_anchor=(1.1, 0.5), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
        frame20 = legend120.get_frame()
        frame20.set_facecolor('lightgray')
        frame20.set_edgecolor('black')
        leg = ax10.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')

        legend21 = ax1.legend(handles=swiftSpts,bbox_to_anchor=(1.1, 0.6), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
        frame21 = legend21.get_frame()
        frame21.set_facecolor('lightgray')
        frame21.set_edgecolor('black')
        leg = ax1.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')

        legend121 = ax11.legend(handles=swiftSpts,bbox_to_anchor=(1.1, 0.5), loc=2, borderaxespad=0.,fontsize=9,title='Swift Data' )
        frame21 = legend121.get_frame()
        frame21.set_facecolor('lightgray')
        frame21.set_edgecolor('black')
        leg = ax11.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')

    ##################################### WAVE GLIDER DATA #####################
    if bool(args.gliderIDs):
        IDdict = {'102740746':  "SV3-130",  # SV3-130 needs double quotes because of -
                  '84929357':   "SV3-153",
                  '1628052144': "SV3-245",
                  '511512553':  "SV3-247"}
        fileLocs = f'{args.base_dir}/csv/WaveGliderPositions.csv'
        fh = open(fileLocs,'w+')

        waveGliderTpts=[]
        waveGliderTptsZ=[]
        waveGliderSpts=[]
        waveGliderSptsZ=[]

        IDs = [item for item in args.gliderIDs.split(',')]
        for ID in IDs:
            print(IDdict[ID])
            dfwaveGlider = pfields.getWaveGlider(args,ID)
            fh.write(f"{IDdict[ID]},{dfwaveGlider['DateTime'].iloc[-1]},{dfwaveGlider['Lat'].iloc[-1]:.3f},{dfwaveGlider['Lon'].iloc[-1]:.3f}\n")
            dfwaveGlider.reset_index(inplace=True)  # used for plotting
            columnsWrite = ['Date','Lat','Lon','Temperature','Salinity','Depth']

            # make a plotting mask
            endPlot = dfwaveGlider['DateTime'].iloc[-1]
            startPlot = endPlot - dt.timedelta(hours = args.hourstoPlot)
            plot = (dfwaveGlider['DateTime']>=startPlot) & (dfwaveGlider['DateTime']<=endPlot) #mask

            if not dfwaveGlider['Lon'].isnull().all():
                waveGliderTlabel=f"{IDdict[ID]}: {dfwaveGlider['Temperature'].iloc[-1]:.1f}{degree}C, {dfwaveGlider['Lon'].iloc[-1]:.2f}W, {dfwaveGlider['Lat'].iloc[-1]:.2f}N"
                waveGliderSlabel=f"{IDdict[ID]}: {dfwaveGlider['Salinity'].iloc[-1]:.2f} {dfwaveGlider['Lon'].iloc[-1]:.2f}W, {dfwaveGlider['Lat'].iloc[-1]:.2f}N"
                if not dfwaveGlider.loc[plot,'Temperature'].isnull().all():
                    waveGliderTpts.append(ax0.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot], dfwaveGlider['index'][plot], dfwaveGlider['Temperature'][plot],
                               cmap=cmap, norm=normsst, marker='D', edgecolor='face',transform=ccrs.PlateCarree(), label=waveGliderTlabel))
                    if args.smallDomain is not None:
                        waveGliderTptsZ.append(ax10.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot], dfwaveGlider['index'][plot], dfwaveGlider['Temperature'][plot],
                                   cmap=cmap, norm=normsst, marker='D', edgecolor='face',transform=ccrs.PlateCarree(), label=waveGliderTlabel))
                else:
                    # waveGliderTlabel=f"{IDdict[ID]}: No temperature data."
                    waveGliderTpts.append(ax1.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot], dfwaveGlider['Temperature'][plot], c=dfwaveGlider['Temperature'][plot],  # <- dummy vars
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=waveGliderTlabel))
                    waveGliderTptsZ.append(ax11.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot], dfwaveGlider['Temperature'][plot], c=dfwaveGlider['Temperature'][plot],
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=waveGliderTlabel))

                if not dfwaveGlider.loc[plot,'Salinity'].isnull().all():
                    waveGliderSpts.append(ax1.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot],dfwaveGlider['index'][plot],dfwaveGlider['Salinity'][plot],
                               cmap=cmap, norm=normsss, marker='D', edgecolor='face',transform=ccrs.PlateCarree(), label=waveGliderSlabel))
                    if args.smallDomain is not None:
                        waveGliderSptsZ.append(ax11.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot],dfwaveGlider['index'][plot],dfwaveGlider['Salinity'][plot],
                                   cmap=cmap, norm=normsss, marker='D', edgecolor='face',transform=ccrs.PlateCarree(), label=waveGliderSlabel))
                else:
                    waveGliderSpts.append(ax1.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot], dfwaveGlider['Salinity'][plot], c=dfwaveGlider['Salinity'][plot],  # <- dummy vars
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=waveGliderSlabel))
                    waveGliderSptsZ.append(ax11.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'][plot], dfwaveGlider['Salinity'][plot], c=dfwaveGlider['Salinity'][plot],
                                cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                                edgecolor='face', label=waveGliderSlabel))
                    # waveGliderSlabel=f"{IDdict[ID]}: No salinity data."
            else:
                # dfwaveGlider['Temperature'] = np.nan
                waveGliderTlabel = f"{IDdict[ID]}: no location data"
                waveGliderTpts.append(ax1.scatter(dfwaveGlider['Lon'][plot], dfwaveGlider['Lat'], dfwaveGlider['Temperature'], c=dfwaveGlider['Temperature'],  # <- dummy vars
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=waveGliderTlabel))
                waveGliderTptsZ.append(ax11.scatter(dfwaveGlider['Lon'], dfwaveGlider['Lat'], dfwaveGlider['Temperature'], c=dfwaveGlider['Temperature'],
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=waveGliderTlabel))

                # dfwaveGlider['Salinity'] = np.nan
                waveGliderSlabel = f"{IDdict[ID]}: no location data"
                waveGliderSpts.append(ax1.scatter(dfwaveGlider['Lon'], dfwaveGlider['Lat'], dfwaveGlider['Salinity'], c=dfwaveGlider['Salinity'],  # <- dummy vars
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=waveGliderSlabel))
                waveGliderSptsZ.append(ax11.scatter(dfwaveGlider['Lon'], dfwaveGlider['Lat'], dfwaveGlider['Salinity'], c=dfwaveGlider['Salinity'],
                            cmap=cmap, norm=normsst, transform=ccrs.PlateCarree(),
                            edgecolor='face', label=waveGliderSlabel))

            dfwaveGlider = dfwaveGlider[columnsWrite]
            gliderFile = f'{args.base_dir}/csv/WaveGlider_{IDdict[ID]}.csv'
            dfwaveGlider.to_csv(gliderFile,float_format='%.3f',index=False)

        fh.close()
        # work the legends on the plots
        legend30 = ax0.legend(handles=waveGliderTpts,bbox_to_anchor=(1.1, 0.3), loc=2, borderaxespad=0.,fontsize=9,title='WaveGlider Data' )
        frame30 = legend30.get_frame()
        frame30.set_facecolor('lightgray')
        frame30.set_edgecolor('black')
        leg = ax0.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')

        legend130 = ax10.legend(handles=waveGliderTpts,bbox_to_anchor=(1.1, 0.2), loc=2, borderaxespad=0.,fontsize=9,title='WaveGlider Data' )
        frame30 = legend130.get_frame()
        frame30.set_facecolor('lightgray')
        frame30.set_edgecolor('black')
        leg = ax10.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')

        legend31 = ax1.legend(handles=waveGliderSpts,bbox_to_anchor=(1.1, 0.3), loc=2, borderaxespad=0.,fontsize=9,title='WaveGlider Data' )
        frame31 = legend31.get_frame()
        frame31.set_facecolor('lightgray')
        frame31.set_edgecolor('black')
        leg = ax1.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')

        legend131 = ax11.legend(handles=waveGliderSpts,bbox_to_anchor=(1.1, 0.2), loc=2, borderaxespad=0.,fontsize=9,title='WaveGlider Data' )
        frame31 = legend131.get_frame()
        frame31.set_facecolor('lightgray')
        frame31.set_edgecolor('black')
        leg = ax11.get_legend()
        for ii in range(len(IDs)):
            leg.legendHandles[ii].set_color('k')


    # need to re-apply legends as it only does the lastest
    if bool(args.buoyIDs):
        ax0.add_artist(legend10) # UTO legend: need to add legends back, as by default only last legend is shown
        ax10.add_artist(legend110)
    if bool(args.swiftIDs):
        ax0.add_artist(legend20) # swift legend
        ax10.add_artist(legend120)
    if bool(args.gliderIDs):
        ax0.add_artist(legend30) # waveglider legend
        ax10.add_artist(legend130)
    fig0.savefig(figstr0)
    fig10.savefig(figstr10)


    if bool(args.buoyIDs):
        ax1.add_artist(legend11) # UTO legend
        ax11.add_artist(legend111)
    if bool(args.swiftIDs):
        ax1.add_artist(legend21) # swift legend
        ax11.add_artist(legend121)
    if bool(args.gliderIDs):
        ax1.add_artist(legend31) # swift legend
        ax11.add_artist(legend131)
    fig1.savefig(figstr1)
    fig11.savefig(figstr11)

#    couldn't get this to work: lftp -u sassie -e 'mirror --only-newer /local/path/to/data /FTP/insitu/data' sftp://ftp.polarscience.org
    # send figures to the ftp site
    # os.system('''/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/insitu/images/; mmv -O old *.png; bye"''')
    # os.system('''/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/insitu/images/; ls -l; bye"''')
    os.system(f'/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/insitu/images/; put {figstr0};put {figstr10};put {figstr1};put {figstr11}; bye"')

    # send the data files (csv) to the ftp site
    os.system(f'/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/insitu/data/; mput {args.base_dir}/csv/*.csv; bye"')
    os.system(f'/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/; mput {args.base_dir}/csv/WaveGliderPositions.csv; bye"')
    # os.system(f'/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/insitu/data/; mput Swift*.csv; bye"')
    # os.system(f'/usr/local/bin/lftp sftp://sassie@ftp.polarscience.org/ --password 2Icy2Fresh! -e "cd /FTP/insitu/data/; mput Wave*.csv; bye"')
    # os.system(f'/usr/local/bin/lftp -u sassie --password 2Icy2Fresh! -e "mirror --only-newer /Users/suzanne/SASSIE/csv/ /FTP/insitu/data" sftp://ftp.polarscience.org')
        # plt.show(block=False)
        # plt.pause(0.001)
        # input('Press enter to close figures.')
        # plt.close('all')
        # exit(-1)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,  fromfile_prefix_chars='@')
    parser.add_argument('--base_dir', type=str, default=os.getcwd(), help='root directory')
    # parser.add_argument('--minLon', type=int, help='minimum Longitude for plot domain')  # cant seem to pass list of ints...
    # parser.add_argument('--maxLon', type=int, help='maximum Longitude for plot domain')
    # parser.add_argument('--minLat', type=int, help='minimum Latitude for plot domain')
    # parser.add_argument('--maxLat', type=int, help='maximum Latitude for plot domain')
    # parser.add_argument('--strdate', type=str, default='None', help='date for which to get data, in yyyymmdd format')
    parser.add_argument('--mapDomain', type=str, help='list of lons/lats: W, E, S, N')
    # parser.add_argument('--shipLocation', type=str, help='lon, lat of Ship')
    parser.add_argument('--smallDomain', type=int, help='+/- m from shipLocation')

    # parser.add_argument('--strdateSSS',type=str,default='None',help='date for which to get , in yyyymmdd format')
    # parser.add_argument('--satelliteICE',type=str,default='g02202', help='Ice Concentration source: \n'
    #                                                       '\t g02202 \n'
    #                                                       '\t nsidc-0081')
    # parser.add_argument('--satelliteSSS',type=str,default='JPL-L3', help='Sea Surface Salinity source: \n'
    #                                                       '\t JPL-L3 \n'
    #                                                       '\t RSS-L3 \n'
    #                                                       '\t JPL-L2B-NRT \n'
    #                                                       '\t JPL-L2B-del \n'  <- we're using this one. Not sure if it's actually 'delayed'
    #                                                       '\t these have not all been tested')
    parser.add_argument('--buoyIDs', type=str, help='list of 15-digit UpTempO buoy IDs. Options: \n'
                                                   '\t 300534060649670 (2021-01) \n'
                                                   '\t 300534062158480 (2021-04)')
    parser.add_argument('--swiftIDs', type=str, help='list of 2-digit Swift float IDs. Options: \n'
                                                   '\t 09 \n'
                                                   '\t 13 \n'
                                                   '\t 15')
    parser.add_argument('--gliderIDs', type=str, help='list of Wave Glider IDs. Options: \n'
                                                   '\t 102740746  (SV3-130) \n'
                                                   '\t 84929357   (SV3-153) \n'
                                                   '\t 1628052144 (SV3-245) \n'
                                                   '\t 511512553  (SV3-247)')
    parser.add_argument('--hourstoPlot', type=int, help='plot data this number of hours leading up to latest')

    args, unknown = parser.parse_known_args()
    print(args)
    # IDs = [item for item in args.swiftIDS.split(',')]  # '09', ,'13','15'  # matlab format with the semi colons
    # print('swift ids',IDs)
    # exit(-1)
    plotSuite(args)
