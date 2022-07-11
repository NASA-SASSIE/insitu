#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2/7/2022

@author: suzanne
"""

import os
import requests
import pandas as pd
import glob as glob

# plotting packages
import matplotlib.pyplot as plt # ver 3.5
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib as mpl  # ver 3.5
import numpy as np
import datetime as dt

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import ProcessFields as pfields
# import IABPplots as iplots

def BeaufortSatelliteMap(args,surface='SST',zoom=False): #strdate=None,surface='SST',extent=[-180,180,50,90],src='JPL-L3'): #nors='n'):

    if args.strdate == None:
        objdate = dt.datetime.now() - dt.timedelta(days=1)
        args.strdate = "%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
    else:
        objdate = dt.datetime.strptime(args.strdate,'%Y%m%d')
    extent=[int(item.strip()) for item in args.mapDomain.split(',')]
    shipLoc = [np.float(item.strip()) for item in args.shipLocation.split(',')]

    if zoom:
        extent[0] = shipLoc[0] - np.float(args.smallDomain)/111000 / np.cos(np.radians(shipLoc[1]))
        extent[1] = shipLoc[0] + np.float(args.smallDomain)/111000 / np.cos(np.radians(shipLoc[1]))
        extent[2] = shipLoc[1] - np.float(args.smallDomain)/111000
        extent[3] = shipLoc[1] + np.float(args.smallDomain)/111000
    
    # cmap = plt.get_cmap('turbo')
    cmap = plt.cm.turbo
    normsst = colors.BoundaryNorm(np.arange(-2,6,0.5),cmap.N)
    normsss = colors.BoundaryNorm(np.arange(22,31,0.25),cmap.N)

    # # establish contour levels and colors
    # sstcolors=['k','darkslateblue','blue','deepskyblue','cyan','limegreen','lime','yellow','darkorange','orangered','red','saddlebrown']
    # sstbounds=[-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,2.0,3.0,4.0,5.0]   # SST

    # ssscolors=['k','darkslateblue','blue','deepskyblue','cyan','limegreen','lime','yellow','darkorange','orangered','red','saddlebrown']
    # sssbounds=list(np.arange(22,30))   # salinity
    # sssbounds.insert(0,0)
    # ssscmap = mpl.colors.ListedColormap(ssscolors[1:])

    icecolors=['0.4','0.5','0.6','0.725','0.85','1.0']
    icelevels=[0.2,0.3,0.4,0.5,0.75]
    
    fig1, ax1 = plt.subplots(1,1, figsize=(10,10))
    if not zoom:
        ax1 = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo(central_longitude=-140))
    else:
        ax1 = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo(central_longitude=shipLoc[0]))
       
    ax1.set_position([0.1,0.22,0.55,0.7])  # position on page

    gl = ax1.gridlines(draw_labels=True, dms=True, 
                       x_inline=False, y_inline=False) #, alpha=0.3) crs=ccrs.PlateCarree(), 
    
    gl.xlocator = mticker.FixedLocator(np.arange(-180,180,20))
    gl.xformatter = LONGITUDE_FORMATTER
    
    gl.ylocator = mticker.FixedLocator(np.arange(68,82,2))
    gl.yformatter = LATITUDE_FORMATTER
    
    ax1.add_feature(cfeature.LAND,facecolor='gray')
    ax1.coastlines(resolution='10m',linewidth=0.5,color='lightgray')
    ax1.set_extent(extent,crs=ccrs.PlateCarree())
    ax1.yaxis.tick_left()
    bbox_ax1 = ax1.get_position()

    # get satellite ice concentration data
    # if objdate < dt.datetime.today() - dt.timedelta(days=30):
    #     print('cannot use nrt nsidc-0081')
    #     exit(-1)
    icedate, ice, icexx, iceyy = pfields.getICE(args,nors='n') 
    if ice is not None:
        print('Date of ice map',icedate,ice.shape)
    # get satellite sst data
    if surface == 'SST':# and not zoom:
        sstdate, sst, sstlon, sstlat = pfields.getSST(args.strdate)
        if sst is not None:
            print('Date of sst map',sstdate,sst.shape)    
            # SST filled contours
            # cbsst = ax1.contourf(sstlon, sstlat, sst, transform=ccrs.PlateCarree(),levels=sstbounds, 
            #                       colors=sstcolors, extend='both')
            cbsst = ax1.pcolormesh(sstlon, sstlat, sst, vmin=-2.0, vmax=5.0, transform=ccrs.PlateCarree(), 
                                   cmap=cmap, norm=normsst)
            ax1.gridlines(crs=ccrs.PlateCarree(),xlocs=np.arange(-180,180,45),ylocs=np.arange(70,80,2), color='gray',x_inline=False)
        else:
            ax1.text(-158,67,'SST data unavailable',color='k',fontsize=10,transform=ccrs.PlateCarree())
        # colorbar for SST
        cb_sst_ax = fig1.add_axes([bbox_ax1.x0, bbox_ax1.y0*.67, bbox_ax1.x1-bbox_ax1.x0, 0.02]) # left, bot, width, height
        cbs = plt.colorbar(cbsst, cax=cb_sst_ax, orientation='horizontal')
        cbs.set_label(label=f'NOAA SST, degC on {sstdate}',fontsize=9) #, position=[bbox_ax1.x1+bbox_ax1.x0,0.02])
        cbs.ax.tick_params(labelsize=9)

    if surface == 'SSS':# and not zoom:
        sssdate, ds = pfields.getSSS(args) 
        if ds is not None:
            if 'L3' in args.satelliteSSS:
                # numssscolors = 40
                # cNorm = colors.Normalize(vmin=22,vmax=30)
                # scalarMap = cm.ScalarMappable(norm=cNorm,cmap=cmap)
                what = ds.smap_sss
                what.where(what>=22)
                cbsss = ax1.pcolormesh(ds.longitude, ds.latitude, ds.smap_sss.where((ds.smap_sss>=22) & (ds.smap_sss<=30)), vmin=22, vmax=30, 
                                        transform=ccrs.PlateCarree(),cmap=cmap, norm=normsss)
                ax1.gridlines(crs=ccrs.PlateCarree(),xlocs=np.arange(-180,180,45),ylocs=np.arange(70,80,2), color='gray',x_inline=False)
            elif 'L2B' in args.satelliteSSS:
                print('you are in L2B')
                cbsss = ax1.scatter(ds['longitude'], ds['latitude'], 15, ds['smap_sss'], vmin=22, vmax=30, 
                                    cmap=cmap, norm=normsss, transform=ccrs.PlateCarree())
                ax1.gridlines(crs=ccrs.PlateCarree(),xlocs=np.arange(-180,180,45),ylocs=np.arange(70,80,2), color='gray',x_inline=False)
        # colorbar for SSS
        cb_sss_ax = fig1.add_axes([bbox_ax1.x0, bbox_ax1.y0*.67, bbox_ax1.x1-bbox_ax1.x0, 0.02]) # left, bot, width, height
        cbs = plt.colorbar(cbsss, cax=cb_sss_ax, orientation='horizontal')
        cbs.set_label(label=f'{args.satelliteSSS} SMAP Salinity on {sssdate}',fontsize=9) #, position=[bbox_ax1.x1+bbox_ax1.x0,0.02])
        cbs.ax.tick_params(labelsize=9)

    # ICE filled contour
    if ice is not None:
        kw = dict(central_latitude=90, central_longitude=-45, true_scale_latitude=70)
        cbice = ax1.contourf(icexx,iceyy,ice, colors=icecolors, levels=icelevels, vmin=0, vmax=0.9, 
                             extend='both', transform=ccrs.Stereographic(**kw))   #use either colors or cmap
        # colorbar for ICE
        cb_ice_ax = fig1.add_axes([bbox_ax1.x0, bbox_ax1.y0*.33, bbox_ax1.x1-bbox_ax1.x0, 0.02])
        # cb_ice_ax = fig1.add_axes([bbox_ax1.x1+0.14, bbox_ax1.y0, 0.02, bbox_ax1.y1-bbox_ax1.y0])
        cbi = plt.colorbar(cbice, cax=cb_ice_ax, orientation='horizontal')
        cbi.set_label(label=f'NSIDC Ice Concentration on {icedate}',fontsize=9)
        cbi.ax.tick_params(labelsize=9)
    else:
        ax1.text(-160,69,'ICE data unavailable',color='k',fontsize=10,transform=ccrs.PlateCarree())
    
    ax1.plot(shipLoc[0],shipLoc[1],'r*',markersize=10, transform=ccrs.PlateCarree())
    # if zoom:
    #     ax1.set_aspect('auto')
    
    if surface == 'SST':
        if not zoom:
            ax1.set_title(f'SASSIE Temperature data {objdate.month}/{objdate.day}/{objdate.year}',fontsize=16)
            figstr=f'{args.base_dir}/figs/sassie_Temp_{objdate.month}-{objdate.day}-{objdate.year}.png'
        else:
            ax1.set_title(f'SASSIE Temperature data {objdate.month}/{objdate.day}/{objdate.year} insitu',fontsize=16)
            gl = ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            gl.xlocator = mticker.FixedLocator(np.arange(-180,180,1))
            gl.xformatter = LONGITUDE_FORMATTER            
            gl.ylocator = mticker.FixedLocator(np.arange(68,82,0.25))
            gl.yformatter = LATITUDE_FORMATTER
            figstr=f'{args.base_dir}/figs/sassie_Temp_{objdate.month}-{objdate.day}-{objdate.year}Zoom.png'
    if surface == 'SSS':
        if not zoom:
            ax1.set_title(f'SASSIE Salinity data {objdate.month}/{objdate.day}/{objdate.year}',fontsize=16)
            figstr=f'{args.base_dir}/figs/sassie_Sali_{args.satelliteSSS}_{objdate.month}-{objdate.day}-{objdate.year}.png'
        else:
            ax1.set_title(f'SASSIE Salinity data {objdate.month}/{objdate.day}/{objdate.year} insitu',fontsize=16)
            gl = ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False) #, alpha=0.3) crs=ccrs.PlateCarree(), 
            gl.xlocator = mticker.FixedLocator(np.arange(-180,180,1))
            gl.xformatter = LONGITUDE_FORMATTER            
            gl.ylocator = mticker.FixedLocator(np.arange(68,82,0.25))
            gl.yformatter = LATITUDE_FORMATTER
            figstr=f'{args.base_dir}/figs/sassie_Sali_{args.satelliteSSS}_{objdate.month}-{objdate.day}-{objdate.year}Zoom.png'
    # establish labels dictionary
    newline = '\n'

    outlabs = {0: {'name':'Utqiagvik',         'lon':-156.8, 'lat':70.3, 'rot1':0,   'col':'k'},
                # 1: {'name':f'Banks{newline}Island',  'lon':-121,   'lat':73,   'rot1':0,   'col':'k'},
                # 2: {'name':'Inuvik',        'lon':-133.7, 'lat':68.4, 'rot1':0,   'col':'k'},
               # 3: {'name':'85N',       'lon': 162, 'lat':86, 'rot1':-25, 'col':'dimgrey'},
               # 4: {'name':'80N',       'lon': 162, 'lat':81, 'rot1':-25, 'col':'dimgrey'},
               # 5: {'name':'75N',       'lon': 162, 'lat':76, 'rot1':-25, 'col':'dimgrey'},
               # 6: {'name':'90E',       'lon':  90, 'lat':84, 'rot1':0,   'col':'dimgrey'},
               # 7: {'name':'135E',      'lon': 141, 'lat':84, 'rot1':45,  'col':'dimgrey'},
               # 8: {'name':'180E',      'lon': 188, 'lat':83, 'rot1':90,  'col':'dimgrey'},
               # 9: {'name':'135W',      'lon':-124, 'lat':82, 'rot1':-45, 'col':'dimgrey'},
               # 10:{'name':'90W',       'lon': -90, 'lat':82, 'rot1':0,   'col':'dimgrey'}
              }    
    for i,o in enumerate(outlabs):
        if i <= 2: fs=14
        else: fs=11
        ax1.text(outlabs[i]['lon'],outlabs[i]['lat'],outlabs[i]['name'],rotation=outlabs[i]['rot1'],
                color=outlabs[i]['col'],fontsize=fs,transform=ccrs.PlateCarree())

    return ax1, fig1, figstr

# def smapBeaufort(strdate=None,level='L3',extent=None):
#     if strdate == None:
#         objdate = dt.datetime.now() - dt.timedelta(days=1)
#         strdate = "%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
#     else:
#         objdate = dt.datetime.strptime(strdate,'%Y%m%d')
    
#     sssdate, ds = pfields.getSSS(strdate)  # already been trimmed lat/lon
#     print(ds.smap_sss.shape)
#     print(ds.longitude.shape)
#     print(ds.latitude.shape)
#     exit(-1)
#     # plot
# map_proj = ccrs.NorthPolarStereo(central_longitude=-150)
# fig = plt.figure(figsize=(5,5))
# ax1 = fig.add_subplot(projection=map_proj)


# sss_L3_img = ax1.pcolormesh(ds_smap_L3.longitude, ds_smap_L3.latitude, ds_smap_L3.smap_sss, 
#                         vmin=20, vmax=30,  # Set max and min values for plotting
#                         cmap='viridis', shading='auto',   # shading='auto' to avoid warning
#                         transform=ccrs.PlateCarree())  # coords are lat,lon but map if NPS 


# #  ----- map stuff
# ax1.coastlines(color='none')  # coastline
# ax1.set_extent([-170, -130,68, 80], crs=ccrs.PlateCarree())
# gl = ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha=0.3)
# ax1.set_title('SMAP L3 SSS')
# plt.colorbar(sss_L3_img, ax=ax1, orientation="horizontal", pad=0.05).set_label('psu')
# # land color
# ax1.add_feature(cartopy.feature.LAND , facecolor=(.7,.7,.7))
# ax1.coastlines('10m')
# ax1.add_feature(cartopy.feature.RIVERS)
# # ticks
# gl.ylocator = mticker.FixedLocator([68, 70, 72, 74, 76, 78])
# gl.xlocator = mticker.FixedLocator([-170, -160, -150, -140, -130])
# gl.top_labels = False
# gl.bottom_labels = True
    
#     pass




