#!/usr/bin/python
import numpy as np
import pandas as pd
import os, shutil, glob
import datetime as dt
import netCDF4 as nc
import urllib
from pathlib import Path
import xarray as xr
import h5netcdf # needed for hdf5 format netcdf files

from urllib.parse import urlparse
from urllib.request import urlopen, Request, build_opener, HTTPCookieProcessor
from urllib.error import HTTPError, URLError
import requests
from nsidc_download_0081_v02 import nsidc_download_0081_v02
import UpTempO_HeaderCodes as HC
import UpTempO_BuoyMaster as BM
import matlab.engine
import scipy.io as sio
from itertools import chain
import matplotlib.pyplot as plt # ver 3.5
from itertools import chain
import subprocess
import json
from scipy import interpolate

#======= Get ICE concentration map from Bremen ==============
def getICE(args,nors='n'):

    args.strdate = None
    if args.strdate is None:
        objdate = dt.datetime.now() - dt.timedelta(days=1)
        args.strdate = "%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
    else:
        objdate = dt.datetime.strptime(args.strdate,'%Y%m%d')

    year="%d" % (objdate.year)   # this info for url sub directories
    noFile = True
    numDaysBack = 0
    strdate = args.strdate

    while noFile:
        if numDaysBack ==2:
            ice,icexx,iceyy = None,None,None
            break

        thefile=f'asi-AMSR2-n6250-{strdate}-v5.4.nc'
        theurl=(f'https://seaice.uni-bremen.de/data/amsr2/asi_daygrid_swath/n6250/netcdf/{year}/{thefile}')

        ncpath=f'{args.base_dir}/SatelliteFields/Bremen_SIC/{thefile}'

        request = urllib.request.Request(theurl)
        request.get_method = lambda: 'HEAD'
        urllib.request.urlopen(request)
        urllib.request.urlretrieve(theurl,ncpath)

        if os.path.exists(ncpath):
            ncdata=nc.Dataset(ncpath)
            ice=np.squeeze(ncdata['z'])
            if np.nanmax(ice) >= 100.:
                ice /= 100.
            ice[ice==0.] = np.nan
            # sst[sst<-900] = np.nan  # set invalids=-999 to nan
            y=ncdata['y'][:]
            x=ncdata['x'][:]
            icexx, iceyy = np.meshgrid(x,y)
            noFile = False
        else:
            objdate = objdate - dt.timedelta(days=1)
            strdate = "%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
            numDaysBack += 1


    return strdate,ice,icexx,iceyy


#======= Get NSIDC ICE concentration map ============== NOT USING THIS PRODUCT
def getICE_nsidc(args,nors='n'):

    # hardwiring 'today' to check working with 0081 product
    args.strdate = None
    if args.strdate is None:
        objdate = dt.datetime.now() - dt.timedelta(days=1)
        args.strdate = "%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
    else:
        objdate = dt.datetime.strptime(args.strdate,'%Y%m%d')

    print('line 35 in pfields: ICE', args.satelliteICE, args.strdate)
    #
    if args.satelliteICE == 'nsidc-0081' : #objdate >= dt.datetime(2021, 11, 1):  # get nsidc-0081
        noFile = True
        numDaysBack = 0

        while noFile:
            print('ICE',args.strdate)
            if numDaysBack == 3:
                ice, icexx, iceyy = None, None, None
                break
            nsidc_download_0081_v02(args.strdate)
            icepath = f'{args.base_dir}/SatelliteFields/NSIDC_ICE/{args.satelliteICE}'
            icefile = f'NSIDC0081_SEAICE_PS_N25km_{args.strdate}_v2.0.nc'

            if os.path.exists(f'{icepath}/{icefile}'):
                print('Ice map avaiable for', args.strdate)
                print()
                ncdata=nc.Dataset(f'{icepath}/{icefile}')
                ice=np.squeeze(ncdata['F18_ICECON'])  # np.squeeze converts nc dataset to np array
                ice[ice==251] = 1.  # fill pole_hole_mask, flags are not scaled
                ice[ice==254] = np.nan  # land
                ice[ice==253] = np.nan  # coast
                ice[ice==0] = np.nan  # no ice
                ice /= 100.
                print('line 114 in pfields',np.nanmax(ice.ravel()))
                y=ncdata['y'][:]
                x=ncdata['x'][:]
                icexx, iceyy = np.meshgrid(x,y)
                noFile = False
            else:
                objdate = objdate - dt.timedelta(days=1)
                strdate = "%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
                numDaysBack += 1
    elif args.satelliteICE == 'g02202': #objdate<=dt.datetime(2021,5,31):   # get G02022
        icepath = f'{args.base_dir}/SatelliteFields/NSIDC_ICE/{args.satelliteICE}'
        icefile = f'seaice_conc_daily_nh_{objdate.year}{objdate.month:02}{objdate.day:02}_f17_v04r00.nc'
        print(f'{icepath}/{icefile}')
        if os.path.exists(f'{icepath}/{icefile}'):
            ncdata=nc.Dataset(f'{icepath}/{icefile}')
            ice=np.squeeze(ncdata['cdr_seaice_conc'])  # np.squeeze converts nc dataset to np array
            ice[ice==251] = 1.  # fill pole_hole_mask, flags are not scaled
            ice[ice==254] = np.nan  # land
            ice[ice==253] = np.nan  # coast
            ice[ice==0] = np.nan  # no ice
            y=ncdata['ygrid'][:]
            x=ncdata['xgrid'][:]
            icexx, iceyy = np.meshgrid(x,y)
    else:
        ice, icexx, iceyy = None, None, None

    return args.strdate,ice,icexx,iceyy


#========= Get SST satellite map ==================
def getSST(args):

    # if args.strdate is None:
    objdate = dt.datetime.now() - dt.timedelta(days=1)
    strdate="%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
    # else:
    #     objdate = dt.datetime.strptime(strdate,'%Y%m%d')
    year="%d" % (objdate.year)   # this info for url sub directories
    month="%.2d" % (objdate.month)
    print(strdate,'in get sst',year,month)
    # SST is available ~9am Pacific

    noFile = True
    numDaysBack = 0

    while noFile:
        if numDaysBack == 2:
            sst,sstlon,sstlat = None,None,None
            break
        thefile='oisst-avhrr-v02r01.'+strdate+'_preliminary.nc'
        print(thefile)
        theurl=('https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/'+year+month+'/'+thefile)
        print('theurl',theurl)
        print()

        ncpath=f'{args.base_dir}/SatelliteFields/NOAA_SST/{thefile}'

        # if not os.path.exists(theurl):
        #     theurl = theurl.replace('_preliminary','')
        # print('theurl',theurl)

        request = urllib.request.Request(theurl)
        request.get_method = lambda: 'HEAD'
        try:
            urllib.request.urlopen(request)
            urllib.request.urlretrieve(theurl,ncpath)

            if os.path.exists(ncpath):
                ncdata=nc.Dataset(ncpath)
                sst=np.squeeze(ncdata['sst'])
                sst[sst<-900] = np.nan  # set invalids=-999 to nan
                lat=ncdata['lat'][:]
                lon=ncdata['lon'][:]
                sstlon, sstlat = np.meshgrid(lon,lat)
                break
        except:
            objdate = objdate - dt.timedelta(days=1)
            strdate = "%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
            numDaysBack += 1

    return objdate, sst, sstlon, sstlat

#========= Get SSS satellite map ================== ONLY USING JPL-L2B-del
def datasetInfo(src,date_range,bounding_box):

    src = src.strip()
    datadict = {'JPL-L3':      {'contentID':'C2208422957-POCLOUD'},
                'RSS-L3':      {'contentID':'C1940468263-POCLOUD'},
                'JPL-L2B-NRT': {'contentID':'C2208418228-POCLOUD'},  # v5.0
                'JPL-L2B-del': {'contentID':'C2208420167-POCLOUD'},  # https://podaac.jpl.nasa.gov/dataset/SMAP_JPL_L2B_SSS_CAP_V5
                }

    granule_url = 'https://cmr.earthdata.nasa.gov/search/granules'

    # use "requests" to find the granules in our date range & domain
    response = requests.get(granule_url,
                            params={
                                'concept_id': datadict[src]['contentID'],
                                'temporal': date_range,
                                'bounding_box': bounding_box,
                                'page_size': 200,
                                },
                            headers={
                               'Accept': 'application/json'
                               }
                      )
    print(response.headers['CMR-Hits'])  # number of granules identified. one day gives 9 granules because it's 8-day running-mean data.
    return response

def getSSS(args):

    # args.strdate = '20220510'
    extent=[int(item.strip()) for item in args.mapDomain.split(',')]
    bounding_box = (f'{extent[0]},{extent[2]},{extent[1]},{extent[3]}')

    if args.strdate is None:
        objdate = dt.datetime.now() - dt.timedelta(days=1)
        args.strdate="%d%.2d%.2d" % (objdate.year,objdate.month,objdate.day)
    else:
        objdate = dt.datetime.strptime(args.strdate,'%Y%m%d')

    savedir = f'{args.base_dir}/SatelliteFields/JPL_SMAP'
    os.makedirs(savedir, exist_ok=True)

    start_date = f'{objdate.year}-{objdate.month:02d}-{objdate.day:02}T00:00:00' #needed for download
    end_date = f'{objdate.year}-{objdate.month:02d}-{objdate.day:02}T23:59:59'
    # convert our start/end dates to datetime format:
    start_date_dt = dt.datetime.strptime(start_date,'%Y-%m-%dT%H:%M:%S') # needed for going back a day
    end_date_dt = dt.datetime.strptime(end_date,'%Y-%m-%dT%H:%M:%S')
    sassiedate = args.strdate

    noFile = True
    numDaysBack = 0
    maxDaysBack = 5

    while noFile:
        if numDaysBack == maxDaysBack:
            sassiedate, ds_smap = start_date, None
            break

        date_range =  start_date + "/" + end_date
        response = datasetInfo('JPL-L2B-del',date_range,bounding_box)

        # download the granules
        granules = response.json()['feed']['entry']

        # => need this to authenticate with Drive
        from requests.auth import HTTPBasicAuth
        granules = response.json()['feed']['entry']

        if bool(granules) == False:   #empty list
            start_date_dt -= dt.timedelta(days=1)
            end_date_dt   -= dt.timedelta(days=1)
            start_date = f'{start_date_dt.year}-{start_date_dt.month:02d}-{start_date_dt.day:02}T00:00:00'
            end_date   = f'{end_date_dt.year}-{end_date_dt.month:02d}-{end_date_dt.day:02}T23:59:59'
            sassiedate = f'{start_date_dt.year}{start_date_dt.month:02d}{start_date_dt.day:02}'
            numDaysBack += 1
            if numDaysBack == maxDaysBack:
                break

        else:
            for granule in granules:
                # loop through the granule links to find which one is the actual file:
                # - identidfied as having 'href' end with '.h5'
                for gl in granule['links']:
                    url = gl['href']
                    if url.endswith('.h5'):
                        # this is the file to download!
                        # the actual filename (without the full path) is buried in the granule link 'title' and preceded with the word 'Download ':
                        # (could also use   filename = granule['title'] + '.h5' )
                        filename = gl['title'].strip('Download ')
                        break
                # use 'requests' to get the file:
                #   Password here is the 'Drive API Password' https://podaac-tools.jpl.nasa.gov/drive/
                #   user here is the Earthdata Login username
                r = requests.get(url, auth=HTTPBasicAuth('YOUR_USERNAME', 'YOUR_PASSWORD'))
                outfile = Path(savedir, filename)
                if not outfile.exists():
                    # print(url, filename)
                    with open(outfile, 'wb') as f:
                        f.write(r.content)
                else:
                    print(filename,'already exists.')
                noFile = False
        # load each file, select data in our domain, and combine into one array called smap_L2
        # this is very klugey! No doubt many cleaner methods exist, but this seems to work OK
        # - code modified from https://stackoverflow.com/questions/59288473/how-do-i-combine-multiple-datasets-h5-files-with-different-dimensions-sizes-i

        # create a new xarray dataset to store the extracted data along the dimension 'points'
        ds_smap = xr.Dataset(
            dict(
                longitude = xr.DataArray([],dims='points'),
                latitude = xr.DataArray([],dims='points'),
                smap_sss = xr.DataArray([],dims='points') #quality_flag = xr.DataArray([],dims='points')
            )
        )

        # the files are listed in "granules" - run through each granule, load, and extract data in our domain
        for granule in granules:
            # local file:
            filename = f'{savedir}/{granule["title"]}.h5'
            # only try to open the local file if it actually exists:
            if Path(filename).exists():
                # load the file into xarray
                ds = xr.open_mfdataset(filename)  # arrays are 76 x ...
                # extract lon, lat, and sss; reshape so they can later be combined more easily
                lon = np.reshape(ds.lon.values,[1,-1])
                lat = np.reshape(ds.lat.values,[1,-1])
                sss = np.reshape(ds.smap_sss.values,[1,-1])
                qflag = np.reshape(ds.quality_flag.values,[1,-1])
                lon=lon[~np.isnan(sss)]
                lat=lat[~np.isnan(sss)]
                qflag=qflag[~np.isnan(sss)]
                sss=sss[~np.isnan(sss)]
                for ii, qf in enumerate(qflag):
                    qbits=bin(np.int(qf))[2:].zfill(16)[::-1]
                    if qbits[5]==1 or qbits[7]==1 or qbits[8]==1:
                        lon[ii]=np.nan
                        lat[ii]=np.nan
                        sss[ii]=np.nan
                        qflag[ii]=np.nan

                lon=lon[~np.isnan(sss)]
                lat=lat[~np.isnan(sss)]
                qflag=qflag[~np.isnan(sss)]
                sss=sss[~np.isnan(sss)]
                # # index of points in our domain:
                # i = ((lon > extent[0]) & (lon < extent[1]) & (lat > extent[2]) & (lat < extent[3]))
                # sss = sss[i]
                # lon = lon[i]
                # lat = lat[i]

                if len(sss)>0:
                    # create a temperary dataset for this file containing
                    this_file = xr.Dataset(
                        dict(
                            longitude = xr.DataArray(lon,dims='points'),
                            latitude = xr.DataArray(lat,dims='points'),
                            smap_sss = xr.DataArray(sss,dims='points')
                        )
                    )
                    # add data from this file to 'smap_L2'
                    ds_smap = xr.concat([ds_smap,this_file], dim='points')
    # print(ds_smap.keys())
    # fig,ax = plt.subplots(1,1)
    # ax.plot(ds_smap['longitude'],ds_smap['latitude'],'.')
    # ax.set_title('NRT SSS locations')
    # plt.show()
    # exit(-1)
    return start_date_dt, ds_smap

# ========= define username and password for apl.pacificgyre.com ==================
def userCred(user):

    psswd={'pscapluw':'microstar'}
    return psswd[user]

def getPGbuoy(args,bid,user,strdate=None):

    binf = BM.BuoyMaster(bid)
    psswd=userCred(user)

    today=dt.datetime.now()
    # enddate="%.2d/%.2d/%d" % (today.month,today.day,today.year)
    downloaddate = f'{today.year}{today.month:02}{today.day:02}'
    # downloaddate =f'{strdate[:4]}{strdate[4:6]}{strdate[6:]}'
    enddate = "%.2d/%.2d/%d" % (today.month,today.day,today.year)
    # objdate = dt.datetime.strptime(strdate,'%Y%m%d')
    twoDayPrev = today - dt.timedelta(hours=args.hourstoPlot) #dt.timedelta(days=800)
    print(twoDayPrev)
    startdate = "%.2d/%.2d/%d" % (twoDayPrev.month,twoDayPrev.day,twoDayPrev.year)
    print('start',startdate)
    print('end',enddate)

    strcommand='http://api.pacificgyre.com/api2/getData.aspx?userName='
    strcommand+=user+'&password='+psswd+'&startDate='+startdate
    strcommand+='&endDate='+enddate+'&commIDs='+bid+'&fileFormat=CSV'
    fid=urllib.request.urlopen(strcommand)  # downloads file with minimal header if bid file not unavailable
    data=fid.read()
    fid.close()
    data=str(data,'utf-8')
    opw=open(f'{args.base_dir}/BuoyData/UTO_{bid}_{downloaddate}.csv','w')
    opw.write(data)
    opw.close()
    df = pd.read_csv(f'{args.base_dir}/BuoyData/UTO_{bid}_{downloaddate}.csv',parse_dates=['DeviceDateTime'])

    if len(df)==0:  # downloaded an empty dataframe, no data for the bids
        df = pd.DataFrame(columns=['Date','Lat','Lon','Temperature','Salinity'])
    else:
        # pick out the pertinent columns
        sas = ['date','lat','lon','sst','temp','sal','press','depth']
        sascol=[s for s in df.columns.values for x in sas if x in s.lower()]
        df = df[sascol]

        # rename columns so they are similar between buoys
        sascolumns = HC.PG_HeaderCodes(sascol)
        df.rename(columns=(dict(zip(sascol,sascolumns.keys()))),inplace=True)
        # add date object column for sorting, we don't need this because we parsed dates above ?
        # df['Date'] = pd.to_datetime(df['Date'],format='%m/%d/%y %H:%M')
        df.sort_values(by='Date', inplace=True)
        df.reset_index(inplace=True)

        # find the shallowest temperature measurement
        try:
            tdepths = binf['tdepths']
            tcols = [col for col in df.columns if col.startswith('T') and col != 'Thull']
            tlist = list(zip(tdepths,tcols))
        except:
            tlist = []
        try:
            CTDtdepths = binf['CTDtdepths']
            tcols = [col for col in df.columns if col.startswith('CTD-T')]
            tlist.extend(list(zip(CTDtdepths,tcols)))
        except:
            tlist.extend([])
        try:
            HULLtdepths = binf['HULLtdepths']
            tcols = [col for col in df.columns if col.startswith('Thull')]
            tlist.extend(list(zip(HULLtdepths,tcols)))
        except:
            tlist.extend([])
        tlist.sort(key=lambda y: y[0])
        print('tlist line 439',tlist)
        # make sure there are valid data in the col
        for ii,item in enumerate(tlist):
            print(item[1])
            if not df[item[1]].isnull().all():
                tcol = item[1]
                print('line 445 tcol',tcol)
                break
        # exit(-1)
        # find the shallowest salinity measurement,     WITH DATA IN IT
        try:
            sdepths = binf['sdepths']
            scols = [col for col in df.columns if col.startswith('S') and col != 'Shull']
            slist = list(zip(sdepths,scols))
        except:
            slist=[]
        try:
            CTDsdepths = binf['CTDsdepths']
            scols = [col for col in df.columns if col.startswith('CTD-S')]
            slist.extend(list(zip(CTDsdepths,scols)))
        except:
            slist.extend([])
        try:
            HULLsdepths = binf['HULLsdepths']
            scols = [col for col in df.columns if col.startswith('Shull')]
            slist.extend(list(zip(HULLsdepths,scols)))
        except:
            slist.extend([])
        slist.sort(key=lambda y: y[0])
        print('slist line 469',slist)

        # keep only the shallowest, where the salinity measurement must be shallower than 5m
        if not slist:
            df = df[['index','Date','Lat','Lon',tcol]]
        else:
            for ii,item in enumerate(slist):
                print(item,item[0],item[1])
                if not df[item[1]].isnull().all():
                    scol = item[1]
                    print('scol',scol)
                    break
            if  slist[ii][0] < 5:
                df = df[['index','Date','Lat','Lon',tcol,scol]]
            else:
                df = df[['index','Date','Lat','Lon',tcol]]
        # formatting
        df['Date'] = df['Date'].dt.round('1s')
        df['Lat'] = df['Lat'].map('{:.03f}'.format).astype(float)
        df['Lon'] = df['Lon'].map('{:.03f}'.format).astype(float)
        df[tcol] = df[tcol].map('{:.03f}'.format).astype(float)
        if slist:
            df[scol] = df[scol].map('{:.03f}'.format).astype(float)

    return df

def getSWIFT(args,ID,starttime,endtime,eng):
    # eng = matlab.engine.start_matlab()
    # eng.addpath('../swift_telemetry')
    swiftpath = f"'{args.base_dir}/swift_telemetry'"  #need single quotes for space in Google Drive dir name
    print('line 625 in pfields, swiftpath',swiftpath,ID,starttime,endtime)

    allbatterylevels, lasttime, lastlat, lastlon = eng.pullSWIFTtelemetry(swiftpath,ID,starttime,endtime,nargout=4)
    swiftpath = swiftpath.strip("'")
    print('swiftpath',swiftpath,allbatterylevels, lasttime, lastlat, lastlon)

    if endtime=='':
        swiftfile = f'buoy-SWIFT {ID}-start-{starttime}-end-None.mat'
    filecsv = f'{swiftpath}/csv/{os.path.splitext(swiftfile)[0]}.csv'
    try:
        swift_struct = sio.loadmat(f'{swiftpath}/{swiftfile}')
        SWIFT = swift_struct['SWIFT']

        buoyid = np.array([item for item in chain(*SWIFT[0,:]['ID'])])
        buoyid = np.array([item.split(' ')[-1] for item in buoyid])
        # buoydate = np.array([jtem for jtem in chain(*[item.tolist() for item in chain(*SWIFT[0,:]['date'])])])
        # print('min buoydate',np.min(np.int(buoydate)))

        CTdepth = np.array([jtem for jtem in chain(*[item.tolist() for item in chain(*SWIFT[0,:]['CTdepth'])])])
        # CTdepth[buoydate==0]
        unique, counts = np.unique(CTdepth, return_counts=True)
        baddepth = unique[counts<10]
        ndepths = len(unique[counts>10])

        # variables come in as list of arrays of nested lists, where each array has one value.
        lon = np.array([jtem for jtem in chain(*[item.tolist() for item in chain(*SWIFT[0,:]['lon'])])])
        lat = np.array([jtem for jtem in chain(*[item.tolist() for item in chain(*SWIFT[0,:]['lat'])])])

        salinity = np.array([jtem for jtem in chain(*[item.tolist() for item in chain(*SWIFT[0,:]['salinity'])])])
        salinity = salinity[CTdepth!=baddepth]
        salinity = np.reshape(salinity,(-1,ndepths))
        salinity[salinity<22] = np.nan

        watertemp = np.array([jtem for jtem in chain(*[item.tolist() for item in chain(*SWIFT[0,:]['watertemp'])])])
        watertemp = watertemp[CTdepth!=baddepth]
        watertemp = np.reshape(watertemp,(-1,ndepths))

        CTdepth = CTdepth[CTdepth!=baddepth]
        CTdepth = np.reshape(CTdepth,(-1,ndepths))

        time = np.array([jtem for jtem in chain(*[item.tolist() for item in chain(*SWIFT[0,:]['time'])])])
        # # # to check in matlab datevec(SWIFT(end).time)
        datetime = [dt.datetime.fromordinal(int(t)) + dt.timedelta(t%1) - dt.timedelta(days=366) for t in time]

        columns = ['DateTime','Lat','Lon','BuoyID']
        [columns.append(f'CTdepth-{ii}') for ii in range(ndepths)]
        [columns.append(f'Salinity-{ii}') for ii in range(ndepths)]
        [columns.append(f'WaterTemp-{ii}') for ii in range(ndepths)]
        # print(columns)

        # create dataFrame
        dfSwift = pd.DataFrame(columns=columns)

        # fill dataFrame
        dfSwift['DateTime'] = datetime
        dfSwift['DateTime'] = dfSwift['DateTime'].dt.round('1s')
        # dfSwift['DateTimeStr'] = dfSwift['DateTime'].dt.strftime('%d/%m/%Y %H:%M:%S')
        dfSwift['Lat'] = lat
        dfSwift['Lon'] = lon
        dfSwift['BuoyID'] = buoyid
        # reduce to three decimals
        dfSwift['Lat'] = dfSwift['Lat'].map('{:.03f}'.format).astype(float)
        dfSwift['Lon'] = dfSwift['Lon'].map('{:.03f}'.format).astype(float)
        for ii in range(ndepths):
            dfSwift[f'CTdepth-{ii}'] = CTdepth[:,ii]
            dfSwift[f'Salinity-{ii}'] = salinity[:,ii]
            dfSwift[f'WaterTemp-{ii}'] = watertemp[:,ii]
            # reduce to three decimals
            dfSwift[f'Salinity-{ii}'] = dfSwift[f'Salinity-{ii}'].map('{:.03f}'.format).astype(float)
            dfSwift[f'WaterTemp-{ii}'] = dfSwift[f'WaterTemp-{ii}'].map('{:.03f}'.format).astype(float)
    except:
        dfSwift = pd.DataFrame(columns=['DateTime','Lat','Lon','WaterTemp-0','Salinity-0'])

    return dfSwift

def getWaveGlider(args,ID):
    # login information
    wsdl            = "https://dataservice.wgms.com/WDS/WGMSData.asmx?wsdl" # Data Service WSDL URL
    org             = "apl-uw"                                               # WGMS Org
    user            = "sdickins"                                      # WGMS Org user
    password        = "$ASSIE"                                          # WGMS Org Password
    resultFormat    = 1                                                     # Default Option

    # date time strings
    now = dt.datetime.now()
    then = now - dt.timedelta(hours=args.hourstoPlot) # 48)
    endDate = now.strftime('%Y-%m-%dT%H:%M:%S')
    startDate = then.strftime('%Y-%m-%dT%H:%M:%S')

    for reportName in ['"AanderraaCT Sensor"','"GPS Waves Sensor Data"']:
        print(reportName)
        print(ID)

        cmd = f'python {args.base_dir}/pyfiles/DataService.py --getReportData --vehicles {np.int(ID)} --startDate {startDate}Z --endDate {endDate}Z --reportName {reportName}'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # waits until p ends and saves output and errors if needed
        out, err = p.communicate()  # both strings

        if str(err, 'utf-8') != '':
            print(f'Error in WaveGlider data retrieval, {dt.datetime.now()} Z: {str(err, "utf-8")}')
            dfWaveGlider = None
            break

        if "AanderraaCT Sensor" in reportName:
            dataAS = list(eval(out))   # a string that is a list of dicts
            dfAS = pd.DataFrame([item['timeStamp'],item['temperature (C)'],item['salinity (PSU)']] for item in dataAS)
            dfAS.columns=['TimeStamp','Temperature','Salinity']
            dfAS['DateTime'] = [dt.datetime.strptime(item,'%Y-%m-%dT%H:%M:%S') for item in dfAS['TimeStamp']]
            # print(dfAS.head())
        elif "Telemetry 6 Report" in reportName:
            dataT6 = json.loads(out)   # a string that json converts to dictionary
            dfT6 = pd.DataFrame([item['gliderTimeStamp'],item['latitude'],item['longitude']] for item in dataT6['recordData'])
            dfT6.columns=['GliderTimeStamp','Latitude','Longitude']
            dfT6['GliderDateTime'] = [dt.datetime.strptime(item,'%Y-%m-%dT%H:%M:%S') for item in dfT6['GliderTimeStamp']]
            # print(dfT6.head())
        elif "GPS Waves Sensor Data" in reportName:
            dataGPS = list(eval(out))   # a string that is a list of dicts
            dfGPS = pd.DataFrame([item['timeStamp'],item['latitude'],item['longitude']] for item in dataGPS)
            dfGPS.columns=['GPSTimeStamp','Latitude','Longitude']
            dfGPS['GPSDateTime'] = [dt.datetime.strptime(item,'%Y-%m-%dT%H:%M:%S') for item in dfGPS['GPSTimeStamp']]
            dfGPS['Latitude'] = dfGPS['Latitude'].astype(float)   # comes as type object, eye roll
            dfGPS['Longitude'] = dfGPS['Longitude'].astype(float)
            # print(dfGPS.head())

    # convert panda series to something scipy.interpolate can use
    secondsSinceAS = (dfAS['DateTime'] - dt.datetime(2022,8,1)).to_numpy() / np.timedelta64(1,'s')  #1e9)  # seconds since
    secondsSinceGPS = (dfGPS['GPSDateTime'] - dt.datetime(2022,8,1)).to_numpy() / np.timedelta64(1,'s')

    fi = interpolate.interp1d(secondsSinceAS, dfAS['Temperature'], fill_value='extrapolate')
    dfGPS['Temperature'] = fi(secondsSinceGPS)
    fi = interpolate.interp1d(secondsSinceAS, dfAS['Salinity'], fill_value='extrapolate')
    dfGPS['Salinity'] = fi(secondsSinceGPS)

    # print(dfGPS['GPSDateTime'].dtype)
    # print(dfGPS['Latitude'].dtype)
    # print(dfGPS['Longitude'].dtype)
    # print(dfGPS['Temperature'].dtype)
    # print(dfGPS['Salinity'].dtype)

    dfWaveGlider = pd.DataFrame()
    # dfWaveGlider['index']
    dfWaveGlider['Date'] = dfGPS['GPSTimeStamp']
    dfWaveGlider['Lat'] = dfGPS['Latitude'].map('{:.03f}'.format).astype(float)
    dfWaveGlider['Lon'] = dfGPS['Longitude'].map('{:.03f}'.format).astype(float)
    dfWaveGlider['Temperature'] = dfGPS['Temperature'].map('{:.03f}'.format).astype(float)
    dfWaveGlider['Salinity'] = dfGPS['Salinity'].map('{:.03f}'.format).astype(float)
    # print(dfWaveGlider.head())

    # fig,ax = plt.subplots(1,1,figsize=(6,3))
    # ax.plot(dfAS['DateTime'],dfAS['Temperature'],'r.')
    # ax.plot(dfGPS['GPSDateTime'],dfGPS['Temperature'],'b.')
    # ax.set_title('ACT temperature (b), interpolated to GPS (r)')
    #
    # fig,ax = plt.subplots(1,1,figsize=(6,3))
    # ax.plot(dfAS['DateTime'],dfAS['Salinity'],'r.')
    # ax.plot(dfGPS['GPSDateTime'],dfGPS['Salinity'],'b.')
    # ax.set_title('ACT salinity (b), interpolated to GPS (r)')

    return dfWaveGlider
