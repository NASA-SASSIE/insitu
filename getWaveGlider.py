#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

# Built-in/Generic Imports
import os, ast
# import sys
# import time
import subprocess
import csv
import datetime as dt
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

# Constants. TODO: move to a config file
wsdl            = "https://dataservice.wgms.com/WDS/WGMSData.asmx?wsdl" # Data Service WSDL URL
org             = "apl-uw"                                               # WGMS Org
user            = "sdickins"                                      # WGMS Org user
password        = "$ASSIE"                                          # WGMS Org Password
resultFormat    = 1                                                     # Default Option

# arguments for the DataService command
now = dt.datetime.now()
endDate = now.strftime('%Y-%m-%dT%H:%M:%S')
ID = '84929357'

# if want data from beginning
if ID == '102740746' or ID == '84929357':   # pull all the data, only plot the last args.hourstoPlot
    startDate = '2022-08-14T19:00:00'
elif ID == '1628052144' or ID == '511512553':
    startDate = '2022-08-12T22:00:00'
# print('line 589 in pfields',ID,startDate,endDate)

# reportName = '"GPS Waves Sensor Data"'
# reportName = '"Telemetry 6 Report"'
# reportName = '"AanderraaCT Sensor"'

# paragraph to test string conversion to something we can use
# cmd = f'python DataService.py --getReportData --vehicles {vid} --startDate {startDate}Z --endDate {endDate}Z --reportName {reportName}'
# p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# # waits until p ends and saves output and errors if needed
# out, err = p.communicate()
# print(list(eval(out)))
# exit(-1)

for reportName in ['"AanderraaCT Sensor"','"GPS Waves Sensor Data"']:
    print(reportName)

    cmd = f'/Users/suzanne/opt/miniconda3/bin/python /Users/suzanne/SASSIE/pyfiles/DataService.py --getReportData --vehicles {ID} --startDate {startDate}Z --endDate {endDate}Z --reportName {reportName}'
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # waits until p ends and saves output and errors if needed
    out, err = p.communicate()  # both strings

    # if str(err, 'utf-8') != '':
    #     print(time.asctime(time.gmtime()) + ' Z: ' + str(err, 'utf-8'))
    #     exit(-1)
    if "AanderraaCT Sensor" in reportName:
        if str(err, 'utf-8') == '' and str(out,'utf-8') != '':
            dataAS = list(eval(out))   # a string that is a list of dicts
            # get the data out individually for each depth by looking at longitude (longitude=0, depth=0.25, longitude=2, depth=1)
            time0 = [item['timeStamp'] for item in dataAS if item['longitude']=='0']
            time2 = [item['timeStamp'] for item in dataAS if item['longitude']=='2']
            temp0 = [item['temperature (C)'] for item in dataAS if item['longitude']=='0']
            temp2 = [item['temperature (C)'] for item in dataAS if item['longitude']=='2']
            sali0 = [item['salinity (PSU)'] for item in dataAS if item['longitude']=='0']
            sali2 = [item['salinity (PSU)'] for item in dataAS if item['longitude']=='2']
            lon0 = [item['longitude'] for item in dataAS if item['longitude']=='0']
            lon2 = [item['longitude'] for item in dataAS if item['longitude']=='2']

            # put all lists into a dictionary as arrays
            d = {'timeStamp0': np.array(time0),
                 'timeStamp2': np.array(time2),
                 'temperature0': np.array(temp0),
                 'temperature2': np.array(temp2),
                 'salinity0': np.array(sali0),
                 'salinity2': np.array(sali2),
                 'longitude0': np.array(lon0),
                 'longitude2': np.array(lon2)
                 }
            # make a dataframe
            dfAS = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in d.items() ]))
            # if there are invalids in the shallow column (ending in 0), fill with data from deeper column (ending in 2)
            dfAS['timeStamp0'].fillna(dfAS['timeStamp2'],inplace=True)
            dfAS['temperature0'].fillna(dfAS['temperature2'],inplace=True)
            dfAS['salinity0'].fillna(dfAS['salinity2'],inplace=True)
            dfAS['longitude0'].fillna(dfAS['longitude2'],inplace=True)
            # drop 'deeper' columns
            dfAS.drop(columns=[item for item in dfAS.columns if item.endswith('2')],inplace=True)
            # rename the columns
            dfAS.rename(columns={item:f'{item[0].upper()}{item[1:-1]}' for item in dfAS.columns},inplace=True)
            # make a depth column from the longitude
            dfAS.loc[dfAS['Longitude']=='0','Longitude'] = 0.25
            dfAS.loc[dfAS['Longitude']=='2','Longitude'] = 1
            dfAS.rename(columns={'Longitude':'Depth'},inplace=True)
            # make a datetime object col for plotting
            dfAS['DateTime'] = [dt.datetime.strptime(item,'%Y-%m-%dT%H:%M:%S') for item in dfAS['TimeStamp']]
            # sort by date
            dfAS.sort_values(by='DateTime',inplace=True)
            print(dfAS.tail())
        else:
            print(f'Error in AanderraaCT data retrieval, {dt.datetime.now()} Z: {str(err, "utf-8")}')
            print(err)
            dfAS = pd.DataFrame(columns=['TimeStamp','Temperature','Salinity','Depth','DateTime'])

    elif "Telemetry 6 Report" in reportName:
        if str(err, 'utf-8') == '':
            dataT6 = json.loads(out)   # a string that json converts to dictionary
            dfT6 = pd.DataFrame([item['gliderTimeStamp'],item['latitude'],item['longitude']] for item in dataT6['recordData'])
            dfT6.columns=['GliderTimeStamp','Latitude','Longitude']
            dfT6['GliderDateTime'] = [dt.datetime.strptime(item,'%Y-%m-%dT%H:%M:%S') for item in dfT6['GliderTimeStamp']]
        else:
            print(f'Error in Telemetry 6 data retrieval, {dt.datetime.now()} Z: {str(err, "utf-8")}')
            print(err)
            dfT6 = pd.DataFrame(columns=['GliderTimeStamp','Latitude','Longitude'])
    elif "GPS Waves Sensor Data" in reportName:
        if str(err, 'utf-8') == '' and str(out,'utf-8') != '':
            dataGPS = list(eval(out))   # a string that is a list of dicts
            dfGPS = pd.DataFrame([item['timeStamp'],item['latitude'],item['longitude']] for item in dataGPS)
            dfGPS.columns=['GPSTimeStamp','Latitude','Longitude']
            # print(dfGPS['GPSTimeStamp'].tail(20))
            dfGPS['GPSDateTime'] = [dt.datetime.strptime(item,'%Y-%m-%dT%H:%M:%S') for item in dfGPS['GPSTimeStamp']]
            dfGPS['Latitude'] = dfGPS['Latitude'].astype(float)   # comes as type object, eye roll
            dfGPS['Longitude'] = dfGPS['Longitude'].astype(float)
        else:
            print(f'Error in GPS Waves data retrieval, {dt.datetime.now()} Z: {str(err, "utf-8")}')
            print(err)
            dfGPS = pd.DataFrame(columns=['GPSTimeStamp','Latitude','Longitude'])

    # convert panda series to something scipy.interpolate can use
    try:
        secondsSinceAS = (dfAS['DateTime'] - dt.datetime(2022,8,1)).to_numpy() / np.timedelta64(1,'s')  # / 1e9)
    except:
        secondsSinceAS = None
    try:
        secondsSinceGPS = (dfGPS['GPSDateTime'] - dt.datetime(2022,8,1)).to_numpy() / np.timedelta64(1,'s')
    except:
        secondsSinceGPS = None

    if secondsSinceAS is not None and secondsSinceGPS is not None:
        fi = interpolate.interp1d(secondsSinceGPS, dfGPS['Latitude'], fill_value='extrapolate')
        dfAS['Latitude'] = fi(secondsSinceAS)
        fi = interpolate.interp1d(secondsSinceGPS, dfGPS['Longitude'], fill_value='extrapolate')
        dfAS['Longitude'] = fi(secondsSinceAS)
print(dfAS.head())
fig,ax = plt.subplots(1,1,figsize=(6,3))
ax.plot(dfAS['DateTime'],dfAS['Longitude'],'r.')
ax.plot(dfGPS['GPSDateTime'],dfGPS['Longitude'],'b.')
ax.set_title('GPS Lon(b), interpolated to AS (r)')
fig.savefig('/Users/suzanne/SASSIE/waveGlider/iLon.png')

fig,ax = plt.subplots(1,1,figsize=(6,3))
ax.plot(dfAS['DateTime'],dfAS['Latitude'],'r.')
ax.plot(dfGPS['GPSDateTime'],dfGPS['Latitude'],'b.')
ax.set_title('GPS Lat (b), interpolated to AS (r)')
fig.savefig('/Users/suzanne/SASSIE/waveGlider/iLat.png')
# plt.show()
# exit(-1)
