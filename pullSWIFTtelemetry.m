% function pullSWIFTtelemetry( IDs, starttime, endtime );
function [ allbatterylevels lasttime lastlat lastlon] = pullSWIFTtelemetry(rootpath, IDs, starttime, endtime );
% pull and compile sbd files from swiftserver for multiple SWIFTs
% use for quasi-realtime situtational awareness during a mission / project
%
%   [ allbatterylevels lasttime lastlat lastlon ] = pullSWIFTtelemetry( IDs, starttime, endtime );
%
%   where inputs are
%   path for swift .m files, .csv dir, .mat files...
%   strings arrays of two-digit SWIFT IDs (e.g., ['16';'17';'18';]
%   startime string (e.g., '2018-09-12T22:00:00')
%   and end time string, which can be empty to pull until present time
%
%   outputs are the minimum battery reading and the lasttime, lastlat, lastlon of telemetry
%
% note that this puts everything in the current working directory
%
%   J. Thomson, 9/2018
%           5/2021 make compatible with microSWIFTs (three digit IDs)

%
%addpath ~/git_repos/SWIFT-codes/GeneralTools
global battery
rootpath

options = weboptions('Timeout', 60);
% starttime
% endtime
% IDs
% pause

IDs = char(IDs);  % do this for python
% starttime = '2018-09-12T22:00:00'
% endtime = [] %'2018-10-12T22:00:00' 
% IDs = ['09';'12';'13';'15';'16';'17']
addpath(rootpath)

try 
  eval(['cd ' rootpath]) % use this with plotSuite.py
catch
  cd(rootpath) % use this with testMatlab.m
end

if  size(IDs,2) == 2 && ischar(IDs), % enforce two digit strings for SWIFT IDs
    
    SWIFTtype = 'v3v4s'
    allbatterylevels = NaN(1,size(IDs,1)); % initialize battery array
    lasttime =  NaN(1,size(IDs,1)); % initialize time array
    lastlat = NaN(1,size(IDs,1)); % initialize
    lastlon = NaN(1,size(IDs,1)); % initialize
    
elseif  size(IDs,2) == 3 && ischar(IDs), % enforce three digit strings for microSWIFT IDs
    
    SWIFTtype = 'micro'
    allbatterylevels = NaN(1,size(IDs,1)); % initialize battery array
    lasttime =  NaN(1,size(IDs,1)); % initialize time array
    lastlat = NaN(1,size(IDs,1)); % initialize
    lastlon = NaN(1,size(IDs,1)); % initialize
    
else
    
    allbatterylevels = NaN(1,size(IDs,1)); % initialize battery array
    lasttime =  NaN(1,size(IDs,1)); % initialize time array
    lastlat = NaN(1,size(IDs,1)); % initialize
    lastlon = NaN(1,size(IDs,1)); % initialize
    
    disp('invalid IDs')
    
  %  return
    
end

%%
% loop thru pulling SBDs for each SWIFT ID

for si=1:size(IDs,1),

    baseurl = 'http://swiftserver.apl.washington.edu/services/buoy?action=get_data&buoy_name='
%     
    
    if SWIFTtype =='v3v4s'
        buoy = ['SWIFT%20' IDs(si,:) ];
    elseif SWIFTtype =='micro'
        buoy = ['microSWIFT%20' IDs(si,:) ];
    else
    end
    search_url = sprintf('%s%s&start=%s&end=%s&format=zip',baseurl,buoy,char(starttime),char(endtime));
 %   filename = sprintf('%s/SWIFT%s.zip',pathData,IDs(si,:));
    filename = sprintf('SWIFT%s.zip',IDs(si,:));
    
    out = websave(filename,search_url,options);
    
  % need to NOT unzip when empty
    eval(['!unzip -q -o SWIFT' IDs(si,:) '.zip ']); %-d ' pathData]);
    delete '*.mat'
    expanded = dir(['*SWIFT ' IDs(si,:) '*'])
    dirFlag = expanded.isdir;

    expanded(dirFlag).name
    if expanded(dirFlag).isdir == true %length(expanded)==1  ~endsWith(expanded(1),'.mat') && 
        cd([expanded(dirFlag).name])
        % run compile SBD, which calls readSWIFT_SBD for each file
        % use a temp file work around the clear all in compile
        save temp si IDs starttime endtime options allbatterylevels lasttime lastlon lastlat SWIFTtype
        % compile is working in buoy
        compileSWIFT_SBDservertelemetry
        load temp
        wd = pwd;
        wdi = find(wd == '/',1,'last');
        wd = wd((wdi+1):end);
        if ~isempty(SWIFT),
            allbatterylevels(si) = battery(end);
            lasttime(si) = max([SWIFT.time]);
            lastlat(si) = SWIFT(end).lat;
            lastlon(si) = SWIFT(end).lon;
        else
            eval(['!rm '  wd '.mat']), % remove file if no results
            allbatterylevels(si) = -999.;
            lasttime(si) = 0.;
            lastlat(si) = -999.;
            lastlon(si) = -999.;
        end
        
        % keep going to next SWIFT
        cd('../')
        
    else
        disp(['multiple expanded directories for SWIFT ' IDs(si,:)])
        allbatterylevels(si) = -999;
        lasttime(si) = 0;
        lastlat(si) = -999;
        lastlon(si) = -999;
    end
    
end
%%
% combine the resulting mat files in the top level directory and make map plots
eval(['!cp *SWIFT*/*SWIFT*start*.mat ./'])
% tempfig = mapSWIFT('watertemp');
% salinityfig = mapSWIFT('salinity');
% wavefig = mapSWIFT('sigwaveheight');

% clc,
% spaces(1:size(IDs,1)) = ' ';
% commas(1:size(IDs,1)) = ',';
% Vs(1:size(IDs,1)) = 'V';
% 
% [IDs commas' spaces' datestr(lasttime) commas' spaces' num2str(lastlat',6) commas' spaces' ...
%     num2str(lastlon',6) commas' spaces' num2str(allbatterylevels',3) Vs']

% %%
% clear
% load 'buoy-SWIFT 17-start-2018-09-12T22:00:00-end-None.mat'
% 
% for ii = 1:length(SWIFT)
%   lon(ii)=SWIFT(ii).lon;
%   lat(ii)=SWIFT(ii).lat;
%   sali(ii,:)=SWIFT(ii).salinity;
%   temp(ii,:)=SWIFT(ii).watertemp;
% end
% 
% whos
% 



















