% wave attenuation during SASSIE play 1
% 
% J. Thomson, Apr 2023

clear all, close all

datapath =  '/Users/jthomson/Dropbox/Projects/SASSIE/SASSIE_data/preliminary/L1' ;

Hs = [];
time = [];
lat = [];
lon = [];

counter = 0;

SWIFTfiles = dir([ datapath '/SWIFT_L1/SWIFTplay1_L1/*.mat' ] );

for i = 1:length(SWIFTfiles)
    load([ datapath '/SWIFT_L1/SWIFTplay1_L1/' SWIFTfiles(i).name ] )
   Hs(counter + [1:length(SWIFT)] ) = [SWIFT.sigwaveheight];
   time(counter + [1:length(SWIFT)] ) = [SWIFT.time];
   lat(counter + [1:length(SWIFT)] ) = [SWIFT.lat];
   lon(counter + [1:length(SWIFT)] ) = [SWIFT.lon];
   counter = length(Hs);
end

WGfiles = dir([ datapath '/WaveGliders_L1/*.mat' ] );

for i = 1:length(WGfiles)
    load([ datapath '/WaveGliders_L1/' WGfiles(i).name ] )
   Hs(counter + [1:length(SV3)] ) = [SV3.sigwaveheight];
   time(counter + [1:length(SV3)] ) = [SV3.time];
   lat(counter + [1:length(SV3)] ) = [SV3.lat];
   lon(counter + [1:length(SV3)] ) = [SV3.lon];
   counter = length(Hs);
end



%% trim time series for a window with stationarity 
maxtime = datenum(2022,9,10,12,0,0); %datenum(2022,9,12,0,0,0);
mintime = datenum(2022,9,9,21,0,0); %datenum(2022,9,9,21,0,0);

figure(1), clf
plot(time,Hs,'x'), hold on
datetick

timetrim = time > maxtime | time < mintime;

Hs(timetrim) = [];
lat(timetrim) = [];
lon(timetrim) = [];
time(timetrim) = [];


figure(1),
plot(time,Hs,'rx'), hold on
datetick
ylabel('Hs [m]')

%% attenuation 

iceedgelat = 72.4;
x = deg2km(lat - iceedgelat) * 1000;
H0 = nanmean(Hs(x<0));
alpha = regress( ( log(H0) - log(Hs(x>0)) )' , x(x>0)' ) 


figure(2), clf
%plot(lat,Hs,'x')
plot([0 x(x>0)],[H0 H0*exp(-alpha*x(x>0))],'r.','linewidth',2), hold on
plot(x,Hs,'kx','linewidth',2), hold on
xlabel('distance to ice edge, x [m]')
ylabel('Hs [m]')

legend(['\alpha = ' num2str(alpha)])
set(gca,'fontsize',18,'fontweight','demi')
print -dpng SASSIE_play1_waveiceattenuation.png
