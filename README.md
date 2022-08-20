## Instructions for set up and how to run SASSIE insitu data code. 

Main code is plotSuite.py, which collects, plots and writes data to .csv files.

Clone this repository. 

Edit the input_args_SASSIE.txt file, changing the --base_dir to your directory. 
You’ll need subdirectories under the base_dir: 
BuoyData/, SatelliteFields/, csv/, figs/, pyfiles/, swift_telemetry/, waveGlider/
	under SatelliteFields/ you’ll need Bremen_SIC/, NOAA_SST/, JPL_SMAP/
	
WaveGlider server relies on the time zone of the local computer. Change the time zone
on your computer to UTC (British Summer) for the duration of the SASSIE experiment.

A cron job runs jobs at specified times. The cron command has five arguments at the beginning:
1. which minutes to run on    * means all
2. which hours
3  which days
4. which months
5. which days of week

Cron jobs are run at local times.

examples:
0 */3 * * *
would 3:00, 6:00, 9:00, ... every day, every month, every day of weel

For debugging and testing I would use something like
*/5 * * * *        which runs every 5 mins on the 5s

Steps to setup a cron job. We will set it up to run the command every three hours. 

In your terminal window:

1. type crontab -e		(opens a vi editor window)
2. type i	    		  (changes editor to ‘insert’ or edit mode)
3. copy/paste your version of the following line into the vi editor window:
0 */3 * * * /Users/suzanne/opt/miniconda3/bin/python /Users/suzanne/SASSIE/pyfiles/plotSuite.py @/Users/suzanne/SASSIE/pyfiles/input_args_SASSIE.txt 2>&1 | mail -s "plotSuite" sdickins@uw.edu

4. press escape key	(gets the editor out of edit mode)
5. type :wq		(saves the file)
6. press enter		(closes the vi window)

To stop a cron job, remove it with crontab -r
To view a cron job, list it with crontab -l


 More explanation of crontab line:
* '0'	run job first minute of the hour
* '*/3'	run job every third hour
* '*'	run every day
* '*'	run every month
* '*'	run every day of the week
* the path to python	(can find with ‘which python’ in terminal window)
* the command, here, the python code and the arguments file (which follows an @)
* 2>&1	this has something to do with debug comments, I think. 
* | mail 	(pipe to mail), sends an email after running
* -s “plotSuite”	let’s you put in whatever subject you want for the email
* lastly 		your email address.


