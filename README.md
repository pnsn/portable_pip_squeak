# portable_pip_squeak
Light weight station quality assessment tool to tell if site is good for 3,4, or 6 channel site

# Overview

This is similar to <a href="https://github.com/pnsn/station_metrics">station_metrics</a>, but is tuned for a quick and dirty assessment of an existing or candidate site.  It is similar to eew_stationreport in that project which is informally known as "pip_squeak".  eew_stationreport is used to assess new seismic stations before acceptance into ShakeAlert.  portable_pip_squeak is similar in that uses many of the same metrics and thresholds, but includes others and the output is intentionally subjective for use by field techs.

# Installation

Built for python 3 with the following major package dependencies:
<a href="http://www.numpy.org/">numpy</a>, <a href="https://github.com/obspy/obspy/wiki">obspy</a>, <a href="https://matplotlib.org">matplotlib</a>.

These can easily be added via command line using <a href="https://www.anaconda.com/">Anaconda</a>.  It is recommended you not use your system python, rather use a virtual environment.

```>> conda create -n pps python=3.7
>> conda activate pps
>> conda install -c anaconda numpy
>> conda install -c conda-forge matplotlib
>> conda install obspy
```
# Files

- *portable_pip_squeak.py*               This is the main ObsPy based python script.
- *config.portable_pip_sqeak*            Parameters are stored in here, defaulted to decent values for PNSN's purposes.
- *parse_and_validate_args_portable.py*  Just a bunch of validating that the input arguments are kosher.

# Speed

This should take around 5 sec for one hour of data.  Plotting adds maybe a minute.

# Examples of running portable_pip_squeak.py

Download 2.5 hours of the braodband (HH? channels) data from UW.TKEY, a 6 channel site:
```
>> ./portable_pip_squeak.py -N UW -S TKEY -C HH -s 2020-04-04T22:00:00 -d 2.5
``` 

Download 1 hours of data from an existing station NC.LMC from a datacenter other than IRIS.  Also make a plot:
```
>> ./portable_pip_squeak.py -N NC -S LMC -C HN -s 2020-04-04T22:00:00 -d 1 -dc NCEDC -p
``` 

Use already existing .SAC files.  Needs 3 SAC files of equal npts as well as gain factor/Stage 0 sensitivity:
```
>> ls MySACfile_Z.SAC MySACfile_N.SAC MySACfile_E.SAC > MyInfile.txt 
>> ./portable_pip_squeak.py -i MyInfile.txt -g 0.102
``` 

Use already existing .mseed file(s).  Needs 3 components of equal npts as either one or three .mseed files as well as gain factor/Stage 0 sensitivity; also make a plot:
```
>> ls My_UW_portable_BB.mseed > MyInfile.txt 
>> ./portable_pip_squeak.py -i MyInfile.txt -g 9.77e8 -p
``` 

# Output columns

- *NoiseFloor*  To approximate the median envelope amplitude quickly, half range of the 2nd to 
    98th percentile amplitudes are used using Acceleration.
- *SAspikes*    This is the same as the eew_stationreport "Spikes" metric and threshold.
- *SArms*       This is the same as the eew_stationreport "RMS" metric and threshold.
- *Big spikes"  Talk about FinDe.
high freq
good for 3/4/5

```
./portable_pip_squeak.py -N UW -S PIER -C HN -s 2020-04-04T22:00:00 -d 1 -dc IRIS -p
UW.PIER  From: 2020-04-04T22:00:00.000000Z to 2020-04-04T23:00:00.000000Z (1.0 hrs)
                           NoiseFloor 	 SAspikes 	 SArms 	 Big spikes
HNE
Without high frequencies:  LOUD 	 FAIL 	 	 FAIL 	 VERY BAD
With high frequencies:     LOUD 	 FAIL 	 	 FAIL 	 VERY BAD
Spikes only at high freq:  No
HNN
Without high frequencies:  LOUD 	 FAIL 	 	 FAIL 	 VERY BAD
With high frequencies:     LOUD 	 FAIL 	 	 FAIL 	 VERY BAD
Spikes only at high freq:  No
HNZ
Without high frequencies:  LOUD 	 FAIL 	 	 FAIL 	 VERY BAD
With high frequencies:     LOUD 	 FAIL 	 	 FAIL 	 VERY BAD
Spikes only at high freq:  No
Good for:  3 channel site-  No   4 channel site-  No   6 channel site-  No
```

```
./portable_pip_squeak.py -N UW -S ALKI -C HN -s 2020-04-04T22:00:00 -d 1 -dc IRIS
UW.ALKI  From: 2020-04-04T22:00:00.000000Z to 2020-04-04T23:00:00.000000Z (1.0 hrs)
                           NoiseFloor 	 SAspikes 	 SArms 	 Big spikes
HNE
Without high frequencies:  Quiet 	 Pass 	 	 Pass 	 None
With high frequencies:     Quiet 	 Pass 	 	 Pass 	 None
Spikes only at high freq:  No
HNN
Without high frequencies:  Quiet 	 Pass 	 	 Pass 	 None
With high frequencies:     Quiet 	 Pass 	 	 Pass 	 None
Spikes only at high freq:  No
HNZ
Without high frequencies:  Quiet 	 Pass 	 	 Pass 	 None
With high frequencies:     Quiet 	 Pass 	 	 Pass 	 None
Spikes only at high freq:  No
Good for:  3 channel site-  Yes   4 channel site-  Yes   6 channel site-  Yes
```
