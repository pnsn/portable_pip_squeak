# portable_pip_squeak
Light weight station quality assessment tool to tell if site is appropriate for a 3,4, or 6 channel station.

# Overview

This is similar to <a href="https://github.com/pnsn/station_metrics">station_metrics</a>, but is tuned for a quick and dirty assessment of an existing or candidate site.  It is similar to eew_stationreport in that project which is informally known as "pip_squeak".  eew_stationreport is used to assess new seismic stations before acceptance into ShakeAlert.  portable_pip_squeak is similar in that uses many of the same metrics and thresholds, but includes others and the output is intentionally subjective for use by field techs.

# Installation

Built for python 3 with the following major package dependencies:
<a href="http://www.numpy.org/">numpy</a>, <a href="https://github.com/obspy/obspy/wiki">obspy</a>, <a href="https://matplotlib.org">matplotlib</a>.

These can easily be added via command line using <a href="https://www.anaconda.com/">Anaconda</a>.  It is recommended you not use your system python, rather use a virtual environment.

```
>> conda create -n pps
>> conda activate pps
>> conda install -c anaconda numpy
>> conda install -c conda-forge matplotlib
>> conda install obspy
```
# Files

- *portable_pip_squeak.py*               This is the main ObsPy based python script.
- *config.portable_pip_sqeak*            Parameters are stored in here, defaulted to decent values for PNSN's purposes.
- *parse_and_validate_args_portable.py*  Just a bunch of validation of the input arguments.

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

Use already existing .SAC files.  Needs 3 SAC files of equal length as well as gain factor/Stage 0 sensitivity:
```
>> ls MySACfile_Z.SAC MySACfile_N.SAC MySACfile_E.SAC > MyInfile.txt 
>> ./portable_pip_squeak.py -i MyInfile.txt -g 0.102
``` 

Use already existing .mseed file(s).  Needs 3 traces (usualy Z, N, E) of equal length as either one or three .mseed files as well as gain factor/Stage 0 sensitivity; also make a plot:
```
>> ls My_UW_portable_BB.mseed > MyInfile.txt 
>> ./portable_pip_squeak.py -i MyInfile.txt -g 9.77e8 -p
``` 

# Simple shell script to run assessment on a bunch of stations.

```
#!/bin/sh

./portable_pip_squeak.py -N UW -S TKEY -C HH -s 2020-04-04T22:00:00 -d 2.5
./portable_pip_squeak.py -N UW -S TKEY -C EN -s 2020-04-04T22:00:00 -d 2.5
./portable_pip_squeak.py -N UW -S ALKI -C HN -s 2020-04-04T22:00:00 -d 2.5
./portable_pip_squeak.py -N UW -S PIER -C HN -s 2020-04-04T22:00:00 -d 2.5
```
- Make sure you make the script executable via:  chmod +x myshellscript.sh
- Execute it, dump the output into a text file, and put it in the background so you can log off and see results tomorrow:
nohup ./myshellscript.sh & > myoutpout.txt

# Output columns

- *NoiseFloor*  To approximate the median envelope amplitude quickly, half range of the 2nd to 
    98th percentile amplitudes are used using acceleration. <a href="https://github.com/pnsn/station_metrics/tree/master/station_metrics/metrics#noise-floor">details</a>
- *SAspikes*    This is the same as the <a href="https://github.com/pnsn/station_metrics/tree/master/station_metrics/metrics">eew_stationreport</a> "Spikes" metric and threshold. <a href="https://github.com/pnsn/station_metrics/tree/master/station_metrics/metrics#metrics-counting-spikestriggers">details</a>
- *SArms*       This is the same as the <a href="https://github.com/pnsn/station_metrics/tree/master/station_metrics/metrics">eew_stationreport</a> "RMS" metric and threshold. <a href="https://github.com/pnsn/station_metrics/tree/master/station_metrics/metrics#rms-noise">details</a>
- *Big spikes*  This uses occurrences of amplitudes greater than 1 or 2 cm/s^2, which is enough for a human to feel, to mimic how problematic the station would be in the ShakeAlert FinDer algorithm.

*Frequencies* A lot of spikes occur only at high frequencies that can be ignored for some earthquake applications, possibly ShakeAlert in the future.  All of the above measures are calculated twice.  Once using the current filtering in the ShakeAlert EPIC waveform processor, a highpass at 0.075 Hz, and again using a bandpass filter of 0.075 - 15 Hz.  The *Spikes at only high freq* is based on the ratio of *SAspikes* using the highpass to *SAspikes* using the bandpass.

*Good for:  3/4/6 channel site*  These are subjective (and likely to be updated) guidelines for what type of station a site might be good for using a combination of *Big Spikes* and the noise floor for the 3 and 4 channel assessment.  The 6 channel assessment additionaly considers *SAspikes* and *SArms*.

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

# Output plots:
Adding -p will spit out these type plots which show highpass filtered acceleration data along with colored lines demarking threshold amplitudes.  These plots are limited to at most 1.5 hours, so if the input duration or files are 4 hours long, only a time window of 1.5 hours in the center of the analyzed seismograms will be plotted.

In these examples, UW.PIER, a notoriously noisy station, clearly does not pass ShakeAlert station acceptance and would be a nuisance to FinDer.  UW.ALKI is a much quieter station and visibly passes fairly easily. 

<img src="https://github.com/pnsn/portable_pip_squeak/blob/master/station_report.UW.PIER.png" width=550 alt="Metric: Noise Floor" />

<img src="https://github.com/pnsn/portable_pip_squeak/blob/master/station_report.UW.ALKI.png" width=550 alt="Metric: Noise Floor" />
