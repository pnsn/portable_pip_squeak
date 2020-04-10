#!/home/ahutko/anaconda3/bin/python

# Quick and dirty rejiggering of 
# https://github.com/pnsn/station_metrics/blob/master/eew_stationreport
# for CSN data in SAC files.  Feb 2020
# Alex H.

import obspy
from obspy import read
from obspy import UTCDateTime
from obspy.signal.trigger import classic_sta_lta
from obspy.clients.fdsn import Client
import numpy as np
import datetime
import timeit
from noise_metrics import *
from parse_and_validate_args import validate_args_and_get_times, parse_args
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.dates import DateFormatter

try:
    from configparser import SafeConfigParser
except:
    from ConfigParser import SafeConfigParser

#---- set up some parameters 

# config file
configfile = 'config.CSN'

# seconds to make amplitude measurement following stalta trig in count_peaks_stalta_new 
twin = 4.0 

# minimum time (sec) bw peaks in count_triggers_FinDer 
TbwPeaksFD = 30. 

#----- command-line arguments specify what data to get and where
#      Note: plotting can slow it down so be careful.
args = parse_args()
[starttime, endtime, network, station, chantype, datacenter, infile, sensitivity ] = validate_args_and_get_times(args)
iplot = args.iplot

#---- This version is simplified from the pip_squeak version to not worry about RESP.

def raw_trace_to_ground_motion_filtered(TraceOrig,AccVelDisp,f1,f2):
    '''Input is an ObsPy trace that will get demeaned, converted to ACC/VEL/DISP, and then filtered.
      set f1 = 0 to apply a lowpass at f2 Hz
      set f2 = 0 to apply a highpass at f1 Hz
    '''
    if ( "cc" in AccVelDisp or "CC" in AccVelDisp ):
        GM = "ACC"
    elif ( "V" in AccVelDisp or "v" in AccVelDisp ):
        GM = "VEL"
    else:
        GM = "DISP"
    trace = TraceOrig.copy()
    trace.detrend(type='demean')
    if ( f1 == 0 ):
        trace.filter("lowpass", freq = f2)
    elif ( f2 == 0 ):
        trace.filter("highpass", freq = f1)
    else:
        trace.filter("bandpass",freqmin=f1,freqmax=f2,corners=2,zerophase=False)
    if ( trace.stats.channel[1:2] == "N" ):
        GMnative = "ACC"
    else:
        GMnative = "VEL"
    if ( GMnative == GM ):
        idonothing = 1
    elif ( GMnative == "ACC" and GM == "VEL" and trace.stats.npts > 2 ):
        trace.integrate(method='cumtrapz')
    elif ( GMnative == "ACC" and GM == "DISP" ):
        trace.integrate(method='cumtrapz')
        trace.integrate(method='cumtrapz')
    elif ( GMnative == "VEL" and GM == "ACC" and trace.stats.npts > 2 ):
        trace.differentiate(method='gradient')
    trace.detrend(type='demean')              #--- of order 0.002 sec/hr trace

    return trace

#---- Read Processing parameters from a config file

parser = SafeConfigParser()
parser.read(configfile)

tpadding = float(parser.get('SectionOne','tpadding'))
sta = float(parser.get('SectionOne','sta'))
lta = float(parser.get('SectionOne','lta'))
tbuffer = sta + lta
freq1 = float(parser.get('SectionOne','freq1'))
freq2 = float(parser.get('SectionOne','freq2'))
freqHP = float(parser.get('SectionOne','freqHP'))
mpd = float(parser.get('SectionOne','mpd'))
minSTALTAlen = float(parser.get('SectionOne','minSTALTAlen'))
RMSlen = float(parser.get('SectionOne','RMSlen'))
FDSNtimeout = float(parser.get('SectionOne','FDSNtimeout'))

#---- If infile exists, go through the lines of sac/mseed files and combine
#        into one stream.  If a Net Stat Chantype given, go download data.

if ( infile is not None ):
    try:
        f = open(infile)
        lines = f.read().splitlines()
        firstline = True
        for line in lines:
            if ( firstline == True ):
                st = read('/home/ahutko/TEMP/NICHOLAS_TEMP_TEST/KOPA/kopa/KOPA/UW.KOPA..HHE.D.2020.057.232452.SAC')
                st = read(line)
                firstline = False
            else:
                st = st + read(line)
    except:
        print("Error:  couldn't read the SAC or mseed file: " + line + " in: " + infile)
        exit()
else:
    # Set up the FDSN client
    try:
        client = Client(datacenter,timeout = FDSNtimeout)
    except Exception as error:
        print ("Error: {}, Failed to connect to FDSN client. Used a timeout of {}".format(error, FDSNtimeout))
    # Download data
    try:
        channels =  chantype + "*" 
        t1 = UTCDateTime(starttime)
        t2 = UTCDateTime(endtime)
        st = client.get_waveforms(network,station,"*",channels,t1,t2,longestonly=True)
        if ( len(st) < 3 ):
            print("Error:  Not all channels are available during this period " + st)
        inv = client.get_stations(network=network, station=station, channel=channels, starttime=t1, endtime=t2, level='response')
        for i in range(0,len(inv[0][0])):
            chan = inv[0][0][i]
            if ( chan.code[2:3] == 'Z' ):
                sensitivity = chan.response.instrument_sensitivity.value
    except:
        print("Error:  Problem downloading data " + network + "." + station + "." 
              + channels + " " + str(t1) + " to " + str(t2) )

GoodFor3Ch = 'Yes'
GoodFor4Ch = 'Yes'
GoodFor6Ch = 'Yes'
n = 0

fig = plt.figure(figsize=(8,8),dpi=80)
if ( iplot == True ):
    gs1 = gridspec.GridSpec(1,1)
    gs1.update(left=0.13,right=0.91,bottom=0.85,top=0.95,wspace=0.01)
    ax1 = plt.subplot(gs1[:,:])
    gs2 = gridspec.GridSpec(1,1)
    gs2.update(left=0.13,right=0.91,bottom=0.7,top=0.8,wspace=0.01)
    ax2 = plt.subplot(gs2[:,:])
    gs3 = gridspec.GridSpec(1,1)
    gs3.update(left=0.13,right=0.91,bottom=0.55,top=0.65,wspace=0.01)
    ax3 = plt.subplot(gs3[:,:])
    gs4 = gridspec.GridSpec(1,1)
    gs4.update(left=0.13,right=0.91,bottom=0.4,top=0.5,wspace=0.01)
    ax4 = plt.subplot(gs4[:,:])
    gs5 = gridspec.GridSpec(1,1)
    gs5.update(left=0.13,right=0.91,bottom=0.25,top=0.35,wspace=0.01)
    ax5 = plt.subplot(gs5[:,:])
    gs6 = gridspec.GridSpec(1,1)
    gs6.update(left=0.13,right=0.91,bottom=0.1,top=0.2,wspace=0.01)
    ax6 = plt.subplot(gs6[:,:])
else:
    iplot = False

for tr in st:
    try:
        start = timeit.default_timer()
        dt = tr.stats.delta
        Tlen = tr.stats.npts * tr.stats.delta
        net = tr.stats.network
        stat = tr.stats.station
        chan = tr.stats.channel
        loc = tr.stats.location

        #--- get rid of long period noise and correct for sensitivity
        tr.filter("highpass", freq = 0.02)
        tr.taper(0.1,max_length = 100)
        tr.detrend(type='demean')
        tr.data = tr.data / sensitivity

        #--- get a numpy array for the STA/LTA function used for trig detection in Elarms
        trVelTrig = raw_trace_to_ground_motion_filtered(tr,"Vel",freqHP,0)
        datastalta = classic_sta_lta(trVelTrig,int(sta/dt),int(lta/dt))  # 0.13 s/hr long trace

        #--- get filtered traces using either a BP or HP filter
        trAccBP = raw_trace_to_ground_motion_filtered(tr,"Acc",freq1,freq2)
        trVelBP = raw_trace_to_ground_motion_filtered(tr,"Vel",freq1,freq2)
        trDisBP = raw_trace_to_ground_motion_filtered(tr,"Dis",freq1,freq2)
        trAccHP = raw_trace_to_ground_motion_filtered(tr,"Acc",freq1,0)
        trVelHP = raw_trace_to_ground_motion_filtered(tr,"Vel",freq1,0)
        trDisHP = raw_trace_to_ground_motion_filtered(tr,"Dis",freq1,0)

        dataAccBP = trAccBP.data
        dataVelBP = trVelBP.data
        dataDisBP = trDisBP.data
        dataAccHP = trAccHP.data
        dataVelHP = trVelHP.data
        dataDisHP = trDisHP.data

        #--- calculate your metrics
        AccMaxBP = max(abs(dataAccBP))
        AccMaxHP = max(abs(dataAccHP))
        AccNoiseFloorBP = noise_floor(dataAccBP)  # 0.01 s/hr long trace
        AccNoiseFloorHP = noise_floor(dataAccHP)

        #--- These is the first step Elarms takes for trigger detection.
        #    It's a simple sta/lta > 20 on the 3Hz highpass Vel trace
        #    and amplitude > threshold within twin sec of the stalta trig.
        #    Current ShakeAlert station acceptance is < 1/hr above 0.34 cm/s^2
        snr20_0p034cm_BP = count_peaks_stalta_new(dataAccBP,datastalta,sta,lta,mpd,20.,dt,twin,0.0034) *Tlen/3600 # 0.002 s/hr long trace
        snr20_0p034cm_HP = count_peaks_stalta_new(dataAccHP,datastalta,sta,lta,mpd,20.,dt,twin,0.0034) *Tlen/3600 

        #--- This mimics Elarms trigger counting and includes Elarms amp thresholds for Vel,Acc,Dis.
        NTrigElarmS_BP,Nboxcars_BP = count_peaks_stalta_Elarms(dataAccBP,dataVelBP,dataDisBP,datastalta,sta,lta,mpd,20,dt,twin,tr.stats.channel) # 0.003 s/hr long trace
        NTrigElarmS_HP,Nboxcars_HP = count_peaks_stalta_Elarms(dataAccHP,dataVelHP,dataDisHP,datastalta,sta,lta,mpd,20,dt,twin,tr.stats.channel)
        NTrigElarmS_BP = NTrigElarmS_BP *Tlen/3600
        NTrigElarmS_HP = NTrigElarmS_HP *Tlen/3600

        #--- "Pump noise" metric gets duration (sec) that RMS > threshold in m/s^2.
        #    Current ShakeAlert station acceptance is < 60s/hr above 0.07cm/s^2.
        RMSduration_0p07cm_BP = duration_exceed_RMS(dataAccBP,0.0007,RMSlen,dt) *Tlen/3600 # 0.015 s/hr long trace
        RMSduration_0p07cm_HP = duration_exceed_RMS(dataAccHP,0.0007,RMSlen,dt) *Tlen/3600

        #--- simulated FinDer triggers.  Threshold in m/s^2.
        finder_1cm_BP = count_triggers_FinDer(dataAccBP,TbwPeaksFD, 0.01, dt) *Tlen/3600 # 0.002 s/hr long trace
        finder_1cm_HP = count_triggers_FinDer(dataAccHP,TbwPeaksFD, 0.01, dt) *Tlen/3600
        finder_2cm_BP = count_triggers_FinDer(dataAccBP,TbwPeaksFD, 0.02, dt) *Tlen/3600
        finder_2cm_HP = count_triggers_FinDer(dataAccHP,TbwPeaksFD, 0.02, dt) *Tlen/3600

        #--- Test measurements against subjective thresholds
        if ( AccNoiseFloorBP < 0.0007 ):
            NoiseFloorBP = 'Quiet'
        elif ( AccNoiseFloorBP > 10 * 0.0007 ):
            NoiseFloorBP = 'VERY LOUD'
        else:
            NoiseFloorBP = 'LOUD'

        if ( AccNoiseFloorHP < 0.0007 ):
            NoiseFloorHP = 'Quiet'
        elif ( AccNoiseFloorHP > 10 * 0.0007 ):
            NoiseFloorHP = 'VERY LOUD'
        else:
            NoiseFloorHP = 'LOUD'

        if ( snr20_0p034cm_BP > 1.0 ):
            ShakeAlertSpikes = 'FAIL'
        else:
            ShakeAlertSpikes = 'Pass'

        if ( snr20_0p034cm_HP > 1.0 ):
            ShakeAlertSpikesHP = 'FAIL'
        else:
            ShakeAlertSpikesHP = 'Pass'

        if ( RMSduration_0p07cm_BP > 60. ):
            ShakeAlertRMS = 'FAIL'
        else:
            ShakeAlertRMS = 'Pass'

        if ( RMSduration_0p07cm_HP > 60. ):
            ShakeAlertRMSHP = 'FAIL'
        else:
            ShakeAlertRMSHP = 'Pass'

        SpikeAmpBP = 'None'
        if ( finder_1cm_BP > 0 and finder_1cm_BP <= 1. ):
            SpikeAmpBP = 'Not many'
        elif ( finder_1cm_BP > 1. ):
            SpikeAmpBP = 'BAD'
        if ( finder_2cm_BP > 1. ):
            SpikeAmpBP = 'VERY BAD'

        SpikeAmpHP = 'None'
        if ( finder_1cm_HP > 0 and finder_1cm_HP <= 1. ):
            SpikeAmpHP = 'Not many'
        elif ( finder_1cm_HP > 1. ):
            SpikeAmpHP = 'BAD'
        if ( finder_2cm_HP > 1. ):
            SpikeAmpHP = 'VERY BAD'

        MostSpikesHighFreqOnly = 'No'
        if ( snr20_0p034cm_BP == 0 ):
           if ( snr20_0p034cm_HP > 0 ):
               MostSpikesHighFreqOnly = 'Yes'
        elif ( (snr20_0p034cm_HP / snr20_0p034cm_BP) > 2 ):
            MostSpikesHighFreqOnly = 'Yes'

        if ( finder_2cm_BP > 2 or AccNoiseFloorBP > ( 10 * 0.0007) ):
            GoodFor3Ch = 'No'
        if ( finder_1cm_BP > 1 or AccNoiseFloorBP > 0.0007 ):
            GoodFor4Ch = 'No'
        if ( finder_1cm_HP > 0.5 or AccNoiseFloorBP > (0.0007 / 10.) or ShakeAlertRMS == 'FAIL' or ShakeAlertSpikes == 'FAIL' ):
            GoodFor6Ch = 'No'

        #--- print out the results
        if ( n == 0 ):
            print(tr.stats.network + "." + tr.stats.station + "  From: " + str(tr.stats.starttime) + " to " + str(tr.stats.endtime) + " (" + str(int(1000*Tlen/3600.)/1000) +  " hrs)" )
            print('                          ','NoiseFloor',"\t","SAspikes","\t","SArms","\t","Big spikes")
        print(tr.stats.channel)
        print('Without high frequencies: ',NoiseFloorBP, "\t", ShakeAlertSpikes, "\t","\t", ShakeAlertRMS, "\t", SpikeAmpBP )
        print('With high frequencies:    ',NoiseFloorHP, "\t", ShakeAlertSpikesHP, "\t","\t", ShakeAlertRMSHP, "\t", SpikeAmpHP )
        print('Spikes only at high freq:  ' + MostSpikesHighFreqOnly )

        if ( iplot == True ):
            t = tr.stats.starttime
            timearray = np.arange(tr.stats.npts)*dt
            x =  np.array([datetime.datetime(t.year,t.month,t.day,t.hour,t.minute,t.second,0) + datetime.timedelta(seconds=jj) for jj in timearray ])
            chan = tr.stats.channel
            if ( chan[2:3] == 'Z' ):
               axA = ax1
               axB = ax4
            elif ( chan[2:3] == 'N' ):
               axA = ax2
               axB = ax5
            elif ( chan[2:3] == 'E' ):
               axA = ax3
               axB = ax6
            if ( len(x) > (5401/dt) ):
                iplot1 = int(len(x)/2) - int(2700/dt)
                iplot2 = int(len(x)/2) + int(2700/dt)
            else:
                iplot1 = 0
                iplot2 = len(x)-1
            xtextR = x[iplot2] + datetime.timedelta(seconds=int(iplot2-iplot1)*dt*0.06)
            xtextL = x[iplot1] - datetime.timedelta(seconds=int(iplot2-iplot1)*dt*0.17)
            axA.plot(x[iplot1:iplot2],trAccHP.data[iplot1:iplot2]*100,color='k')
            axA.plot([x[iplot1],x[iplot2]],[0.34,0.34],color='r')
            axA.set_xticklabels([])
            axA.text(xtextR,0,chan,fontsize=12)
            axB.plot(x[iplot1:iplot2],trAccHP.data[iplot1:iplot2]*100,color='k')
            axB.plot([x[iplot1],x[iplot2]],[2,2],color='green')
            axB.set_xticklabels([])
            axB.text(xtextR,0,chan,fontsize=12)
            if ( chan[2:3] == 'N'):
                axA.text(xtextL,0.,'cm/s^2',backgroundcolor='w',fontsize=12,rotation=90)
                axB.text(xtextL,0.,'cm/s^2',backgroundcolor='w',fontsize=12,rotation=90)
            if ( chan[2:3] == 'E' ):
                miny = min(trAccHP.data[iplot1:iplot2]*100)
                maxy = max(2,max(trAccHP.data[iplot1:iplot2]*100))
                miny = miny - ((maxy-miny)*0.7)
        n = n + 1
    except:
        print('failed to read SAC file: ',line)

if ( iplot == True ):
    ax6.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ilabel1 = iplot2 - int((iplot2-iplot1)*0.99)
    ilabel2 = iplot2 - int((iplot2-iplot1)*0.35)
    ylabel = min(-1,miny)
    ax6.text(x[ilabel1],ylabel,'ShakeAlert acceptance spike threshold',color='r',fontsize=12)
    ax6.text(x[ilabel2],ylabel,'FinDer trigger threshold',color='green',fontsize=12)
    title = tr.stats.network + "." + tr.stats.station + "   " + str(t.year) + "/" + str(t.month) + "/" + str(t.day) + "   (0.075 Hz highpass)"
    fig.suptitle(title,size=15)
    figname = "station_report." + tr.stats.network + "." + tr.stats.station + ".png"
    plt.savefig(figname)
    plt.close("all")

print("Good for:  3 channel site- ", GoodFor3Ch, "  4 channel site- ", GoodFor3Ch, "  6 channel site- ", GoodFor6Ch )

