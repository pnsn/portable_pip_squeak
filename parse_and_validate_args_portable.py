#!/usr/bin/env python

import argparse
import sys
from uuid import uuid4
import datetime

def parse_args():
    """
    Parse input arguments and make an extended help menu.
    If Loc = "--", change it to a token string and then back to a string.
    """
    sentinel_dict = {}

    def _preprocess_sysargv(argv):
        inputs = []
        for arg in argv[1:]:
            # handles case where values contain --, otherwise they will
            # be interpreted as arguments.
            if '--,' in arg or ',--' in arg or arg == '--':
                sentinel = uuid4().hex
                key = '%s' % sentinel
                sentinel_dict[key] = arg
                inputs.append(sentinel)
            else:
                inputs.append(arg)
        return inputs

    def _postprocess_sysargv(v):
        if v in sentinel_dict:
            return sentinel_dict.get(v)
        else:
            return v

    #----- read input arguments
    for i, arg in enumerate(sys.argv):
       if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
    parser = argparse.ArgumentParser()

    parser.add_argument('-u', action='store_true', dest='helpmenu',help='extended HELP MENU with examples')
    parser.add_argument('-i','--infile',action='store', dest='infile',help='name of file with SAC or mseed file(s)')
    parser.add_argument('-g','--gain',action='store', dest='sensitivity',help='Stage 0 sensitivity')
    parser.add_argument('-N','--net', action='store', dest='network',help='network')
    parser.add_argument('-S','--sta', action='store', dest='station',help='station')
    parser.add_argument('-C','--cha', action='store', dest='chantype',help='chantype')
    parser.add_argument('-s','--start', action='store', dest='startstring',help='start time YYYY-MM-DDTHH:MM:SS')
    parser.add_argument('-e','--end', action='store', dest='endstring',help='end time YYYY-MM-DDTHH:MM:SS')
    parser.add_argument('-d','--duration', action='store', dest='durationinhours',help='duration in hours')
    parser.add_argument('-dc','--dc','--datacenter', action='store', dest='datacenter',default='IRIS',help='FDSN data center (e.g. IRIS, SCEDC, NCEDC)')
    parser.add_argument('-p','--plot',action='store_true',dest='iplot',help='make plots of each hourly trace (NOTE: can be slow)')

    helpextended = parser.parse_args(_preprocess_sysargv(sys.argv)).helpmenu
    if ( helpextended is True  ):
        print ('')
        print ('portable_pip_squeak: assess a station either using local data or to be downloaded')
        print ('')
        print ('Usage: portable_pip_squeak.py [options]')
        print ('')
        print ('EXAMPLES:')
        print ('portable_pip_squeak.py --infile my_SAC_files.txt')
        print ('portable_pip_squeak.py -N UW -S TKEY -C HH -s 2018-01-01T00:00:00 -d 2 -p')
        print ('portable_pip_squeak.py -N CI -S LEO -C HN -s 2020-01-01T00:00:00 -d 24 -dc SCEDC')
        print ('')
        print ('Inputs if supplying your own data:')
        print ('  -i,  --infile         Name of text file with SAC/mseed file(s) of 3 (Z,N,E) traces.')
        print ('  -g,  --gain           Gain or Stage 0 sensitivity')
        print (' ')
        print ('Inputs if downloading data:')
        print ('  -s,  --starttime      Trace start time (YYYY-MM-DD,HH:MM:SS)')
        print ('')
        print ('  One of these:')
        print ('  -e,  --endtime        Trace end time (YYYY-MM-DD,HH:MM:SS)')
        print ('  -d,  --duration       Duration in hours from starttime')
        print ('                      Note: if duration is neg, starttime becomes endtime')
        print ('  N, S, C and a datacenter if other than IRIS')
        print ('  -N,  --net            Network code')
        print ('  -S,  --sta            Station code')
        print ('  -C,  --cha            Channel type, e.g. EN or HH')
        print ('  -dc, --datacenter     Name of FDSN data center if not IRIS, e.g. SCEDC, NCEDC')
        print (' ')
        print ('Optional flags:')
        print ('-P,  --plot           Flag to make a figure for each hour.  Note: can be slow.')
        print ('-u                    Print this extended help menu')
        print ('')


    return parser.parse_args(_preprocess_sysargv(sys.argv))


def validate_args_and_get_times(args):
    """
    Validate arguments and pass back start/end times and durations.
    """
    infile = args.infile
    sensitivity = args.sensitivity
    network = args.network
    station = args.station
    chantype = args.chantype
    startstring = args.startstring
    endstring = args.endstring
    durationinhours = args.durationinhours
    datacenter = args.datacenter

    #----- make sure there is either an entire SNCL or an input file of SNCLs
    if ( infile is not None ):
        try:
#            itest = float(sensitivity) + 1
            sensitivity = float(sensitivity)
            starttime = None
            endtime = None
            network = None
            station = None
            chantype = None
            datacenter = None
        except:
            print("Error:  need a valid --gain value for the Stage 0 sensitivity")
            exit()
    elif ( network is None or station is None or chantype is None ):
        print ("Error:  need either network + station + chantype OR a file with SAC/mseed files")
        exit()
    else:
        infile = None
        sensitivity = None
        #----- validate start/end time
        if ( startstring is None and ( endstring is not None or durationinhours 
             is not None ) ):
            print ("Error: need a starttime and [endtime or duraiton]")
            exit()
        else:
            try:
                starttime = datetime.datetime.strptime(startstring, "%Y-%m-%dT%H:%M:%S")
            except:
                print ("Error: invalid starttime.  Use format: 2018-01-31T23:59:01  or  2018-1-31T23:59:1")
                exit()
            if ( durationinhours is not None ):
                try:
                    if ( float(durationinhours) > 0 ):
                        durationinhours = float(durationinhours) 
                        durationinsec = durationinhours * 3600.
                        endtime = starttime + datetime.timedelta(seconds=durationinsec)
                    else:
                        durationinhours = abs(float(durationinhours))
                        durationinsec = durationinhours * 3600.
                        endtime = starttime
                        starttime = endtime - datetime.timedelta(seconds=durationinsec)
                except:
                    print ("Error: invalid durationinhours; requires integer")
                    exit()
            else:
                try:
                    endtime = datetime.datetime.strptime(endstring, "%Y-%m-%dT%H:%M:%S")
                    durationinsec = (endtime - starttime).total_seconds()
                    durationinhours = durationinsec / 3600.
                except:
                    print ("Error: invalid endtime.  Use format: 2018-01-31T23:59:01  or  2018-1-31T23:59:1")
                    exit()

    #----- validate chantype
    if ( chantype is not None ):
        try:
            chantype = chantype[0:2].upper()
            if ( chantype not in [ 'EN', 'HN', 'HH', 'BH', 'LH', 'EH', 'SH' ] ):
                print ("Error:  need a valid chantype, e.g. BH or EN")
                exit()
        except:
            print ("Error:  need a valid chantype, e.g. BH or EN")
            exit()

    #----- half-validate network
    if ( network is not None ):
        try:
            network = network.upper()
            if ( len(network) < 1 or len(network) > 2 ):
                print ("Error:  need a valid network, e.g. UW")
                exit()
        except:
            print ("Error:  need a valid network, e.g. UW")
            exit()

    #----- half-validate station
    if ( station is not None ):
        try:
            station = station.upper()
            if ( len(station) < 2 or len(station) > 6 ):
                print ("Error:  need a valid station, e.g. ALKI")
                exit()
        except:
            print ("Error:  need a valid station, e.g. ALKI")
            exit()

    #----- validate the FDSN data center against a list
    if (datacenter is not None):
        if (datacenter not in [ "IRIS", "NCEDC", "SCEDC", "BGR", "EMSC", "ETH", "GEONET", "GFZ", "ICGC", "INGV", "IPGP", "ISC", "KOERI", "LMU", "NIEP", "NOA", "ODC", "ORFEUS", "RESIF", "TEXNET", "USGS", "USP" ]):
            print ("Error: need a valid FDSN data center")
            exit()

    return [starttime, endtime, network, station, chantype, datacenter, infile, sensitivity ]


