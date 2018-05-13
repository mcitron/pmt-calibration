#! /usr/bin/env python

## take a binary file output from DRS scope
## create a ROOT file with time and voltage arrays
##
## note that time values are different for every event since
## DRS binning is uneven, and the trigger can occur in any bin

import struct
import datetime as dt
import json
import numpy as np
import ROOT as r

# name = "r878b_1450V_2p7V_20ns_300Hz_100000evnts_v4"
# name = "no_pmt_fixtrig_300Hz_20000"
name = "afterpulses/r878_1450V_2p7V_13ns_300Hz_0p6GHz_200000evnts"

indir =  "/nfs-7/userdata/bemarsh/milliqan/pmt_calib/scope_output"
outdir = "/nfs-7/userdata/bemarsh/milliqan/pmt_calib/processed"

fin = indir+"/{0}.dat".format(name)
fout = r.TFile(outdir+"/{0}.root".format(name), "RECREATE")
ts = np.zeros(1024, dtype='float')
vs = np.zeros(1024, dtype='float')
t = r.TTree("Events","Events")
t.Branch("times", ts, 'times[1024]/D')
t.Branch("voltages", vs, 'voltages[1024]/D')

READ_CHN = 2
N_BINS = 1024 # number of timing bins per channel
POLARITY = -1 # use -1 to invert. peaks must be positive

def getStr(fid, length):
    data = fid.read(length)
    if len(data)==0:
        return None
    res = struct.unpack("c"*length, data)
    return "".join(res)

def getShort(fid, num=1):
    data = fid.read(2*num)
    if len(data)==0:
        return None
    res = struct.unpack("H"*num, data)
    return res[0] if num==1 else res

def getFloat(fid, num=1):
    data = fid.read(4*num)
    if len(data)==0:
        return None
    res = struct.unpack("f"*num, data)
    return res[0] if num==1 else res

def getInt(fid, num=1):
    data = fid.read(4*num)
    if len(data)==0:
        return None
    res = struct.unpack("I"*num, data)
    return res[0] if num==1 else res

fid = open(fin,'rb')

# make sure file header is correct
fhdr = getStr(fid, 4)
if fhdr != "DRS2":
    print "ERROR: unrecognized header "+fhdr
    exit(1)

# make sure timing header is correct
thdr = getStr(fid, 4)
if thdr != "TIME":
    print "ERROR: unrecognized time header "+thdr
    exit(1)

# get the boards in file
n_boards = 0
channels = []
board_ids = []
bin_widths = []
while True:
    bhdr = getStr(fid, 2)
    if bhdr != "B#":
        fid.seek(-2, 1)
        break
    board_ids.append( getShort(fid) )
    n_boards += 1
    bin_widths.append([])
    print "Found Board #"+str(board_ids[-1])

    while True:
        chdr = getStr(fid, 3)
        if chdr != "C00":
            fid.seek(-3, 1)
            break
        cnum = int(getStr(fid, 1))
        print "  Found channel #"+str(cnum)
        channels.append(cnum)
        bin_widths[n_boards-1].append(getFloat(fid, N_BINS))

    if len(bin_widths[n_boards-1])==0:
        print "ERROR: Board #{0} doesn't have any channels!".format(board_ids[-1])

if n_boards==0:
    print "ERROR: didn't find any valid boards!"
    exit(1)

if n_boards > 1:
    print "ERROR: only support one board. Found {0} in file.".format(n_boards)
    exit(1)

bin_widths = bin_widths[0]
n_chan = len(bin_widths)
rates = []

n_evt = 0
while True:
    ehdr = getStr(fid, 4)
    if ehdr == None:
        break
    if ehdr != "EHDR":
        raise Exception("Bad event header!")

    n_evt += 1
    # print "Found Event #"+str(n_evt)
    serial = getInt(fid)
    # print "  Serial #"+str(serial)
    date = getShort(fid, 7)
    date = dt.datetime(*date[:6], microsecond=1000*date[6])
    # print "  Date: "+str(date)
    rangeCtr = getShort(fid)
    # print "  Range Center: "+str(rangeCtr)
    getStr(fid, 2)
    b_num = getShort(fid)
    getStr(fid, 2)
    trig_cell = getShort(fid)
    # print "  Trigger Cell: "+str(trig_cell)

    for ichn in range(n_chan):
        chdr = getStr(fid, 4)
        if chdr != "C00"+str(channels[ichn]):
            print "ERROR: bad event data!"
            exit(1)

        scaler = getInt(fid)
        voltages = np.array(getShort(fid, N_BINS))
        if READ_CHN != channels[ichn]:
            continue
        voltages = voltages/65535. * 1000 - 500 + rangeCtr
        times = np.roll(np.array(bin_widths[ichn]), N_BINS-trig_cell)
        times = np.cumsum(times)
        times = np.append([0], times[:-1])
        rates.append((times[-1]-times[0])/(times.size-1))
        np.copyto(ts, times)
        np.copyto(vs, voltages)
        
        t.Fill()

print "Measured sampling rate: {0:.2f} GHz".format(1.0/np.mean(rates))
# t.Write()
fout.Close()




