#! /usr/bin/env/python

## take output root file from processBinary.py and
## add higher level information like pulse time, area, etc

import sys
import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle

tstart = 140
tend   = 290

# tstart = 200
# tend   = 310

# tstart = 280
# tend   = 390

doTiming = True
waveForm = False
if doTiming:
    template = pickle.load(open("peak_templates/template_peak_70_75_0p7GHz.pkl", 'rb'))[::2]
    template /= np.sum(template)

f = r.TFile(sys.argv[1], "UPDATE")
t = f.Get("Events")

times = np.zeros(1024, dtype=float)
voltages = np.zeros(1024, dtype=float)
area = np.array([0], dtype=float)
offset = np.array([0], dtype=float)
noise = np.array([0], dtype=float)
smoothed_max = np.array([0], dtype=float)
tmax = np.array([0], dtype=float)
thalfmax = np.array([0], dtype=float)


t.SetBranchStatus("*", 0)
t.SetBranchStatus("times", 1)
t.SetBranchStatus("voltages", 1)
nt = t.CloneTree()

nt.SetBranchAddress("times",times)
nt.SetBranchAddress("voltages",voltages)
b_area = nt.Branch("area", area, "area/D")
b_offset = nt.Branch("offset", offset, "offset/D")
b_noise = nt.Branch("noise", noise, "noise/D")
b_smoothed_max = nt.Branch("smoothed_max", smoothed_max, "smoothed_max/D")
b_tmax = nt.Branch("tmax", tmax, "tmax/D")
b_thalfmax = nt.Branch("thalfmax", thalfmax, "thalfmax/D")

Nevt = nt.GetEntries()
for ievt in range(Nevt):
# for ievt in range(10000):
    nt.GetEntry(ievt)

    if ievt%10000==0:
        print "iEvt:", ievt

    vs = -voltages

    # for j in range(1,vs.size):
    #     if abs(vs[j] - vs[j-1]) > 10:
    #         vs[j] = vs[j-1]

    istart = np.argmax(times>tstart+times[0])
    iend = np.argmax(times>tend+times[0])

    # offset[0] = np.median(vs)
    offset[0] = np.mean(vs[30:istart*3/4])
    noise[0] = 0.5*(np.percentile(vs[:istart*3/4],95) - np.percentile(vs[:istart*3/4],5))
    
    vs -= offset[0]

    area[0] = np.trapz(vs[istart:iend], times[istart:iend])

    if not doTiming:
        smoothed_max[0] = -999
        tmax[0] = -999
        thalfmax[0] = -999
    else:
        convolved = np.convolve(vs, template[::-1], mode='valid')
        imax = np.argmax(template)
        convolved_time = times[imax:imax+convolved.size]
        icstart = np.argmax(convolved_time>tstart + times[0])
        icend = np.argmax(convolved_time>tend + times[0])
        icmax = np.argmax(convolved[icstart:icend]) + icstart
        cmax = convolved[icmax]
        ihm = icmax
        while ihm > icstart and convolved[ihm] > cmax/2:
            ihm -= 1
        if cmax < 0.5 or convolved[ihm] > cmax/2:
            thalfmax[0] = -999
        else:
            thalfmax[0] = convolved_time[ihm] + (convolved_time[ihm+1]-convolved_time[ihm])/(convolved[ihm+1]-convolved[ihm]) * (cmax/2 - convolved[ihm])

        smoothed_max[0] = cmax
        tmax[0] = convolved_time[icmax]

    if waveForm:
        plt.clf()
        plt.plot(times,vs,color="red")
        plt.axhspan(offset[0]-noise[0],offset[0]+noise[0],color='red',alpha=0.3)
        plt.axvline(tstart,ls="dashed",color="black")
        plt.axvline(tend,ls="dashed",color="black")
        plt.axhline(offset[0],ls="dashed",color="black")
        plt.axhline(offset[0]-noise[0],ls="dashed",color="black",lw=1.2)
        plt.axhline(offset[0]+noise[0],ls="dashed",color="black",lw=1.2)
        plt.xlim(0,1024)
        plt.xlabel('Time (ns)')
        plt.ylim(-10,20)
        plt.ylabel('Vout (mV)')
        peaks = findAPs(times,vs,offset[0],noise[0],template,tstart,tend,plot=False,outName=outdir+"/{0:05d}.png".format(i))
        shift=0
        for n in range(len(peaks)-1):
            if peaks[n]>280:
                peaks=peaks[n:]
                break
        for n in range(len(peaks)):
            for x in range(peaks[n]-20,peaks[n]+20):
                if x>1023:
                    x=1023
                if vs[x]>vs[peaks[n]]:
                    shift+=1
                    peaks[n]=x
            if n-shift>len(peaks):
                if vs[n] < noise[0]+offset[0]:
                    peaks=np.delete(peaks,n)
        plt.grid(True)
        plt.title("waveform "+str(i))
        plt.plot(times[peaks],vs[peaks],"x",markersize=10,markeredgewidth=3,color='Blue',alpha=0.5)
        for peak in peaks:
            n_APs[0] +=1 

        for j in range(len(peaks)):
            end1=peaks[j]
            end2=peaks[j]
            while vs[end1]>offset[0]:
                end1-=1
            while vs[end2]>offset[0] and end2<1023:
                end2+=1
            plt.axvspan(times[end1],times[end2],ymin=0,ymax=vs[xv],color="grey",alpha=0.5)
            AP_time[j] = times[peaks[j]]
            AP_area[j] = np.trapz(vs[end1:end2] - offset[0], times[end1:end2])
        plt.savefig(outdir+"wavform "+str(i)+'.png',dpi=90)
    b_area.Fill()
    b_offset.Fill()
    b_noise.Fill()
    b_smoothed_max.Fill()
    b_tmax.Fill()
    b_thalfmax.Fill()

nt.Write("Events",r.TObject.kWriteDelete)
f.Close()
