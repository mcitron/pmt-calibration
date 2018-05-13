import time
import sys, os
import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle

def findAPs(time, voltage, offset, noise, template, main_peak_start, main_peak_end, plot=False, outName=None):

    voltage = np.array(voltage)-offset
    convolved = np.convolve(voltage, template[::-1], mode='valid')
    imax = np.argmax(template)
    convolved_time = time[imax:imax+convolved.size]

    max_inds = []
    for i in range(convolved.size):
        if convolved_time[i] > main_peak_end and convolved[i] >= np.amax(convolved[max(0,i-5):min(convolved.size,i+6)]) and convolved[i] > noise/2:
            oi = i + imax
            max_inds.append(oi)

    if plot and outName is not None:
        plt.figure(1, figsize=(12,9))
        plt.clf()
        plt.plot(time,voltage, 'b-')
        plt.plot([main_peak_start]*2, [-5,30], 'k--')
        plt.plot([main_peak_end]*2, [-5,30], 'k--')
        plt.plot([0,time[-1]], [0,0], 'k--')
        plt.plot([0,time[-1]], [noise]*2, 'k--')
        plt.plot([0,time[-1]], [-noise]*2, 'k--')
        plt.plot(convolved_time, convolved, 'r-', linewidth=2)
        plt.plot(time[:template.size], template*75, 'g-', linewidth=2)
        for oi in max_inds:
            plt.plot([time[oi]]*2, [0, 20], 'k:')
        plt.gca().set_ylim(-5, 30)
        plt.gca().set_xlim(0, time[-1])
        plt.xlabel("Time (ns)")
        plt.ylabel("Voltage (mV)")
        plt.savefig(outName)

    return max_inds

if __name__=="__main__":

    template = pickle.load(open("../../pmt_calib/peak_templates/template_peak_70_75_1GHz.pkl", 'rb'))[::2]
    template /= np.sum(template)
    
    tstart = 200
    tend = 310

    fname = "/nfs-7/userdata/bemarsh/milliqan/pmt_calib/processed/afterpulses/test_postprocessed.root"
    f = r.TFile(fname)
    t = f.Get("Events")

    bn = fname.split("/")[-1].split(".")[0]
    outdir = "/home/users/bemarsh/public_html/milliqan/aps/" + bn
    os.system("mkdir -p " + outdir)
    os.system("cp ~/scripts/index.php " + outdir)

    start = time.time()

    # for i in range(200):
    for i in range(t.GetEntries()):
        t.GetEntry(i)

        # print i
        
        times = np.array(list(t.times))
        voltage = -np.array(list(t.voltages))
                            
        max_inds = findAPs(times, voltage, t.offset, t.noise, template, tstart, tend, plot=False, outName=outdir+"/{0:05d}.png".format(i))
        
        if i%10==0:
            print i, [mi for mi in max_inds]


    end = time.time()

    print "Time elapsed:", end-start
