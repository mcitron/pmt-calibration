import ROOT as r
import numpy as np
import cPickle as pickle

model_defs = {
    "r878": {
        "template": "assets/template_peak_70_75_2GHz.pkl",
        "input_rate": 2.5,
        }
    }

def GenerateSignal(model="r878", output_samp_freq=2.0, noise=None, extend=None, verbose=False):
    # model defines which PMT to use (see model_defs above)
    # output_samp_freq is the sampling frequency of the ouput waveform
    # (in GHz, so 2.0 means time bin spacing of 0.5 ns)
    # if noise is not None, add random gaussian noise of rms 'noise'
    # if extend is not None, extend each side of output waveform by 'extend' ns of zero voltage

    # Set seed
    r.gRandom.SetSeed(0);

    if model not in model_defs:
        raise Exception("No model definition for PMT model {0}!".format(model))

    template_file = model_defs[model]["template"]
    input_rate = model_defs[model]["input_rate"]
    dt = 1.0 / input_rate

    # Load template
    with open(template_file, "rb") as pin:
        # Get template y-array (millivolts)
        template_y = pickle.load(pin)
        # Construct template x-array (ns)
        template_x = np.arange(0, dt * template_y.size, dt)

        
    # re-sample with new time values
    new_dt = 1.0 / output_samp_freq
    times = np.arange(0, template_x[-1] + 1e-9, new_dt)
    voltages = []
    curj = 0
    for i in range(times.size):
        while curj < template_x.size and template_x[curj] < times[i]:
            curj += 1
        if curj==0 or curj==template_x.size:
            voltages.append(0.0)
        else:
            voltages.append( template_y[curj-1] + (template_y[curj] - template_y[curj-1])/(template_x[curj] - template_x[curj-1]) * (times[i] - template_x[curj-1]) )
    voltages = np.array(voltages)
        
    prenorm_area = np.trapz(voltages, times)

    # Load areas
    areas = r.TFile.Open("assets/r878_areas.root")
    # SPE area dist
    spe_areas = (areas.Get("ht")).Clone("spe")
    for b in range(0, spe_areas.GetSize()):
        # Kill any negative areas
        if spe_areas.GetBinContent(b) < 0: spe_areas.SetBinContent(b, 0)

    # Get random area        
    area = spe_areas.GetRandom()
    # Scale template y-values
    voltages *= area/prenorm_area
    if verbose: 
        postnorm_area = np.trapz(voltages, times)
        print("SPE Area: {0}, Original Template Area: {1}, Scaled Template Area: {2}".format(area, prenorm_area, postnorm_area))

    if extend is not None:
        nsamp = int(extend/new_dt)
        times = np.arange(0, new_dt*(times.size + 2*nsamp), new_dt)
        voltages = np.append(voltages, np.zeros(nsamp))
        voltages = np.append(np.zeros(nsamp), voltages)
        
    if noise is not None:
        voltages += np.random.normal(loc=0.0, scale=noise, size=voltages.size)

    return times, voltages

if __name__ == "__main__":
    times, voltages = GenerateSignal(output_samp_freq=2.0, noise=0.38, extend=150, verbose=True)

    import matplotlib.pyplot as plt
    plt.plot(times, voltages)
    plt.gca().set_xlim(0,400)
    plt.show()

