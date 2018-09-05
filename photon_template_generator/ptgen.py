import ROOT as r
import numpy as np
import cPickle as pickle


def ptgen(verbose=False):

    # Set seed
    r.gRandom.SetSeed(0);

    # Load template
    with open("assets/template_peak_70_75_2GHz.pkl", "rb") as pin:
        # Get template y-array (Volts)
        template_y = pickle.load(pin)
        # Construct template x-array (ns)
        template_x = []
        x = 0.0
        for y in enumerate(template_y):
            template_x.append(x)
            x += 0.5
        template_x = np.array(template_x)
        template_area = np.trapz(template_y, template_x)

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
    template_y *= area/template_area
    if verbose: 
        new_template_area = np.trapz(template_y, template_x)
        print("SPE Area: {0}, Original Template Area: {1}, Scaled Template Area: {2}".format(area, template_area, new_template_area))

    return template_x, template_y

if __name__ == "__main__":
    ptgen(verbose=True)
