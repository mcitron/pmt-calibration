// Root Includes
#include "TTree.h"

#include "SPEGen.h"

// Library Includes
#include <iostream>
#include <cmath>

// Generate SPE Waveform
TH1D SPE::Generate() {

    // Load template
    double dt = 1.0 / input_sample_freq;
    TH1D* h_template = (TH1D*)template_file->Get("temp")->Clone("h_template");
    int length = h_template->GetSize();
    double *template_y = new double[length];
    double *template_x = new double[length];
    for (unsigned int i = 0; i < length; i++) {
        template_y[i] = h_template->GetBinContent(i);
        template_x[i] = dt * i;
    }

    // Resample with new time values
    double new_dt = 1.0 / output_sample_freq;
    int new_length = (int)std::round(template_x[length - 1] / new_dt);
    // Fill times
    double *times = new double[new_length];
    for (unsigned int i = 0; i < new_length; i++) {
        times[i] = new_dt * i;
    }
    // Fill voltages
    TH1D h_volts = TH1D("spe", "Generated Signal", new_length, 0, times[new_length-1] + times[1]-times[0]);
    int j = 0; // index of template_x, template_y
    for (unsigned int i = 0; i < new_length; i++) {
        while (j < length && template_x[j] < times[i]) {
            j++;
        }
        if (j == 0 || j == length) {
            h_volts.SetBinContent(i+1, 0);
            h_volts.SetBinError(i+1, 0);
        }
        else {
            h_volts.SetBinContent(i+1, template_y[j-1] + (template_y[j] - template_y[j-1])/(template_x[j] - template_x[j-1]) * (times[i] - template_x[j-1]));
            h_volts.SetBinError(i+1, 0);
        }
    }
    // Get prenormalized area
    double prenorm_area = h_volts.Integral();

    // Load areas
    TH1D* h_areas = (TH1D*)areas_file->Get("ht");
    for (unsigned int i = 0; i < h_areas->GetSize(); i++) {
        if (h_areas->GetBinContent(i) < 0) {
            h_areas->SetBinContent(i, 0);
        }
    }
    // Get random area from distribution
    double area = h_areas->GetRandom();
    // Scale h_volts
    h_volts.Scale(area / prenorm_area);
    if (verbose) {
        double postnorm_area = h_volts.Integral();
        std::cout << "SPE Area: " << area << std::endl 
                  << "Original Template Area: " << prenorm_area << std::endl 
                  << "Scaled Template Area: " << postnorm_area << std::endl;
    }

    delete [] template_x;
    delete [] template_y;
    delete [] times;

    return h_volts;
}

void SPEGen() {
    gRandom->SetSeed(0);

    SPE spe = SPE();
    spe.SetVerbose(true);
    spe.SetOutputSampleFreq(1.5); // generate at 1.5 GHz
    TH1D out = spe.Generate();
    out.Draw();
    return;
}
