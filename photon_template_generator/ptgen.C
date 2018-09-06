// Root Includes
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

// Library Includes
#include <iostream>
#include <cmath>

class SPE {

    /* Generate a SPE signal
     *
     * Member Variable Definitions:
     *   areas_file: TFile containing area distribution histogram
     *     - can be changed by specifying another path
     *   template_file: TFile containing signal template
     *     - can be changed by specifying another path
     *   model: defines which PMT to use (defined above)
     *   output_sample_freq: sampling frequency of the output waveform (GHz)
     *   verbose: toggle printout of computed areas
     */

    TFile* areas_file;
    TFile* template_file;
    double input_rate;
    double output_sample_freq;
    bool verbose;

    public:
        // Constructors
        SPE();
        SPE(TString, TString, double, double, bool);
        // Setters
        void SetAreasFile(TString);
        void SetTemplateFile(TString);
        void SetInputRate(double);
        void SetOutputSampleFreq(double);
        void SetVerbose(bool);
        // Utility
        TH1D* Render();

};
// Default Constructor
SPE::SPE() {
    areas_file = new TFile("assets/r878_areas.root"); 
    template_file = new TFile("assets/template_peak_70_75_2GHz.root");
    input_rate = 2.5;
    output_sample_freq = 2.0;
    verbose = false;
}
// Overload Constructor
SPE::SPE(TString areas_path, TString template_path, double new_input_rate, double new_output_sample_freq, bool is_verbose) {
    areas_file = new TFile(areas_path); 
    template_file = new TFile(template_path);
    input_rate = new_input_rate;
    output_sample_freq = new_output_sample_freq;
    verbose = is_verbose;
}
// Setters
void SPE::SetAreasFile(TString areas_path) {
    areas_file = new TFile(areas_path); 
    return;
}
void SPE::SetTemplateFile(TString template_path) {
    template_file = new TFile(template_path);
    return;
}
void SPE::SetInputRate(double new_input_rate) {
    input_rate = new_input_rate;
    return;
}
void SPE::SetOutputSampleFreq(double new_output_sample_freq) {
    output_sample_freq = new_output_sample_freq;
    return;
}
void SPE::SetVerbose(bool is_verbose) {
    verbose = is_verbose;
    return;
}
// Render SPE Waveform
TH1D* SPE::Render() {
    // Set seed
    gRandom->SetSeed(0);

    // Load template
    double dt = 1.0 / input_rate;
    TH1D* h_template = (TH1D*)template_file->Get("temp")->Clone("h_template");
    int length = h_template->GetSize();
    double template_y[length];
    double template_x[length];
    for (unsigned int i = 0; i < length; i++) {
        template_y[i] = h_template->GetBinContent(i);
        template_x[i] = dt * i;
    }

    // Resample with new time values
    double new_dt = 1.0 / output_sample_freq;
    int new_length = (int)std::round(template_x[length - 1] / new_dt);
    // Fill times
    double times[new_length];
    for (unsigned int i = 0; i < new_length; i++) {
        times[i] = new_dt * i;
    }
    // Fill voltages
    TH1D* h_volts = new TH1D("spe", "Generated Signal", new_length, 0, times[new_length-1]);
    int j = 0; // index of template_x, template_y
    for (unsigned int i = 0; i < new_length; i++) {
        while (j < length && template_x[j] < times[i]) {
            j++;
        }
        if (j == 0 || j == length) {
            h_volts->SetBinContent(i, 0);
            h_volts->SetBinError(i, 0);
        }
        else {
            h_volts->SetBinContent(i, template_y[j-1] + (template_y[j] - template_y[j-1])/(template_x[j] - template_x[j-1]) * (times[i] - template_x[j-1]));
            h_volts->SetBinError(i, 0);
        }
    }
    // Get prenormalized area
    double prenorm_area = h_volts->Integral();

    // Load areas
    TH1D* h_areas = (TH1D*)areas_file->Get("ht")->Clone("h_areas");
    for (unsigned int i = 0; i < h_areas->GetSize(); i++) {
        if (h_areas->GetBinContent(i) < 0) {
            h_areas->SetBinContent(i, 0);
        }
    }
    // Get random area from distribution
    double area = h_areas->GetRandom();
    // Scale h_volts
    h_volts->Scale(area / prenorm_area);
    if (verbose) {
        double postnorm_area = h_volts->Integral();
        std::cout << "SPE Area: " << area << std::endl 
                  << "Original Template Area: " << prenorm_area << std::endl 
                  << "Scaled Template Area: " << postnorm_area << std::endl;
    }

    return h_volts;
}

void ptgen() {
    SPE* spe = new SPE();
    spe->SetVerbose(true);
    TH1D* out = spe->Render();
    out->Draw();
    return;
}
