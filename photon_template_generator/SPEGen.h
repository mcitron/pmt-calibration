#include "TFile.h"
#include "TString.h"

class SPE {

    /* Generate a SPE signal
     *
     * Member Variable Definitions:
     *   areas_file: TFile containing area distribution histogram
     *     - can be changed by specifying another path
     *   template_file: TFile containing signal template
     *     - can be changed by specifying another path
     *   output_sample_freq: sampling frequency of the output waveform (GHz)
     *   verbose: toggle printout of computed areas
     */

    TFile* areas_file;
    TFile* template_file;
    double input_sample_freq;
    double output_sample_freq;
    bool verbose;

    public:
        // Constructors
        SPE();
       ~SPE();
        SPE(TString, TString, double, double, bool);
        // Setters
        void SetAreasFile(TString);
        void SetTemplateFile(TString);
        void SetInputSampleFreq(double);
        void SetOutputSampleFreq(double);
        void SetVerbose(bool);
        // Utility
        TH1D Generate();

};
// Default Constructor
SPE::SPE() {
    areas_file = new TFile("assets/r878_areas.root"); 
    template_file = new TFile("assets/template_peak_70_75_2GHz.root");
    input_sample_freq = 2.5;
    output_sample_freq = 2.0;
    verbose = false;
}
// Overload Constructor
SPE::SPE(TString areas_path, TString template_path, double new_input_rate, double new_output_sample_freq, bool is_verbose) {
    areas_file = new TFile(areas_path); 
    template_file = new TFile(template_path);
    input_sample_freq = new_input_rate;
    output_sample_freq = new_output_sample_freq;
    verbose = is_verbose;
}

SPE::~SPE(){
    if(areas_file != NULL){
        areas_file->Close();
        delete areas_file;
    }
    if(template_file != NULL){
        template_file->Close();
        delete template_file;
    }
}

// Setters
void SPE::SetAreasFile(TString areas_path) {
    if(areas_file != NULL){
        areas_file->Close();
        delete areas_file;
    }
    areas_file = new TFile(areas_path); 
}
void SPE::SetTemplateFile(TString template_path) {
    if(template_file != NULL){
        template_file->Close();
        delete template_file;
    }
    template_file = new TFile(template_path);
}
void SPE::SetInputSampleFreq(double new_input_rate) {
    input_sample_freq = new_input_rate;    
}
void SPE::SetOutputSampleFreq(double new_output_sample_freq) {
    output_sample_freq = new_output_sample_freq;    
}
void SPE::SetVerbose(bool is_verbose) {
    verbose = is_verbose;    
}
