#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TMath.h"
#include "TColor.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "Drawing_Helper_Functions.h"
#include "BDT_Helper_Functions.h"

//------------------------------------------------------------------------------------------------------------------------------------------

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/fullreco_bdtresponse.root");
TTree *t1 = (TTree*)f1->Get("EventSelection");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/1.root");
TTree *t2 = (TTree*)f2->Get("EventSelection");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/2.root");
TTree *t3 = (TTree*)f3->Get("EventSelection");

TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/3.root");
TTree *t4 = (TTree*)f4->Get("EventSelection");

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareDistributions(TTree* pTree, int targetMultiplicity, std::string variableName, std::string typeName, std::vector<std::pair<std::string, int>> filterValuesSignal, std::vector<std::pair<std::string, int>> filterValuesBackground, std::vector<std::pair<std::string, float>> &bdtVector, int nBins, float lowerBound, float upperBound, int &counter)
{
    TH1F* pSignalHistogram = new TH1F("pSignalHistogram", "Signal Histogram", nBins, lowerBound, upperBound); 
    TH1F* pBackgroundHistogram = new TH1F("pBackgroundHistogram", "Background Histogram", nBins, lowerBound, upperBound); 

    filterValuesSignal.push_back(std::make_pair("RecoNeutrinoNumberAssociatedParticles", targetMultiplicity));
    filterValuesBackground.push_back(std::make_pair("RecoNeutrinoNumberAssociatedParticles", targetMultiplicity));

    CreateFilteredHistogram(pTree, pSignalHistogram, variableName, filterValuesSignal, bdtVector, typeName, counter, true);
    CreateFilteredHistogram(pTree, pBackgroundHistogram, variableName, filterValuesBackground, bdtVector, typeName, counter, true);

    Draw(pSignalHistogram, pBackgroundHistogram, kBlue, kRed, false, variableName + "_Distribution_SB_Comparison", variableName, "Fraction of Events", variableName + " Normalised Variable Distribution Comparison for Signal and Background", "Signal", "Background");

    delete pSignalHistogram;
    delete pBackgroundHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Signal_Background_Comparison(void)
{
    int targetMultiplicity(2);

    std::vector<std::string> variableNamesVector = {"ShortestPfoMinChiSquaredPerHit3D", "ShortestPfoLength"}; 
    std::vector<std::string> typeNamesVector = {"float", "float"}; 
    std::vector<int> numberBinsVector = {101, 101}; 
    std::vector<float> lowerBoundsVector = {-1.0, -1.0}; 
    std::vector<float> upperBoundsVector = {100.0, 100.0}; 
    std::vector<std::pair<std::string, int>> filterValuesSignal = {std::make_pair("ShortestPfoMCPDG", 13)};
    std::vector<std::pair<std::string, int>> filterValuesBackground = {std::make_pair("ShortestPfoMCPDG", 2212)};
    //std::vector<std::pair<std::string, float>> bdtVector = {std::make_pair("TMVA_BDT_AllVariables_N2_Contained", -0.05)};
    std::vector<std::pair<std::string, float>> bdtVector = {};
    int counter(0);

    for (int i = 0; i < static_cast<int>(variableNamesVector.size()); ++i)
        CompareDistributions(t1, targetMultiplicity, variableNamesVector.at(i), typeNamesVector.at(i), filterValuesSignal, filterValuesBackground, bdtVector, numberBinsVector.at(i), lowerBoundsVector.at(i), upperBoundsVector.at(i), counter);
}

//------------------------------------------------------------------------------------------------------------------------------------------
