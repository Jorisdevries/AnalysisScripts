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

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/data_mc_comparison/mcc_87_tiny_bdtresponse.root");
TTree *t1 = (TTree*)f1->Get("EventSelection");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/data_mc_comparison/data_tiny_bdtresponse.root");
TTree *t2 = (TTree*)f2->Get("EventSelection");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/data_mc_comparison/cr_data_small.root");
TTree *t3 = (TTree*)f3->Get("EventSelection");

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareDistributions(TTree* pMCTree, TTree* pBNBDataTree, TTree* pEXTBNBDataTree, int targetMultiplicity, std::string variableName, std::string variableToken, std::string typeName, std::vector<std::pair<std::string, int>> filterValues, std::vector<std::pair<std::string, float>> &bdtVector, int nBins, float lowerBound, float upperBound, int &counter)
{
    TH1F* pMCHistogram = new TH1F("pMCHistogram", "MC Histogram", nBins, lowerBound, upperBound); 
    TH1F* pBNBDataHistogram = new TH1F("pBNBDataHistogram", "BNB Data Histogam", nBins, lowerBound, upperBound); 
    TH1F* pEXTBNBDataHistogram = new TH1F("pEXTBNBDataHistogram", "EXTBNB Data Histogam", nBins, lowerBound, upperBound); 

    filterValues.push_back(std::make_pair("RecoNeutrinoMultiplicity", targetMultiplicity));

    CreateFilteredHistogram(pMCTree, pMCHistogram, variableName, filterValues, bdtVector, typeName, counter, true);
    CreateFilteredHistogram(pBNBDataTree, pBNBDataHistogram, variableName, filterValues, bdtVector, typeName, counter, true);
    CreateFilteredHistogram(pEXTBNBDataTree, pEXTBNBDataHistogram, variableName, filterValues, bdtVector, typeName, counter, true);

    pMCHistogram->Scale(1/(pMCHistogram->GetEntries()));
    pBNBDataHistogram->Scale(1/(pBNBDataHistogram->GetEntries()));
    pEXTBNBDataHistogram->Scale(1/(pEXTBNBDataHistogram->GetEntries()));

    TH1F *pNormalisedDataHistogram = (TH1F*) pBNBDataHistogram->Clone();
    pNormalisedDataHistogram->Add(pEXTBNBDataHistogram, -1);

    Draw(pMCHistogram, pNormalisedDataHistogram, kBlue, kRed, false, "N" + std::to_string(targetMultiplicity) + "_" + variableName + "_BDT_" + (bdtVector.size() != 0 ? bdtVector.front().first : "") + "_Distribution_MC_Data_Comparison", variableToken, "Fraction of Events", "N=" + std::to_string(targetMultiplicity) + " " + variableToken + " Normalised Variable Distribution Comparison for MC and Data" + (bdtVector.size() != 0 ? " Post BDT" : ""), "MC", "Data", "DataMCComparison");

    delete pMCHistogram;
    delete pBNBDataHistogram;
    delete pEXTBNBDataHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Data_MC_Comparison(void)
{
    int targetMultiplicity(1);
    std::string prefix("Longest");

    std::vector<std::string> variableNamesVector = {"NeutrinoMomentumX", "NeutrinoMomentumY", "NeutrinoMomentumZ", "OpeningAngle", prefix + "PfoCharge", prefix + "PfoLength", prefix + "PfoTrackProbability", prefix + "PfoPolarAngle", prefix + "PfoCorrectedMinX", prefix + "PfoCorrectedMaxX", prefix + "PfoCorrectedMaxY", prefix + "PfoUpDownDeltaChiSquaredPerHit", prefix + "PfoDeltaChiSquaredPerHit", prefix + "PfoMinChiSquaredPerHit", prefix + "PfoFitChargeRangeOverLength", prefix + "PfoMinSumQW", prefix + "PfoRecoDeltaY", prefix + "PfoRecoLowestTenCmTotalCharge", prefix + "PfoMinChiSquaredPerHit3D"}; 
    std::vector<std::string> variableTokensVector = {"|#hat{p}_{#nu} #upoint #hat{x}|", "|#hat{p}_{#nu} #upoint #hat{y}|", "|#hat{p}_{#nu} #upoint #hat{z}|", "Opening Angle #zeta", "Q_{PFO, l}", "L_{PFO, l}", "P_{t, l}", "#theta_{l}", "x_{min, l}", "x_{max, l}", "y_{max, l}", "#Delta#chi^{2}_{DU,l}", "#Delta#chi^{2}_{FB,l}", "#chi^{2}_{min,l}", "Q_{r,l}/L_{l}", "Q_{min,l}", "#Deltay_{l}", "Q_{10cm,l}", "#chi^{2}_{proton}/N"}; 
    std::vector<std::string> typeNamesVector = {"float", "float", "float", "float", "float", "float", "float", "float", "float", "float", "float", "float","float","float","float","float","float","float", "float"}; 
    std::vector<int> numberBinsVector = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50}; 
    std::vector<float> lowerBoundsVector = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 80.0, -5.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
    std::vector<float> upperBoundsVector = {1.2, 1.2, 1.2, 3.5, 3e5, 500.0, 1.2, 3.5, 300.0, 300.0, 150.0, 5.0, 5.0, 5.0, 0.3, 20.0, 250.0, 25e3, 100.0}; 
    std::vector<std::pair<std::string, int>> filterValues = {};
    //std::vector<std::pair<std::string, float>> bdtVector = {std::make_pair("BDT_N1_Contained_CleanedVariables", -0.05)};
    std::vector<std::pair<std::string, float>> bdtVector = {};
    int counter(0);

    for (int i = 0; i < static_cast<int>(variableNamesVector.size()); ++i)
        CompareDistributions(t1, t2, t3, targetMultiplicity, variableNamesVector.at(i), variableTokensVector.at(i), typeNamesVector.at(i), filterValues, bdtVector, numberBinsVector.at(i), lowerBoundsVector.at(i), upperBoundsVector.at(i), counter);

}

//------------------------------------------------------------------------------------------------------------------------------------------
