#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLegend.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>

#include "Drawing_Helper_Functions.h"

//------------------------------------------------------------------------------------------------------------------------------------------

TFile *f = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/bpf_data.root");
TTree *t1 = (TTree*)f->Get("BraggPeakFilterTree");

int nentries = (Int_t)t1->GetEntries();

//------------------------------------------------------------------------------------------------------------------------------------------

float GetSymmetricCut(TH1F* braggPeakNearestRatio, TH1F* nonBraggPeakNearestRatio, bool symmetric)
{
    float nearestSymmetricCut(0.f);
    for (int bin = 1; bin <= braggPeakNearestRatio->GetXaxis()->GetNbins(); bin++)
    {
        float braggBinEntry(braggPeakNearestRatio->GetBinContent(bin));
        float nextBraggBinEntry(braggPeakNearestRatio->GetBinContent(bin+1));
        float nextNextBraggBinEntry(braggPeakNearestRatio->GetBinContent(bin+2));
        float nonBraggBinEntry(nonBraggPeakNearestRatio->GetBinContent(bin));
        float nextNonBraggBinEntry(nonBraggPeakNearestRatio->GetBinContent(bin+1));
        float nextNextNonBraggBinEntry(nonBraggPeakNearestRatio->GetBinContent(bin+2));

        if (braggPeakNearestRatio->GetBinCenter(bin) > 1.4 && nonBraggBinEntry/braggBinEntry > 1.0 && nextNonBraggBinEntry/nextBraggBinEntry > 1.0)
        {
            if (symmetric)
                nearestSymmetricCut = std::abs(1.0 - braggPeakNearestRatio->GetBinCenter(bin));
            else
                nearestSymmetricCut = braggPeakNearestRatio->GetBinCenter(bin);
            break;
        }
    }

    return nearestSymmetricCut;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BraggPeakFilter( void )
{
    TH1F *braggPeakNearestRatio = new TH1F("braggPeakNearestRatio","", 100, 0.0, 3.0);
    TH1F *nonBraggPeakNearestRatio = new TH1F("nonBraggPeakNearestRatio","", 100, 0.0, 3.0);

    TH1F *braggPeakPlusMinusNRatio = new TH1F("braggPeakPlusMinusNRatio","", 100, 0.0, 3.0);
    TH1F *nonBraggPeakPlusMinusNRatio = new TH1F("nonBraggPeakPlusMinusNRatio","", 100, 0.0, 3.0);

    TH1F *braggPeakDistanceFromAverageCharge = new TH1F("braggPeakDistanceFromAverageCharge","", 100, 0.0, 3.0);
    TH1F *nonBraggPeakDistanceFromAverageCharge = new TH1F("nonBraggPeakDistanceFromAverageCharge","", 100, 0.0, 3.0);

    TH1F *braggPeakDistanceFromTrackBody = new TH1F("braggPeakDistanceFromTrackBody","", 100, 0.0, 3.0);
    TH1F *nonBraggPeakDistanceFromTrackBody = new TH1F("nonBraggPeakDistanceFromTrackBody","", 100, 0.0, 3.0);

    std::vector<std::pair<std::string, int>> braggFilterValues = {std::make_pair("IsBraggPeakHit", 1)};
    std::vector<std::pair<std::string, int>> nonBraggFilterValues = {std::make_pair("IsBraggPeakHit", 0)};

    CreateFilteredHistogram(t1, braggPeakNearestRatio, "NearestNeighbourRatio", braggFilterValues, "float");
    CreateFilteredHistogram(t1, nonBraggPeakNearestRatio, "NearestNeighbourRatio", nonBraggFilterValues, "float");
    Draw(braggPeakNearestRatio, nonBraggPeakNearestRatio, kGreen, kRed, true, "NearestRatioComparison", "max #left(#tilde{Q}_{i}/#tilde{Q}_{i #pm 1}#right)", "Fraction of Entries", "Comparison of max #left(#tilde{Q}_{i}/#tilde{Q}_{i #pm 1}#right) for Bragg and Non-Bragg Hits", "Bragg Peak Hits", "Non-Bragg Peak Hits");

    CreateFilteredHistogram(t1, braggPeakPlusMinusNRatio, "ThirdNearestNeighbourRatio", braggFilterValues, "float");
    CreateFilteredHistogram(t1, nonBraggPeakPlusMinusNRatio, "ThirdNearestNeighbourRatio", nonBraggFilterValues, "float");
    Draw(braggPeakPlusMinusNRatio, nonBraggPeakPlusMinusNRatio, kGreen, kRed, true, "PlusMinusNRatioComparison", "max #left(#tilde{Q}_{i}/#tilde{Q}_{i #pm 3}#right)", "Fraction of Entries", "Comparison of max #left(#tilde{Q}_{i}/#tilde{Q}_{i #pm 3}#right) for Bragg and Non-Bragg Hits", "Bragg Peak Hits", "Non-Bragg Peak Hits");

    for (int i = 1; i < braggPeakPlusMinusNRatio->GetNbinsX(); ++i)
    {
        std::cout << i << std::endl;
        std::cout << static_cast<float>(braggPeakPlusMinusNRatio->GetBinContent(i))/nonBraggPeakPlusMinusNRatio->GetBinContent(i) << std::endl;
        std::cout << braggPeakPlusMinusNRatio->GetBinCenter(i) << std::endl;
        std::cout << "-------------------" << std::endl;
    }

    int crossoverBin = 58;
    std::cout << "Integral to cross-over for Bragg hits: " << braggPeakPlusMinusNRatio->Integral(1, crossoverBin) << std::endl; 
    std::cout << "Integral to cross-over for Bragg hits: " << nonBraggPeakPlusMinusNRatio->Integral(1, crossoverBin) << std::endl; 

    CreateFilteredHistogram(t1, braggPeakDistanceFromAverageCharge, "DistanceFromAverageCharge", braggFilterValues, "float");
    CreateFilteredHistogram(t1, nonBraggPeakDistanceFromAverageCharge, "DistanceFromAverageCharge", nonBraggFilterValues, "float");
    Draw(braggPeakDistanceFromAverageCharge, nonBraggPeakDistanceFromAverageCharge, kGreen, kRed, true, "DistanceFromAverageChargeComparison", "#left|#tilde{Q}_{i} - #LT#tilde{Q}#GT#right|", "Fraction of Entries", "Comparison of #left|#tilde{Q}_{i} - #LT#tilde{Q}#GT#right| for Bragg and Non-Bragg Hits", "Bragg Peak Hits", "Non-Bragg Peak Hits");

    CreateFilteredHistogram(t1, braggPeakDistanceFromTrackBody, "DistanceFromTrackBody", braggFilterValues, "float");
    CreateFilteredHistogram(t1, nonBraggPeakDistanceFromTrackBody, "DistanceFromTrackBody", nonBraggFilterValues, "float");
    Draw(braggPeakDistanceFromTrackBody, nonBraggPeakDistanceFromTrackBody, kGreen, kRed, true, "DistanceFromTrackBodyComparison", "#left|#tilde{Q}_{i} - #LT#tilde{Q}#GT#right|", "Fraction of Entries", "Comparison of #left|#tilde{Q}_{i} - #LT#tilde{Q}#GT#right| for Bragg and Non-Bragg Hits", "Bragg Peak Hits", "Non-Bragg Peak Hits");

    /*
    TH1F *c_clone1 = (TH1F*)braggPeakNearestRatio->Clone(); 
    int lower1(c_clone1->GetXaxis()->FindBin(1.0 - nearestSymmetricCut)), upper1(c_clone1->GetXaxis()->FindBin(1.0 + nearestSymmetricCut));
    c_clone1->SetFillColor(kGreen+3);
    c_clone1->SetFillStyle(3003);
    c_clone1->GetXaxis()->SetRange(lower1,upper1);
    c_clone1->Draw("histsame");
    */
}
