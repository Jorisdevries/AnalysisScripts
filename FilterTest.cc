#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRandom3.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>

#include "Drawing_Helper_Functions.h"

//------------------------------------------------------------------------------------------------------------------------------------------

void FilterTest( void )
{
    TFile *f = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/filter_test/filter_test_2.root");
    TTree *t1 = (TTree*)f->Get("FitTree");

    gStyle->SetOptStat(0);

    /*
        float fractionFilteredBelongingToMuon(((float)nMuonHitsFiltered)/nTotalHitsFiltered), fractionFilteredBelongingToNonMuon(((float)nNonMuonHitsFiltered)/nTotalHitsFiltered),
        fractionMuonHitsFiltered(((float)nMuonHitsFiltered)/nMuonHits), fractionNonMuonHitsFiltered(((float)nNonMuonHitsFiltered)/nNonMuonHits), fractionClusterHitsFiltered(((float)nTotalHitsFiltered)/totalNumberHits);
    */
    
    std::vector<std::string> variableNames = {"fractionFilteredBelongingToNonMuon", "fractionNonMuonHitsFiltered", "fractionFilteredBelongingToMuon", "fractionMuonHitsFiltered", "fractionImpureHits"};
    std::vector<std::string> xTitles  = {"Filter Purity P", "Filter Completeness C", "Fraction Cluster Hits Filtered", "Fraction Pure Hits Retained", "Fraction of Cluster Hits Classified Impure"};
    std::vector<std::string> yTitles = {"Fraction of Entries", "Fraction of Entries", "Fraction of Entries", "Fraction of Entries", "Fraction of Entries"};
    std::vector<std::string> plotTitles = {"Normalised Distribution of Filter Purity P", "Normalised Distribution of Filter Completeness C", "Normalised Distribution of Fraction Filtered Hits Belonging to Muon", "Normalised Distribution of Fraction Pure Hits Retained", "Normalised Distribution of Fraction of Cluster Hits Classified Impure"};
    
    for (int i = 1; i < static_cast<int>(variableNames.size()); ++i)
    {
        TH1F *pHistogram = new TH1F(variableNames.at(i).c_str(),"", 100, 0.0, 1.0);
        CreateFloatHistogramNoise(t1, pHistogram, variableNames.at(i));
        Draw(pHistogram, kBlue, true, variableNames.at(i), xTitles.at(i), yTitles.at(i), plotTitles.at(i));

        delete pHistogram;
    }

    //TH1F *pHistogram2 = new TH1F("pHistogram2", "", 100, 0, 400);
    //CreateHistogram(t1, pHistogram2, "nTotalHitsFiltered");
    //Draw(pHistogram2, kBlue, true, "nTotalHitsFiltered", "Total Number Cluster Hits Filtered", "Fraction of Entries", "Normalised Distribution of Total Number Cluster Hits Filtered");
}
