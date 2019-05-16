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
#include "TStyle.h"
#include "TColor.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>

#include "Drawing_Helper_Functions.h"

//------------------------------------------------------------------------------------------------------------------------------------------

TFile *f = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/alldata_final.root");
TTree *t1 = (TTree*)f->Get("DirectionFitTree");
int nentries = (Int_t)t1->GetEntries();

//------------------------------------------------------------------------------------------------------------------------------------------

float CountFractionSignedEntries(TH1F* histogram)
{
    int nPositiveEntries(0), nNegativeEntries(0);

    for (int i = 0; i <= histogram->GetNbinsX() + 1; ++i)
    {
        if (histogram->GetBinCenter(i) > 0)
            nPositiveEntries += histogram->GetBinContent(i);
        else
            nNegativeEntries += histogram->GetBinContent(i);
    }

    std::cout << "Fraction of positive entries: " << static_cast<float>(nPositiveEntries)/histogram->GetEntries() << std::endl;
    std::cout << "Fraction of negative entries: " << static_cast<float>(nNegativeEntries)/histogram->GetEntries() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquared_Modern(void)
{
    TH1F *forwardsChiSquaredPerHit = new TH1F("forwardsChiSquaredPerHit","", 100, 0.0, 10.0);
    TH1F *backwardsChiSquaredPerHit = new TH1F("backwardsChiSquaredPerHit","", 100, 0.0, 10.0);

    TH1F *forwardsDeltaChiSquaredPerHit = new TH1F("forwardsDeltaChiSquaredPerHit","", 100, -50.0, 50.0);
    TH1F *backwardsDeltaChiSquaredPerHit = new TH1F("backwardsDeltaChiSquaredPerHit","", 100, -50.0, 50.0);

    /*
    TH1F *upwardsChiSquaredPerHit = new TH1F("upwardsChiSquaredPerHit","", 75, 0.0, 5.0);
    TH1F *downwardsChiSquaredPerHit = new TH1F("downwardsChiSquaredPerHit","", 75, 0.0, 5.0);

    TH1F *upwardsDeltaChiSquared = new TH1F("upwardsDeltaChiSquared","", 75, -4000.0, 4000.0);
    TH1F *downwardsDeltaChiSquared = new TH1F("downwardsDeltaChiSquared","", 75, -4000.0, 4000.0);

    TH1F *upwardsDeltaChiSquaredPerHit = new TH1F("upwardsDeltaChiSquaredPerHit","", 75, -25.0, 25.0);
    TH1F *downwardsDeltaChiSquaredPerHit = new TH1F("downwardsDeltaChiSquaredPerHit","", 75, -25.0, 25.0);
    */

    std::vector<std::pair<std::string, int>> forwardsFilterValues = {std::make_pair("TrueForwards", 1)};
    std::vector<std::pair<std::string, int>> backwardsFilterValues = {std::make_pair("TrueForwards", 0)};

    CreateFilteredHistogram(t1, forwardsChiSquaredPerHit, "ForwardsChiSquaredPerHit", forwardsFilterValues, "float"); 
    CreateFilteredHistogram(t1, backwardsChiSquaredPerHit, "BackwardsChiSquaredPerHit", backwardsFilterValues, "float"); 

    CountFractionSignedEntries(forwardsChiSquaredPerHit);
    CountFractionSignedEntries(backwardsChiSquaredPerHit);

    //std::cout << "Maximum forwardsChiSquaredPerHit bin center: " << forwardsChiSquaredPerHit->GetXaxis()->GetBinCenter(forwardsChiSquaredPerHit->GetMaximumBin()) << std::endl;
    //std::cout << "Maximum forwardsChiSquaredPerHit bin : " << (forwardsChiSquaredPerHit->GetMaximumBin()) << std::endl;
    //std::cout << "Maximum backwardsChiSquaredPerHit bin center: " << backwardsChiSquaredPerHit->GetXaxis()->GetBinCenter(backwardsChiSquaredPerHit->GetMaximumBin()) << std::endl;

    CreateFilteredHistogram(t1, forwardsDeltaChiSquaredPerHit, "DeltaChiSquaredPerHit", forwardsFilterValues, "float"); 
    CreateFilteredHistogram(t1, backwardsDeltaChiSquaredPerHit, "DeltaChiSquaredPerHit", backwardsFilterValues, "float"); 

    CountFractionSignedEntries(forwardsDeltaChiSquaredPerHit);
    CountFractionSignedEntries(backwardsDeltaChiSquaredPerHit);

    Draw(forwardsChiSquaredPerHit, kBlue, "forwardsChiSquaredPerHit", "#chi^{2}_{F}/N", "Number of Entries", "True Forwards Muons: #chi^{2}_{F}/N");
    Draw(backwardsChiSquaredPerHit, kRed, "backwardsChiSquaredPerHit", "#chi^{2}_{B}/N", "Number of Entries", "True Backwards Muons: #chi^{2}_{B}/N");

    Draw(forwardsDeltaChiSquaredPerHit, kBlue, "forwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N", "Number of Entries", "True Forwards Muons: #Delta#chi^{2}/N");
    Draw(backwardsDeltaChiSquaredPerHit, kRed, "backwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N", "Number of Entries", "True Backwards Muons: #Delta#chi^{2}/N");

    gStyle->SetPalette(55);

    TH2F *forwardsChiSquaredPerHitWithN = new TH2F("forwardsChiSquaredPerHitWithN","", 60, 0.0, 10.0, 60, 0.0, 1000.0);
    CreateFiltered2DHistogram(t1, forwardsChiSquaredPerHitWithN, "ForwardsChiSquaredPerHit", "NumberHits", forwardsFilterValues, "float", "int"); 

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    forwardsChiSquaredPerHitWithN->SetXTitle("#chi^{2}_{F}/N");
    forwardsChiSquaredPerHitWithN->SetYTitle("Numbe of Hits N");
    forwardsChiSquaredPerHitWithN->SetTitle("True Forwards Muons: #chi^{2}_{F}/N with N");
    forwardsChiSquaredPerHitWithN->GetYaxis()->SetTitleOffset(1.3);
    forwardsChiSquaredPerHitWithN->SetStats(kFALSE);
    forwardsChiSquaredPerHitWithN->Draw("COLZ");
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/ChiSquared/Forwards_2D.pdf");
    canvas->SaveAs(savePath.c_str());
    delete canvas;

    TH2F *backwardsChiSquaredPerHitWithN = new TH2F("backwardsChiSquaredPerHitWithN","", 60, 0.0, 10.0, 60, 0.0, 1000.0);
    CreateFiltered2DHistogram(t1, backwardsChiSquaredPerHitWithN, "BackwardsChiSquaredPerHit", "NumberHits", backwardsFilterValues, "float", "int"); 

    TCanvas *canvas1 = new TCanvas("Canvas", "Canvas", 900, 600);
    backwardsChiSquaredPerHitWithN->SetXTitle("#chi^{2}_{B}/N");
    backwardsChiSquaredPerHitWithN->SetYTitle("Numbe of Hits N");
    backwardsChiSquaredPerHitWithN->SetTitle("True Backwards Muons: #chi^{2}_{B}/N with N");
    backwardsChiSquaredPerHitWithN->GetYaxis()->SetTitleOffset(1.3);
    backwardsChiSquaredPerHitWithN->SetStats(kFALSE);
    backwardsChiSquaredPerHitWithN->Draw("COLZ");
    std::string savePath2("/usera/jjd49/pandora_direction/Scripts/Figures/ChiSquared/Backwards_2D.pdf");
    canvas1->SaveAs(savePath2.c_str());
    delete canvas1;

    /*
    Draw(upwardsChiSquaredPerHit, kBlue, "upwardsChiSquaredPerHit", "#chi^{2}_{U}/N", "Number of Entries", "True Upwards Muons: #chi^{2}_{U}/N");
    Draw(upwardsDeltaChiSquared, kBlue, "upwardsDeltaChiSquared", "#Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}", "Number of Entries", "True Upwards Muons: #Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}");
    Draw(upwardsDeltaChiSquaredPerHit, kBlue, "upwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N", "Number of Entries", "True Upwards Muons: #Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N");

    Draw(downwardsChiSquaredPerHit, kRed, "downwardsChiSquaredPerHit", "#chi^{2}_{D}/N", "Number of Entries", "True Downwards Muons: #chi^{2}_{D}/N");
    Draw(downwardsDeltaChiSquared, kRed, "downwardsDeltaChiSquared", "#Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}", "Number of Entries", "True Downwards Muons: #Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}");
    Draw(downwardsDeltaChiSquaredPerHit, kRed, "downwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N", "Number of Entries", "True Downwards Muons: #Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N");

    */

    int     trueDirection;
    float   forwardsChiSquaredPerHitValue;
    float   backwardsChiSquaredPerHitValue;
    int     numberHits;
    int     fileIdentifier;
    int     eventNumber;
    
    t1->SetBranchAddress("TrueForwards",                                             &trueDirection);
    t1->SetBranchAddress("ForwardsChiSquaredPerHit",                                             &forwardsChiSquaredPerHitValue);
    t1->SetBranchAddress("BackwardsChiSquaredPerHit",                                             &backwardsChiSquaredPerHitValue);
    t1->SetBranchAddress("NumberHits",                                             &numberHits);
    t1->SetBranchAddress("FileIdentifier",                                             &fileIdentifier);
    t1->SetBranchAddress("EventNumber",                                             &eventNumber);

    for (int j = 1; j < t1->GetEntries(); j++)
    {   
        t1->GetEntry(j);

        if (trueDirection == 1 && forwardsChiSquaredPerHitValue >= 5.0)
            std::cout << "Bad forwards event: " << fileIdentifier << ":" << eventNumber << " (" << forwardsChiSquaredPerHitValue << ", N: " << numberHits << ")" << std::endl;
        if (trueDirection == 0 && backwardsChiSquaredPerHitValue >= 5.0)
            std::cout << "Bad backwards event: " << fileIdentifier << ":" << eventNumber << " (" << backwardsChiSquaredPerHitValue << ", N: " << numberHits << ")" << std::endl;
    }   

    //--------------------------------------------------------------------------------------------------------------------------------------
}
