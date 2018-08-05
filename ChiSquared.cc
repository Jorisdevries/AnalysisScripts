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
#include "TColor.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH1F* pHistogram, EColor colour, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram->SetXTitle(xTitle.c_str());
    pHistogram->SetYTitle(yTitle.c_str());
    pHistogram->SetTitle(plotTitle.c_str());
    pHistogram->GetYaxis()->SetTitleOffset(1.0);
    pHistogram->SetLineColor(colour);
    pHistogram->Draw();
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/chisquared/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH2F* pHistogram, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram->SetXTitle(xTitle.c_str());
    pHistogram->SetYTitle(yTitle.c_str());
    pHistogram->SetTitle(plotTitle.c_str());
    pHistogram->GetYaxis()->SetTitleOffset(1.0);
    pHistogram->Draw("COLZ");
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/chisquared/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquared(void)
{
    TFile *f = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/catroots/1.root");
    
    TTree *t1 = (TTree*)f->Get("PFO");
    int nentries = (Int_t)t1->GetEntries();
    
    TH1F *forwardsChiSquaredPerHit = new TH1F("forwardsChiSquaredPerHit","", 75, 0.0, 5.0);
    TH1F *backwardsChiSquaredPerHit = new TH1F("backwardsChiSquaredPerHit","", 75, 0.0, 5.0);
    TH1F *upwardsChiSquaredPerHit = new TH1F("upwardsChiSquaredPerHit","", 75, 0.0, 5.0);
    TH1F *downwardsChiSquaredPerHit = new TH1F("downwardsChiSquaredPerHit","", 75, 0.0, 5.0);

    TH1F *forwardsDeltaChiSquared = new TH1F("forwardsDeltaChiSquared","", 75, -4000.0, 4000.0);
    TH1F *backwardsDeltaChiSquared = new TH1F("backwardsDeltaChiSquared","", 75, -4000.0, 4000.0);
    TH1F *upwardsDeltaChiSquared = new TH1F("upwardsDeltaChiSquared","", 75, -4000.0, 4000.0);
    TH1F *downwardsDeltaChiSquared = new TH1F("downwardsDeltaChiSquared","", 75, -4000.0, 4000.0);

    TH1F *forwardsDeltaChiSquaredPerHit = new TH1F("forwardsDeltaChiSquaredPerHit","", 75, -25.0, 25.0);
    TH1F *backwardsDeltaChiSquaredPerHit = new TH1F("backwardsDeltaChiSquaredPerHit","", 75, -25.0, 25.0);
    TH1F *upwardsDeltaChiSquaredPerHit = new TH1F("upwardsDeltaChiSquaredPerHit","", 75, -25.0, 25.0);
    TH1F *downwardsDeltaChiSquaredPerHit = new TH1F("downwardsDeltaChiSquaredPerHit","", 75, -25.0, 25.0);
    
    int   MCDirection;
    int   mcForwards;
    int   mcDownwards;
    float MCPhi;
    float MCTheta;
    float BeginX;
    float BeginY;
    float BeginZ;
    float EndX;
    float EndY;
    float EndZ;
    float Probability;
    int   Hypothesis;
    float DeltaChiSquaredPerHit;
    float MinChiSquaredPerHit;
    float ForwardsChiSquared;
    float BackwardsChiSquared;
    float ForwardsChiSquaredPerHit;
    float BackwardsChiSquaredPerHit;
    float UpDownDeltaChiSquaredPerHit;
    float UpwardsChiSquared;
    float DownwardsChiSquared;
    float UpwardsChiSquaredPerHit;
    float DownwardsChiSquaredPerHit;
    int   NumberHits;
    int   MCIsClusterTwoParticles;
    int   FileIdentifier;
    int   EventNumber;
    
    t1->SetBranchAddress("MCDirection", &MCDirection);
    t1->SetBranchAddress("MCForwards", &mcForwards);
    t1->SetBranchAddress("MCDownwards", &mcDownwards);
    t1->SetBranchAddress("MCPhi", &MCPhi);
    t1->SetBranchAddress("MCTheta", &MCTheta);
    t1->SetBranchAddress("BeginX", &BeginX);
    t1->SetBranchAddress("BeginY", &BeginY);
    t1->SetBranchAddress("BeginZ", &BeginZ);
    t1->SetBranchAddress("EndX", &EndX);
    t1->SetBranchAddress("EndY", &EndY);
    t1->SetBranchAddress("EndZ", &EndZ);
    t1->SetBranchAddress("Probability", &Probability);
    t1->SetBranchAddress("Hypothesis", &Hypothesis);
    t1->SetBranchAddress("MinChiSquaredPerHit", &MinChiSquaredPerHit);
    t1->SetBranchAddress("DeltaChiSquaredPerHit", &DeltaChiSquaredPerHit);
    t1->SetBranchAddress("ForwardsChiSquared", &ForwardsChiSquared);
    t1->SetBranchAddress("BackwardsChiSquared", &BackwardsChiSquared);
    t1->SetBranchAddress("ForwardsChiSquaredPerHit", &ForwardsChiSquaredPerHit);
    t1->SetBranchAddress("BackwardsChiSquaredPerHit", &BackwardsChiSquaredPerHit);
    t1->SetBranchAddress("UpDownDeltaChiSquaredPerHit", &UpDownDeltaChiSquaredPerHit);
    t1->SetBranchAddress("UpwardsChiSquared", &UpwardsChiSquared);
    t1->SetBranchAddress("DownwardsChiSquared", &DownwardsChiSquared);
    t1->SetBranchAddress("UpwardsChiSquaredPerHit", &UpwardsChiSquaredPerHit);
    t1->SetBranchAddress("DownwardsChiSquaredPerHit", &DownwardsChiSquaredPerHit);
    t1->SetBranchAddress("NumberHits", &NumberHits);
    t1->SetBranchAddress("MCIsClusterTwoParticles", &MCIsClusterTwoParticles);
    t1->SetBranchAddress("FileIdentifier", &FileIdentifier);
    t1->SetBranchAddress("EventNumber", &EventNumber);
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::cout << "There are " << nentries << " entries in the tree" << std::endl;
    
    for (int i = 1; i < nentries; i++)
    {
        t1->GetEntry(i);

        if (mcForwards == 1)
        {
            forwardsChiSquaredPerHit->Fill(ForwardsChiSquaredPerHit);
            forwardsDeltaChiSquared->Fill(DeltaChiSquaredPerHit * NumberHits);
            forwardsDeltaChiSquaredPerHit->Fill(DeltaChiSquaredPerHit);
        }
        else
        {
            backwardsChiSquaredPerHit->Fill(BackwardsChiSquaredPerHit);
            backwardsDeltaChiSquared->Fill(DeltaChiSquaredPerHit * NumberHits);
            backwardsDeltaChiSquaredPerHit->Fill(DeltaChiSquaredPerHit);
        }

        if (mcDownwards == 1)
        {
            downwardsChiSquaredPerHit->Fill(DownwardsChiSquaredPerHit);
            downwardsDeltaChiSquared->Fill(UpDownDeltaChiSquaredPerHit * NumberHits);
            downwardsDeltaChiSquaredPerHit->Fill(UpDownDeltaChiSquaredPerHit);
        }
        else
        {
            upwardsChiSquaredPerHit->Fill(UpwardsChiSquaredPerHit);
            upwardsDeltaChiSquared->Fill(UpDownDeltaChiSquaredPerHit * NumberHits);
            upwardsDeltaChiSquaredPerHit->Fill(UpDownDeltaChiSquaredPerHit);
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------------------------------

    Draw(forwardsChiSquaredPerHit, kBlue, "forwardsChiSquaredPerHit", "#chi^{2}_{F}/N", "Number of Entries", "True Forwards Muons: #chi^{2}_{F}/N");
    Draw(forwardsDeltaChiSquared, kBlue, "forwardsDeltaChiSquared", "#Delta#chi^{2}", "Number of Entries", "True Forwards Muons: #Delta#chi^{2}");
    Draw(forwardsDeltaChiSquaredPerHit, kBlue, "forwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N", "Number of Entries", "True Forwards Muons: #Delta#chi^{2}/N");

    Draw(backwardsChiSquaredPerHit, kRed, "backwardsChiSquaredPerHit", "#chi^{2}_{B}/N", "Number of Entries", "True Backwards Muons: #chi^{2}_{B}/N");
    Draw(backwardsDeltaChiSquared, kRed, "backwardsDeltaChiSquared", "#Delta#chi^{2}", "Number of Entries", "True Backwards Muons: #Delta#chi^{2}");
    Draw(backwardsDeltaChiSquaredPerHit, kRed, "backwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N", "Number of Entries", "True Backwards Muons: #Delta#chi^{2}/N");

    Draw(upwardsChiSquaredPerHit, kBlue, "upwardsChiSquaredPerHit", "#chi^{2}_{U}/N", "Number of Entries", "True Upwards Muons: #chi^{2}_{U}/N");
    Draw(upwardsDeltaChiSquared, kBlue, "upwardsDeltaChiSquared", "#Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}", "Number of Entries", "True Upwards Muons: #Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}");
    Draw(upwardsDeltaChiSquaredPerHit, kBlue, "upwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N", "Number of Entries", "True Upwards Muons: #Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N");

    Draw(downwardsChiSquaredPerHit, kRed, "downwardsChiSquaredPerHit", "#chi^{2}_{D}/N", "Number of Entries", "True Downwards Muons: #chi^{2}_{D}/N");
    Draw(downwardsDeltaChiSquared, kRed, "downwardsDeltaChiSquared", "#Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}", "Number of Entries", "True Downwards Muons: #Delta#chi^{2} #equiv #chi^{2}_{D} - #chi^{2}_{U}");
    Draw(downwardsDeltaChiSquaredPerHit, kRed, "downwardsDeltaChiSquaredPerHit", "#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N", "Number of Entries", "True Downwards Muons: #Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N");

    //--------------------------------------------------------------------------------------------------------------------------------------
}
