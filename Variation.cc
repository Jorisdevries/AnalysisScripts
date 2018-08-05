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
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/variation/" + histogramName + ".png");
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
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/variation/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Variation(void)
{
    TFile *f = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/catroots/0.root");
    
    TTree *t1 = (TTree*)f->Get("PFO");
    int nentries = (Int_t)t1->GetEntries();
    
    TH2F *forwardsMinChiSquaredPerHitWithN = new TH2F("forwardsMinChiSquaredPerHitWithN","", 50, -1.0, 10.0, 50, 0, 1100);
    TH2F *backwardsMinChiSquaredPerHitWithN = new TH2F("backwardsMinChiSquaredPerHitWithN","", 50, -1.0, 10.0, 50, 0, 1100);
    TH2F *upwardsMinChiSquaredPerHitWithN = new TH2F("upwardsMinChiSquaredPerHitWithN","", 50, -1.0, 10.0, 50, 0, 1100);
    TH2F *downwardsMinChiSquaredPerHitWithN = new TH2F("downwardsMinChiSquaredPerHitWithN","", 50, -1.0, 10.0, 50, 0, 1100);

    TH2F *forwardsMinChiSquaredPerHitWithLength = new TH2F("forwardsMinChiSquaredPerHitWithLength","", 50, -1.0, 4.0, 50, 0, 600);
    TH2F *backwardsMinChiSquaredPerHitWithLength = new TH2F("backwardsMinChiSquaredPerHitWithLength","", 50, -1.0, 10.0, 50, 0, 1100);
    TH2F *upwardsMinChiSquaredPerHitWithLength = new TH2F("upwardsMinChiSquaredPerHitWithLength","", 50, -1.0, 10.0, 50, 0, 1100);
    TH2F *downwardsMinChiSquaredPerHitWithLength = new TH2F("downwardsMinChiSquaredPerHitWithLength","", 50, -1.0, 10.0, 50, 0, 1100);

    TH2F *forwardsMinChiSquaredPerHitWithPhi = new TH2F("forwardsMinChiSquaredPerHitWithPhi","", 50, -1.0, 4.0, 50, 0, 3.14);
    TH2F *forwardsMinChiSquaredPerHitWithTheta = new TH2F("forwardsMinChiSquaredPerHitWithTheta","", 50, -1.0, 4.0, 50, 0, 3.14);
    
    int   MCDirection;
    float recoLength;
    float MCLength;
    int   MCForwards;
    int   MCDownwards;
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
    float DownwardsChiSquared;
    float UpwardsChiSquared;
    float DownwardsChiSquaredPerHit;
    float UpwardsChiSquaredPerHit;
    int   NumberHits;
    int   MCIsClusterTwoParticles;
    int   FileIdentifier;
    int   EventNumber;
    
    t1->SetBranchAddress("MCDirection", &MCDirection);
    t1->SetBranchAddress("recoLength", &recoLength);
    t1->SetBranchAddress("MCLength", &MCLength);
    t1->SetBranchAddress("MCForwards", &MCForwards);
    t1->SetBranchAddress("MCDownwards", &MCDownwards);
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
    t1->SetBranchAddress("DownwardsChiSquared", &DownwardsChiSquared);
    t1->SetBranchAddress("UpwardsChiSquared", &UpwardsChiSquared);
    t1->SetBranchAddress("DownwardsChiSquaredPerHit", &DownwardsChiSquaredPerHit);
    t1->SetBranchAddress("UpwardsChiSquaredPerHit", &UpwardsChiSquaredPerHit);
    t1->SetBranchAddress("NumberHits", &NumberHits);
    t1->SetBranchAddress("MCIsClusterTwoParticles", &MCIsClusterTwoParticles);
    t1->SetBranchAddress("FileIdentifier", &FileIdentifier);
    t1->SetBranchAddress("EventNumber", &EventNumber);
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::cout << "There are " << nentries << " entries in the tree" << std::endl;
    
    for (int i = 1; i < nentries; i++)
    {
        t1->GetEntry(i);

        if (MCForwards == 1)
        {
            forwardsMinChiSquaredPerHitWithN->Fill(MinChiSquaredPerHit, NumberHits);
            forwardsMinChiSquaredPerHitWithLength->Fill(MinChiSquaredPerHit, recoLength);
            forwardsMinChiSquaredPerHitWithPhi->Fill(MinChiSquaredPerHit, MCPhi);
            forwardsMinChiSquaredPerHitWithTheta->Fill(MinChiSquaredPerHit, MCTheta);
        }
        else
        {
            backwardsMinChiSquaredPerHitWithN->Fill(MinChiSquaredPerHit, NumberHits);
        }

        if (MCDownwards == 1)
        {
            upwardsMinChiSquaredPerHitWithN->Fill(MinChiSquaredPerHit, NumberHits);
        }
        else
        {
            downwardsMinChiSquaredPerHitWithN->Fill(MinChiSquaredPerHit, NumberHits);
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------------------------------

    Draw(forwardsMinChiSquaredPerHitWithN, "forwardsMinChiSquaredPerHitWithN", "#chi^{2}_{Min.}/N", "Number of Track Hits N", "True Forwards Muons: #chi^{2}_{Min.}/N with N");
    Draw(backwardsMinChiSquaredPerHitWithN, "backwardsMinChiSquaredPerHitWithN", "#chi^{2}_{Min.}/N", "Number of Track Hits N", "True Backwards Muons: #chi^{2}_{Min.}/N with N");
    Draw(upwardsMinChiSquaredPerHitWithN, "upwardsMinChiSquaredPerHitWithN", "#chi^{2}_{Min.}/N", "Number of Track Hits N", "True Downwards Muons: #chi^{2}_{Min.}/N with N");
    Draw(downwardsMinChiSquaredPerHitWithN, "downwardsMinChiSquaredPerHitWithN", "#chi^{2}_{Min.}/N", "Number of Track Hits N", "True Downwards Muons: #chi^{2}_{Min.}/N with N");

    Draw(forwardsMinChiSquaredPerHitWithLength, "forwardsMinChiSquaredPerHitWithLength", "#chi^{2}_{Min.}/N", "Track Length L (cm)", "True Forwards Muons: #chi^{2}_{Min.}/N with Track Length");
    Draw(forwardsMinChiSquaredPerHitWithPhi, "forwardsMinChiSquaredPerHitWithPhi", "#chi^{2}_{Min.}/N", "Track Phi L (cm)", "True Forwards Muons: #chi^{2}_{Min.}/N with Track Phi");
    Draw(forwardsMinChiSquaredPerHitWithTheta, "forwardsMinChiSquaredPerHitWithTheta", "#chi^{2}_{Min.}/N", "Track Theta L (cm)", "True Forwards Muons: #chi^{2}_{Min.}/N with Track Theta");

    //--------------------------------------------------------------------------------------------------------------------------------------
}
