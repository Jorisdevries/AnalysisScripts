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
#include "TLegend.h"

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
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/vertex/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH1F* pHistogram1, TH1F* pHistogram2, EColor colour1, EColor colour2, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    //pHistogram1->Scale(1/(pHistogram1->GetEntries()));
    //pHistogram2->Scale(1/(pHistogram2->GetEntries()));

    float largestBinEntry(0.f);
 
    for (int i = 1; i < pHistogram1->GetNbinsX(); i++) 
    {    
        if (pHistogram1->GetBinContent(i) > largestBinEntry) 
            largestBinEntry = pHistogram1->GetBinContent(i);

        if (pHistogram2->GetBinContent(i) > largestBinEntry) 
            largestBinEntry = pHistogram2->GetBinContent(i);
    }  

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram1->SetXTitle(xTitle.c_str());
    pHistogram1->SetYTitle(yTitle.c_str());
    pHistogram1->SetTitle(plotTitle.c_str());
    pHistogram1->GetYaxis()->SetTitleOffset(1.3);
    pHistogram1->GetYaxis()->SetRangeUser(0.0, 1.05 * largestBinEntry);
    pHistogram1->SetLineColor(colour1);
    pHistogram1->Draw();
    pHistogram2->SetLineColor(colour2);

    auto legend = new TLegend(0.55,0.75,0.75,0.85);
    legend->SetHeader("Event Type"); 
    legend->AddEntry(pHistogram1, "Unmodified #DeltaR","l");
    legend->AddEntry(pHistogram2, "Modified #DeltaR","l");
    legend->Draw("same");

    pHistogram2->Draw("same");
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/vertex/" + histogramName + ".png");
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
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/vertex/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Splitting_Vertex_Creation(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/saved_results/spacecharge.root");
    TTree *t1 = (TTree*)f1->Get("Splitting");

    int nentries1 = (Int_t)t1->GetEntries();
    
    TH1F *pVertexDRHistogram = new TH1F("pVertexDRHistogram","", 100, -1, 10);
    TH1F *pSCCVertexDRHistogram = new TH1F("pSCCVertexDRHistogram","", 100, -1, 10);

    TH1F *pNewVertexDRHistogram = new TH1F("pNewVertexDRHistogram","", 82, -1, 40);
    TH1F *pSCCNewVertexDRHistogram = new TH1F("pSCCNewertexDRHistogram","", 82, -1, 40);

    TH2F *pOldAndNewDRHistogram = new TH2F("pOldAndNewDRHistogram", "", 100, 0, 100, 100, 0, 100);
    TH2F *pSCCOldAndNewDRHistogram = new TH2F("pSCCOldAndNewDRHistogram", "", 100, 0, 100, 100, 0, 100);

    TH2F *pSCCEffect = new TH2F("pSCCEffect", "", 100, 0, 100, 100, 0, 100);

    TH2F *pPositionShiftX = new TH2F("pPositionShiftX", "", 100, 0, 256, 100, -10, 10);
    TH2F *pPositionShiftY = new TH2F("pPositionShiftY", "", 100, -160, 160, 100, -10, 10);
    TH2F *pPositionShiftZ = new TH2F("pPositionShiftY", "", 100, 0, 1040, 100, -10, 10);

    int   MCIsClusterTwoParticles;
    float VertexDR;
    float SCCVertexDR;
    float AfterSplitVertexDR;
    float SCCAfterSplitVertexDR;
    float VertexX;
    float VertexY;
    float VertexZ;
    float PositionOffsetX;
    float PositionOffsetY;
    float PositionOffsetZ;
    float SCCVertexDRChangeX;
    float SCCVertexDRChangeY;
    float SCCVertexDRChangeZ;
    int   FileIdentifier;
    int   EventNumber;
    
    t1->SetBranchAddress("MCIsClusterTwoParticles", &MCIsClusterTwoParticles);
    t1->SetBranchAddress("VertexDR", &VertexDR);
    t1->SetBranchAddress("SCCVertexDR", &SCCVertexDR);
    t1->SetBranchAddress("AfterSplitVertexDR", &AfterSplitVertexDR);
    t1->SetBranchAddress("SCCAfterSplitVertexDR", &SCCAfterSplitVertexDR);
    t1->SetBranchAddress("VertexX", &VertexX);
    t1->SetBranchAddress("VertexY", &VertexY);
    t1->SetBranchAddress("VertexZ", &VertexZ);
    t1->SetBranchAddress("PositionOffsetX", &PositionOffsetX);
    t1->SetBranchAddress("PositionOffsetY", &PositionOffsetY);
    t1->SetBranchAddress("PositionOffsetZ", &PositionOffsetZ);
    t1->SetBranchAddress("SCCVertexDRChangeX", &SCCVertexDRChangeX);
    t1->SetBranchAddress("SCCVertexDRChangeY", &SCCVertexDRChangeY);
    t1->SetBranchAddress("SCCVertexDRChangeZ", &SCCVertexDRChangeZ);
    t1->SetBranchAddress("FileIdentifier", &FileIdentifier);
    t1->SetBranchAddress("EventNumber", &EventNumber);
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::cout << "There are " << nentries1 << " entries in tree 1" << std::endl;

    for (int i = 0; i < nentries1; i++)
    {
        t1->GetEntry(i);

        pVertexDRHistogram->Fill(VertexDR);
        pSCCVertexDRHistogram->Fill(SCCVertexDR);

        pNewVertexDRHistogram->Fill(AfterSplitVertexDR);
        pSCCNewVertexDRHistogram->Fill(SCCAfterSplitVertexDR);

        pOldAndNewDRHistogram->Fill(VertexDR, AfterSplitVertexDR);
        pSCCOldAndNewDRHistogram->Fill(SCCVertexDR, SCCAfterSplitVertexDR);

        pSCCEffect->Fill(VertexDR, SCCVertexDR);

        pPositionShiftX->Fill(VertexX, PositionOffsetX);
        pPositionShiftY->Fill(VertexY, PositionOffsetY);
        pPositionShiftZ->Fill(VertexZ, PositionOffsetZ);

        //std::cout << VertexDR << "->" << SCCVertexDR << std::endl;
        //std::cout << "--------------------" << std::endl;
    
        //if (VertexDR < 10.0 && AfterSplitVertexDR > 50.0 && MCIsClusterTwoParticles == 1)
        //    std::cout << FileIdentifier << ":" << EventNumber << std::endl;
        
        if (SCCAfterSplitVertexDR - SCCVertexDR > 50.0 && MCIsClusterTwoParticles == 1)
            std::cout << FileIdentifier << ":" << EventNumber << std::endl;
    }
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    Draw(pPositionShiftX, "PositionShiftX", "Before Space Charge Correction Vertex X Position", "Space Charge Correction X Component", "X: Space Charge Correction");
    Draw(pPositionShiftY, "PositionShiftY", "Before Space Charge Correction Vertex Y Position", "Space Charge Correction Y Component", "Y: Space Charge Correction");
    Draw(pPositionShiftZ, "PositionShiftZ", "Before Space Charge Correction Vertex Z Position", "Space Charge Correction Z Component", "Z: Space Charge Correction");
    Draw(pVertexDRHistogram, pSCCVertexDRHistogram, kRed, kBlue, "SpaceChargeCorrection", "Vertex #DeltaR", "Number of Entries", "Vertex #DeltaR (1 #mu 1p CCQEL) Before And After Space Charge Correction");

    Draw(pOldAndNewDRHistogram, "DRChange", "Before Splitting Vertex #DeltaR", "After Splitting Vertex #DeltaR", "Change in Vertex #DeltaR From Splitting");
    Draw(pSCCOldAndNewDRHistogram, "SCCDRChange", "Before Splitting Vertex #DeltaR", "After Splitting Vertex #DeltaR", "Change in Vertex #DeltaR From Splitting (With Space Charge Correction)");

    /*
    Draw(pSCCEffect, "SCCEffect", "Before Space Charge Correction Vertex #DeltaR", "After Space Charge Correction Vertex #DeltaR", "Change in Vertex #DeltaR Due To Space Charge Correction");
    Draw(pVertexDRHistogram, pNewVertexDRHistogram, kRed, kBlue, "VertexComparison", "Vertex #DeltaR", "Number of Entries", "Vertex #DeltaR (1 #mu 1p CCQEL) Before And After Splitting Vertex Creation");
    Draw(pSCCVertexDRHistogram, pSCCNewVertexDRHistogram, kRed, kBlue, "SCCVertexComparison", "Vertex #DeltaR", "Number of Entries", "Vertex #DeltaR (1 #mu 1p CCQEL) Before And After Splitting Vertex Creation (With Space Charge Correction)");
    Draw(pOldAndNewDRHistogram, "DRChange", "Before Splitting Vertex #DeltaR", "After Splitting Vertex #DeltaR", "Change in Vertex #DeltaR From Splitting");
    Draw(pSCCOldAndNewDRHistogram, "SCCDRChange", "Before Splitting Vertex #DeltaR", "After Splitting Vertex #DeltaR", "Change in Vertex #DeltaR From Splitting (With Space Charge Correction)");
    */

    //--------------------------------------------------------------------------------------------------------------------------------------
}
