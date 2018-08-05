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
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/mc_data_comparison/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH1F* pHistogram1, TH1F* pHistogram2, EColor colour1, EColor colour2, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    pHistogram1->Scale(1/(pHistogram1->GetEntries()));
    pHistogram2->Scale(1/(pHistogram2->GetEntries()));

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
    pHistogram1->GetYaxis()->SetRangeUser(0.0, 1.25 * largestBinEntry);
    pHistogram1->SetLineColor(colour1);
    pHistogram1->Draw();
    pHistogram2->SetLineColor(colour2);

    auto legend = new TLegend(0.15,0.75,0.35,0.85);
    legend->SetHeader("Event Type"); 
    legend->AddEntry(pHistogram1, "MC Events","l");
    legend->AddEntry(pHistogram2, "Data Events","l");
    legend->Draw("same");

    pHistogram2->Draw("same");
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/mc_data_comparison/" + histogramName + ".png");
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
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/mc_data_comparison/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cosmic_MC_Data_Comparison(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/saved_results/mc_cosmic.root");
    TFile *f2 = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/saved_results/data.root");
    
    TTree *t1 = (TTree*)f1->Get("PFO");
    TTree *t2 = (TTree*)f2->Get("PFO");

    int nentries1 = (Int_t)t1->GetEntries();
    int nentries2 = (Int_t)t2->GetEntries();
    
    TH1F *pMCNumberHitsHistogram = new TH1F("pMCNumberHitsHistogram","", 100, 0, 500);
    TH1F *pDataNumberHitsHistogram = new TH1F("pDataNumberHitsHistogram","", 100, 0, 500);

    TH1F *pMCRecoLengthHistogram = new TH1F("pMCRecoLengthHistogram","", 100, 0, 500);
    TH1F *pDataRecoLengthHistogram = new TH1F("pDataRecoLengthHistogram","", 100, 0, 500);

    TH1F *pMCChiSquaredMinHistogram = new TH1F("pMCChiSquaredMinHistogram","", 100, 0, 5.0);
    TH1F *pDataChiSquaredMinHistogram = new TH1F("pDataChiSquaredMinHistogram","", 100, 0, 5.0);

    TH1F *pMCBeginXHistogram = new TH1F("pMCBeginXHistogram","", 100, 0, 320);
    TH1F *pDataBeginXHistogram = new TH1F("pDataBeginXHistogram","", 100, 0, 320);

    TH1F *pMCBeginYHistogram = new TH1F("pMCBeginYHistogram","", 100, -256.35/2, 256.35/2);
    TH1F *pDataBeginYHistogram = new TH1F("pDataBeginYHistogram","", 100, -256.35/2, 256.35/2);

    TH1F *pMCBeginZHistogram = new TH1F("pMCBeginZHistogram","", 50, 0, 1037);
    TH1F *pDataBeginZHistogram = new TH1F("pDataBeginZHistogram","", 50, 0, 1037);

    TH1F *pMCRecoPhiHistogram = new TH1F("pMCRecoPhiHistogram","", 100, 0.0, 3.2);
    TH1F *pDataRecoPhiHistogram = new TH1F("pDataRecoPhiHistogram","", 100, 0.0, 3.2);

    TH1F *pMCRecoThetaHistogram = new TH1F("pMCRecoThetaHistogram","", 100, 0.0, 3.2);
    TH1F *pDataRecoThetaHistogram = new TH1F("pDataRecoThetaHistogram","", 100, 0.0, 3.2);
    
    int   MCDownwards;
    int   NeutrinoInduced;
    float recoLength;
    float RecoPhi;
    float RecoTheta;
    float BeginX;
    float BeginY;
    float BeginZ;
    float EndX;
    float EndY;
    float EndZ;
    float MinChiSquaredPerHit;
    float UpDownDeltaChiSquaredPerHit;
    int   NumberHits;
    int   FileIdentifier;
    int   EventNumber;
    
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::cout << "There are " << nentries1 << " entries in tree 1" << std::endl;
    std::cout << "There are " << nentries2 << " entries in tree 2" << std::endl;

    t1->SetBranchAddress("MCDownwards", &MCDownwards);
    t1->SetBranchAddress("NeutrinoInduced", &NeutrinoInduced);
    t1->SetBranchAddress("recoLength", &recoLength);
    t1->SetBranchAddress("RecoPhi", &RecoPhi);
    t1->SetBranchAddress("RecoTheta", &RecoTheta);
    t1->SetBranchAddress("BeginX", &BeginX);
    t1->SetBranchAddress("BeginY", &BeginY);
    t1->SetBranchAddress("BeginZ", &BeginZ);
    t1->SetBranchAddress("EndX", &EndX);
    t1->SetBranchAddress("EndY", &EndY);
    t1->SetBranchAddress("EndZ", &EndZ);
    t1->SetBranchAddress("MinChiSquaredPerHit", &MinChiSquaredPerHit);
    t1->SetBranchAddress("UpDownDeltaChiSquaredPerHit", &UpDownDeltaChiSquaredPerHit);
    t1->SetBranchAddress("NumberHits", &NumberHits);
    t1->SetBranchAddress("FileIdentifier", &FileIdentifier);
    t1->SetBranchAddress("EventNumber", &EventNumber);
    
    for (int i = 0; i < nentries1; i++)
    {
        t1->GetEntry(i);
        pMCNumberHitsHistogram->Fill(NumberHits);
        pMCRecoLengthHistogram->Fill(recoLength);
        pMCChiSquaredMinHistogram->Fill(MinChiSquaredPerHit);
        if (BeginY < EndY)
        {
            pMCBeginXHistogram->Fill(BeginX);
            pMCBeginYHistogram->Fill(BeginY);
            pMCBeginZHistogram->Fill(BeginZ);
        }
        else
        {
            pMCBeginXHistogram->Fill(EndX);
            pMCBeginYHistogram->Fill(EndY);
            pMCBeginZHistogram->Fill(EndZ);
        }
    }

    t2->SetBranchAddress("MCDownwards", &MCDownwards);
    t2->SetBranchAddress("NeutrinoInduced", &NeutrinoInduced);
    t2->SetBranchAddress("recoLength", &recoLength);
    t2->SetBranchAddress("RecoPhi", &RecoPhi);
    t2->SetBranchAddress("RecoTheta", &RecoTheta);
    t2->SetBranchAddress("BeginX", &BeginX);
    t2->SetBranchAddress("BeginY", &BeginY);
    t2->SetBranchAddress("BeginZ", &BeginZ);
    t2->SetBranchAddress("EndX", &EndX);
    t2->SetBranchAddress("EndY", &EndY);
    t2->SetBranchAddress("EndZ", &EndZ);
    t2->SetBranchAddress("MinChiSquaredPerHit", &MinChiSquaredPerHit);
    t2->SetBranchAddress("UpDownDeltaChiSquaredPerHit", &UpDownDeltaChiSquaredPerHit);
    t2->SetBranchAddress("NumberHits", &NumberHits);
    t2->SetBranchAddress("FileIdentifier", &FileIdentifier);
    t2->SetBranchAddress("EventNumber", &EventNumber);

    for (int i = 0; i < nentries2; i++)
    {
        t2->GetEntry(i);
        pDataNumberHitsHistogram->Fill(NumberHits);
        pDataRecoLengthHistogram->Fill(recoLength);
        pDataChiSquaredMinHistogram->Fill(MinChiSquaredPerHit);

        if (BeginY < EndY)
        {
            pDataBeginXHistogram->Fill(BeginX);
            pDataBeginYHistogram->Fill(BeginY);
            pDataBeginZHistogram->Fill(BeginZ);
        }
        else
        {
            pDataBeginXHistogram->Fill(EndX);
            pDataBeginYHistogram->Fill(EndY);
            pDataBeginZHistogram->Fill(EndZ);
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------------------------------

    Draw(pMCNumberHitsHistogram, pDataNumberHitsHistogram, kBlue, kRed, "NumberHits", "Number of PFO Hits N", "Number of Entries", "Number of PFO Hits N for MC and Data Cosmics");
    Draw(pMCRecoLengthHistogram, pDataRecoLengthHistogram, kBlue, kRed, "RecoLength", "Reconstructed PFO Length L", "Number of Entries", "Reconstructed PFO Length L for MC and Data Cosmics");
    Draw(pMCChiSquaredMinHistogram, pDataChiSquaredMinHistogram, kBlue, kRed, "MinChiSquared", "Best Fit #chi^{2}_{min}/N", "Number of Entries", "Best Fit #chi^{2}_{min}/N for MC and Data Cosmics");
    Draw(pMCBeginXHistogram, pDataBeginXHistogram, kBlue, kRed, "BeginX", "X Position of Reconstructed Track Low Y Endpoint (cm)", "Number of Entries", "X Position of Reconstructed Track Low Y Endpoint for MC and Data Cosmics");
    Draw(pMCBeginYHistogram, pDataBeginYHistogram, kBlue, kRed, "BeginY", "Y Position of Reconstructed Track Low Y Endpoint (cm)", "Number of Entries", "Y Position of Reconstructed Track Low Y Endpoint for MC and Data Cosmics");
    Draw(pMCBeginZHistogram, pDataBeginZHistogram, kBlue, kRed, "BeginZ", "Z Position of Reconstructed Track Low Y Endpoint (cm)", "Number of Entries", "Z Position of Reconstructed Track Low Y Endpoint for MC and Data Cosmics");

    //--------------------------------------------------------------------------------------------------------------------------------------
}
