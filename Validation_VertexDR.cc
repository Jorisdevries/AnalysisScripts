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

void Validation_VertexDR(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_direction_results/ccqel_newest.root");
    TTree *t1 = (TTree*)f1->Get("Validation");

    TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_direction_results/1.root");
    TTree *t2 = (TTree*)f2->Get("Validation");

    int nentries1 = (Int_t)t1->GetEntries();
    int nentries2 = (Int_t)t2->GetEntries();
    
    //TH1F *pVertexDRHistogram = new TH1F("pVertexDRHistogram","", 100, 0, 100);

    float targetVertexX1, targetVertexY1, targetVertexZ1;
    float recoVertexX1, recoVertexY1, recoVertexZ1;
    int   FileIdentifier1;
    int   EventNumber1;

    float targetVertexX2, targetVertexY2, targetVertexZ2;
    float recoVertexX2, recoVertexY2, recoVertexZ2;
    int   FileIdentifier2;
    int   EventNumber2;
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::cout << "There are " << nentries1 << " entries in tree 1" << std::endl;

    t1->SetBranchAddress("targetVertexX", &targetVertexX1);
    t1->SetBranchAddress("targetVertexY", &targetVertexY1);
    t1->SetBranchAddress("targetVertexZ", &targetVertexZ1);
    t1->SetBranchAddress("recoVertexX", &recoVertexX1);
    t1->SetBranchAddress("recoVertexY", &recoVertexY1);
    t1->SetBranchAddress("recoVertexZ", &recoVertexZ1);
    t1->SetBranchAddress("fileIdentifier", &FileIdentifier1);
    t1->SetBranchAddress("eventNumber", &EventNumber1);

    t2->SetBranchAddress("targetVertexX", &targetVertexX2);
    t2->SetBranchAddress("targetVertexY", &targetVertexY2);
    t2->SetBranchAddress("targetVertexZ", &targetVertexZ2);
    t2->SetBranchAddress("recoVertexX", &recoVertexX2);
    t2->SetBranchAddress("recoVertexY", &recoVertexY2);
    t2->SetBranchAddress("recoVertexZ", &recoVertexZ2);
    t2->SetBranchAddress("fileIdentifier", &FileIdentifier2);
    t2->SetBranchAddress("eventNumber", &EventNumber2);
    
    //Find events where direction doesn't save the event
    /*
    for (int i = 0; i < nentries1; i++)
    {
        t1->GetEntry(i);

        if (i > 10000)
            break;

        float vertexDR1(std::sqrt(((targetVertexX1 - recoVertexX1) * (targetVertexX1 - recoVertexX1)) + ((targetVertexY1 - recoVertexY1) * (targetVertexY1 - recoVertexY1)) + ((targetVertexZ1 - recoVertexZ1) * (targetVertexZ1 - recoVertexZ1))));

        for (int j = 0; j < nentries2; j++)
        {
            t2->GetEntry(j);

            if (j > 10000)
                break;

            if (FileIdentifier1 != FileIdentifier2 || EventNumber1 != EventNumber2)
                continue;

            float vertexDR2(std::sqrt(((targetVertexX2 - recoVertexX2) * (targetVertexX2 - recoVertexX2)) + ((targetVertexY2 - recoVertexY2) * (targetVertexY2 - recoVertexY2)) + ((targetVertexZ2 - recoVertexZ2) * (targetVertexZ2 - recoVertexZ2))));

            if (vertexDR2 < 5.0 && vertexDR1 > 20.0 && EventNumber1 <= 5)
                std::cout << "DR difference of " << vertexDR1 - vertexDR2 << " in event " << FileIdentifier1 << ":" << EventNumber1 << std::endl;

        }
    }
    */

    //Find events that are made worse by direction
    for (int i = 0; i < nentries1; i++)
    {
        t1->GetEntry(i);

        if (i > 10000)
            break;

        float vertexDR1(std::sqrt(((targetVertexX1 - recoVertexX1) * (targetVertexX1 - recoVertexX1)) + ((targetVertexY1 - recoVertexY1) * (targetVertexY1 - recoVertexY1)) + ((targetVertexZ1 - recoVertexZ1) * (targetVertexZ1 - recoVertexZ1))));

        for (int j = 0; j < nentries2; j++)
        {
            t2->GetEntry(j);

            if (j > 10000)
                break;

            if (FileIdentifier1 != FileIdentifier2 || EventNumber1 != EventNumber2)
                continue;

            float vertexDR2(std::sqrt(((targetVertexX2 - recoVertexX2) * (targetVertexX2 - recoVertexX2)) + ((targetVertexY2 - recoVertexY2) * (targetVertexY2 - recoVertexY2)) + ((targetVertexZ2 - recoVertexZ2) * (targetVertexZ2 - recoVertexZ2))));

            /*
            std::cout << "recoVertexX1: " << recoVertexX1 << std::endl;
            std::cout << "recoVertexY1: " << recoVertexY1 << std::endl;
            std::cout << "recoVertexZ1: " << recoVertexZ1 << std::endl;
            std::cout << "targetVertexX1: " << targetVertexX1 << std::endl;
            std::cout << "targetVertexY1: " << targetVertexY1 << std::endl;
            std::cout << "targetVertexZ1: " << targetVertexZ1 << std::endl;
            std::cout << "-----------------" << std::endl;
            */

            if (vertexDR1 > vertexDR2 + 50)
                std::cout << "DR difference of (dir - nodir) " << vertexDR1 - vertexDR2 << " (" << vertexDR1 << " - " << vertexDR2 << ") in event " << FileIdentifier1 << ":" << EventNumber1 << std::endl;

        }
    }

    //Find events that are present in the no-direction file but not in the direction file 
    /*
    ofstream outfile;
    outfile.open("events_dir.txt");
    int currentFileIdentidier(-1);

    for (int i = 0; i < nentries1; i++)
    {
        t1->GetEntry(i);
    
            outfile << FileIdentifier1 << ":" << EventNumber1 << std::endl;
            currentFileIdentidier = FileIdentifier1;
    }


    ofstream outfile2;
    outfile2.open("events_nodir.txt");
    int currentFileIdentidier2(-1);
    for (int j = 0; j < nentries2; j++)
    {
        t2->GetEntry(j);

            outfile2 << FileIdentifier2 << ":" << EventNumber2 << std::endl;
            currentFileIdentidier2 = FileIdentifier2;
    }
    */

    //--------------------------------------------------------------------------------------------------------------------------------------

    //Draw(pVertexDRHistogram, pSCCVertexDRHistogram, kRed, kBlue, "SpaceChargeCorrection", "Vertex #DeltaR", "Number of Entries", "Vertex #DeltaR (1 #mu CCQEL) Before And After Space Charge Correction");

    //--------------------------------------------------------------------------------------------------------------------------------------
}
