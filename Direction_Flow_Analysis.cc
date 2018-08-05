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
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/directionflow/" + histogramName + ".png");
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
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/directionflow/" + histogramName + ".png");
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
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/vertex/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FindBadEvents(TTree* noDirectionTree, TTree* directionTree)
{
    float noDirectionSCCVertexDR, directionSCCVertexDR;
    noDirectionTree->SetBranchAddress("SCCVertexDR", &noDirectionSCCVertexDR);
    directionTree->SetBranchAddress("SCCVertexDR", &directionSCCVertexDR);

    int noDirectionFileIdentifier, noDirectionEventNumber;
    noDirectionTree->SetBranchAddress("FileIdentifier", &noDirectionFileIdentifier);
    noDirectionTree->SetBranchAddress("EventNumber", &noDirectionEventNumber);

    int directionFileIdentifier, directionEventNumber;
    directionTree->SetBranchAddress("FileIdentifier", &directionFileIdentifier);
    directionTree->SetBranchAddress("EventNumber", &directionEventNumber);

    int noDirectionNMuons, noDirectionNProtons, noDirectionNOthers;
    noDirectionTree->SetBranchAddress("nMuons", &noDirectionNMuons);
    noDirectionTree->SetBranchAddress("nProtons", &noDirectionNProtons);
    noDirectionTree->SetBranchAddress("nOthers", &noDirectionNOthers);

    for (int i = 0; i < noDirectionTree->GetEntries(); ++i)
    {
        noDirectionTree->GetEntry(i);

        for (int j = 0; j < directionTree->GetEntries(); ++j)
        {
            directionTree->GetEntry(j);

            if (noDirectionFileIdentifier != directionFileIdentifier || noDirectionEventNumber != directionEventNumber)
                continue;
        
            //if (directionSCCVertexDR > noDirectionSCCVertexDR)
            if (noDirectionEventNumber <= 5 && noDirectionNMuons == 1 && noDirectionNProtons == 0 && noDirectionNOthers == 0 && directionSCCVertexDR > 50 && noDirectionSCCVertexDR > 50 && directionSCCVertexDR == noDirectionSCCVertexDR)
                std::cout << "Bad event: " << noDirectionFileIdentifier << ":" << noDirectionEventNumber << std::endl; 
        }
    }
       
}
//------------------------------------------------------------------------------------------------------------------------------------------

void CompareEventType(TTree* noDirectionTree, TTree* directionTree, int nMuons, int nProtons, int nOthers, bool enableSCE)
{
    int upperBound(10);
    TH1F *pVertexDRHistogram = new TH1F("pVertexDRHistogram","", 100, -1, upperBound);
    TH1F *pDirectionVertexDRHistogram = new TH1F("pDirectionVertexDRHistogram","", 100, -1, upperBound);

    float noDirectionVertexDR, noDirectionSCCVertexDR;
    int noDirectionNMuons, noDirectionNProtons, noDirectionNOthers;
    float directionVertexDR, directionSCCVertexDR;
    int directionNMuons, directionNProtons, directionNOthers;

    noDirectionTree->SetBranchAddress("VertexDR", &noDirectionVertexDR);
    noDirectionTree->SetBranchAddress("SCCVertexDR", &noDirectionSCCVertexDR);
    noDirectionTree->SetBranchAddress("nMuons", &noDirectionNMuons);
    noDirectionTree->SetBranchAddress("nProtons", &noDirectionNProtons);
    noDirectionTree->SetBranchAddress("nOthers", &noDirectionNOthers);

    directionTree->SetBranchAddress("VertexDR", &directionVertexDR);
    directionTree->SetBranchAddress("SCCVertexDR", &directionSCCVertexDR);
    directionTree->SetBranchAddress("nMuons", &directionNMuons);
    directionTree->SetBranchAddress("nProtons", &directionNProtons);
    directionTree->SetBranchAddress("nOthers", &directionNOthers);

    int unmodifiedOverflow(0), directionFlowOverflow(0);

    for (int i = 0; i < noDirectionTree->GetEntries(); ++i)
    {
        noDirectionTree->GetEntry(i);

        if (!(noDirectionNMuons == nMuons && noDirectionNProtons == nProtons && noDirectionNOthers == nOthers))
            continue;

        if ((!enableSCE && noDirectionVertexDR > upperBound) || (enableSCE && noDirectionSCCVertexDR > upperBound))
            ++unmodifiedOverflow;
        
        if (!enableSCE)
            pVertexDRHistogram->Fill(noDirectionVertexDR);
        else
            pVertexDRHistogram->Fill(noDirectionSCCVertexDR);
    }

    for (int i = 0; i < directionTree->GetEntries(); ++i)
    {
        directionTree->GetEntry(i);

        if (!(directionNMuons == nMuons && directionNProtons == nProtons && directionNOthers == nOthers))
            continue;

        if ((!enableSCE && directionVertexDR > upperBound) || (enableSCE && directionSCCVertexDR > upperBound))
            ++directionFlowOverflow;
        
        if (!enableSCE)
            pDirectionVertexDRHistogram->Fill(directionVertexDR);
        else
            pDirectionVertexDRHistogram->Fill(directionSCCVertexDR);
    }

    std::cout << noDirectionTree->GetEntries() << std::endl;
    std::cout << directionTree->GetEntries() << std::endl;

    std::cout << "Unmodified overflow: " << unmodifiedOverflow << " (" << 100 * (float)unmodifiedOverflow/noDirectionTree->GetEntries() << "%)" << std::endl;
    std::cout << "Overflow with direction flow enabled: " << directionFlowOverflow << " (" << 100 * (float)directionFlowOverflow/directionTree->GetEntries() << "%)" << std::endl;

    std::string stringNMuons(std::to_string(nMuons)), stringNProtons(std::to_string(nProtons)), stringSCE(enableSCE ? "With" : "Without");

    Draw(pVertexDRHistogram, pDirectionVertexDRHistogram, kRed, kBlue, "DirectionFlowComparison_"+stringNMuons+"_"+stringNProtons+"_"+stringSCE+"SCE", "Vertex #DeltaR", "Number of Entries", "Vertex #DeltaR (" + stringNMuons + "#mu " + stringNProtons + "p CCQEL) Direction Flow Comparison " + stringSCE +" SCE Enabled");

    delete pVertexDRHistogram;
    delete pDirectionVertexDRHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Direction_Flow_Analysis(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/unmodified.root");
    TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/temp.root");

    TTree *noDirectionTree = (TTree*)f1->Get("Vertex");
    TTree *directionTree = (TTree*)f2->Get("Vertex");

    FindBadEvents(noDirectionTree, directionTree);
    CompareEventType(noDirectionTree, directionTree, 1, 1, 0, true);
    CompareEventType(noDirectionTree, directionTree, 1, 0, 0, true);

    //--------------------------------------------------------------------------------------------------------------------------------------
}
