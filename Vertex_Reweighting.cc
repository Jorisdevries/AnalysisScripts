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
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/vertex/" + histogramName + ".png");
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
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/vertex/" + histogramName + ".png");
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

void CreateVertexDRHistogramFromTree(TTree* pTree, TH1F* pHistogram, int targetInteractionType)
{
    int interactionType;
    float vertexDR;

    pTree->SetBranchAddress("interactionType", &interactionType);
    pTree->SetBranchAddress("SCCVertexDR", &vertexDR);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);
    
        if (interactionType != targetInteractionType)
            continue;

        pHistogram->Fill(vertexDR);
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareTrees(TTree* pUnmodifiedTree, TTree* pModifiedTree, int interactionType, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    TH1F *pUnmodifiedHistogram = new TH1F("pUnmodifiedHistogram","", 100, 0, 5);
    TH1F *pModifiedHistogram = new TH1F("pModifiedHistogram","", 100, 0, 5);

    CreateVertexDRHistogramFromTree(pUnmodifiedTree, pUnmodifiedHistogram, interactionType);
    CreateVertexDRHistogramFromTree(pModifiedTree, pModifiedHistogram, interactionType);

    histogramName += std::to_string(interactionType).c_str();
    plotTitle += std::to_string(interactionType).c_str();

    Draw(pUnmodifiedHistogram, pModifiedHistogram, kRed, kBlue, histogramName, xTitle, yTitle, plotTitle); 

    delete pUnmodifiedHistogram;
    delete pModifiedHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareTreeEvents(TTree* pUnmodifiedTree, TTree* pModifiedTree, int targetInteractionType)
{
    int interactionType, fileIdentifier, eventNumber;
    float vertexDR, vertexDRMin, mcLength;

    pUnmodifiedTree->SetBranchAddress("interactionType", &interactionType);
    pUnmodifiedTree->SetBranchAddress("SCCVertexDR", &vertexDR);
    pUnmodifiedTree->SetBranchAddress("VertexDRMin", &vertexDRMin);
    pUnmodifiedTree->SetBranchAddress("MCLength", &mcLength);
    pUnmodifiedTree->SetBranchAddress("FileIdentifier", &fileIdentifier);
    pUnmodifiedTree->SetBranchAddress("EventNumber", &eventNumber);

    int modifiedInteractionType, modifiedFileIdentifier, modifiedEventNumber;
    float modifiedVertexDR, modifiedVertexDRMin, modifiedMCLength;

    pModifiedTree->SetBranchAddress("interactionType", &modifiedInteractionType);
    pModifiedTree->SetBranchAddress("SCCVertexDR", &modifiedVertexDR);
    pModifiedTree->SetBranchAddress("VertexDRMin", &modifiedVertexDRMin);
    pModifiedTree->SetBranchAddress("MCLength", &modifiedMCLength);
    pModifiedTree->SetBranchAddress("FileIdentifier", &modifiedFileIdentifier);
    pModifiedTree->SetBranchAddress("EventNumber", &modifiedEventNumber);

    std::cout << pUnmodifiedTree->GetEntries() << std::endl;
    std::cout << pModifiedTree->GetEntries() << std::endl;

    for (int i = 0; i < pUnmodifiedTree->GetEntries(); i++)
    {
        pUnmodifiedTree->GetEntry(i);

        if (interactionType != targetInteractionType)
            continue;

        for (int j = 0; j < pModifiedTree->GetEntries(); j++)
        {
            pModifiedTree->GetEntry(j);

            if (modifiedInteractionType != targetInteractionType)
                continue;

            if (fileIdentifier != modifiedFileIdentifier || eventNumber != modifiedEventNumber)
                continue;

            if (modifiedVertexDR < 5.0 && vertexDR > 5.0)
            {
                //if (eventNumber < 5)
                std::cout << "Event that is improved by modification at " << fileIdentifier << ":" << eventNumber << " (" << vertexDR << "/" << modifiedVertexDR << ")" << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FindDirectionEvents(TTree* pTree, int targetInteractionType)
{
    int interactionType, fileIdentifier, eventNumber;
    float vertexDR, vertexDRMin, mcLength;

    pTree->SetBranchAddress("interactionType", &interactionType);
    pTree->SetBranchAddress("SCCVertexDR", &vertexDR);
    pTree->SetBranchAddress("VertexDRMin", &vertexDRMin);
    pTree->SetBranchAddress("MCLength", &mcLength);
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier);
    pTree->SetBranchAddress("EventNumber", &eventNumber);

    int nDirectionCandidates(0), nChannelEevents(0);

    std::cout << pTree->GetEntries() << std::endl;

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);
    
        if (interactionType != targetInteractionType)
            continue;

        ++nChannelEevents;

        if ((vertexDR - vertexDRMin)/mcLength > 0.8 && (vertexDR - vertexDRMin)/mcLength < 1.2 && vertexDRMin < 5.0)
        {
            if (eventNumber <= 3)
                std::cout << "Candidate direction event at " << fileIdentifier << ":" << eventNumber << " (" << vertexDR << "/" << mcLength << ")" << std::endl;

            ++nDirectionCandidates;
        }

    }

    std::cout << "Percentage of direction candidates: " << 100.0 * ((float)nDirectionCandidates/nChannelEevents) << std::endl;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void Vertex_Reweighting(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_length_drmin/unmodified.root");
    TTree *t1 = (TTree*)f1->Get("Vertex");

    TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_length_drmin/only_80percent_direction.root");
    TTree *t2 = (TTree*)f2->Get("Vertex");

    TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_length_drmin/reweighted.root");
    TTree *t3 = (TTree*)f3->Get("Vertex");

    TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_length_drmin/directionfeature_svm.root");
    TTree *t4 = (TTree*)f4->Get("Vertex");

    TFile *f5 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_length_drmin/cheated.root");
    TTree *t5 = (TTree*)f5->Get("Vertex");

    TFile *f6 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/vertex_length_drmin/z_candidates.root");
    TTree *t6 = (TTree*)f6->Get("Vertex");

    int interactionType(0);

    CompareTreeEvents(t1, t6, interactionType);

    /*
    std::cout << "Unmodified: " << std::endl;
    FindDirectionEvents(t1, interactionType);

    std::cout << "Only good direction: " << std::endl;
    FindDirectionEvents(t2, interactionType);

    std::cout << "Reweighted: " << std::endl;
    FindDirectionEvents(t3, interactionType);

    std::cout << "Direction feature SVM: " << std::endl;
    FindDirectionEvents(t4, interactionType);
    */

    TH1F *pUnmodifiedHistogram = new TH1F("pUnmodifiedHistogram","", 100, 0, 5);
    TH1F *pGoodDirectionOnlyHistogram = new TH1F("pGoodDirectionOnlyHistogram","", 100, 0, 5);
    TH1F *pReweightedHistogram = new TH1F("pReweightedHistogram","", 100, 0, 5);
    TH1F *pDirectionRegionSVMHistogram = new TH1F("pDirectionRegionSVMHistogram","", 100, 0, 5);
    TH1F *pCheatedHistogram = new TH1F("pCheatedHistogram","", 100, 0, 5);
    TH1F *pZCandidatesHistogram = new TH1F("pZCandidatesHistogram","", 100, 0, 5);

    CreateVertexDRHistogramFromTree(t1, pUnmodifiedHistogram, interactionType);
    CreateVertexDRHistogramFromTree(t2, pGoodDirectionOnlyHistogram, interactionType);
    CreateVertexDRHistogramFromTree(t3, pReweightedHistogram, interactionType);
    CreateVertexDRHistogramFromTree(t4, pDirectionRegionSVMHistogram, interactionType);
    CreateVertexDRHistogramFromTree(t5, pCheatedHistogram, interactionType);
    CreateVertexDRHistogramFromTree(t6, pZCandidatesHistogram, interactionType);

    pUnmodifiedHistogram->SetLineColor(kBlack);
    pGoodDirectionOnlyHistogram->SetLineColor(kRed);
    pReweightedHistogram->SetLineColor(kOrange+1);
    pDirectionRegionSVMHistogram->SetLineColor(kBlue);
    pCheatedHistogram->SetLineColor(kGreen+1);
    pZCandidatesHistogram->SetLineColor(kMagenta);

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    std::string title("Vertex #DeltaR for Event Type " + std::to_string(interactionType));

    pUnmodifiedHistogram->SetXTitle("Vertex #DeltaR");
    pUnmodifiedHistogram->SetYTitle("Number of Events");
    pUnmodifiedHistogram->SetTitle(title.c_str());
    pUnmodifiedHistogram->GetYaxis()->SetTitleOffset(1.0);

    pUnmodifiedHistogram->Draw();
    pZCandidatesHistogram->Draw("same");

    //--------------------------------------------------------------------------------------------------------------------------------------
}
