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
#include "TLine.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/probability_curves_comparison/0.root");
TTree *t1 = (TTree*)f1->Get("DirectionFitTree");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/probability_curves_comparison/1.root");
TTree *t2 = (TTree*)f2->Get("DirectionFitTree");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/probability_curves_comparison/2.root");
TTree *t3 = (TTree*)f3->Get("DirectionFitTree");

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateProbabilityPlot(TTree* tree, TGraph* Probability, std::string forwardsName, std::string backwardsName, int &nBins, bool extendTails, bool writePlots, bool saveHistogram = false)
{
    float histogramWidth(50.f);

    TH1F *forwardsDeltaChiSquared = new TH1F("forwardsDeltaChiSquared","", 80, -histogramWidth, histogramWidth);
    TH1F *backwardsDeltaChiSquared = new TH1F("backwardsDeltaChiSquared","", 80, -histogramWidth, histogramWidth);
    
    int         trueParticleDirection;
    int         fileIdentifier;
    int         eventNumber;
    float       forwardsChiSquared;
    float       backwardsChiSquared;
    
    tree->SetBranchAddress("TrueForwards",                                &trueParticleDirection);
    tree->SetBranchAddress("FileIdentifier",                                       &fileIdentifier);
    tree->SetBranchAddress("EventNumber",                                       &eventNumber);
    tree->SetBranchAddress(forwardsName.c_str(),                                   &forwardsChiSquared);
    tree->SetBranchAddress(backwardsName.c_str(),                                  &backwardsChiSquared);
    
    int nentries = (Int_t)tree->GetEntries();
    
    for (int i = 1; i < nentries; i++)
    {
        tree->GetEntry(i);

        float deltaChiSquaredPerHit(forwardsChiSquared - backwardsChiSquared);

        if (trueParticleDirection == 0)
            backwardsDeltaChiSquared->Fill(deltaChiSquaredPerHit);
        
        if (trueParticleDirection == 1)
            forwardsDeltaChiSquared->Fill(deltaChiSquaredPerHit);
    }

    TCanvas *c5 = new TCanvas("forwardsDeltaChiSquared", "forwardsDeltaChiSquared", 900, 600);
    forwardsDeltaChiSquared->SetYTitle("Number of entries");
    forwardsDeltaChiSquared->SetXTitle("#Delta#chi^{2}/N = (#chi^{2}_{F} - #chi^{2}_{B})/N");
    forwardsDeltaChiSquared->SetTitle("True Forwards Muons: #Delta#chi^{2}/N = (#chi^{2}_{F} - #chi^{2}_{B})/N");
    forwardsDeltaChiSquared->GetYaxis()->SetTitleOffset(1.0);
    if (writePlots)
    {
        forwardsDeltaChiSquared->Draw();
        c5->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Probability/forwardsDeltaChiSquared_pc.pdf");
    }
    
    TCanvas *c6 = new TCanvas("backwardsDeltaChiSquared", "backwardsDeltaChiSquared", 900, 600);
    backwardsDeltaChiSquared->SetYTitle("Number of entries");
    backwardsDeltaChiSquared->SetXTitle("#Delta#chi^{2}/N = (#chi^{2}_{F} - #chi^{2}_{B})/N");
    backwardsDeltaChiSquared->SetTitle("True Backwards Muons: #Delta#chi^{2}/N = (#chi^{2}_{F} - #chi^{2}_{B})/N");
    backwardsDeltaChiSquared->GetYaxis()->SetTitleOffset(1.0);
    backwardsDeltaChiSquared->SetLineColor(kRed);
    if (writePlots)
    {
        backwardsDeltaChiSquared->Draw();
        c6->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Probability/backwardsDeltaChiSquared_pc.pdf");
    }

    if (extendTails)
    {
        int forwardsTailStart(0), backwardsTailStart(0);

        for (int i = 1; i <= nBins; i++)
        {
            int forwardsBinContent(forwardsDeltaChiSquared->GetBinContent(i));
            int backwardsBinContent(backwardsDeltaChiSquared->GetBinContent(i));

            if (forwardsDeltaChiSquared->GetBinCenter(i) > 0 && forwardsBinContent > forwardsDeltaChiSquared->GetEntries() * 0.0025)
                forwardsTailStart = i;

            if (backwardsDeltaChiSquared->GetBinCenter(i) < 0 && backwardsBinContent <= backwardsDeltaChiSquared->GetEntries() * 0.0025)
                backwardsTailStart = i;
        }

        ++forwardsTailStart;

        std::cout << "forwardsDeltaChiSquared->GetEntries() * 0.0025: " << forwardsDeltaChiSquared->GetEntries() * 0.0025 << std::endl;
        std::cout << "backwardsDeltaChiSquared->GetEntries() * 0.0025: " << backwardsDeltaChiSquared->GetEntries() * 0.0025 << std::endl;
        std::cout << "forwardsTailStart: " << forwardsTailStart << std::endl;
        std::cout << "backwardsTailStart: " << backwardsTailStart << std::endl;
    
        //float tailDefinition(10.f);
        //float forwardsTailStart(forwardsDeltaChiSquared->FindBin(tailDefinition));

        float forwardsAverageBinContent(0.f);
        int forwardsNBins(0);

        for (int i = 1; i <= nBins; i++)
        {
            if (i >= forwardsTailStart)
            {
                forwardsAverageBinContent += forwardsDeltaChiSquared->GetBinContent(i);
                forwardsNBins++;
            }
        }

        forwardsAverageBinContent /= forwardsNBins;
        //if (forwardsAverageBinContent < 1.0)
        //    forwardsAverageBinContent = 1.f;

        for (int i = 1; i <= nBins; i++)
        {
            if (i >= forwardsTailStart)
                forwardsDeltaChiSquared->SetBinContent(i, forwardsAverageBinContent);
        }

        //float backwardsTailStart(backwardsDeltaChiSquared->FindBin(-tailDefinition));

        float backwardsAverageBinContent(0);
        int backwardsNBins(0);
        for (int i = 1; i <= nBins; i++)
        {
            if (i <= backwardsTailStart)
            {
                backwardsAverageBinContent += backwardsDeltaChiSquared->GetBinContent(i);
                backwardsNBins++;
            }
        }

        backwardsAverageBinContent /= backwardsNBins;
        //if (backwardsAverageBinContent < 1.0)
        //    backwardsAverageBinContent = 1.f;

        for (int i = 1; i <= nBins; i++)
        {
            if (i <= backwardsTailStart)
                backwardsDeltaChiSquared->SetBinContent(i, backwardsAverageBinContent);
        }

        if (writePlots)
        {
            c5->cd();
            c5->Update();
            forwardsDeltaChiSquared->Draw();
            c5->Update();
            c5->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Probability/forwardsDeltaChiSquared_extendedtails_pc.pdf");
        }
        if (writePlots)
        {
            c6->cd();
            c6->Update();
            backwardsDeltaChiSquared->Draw();
            c6->Update();
            c6->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Probability/backwardsDeltaChiSquared_extendedtails_pc.pdf");
        }
    }

    for (int i = 1; i <= nBins; i++)
    {
        //Normalised distributions
        float forwardsBinEntry(forwardsDeltaChiSquared->GetBinContent(i)/forwardsDeltaChiSquared->GetEntries());
        float backwardsBinEntry(backwardsDeltaChiSquared->GetBinContent(i)/backwardsDeltaChiSquared->GetEntries());

        float binCenter(forwardsDeltaChiSquared->GetBinCenter(i));
        float forwardsProbability(forwardsBinEntry/(forwardsBinEntry + backwardsBinEntry));
        Probability->SetPoint(i, binCenter, forwardsProbability);
    }

    delete forwardsDeltaChiSquared;
    delete backwardsDeltaChiSquared;
    delete c5;
    delete c6;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProbabilityCurves(void)
{
    int nBins(80);
    float histogramWidth(50.0f);

    TGraph *Probability = new TGraph(nBins);
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(kRed);
    Probability->SetTitle("Forwards Probability P_{F} Distribution Comparison;#Delta#chi^{2}_{FB}/N = (#chi^{2}_{F} - #chi^{2}_{B})/N;Forwards Probability P_{F}");
    CreateProbabilityPlot(t1, Probability, "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", nBins, true, false);

    TGraph *Probability_TEF = new TGraph(nBins);
    Probability_TEF->SetMarkerStyle(6);
    Probability_TEF->SetMarkerColor(kBlue);
    CreateProbabilityPlot(t2, Probability_TEF, "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", nBins, true, false);

    TGraph *Probability_Splitting = new TGraph(nBins);
    Probability_Splitting->SetMarkerStyle(6);
    Probability_Splitting->SetMarkerColor(kCyan);
    CreateProbabilityPlot(t3, Probability_Splitting, "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", nBins, true, true);

    //--------------------------------------------------------------------------------------------------------------------------------------

    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    Probability->Draw("AP");
    Probability_TEF->Draw("Psame");
    Probability_Splitting->Draw("Psame");

    Probability->GetXaxis()->SetLimits(-histogramWidth, histogramWidth);
    Probability->GetYaxis()->SetRangeUser(0.0, 1.0);
    Probability_Splitting->GetYaxis()->SetRangeUser(0.0, 1.0);
    TLine *line = new TLine(-histogramWidth,1.0,histogramWidth,1.0);
    line->SetLineColor(kBlue);
    //line->Draw();
    
    auto legend = new TLegend(0.55,0.675,0.85,0.875);
    legend->SetHeader("Colour Legend"); 
    legend->AddEntry(Probability,"Hit Filter","p");
    legend->AddEntry(Probability_TEF,"With Bragg Peak Filter","p");
    legend->AddEntry(Probability_Splitting,"With Splitting","p");
    legend->Draw("same");

    /*
    auto legend = new TLegend(0.55,0.65,0.85,0.85);
    legend->SetHeader("Colour Legend"); 
    legend->AddEntry(Probability,"Track End Filter: No Cut","p");
    legend->AddEntry(Probability_TEF,"Track End Filter: Naive Cut","p");
    legend->AddEntry(Probability_Splitting,"Track End Filter: Multi-N Cut","p");
    legend->Draw("same");
    */

    /*
    auto legend = new TLegend(0.55,0.65,0.85,0.85);
    legend->SetHeader("Colour Legend"); 
    legend->AddEntry(Probability,"Splitting: No Cut","p");
    legend->AddEntry(Probability_TEF,"Splitting: P&C Intersection Cut","p");
    legend->AddEntry(Probability_Splitting,"Splitting: 95% Purity Cut","p");
    legend->Draw("same");
    */

    c0->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Probability/Probability_curves_pc.pdf");
}
