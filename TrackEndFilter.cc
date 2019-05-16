#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLegend.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>

//------------------------------------------------------------------------------------------------------------------------------------------

float GetSymmetricCut(TH1F* braggPeakNearestRatio, TH1F* nonBraggPeakNearestRatio, bool symmetric)
{
    float nearestSymmetricCut(0.f);
    for (int bin = 1; bin <= braggPeakNearestRatio->GetXaxis()->GetNbins(); bin++)
    {
        float braggBinEntry(braggPeakNearestRatio->GetBinContent(bin));
        float nextBraggBinEntry(braggPeakNearestRatio->GetBinContent(bin+1));
        float nextNextBraggBinEntry(braggPeakNearestRatio->GetBinContent(bin+2));
        float nonBraggBinEntry(nonBraggPeakNearestRatio->GetBinContent(bin));
        float nextNonBraggBinEntry(nonBraggPeakNearestRatio->GetBinContent(bin+1));
        float nextNextNonBraggBinEntry(nonBraggPeakNearestRatio->GetBinContent(bin+2));

        if (braggPeakNearestRatio->GetBinCenter(bin) > 1.4 && nonBraggBinEntry/braggBinEntry > 1.0 && nextNonBraggBinEntry/nextBraggBinEntry > 1.0)
        {
            if (symmetric)
                nearestSymmetricCut = std::abs(1.0 - braggPeakNearestRatio->GetBinCenter(bin));
            else
                nearestSymmetricCut = braggPeakNearestRatio->GetBinCenter(bin);
            break;
        }
    }

    return nearestSymmetricCut;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackEndFilter( void )

{
    TFile *f = TFile::Open("/usera/jjd49/new_LAr/CondorUtilities/saved_results/trackendfilter/trackendfilter.root");
    
    TTree *t1 = (TTree*)f->Get("Tree");
    int nentries = (Int_t)t1->GetEntries();

    TH1F *braggPeakNearestRatio = new TH1F("braggPeakNearestRatio","", 75, 0.0, 3.0);
    TH1F *nonBraggPeakNearestRatio = new TH1F("nonBraggPeakNearestRatio","", 75, 0.0, 3.0);

    TH1F *braggPeakPlusMinusNRatio = new TH1F("braggPeakPlusMinusNRatio","", 75, 0.0, 3.0);
    TH1F *nonBraggPeakPlusMinusNRatio = new TH1F("nonBraggPeakPlusMinusNRatio","", 75, 0.0, 3.0);

    TH1F *braggPeakDistanceFromBodyQoverW = new TH1F("braggPeakDistanceFromBodyQoverW","", 75, 0.0, 10.0);
    TH1F *nonBraggPeakDistanceFromBodyQoverW = new TH1F("nonBraggPeakDistanceFromBodyQoverW","", 75, 0.0, 10.0);

    int         trueParticleDirection;
    int         nHitsToSkip;
    int         inBraggPeak;
    float       trackEndRange;
    float       nearestRatio;
    float       plusMinusNRatio;
    float       distanceFromBodyQoverW;
    int         fileIdentifier;
    int         eventNumber;
    
    t1->SetBranchAddress("trueParticleDirection",                               &trueParticleDirection);
    t1->SetBranchAddress("nHitsToSkip",                                         &nHitsToSkip);
    t1->SetBranchAddress("trackEndRange",                                       &trackEndRange);
    t1->SetBranchAddress("inBraggPeak",                                         &inBraggPeak);
    t1->SetBranchAddress("nearestRatio",                                        &nearestRatio);
    t1->SetBranchAddress("plusMinusNRatio",                                     &plusMinusNRatio);
    t1->SetBranchAddress("distanceFromBodyQoverW",                              &distanceFromBodyQoverW);
    t1->SetBranchAddress("fileIdentifier",                                      &fileIdentifier);
    t1->SetBranchAddress("eventNumber",                                         &eventNumber);
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::vector<float> braggPeakNearestVector, braggPeakPlusMinusNVector, braggPeakDistanceQWVector;
    int nBraggBackwardsEntries(0), nBraggForwardsEntries(0);
    int nNonBraggBackwardsEntries(0), nNonBraggForwardsEntries(0);

    for (int i = 0; i <= nentries; i++)
    {
        t1->GetEntry(i);

        if (nBraggBackwardsEntries >= 20000 && inBraggPeak == 1 && trueParticleDirection == 0)
            continue;
        if (nBraggForwardsEntries >= 20000 && inBraggPeak == 1 && trueParticleDirection == 1)
            continue;
        if (nNonBraggBackwardsEntries >= 20000 && inBraggPeak == 0 && trueParticleDirection == 0)
            continue;
        if (nNonBraggForwardsEntries >= 20000 && inBraggPeak == 0 && trueParticleDirection == 1)
            continue;

        if (i == 0)
            std::cout << "Number of hits skipped over N = " << nHitsToSkip << std::endl;

        if (inBraggPeak == 0)
        {
            if (trueParticleDirection == 0)
                nNonBraggBackwardsEntries++;
            if (trueParticleDirection == 1)
                nNonBraggForwardsEntries++;

            nonBraggPeakNearestRatio->Fill(nearestRatio);
            nonBraggPeakPlusMinusNRatio->Fill(plusMinusNRatio);
            nonBraggPeakDistanceFromBodyQoverW->Fill(distanceFromBodyQoverW); 
        }
        if (inBraggPeak == 1)
        {
            if (trueParticleDirection == 0)
                nBraggBackwardsEntries++;
            if (trueParticleDirection == 1)
                nBraggForwardsEntries++;

            braggPeakNearestRatio->Fill(nearestRatio);
            braggPeakPlusMinusNRatio->Fill(plusMinusNRatio);
            braggPeakDistanceFromBodyQoverW->Fill(distanceFromBodyQoverW); 

            braggPeakNearestVector.push_back(nearestRatio);
            braggPeakPlusMinusNVector.push_back(plusMinusNRatio);
            braggPeakDistanceQWVector.push_back(distanceFromBodyQoverW);
        }
    }

    std::cout << "Number of forwards Bragg hits: " << nBraggForwardsEntries << std::endl;
    std::cout << "Number of backwards Bragg hits: " << nBraggBackwardsEntries << std::endl;
    std::cout << "Number of forwards NonBragg hits: " << nNonBraggForwardsEntries << std::endl;
    std::cout << "Number of backwards NonBragg hits: " << nNonBraggBackwardsEntries << std::endl;

    Double_t scale1 = 1/nonBraggPeakNearestRatio->Integral();
    nonBraggPeakNearestRatio->Scale(scale1);
    Double_t scale2 = 1/braggPeakNearestRatio->Integral();
    braggPeakNearestRatio->Scale(scale2);

    Double_t scale3 = 1/nonBraggPeakPlusMinusNRatio->Integral();
    nonBraggPeakPlusMinusNRatio->Scale(scale3);
    Double_t scale4 = 1/braggPeakPlusMinusNRatio->Integral();
    braggPeakPlusMinusNRatio->Scale(scale4);

    Double_t scale5 = 1/nonBraggPeakDistanceFromBodyQoverW->Integral();
    nonBraggPeakDistanceFromBodyQoverW->Scale(scale5);
    Double_t scale6 = 1/braggPeakDistanceFromBodyQoverW->Integral();
    braggPeakDistanceFromBodyQoverW->Scale(scale6);

    std::sort(braggPeakNearestVector.begin(), braggPeakNearestVector.end());
    std::sort(braggPeakPlusMinusNVector.begin(), braggPeakPlusMinusNVector.end());
    std::sort(braggPeakDistanceQWVector.begin(), braggPeakDistanceQWVector.end());

    //std::cout << "Bragg peak 0.05 cut value: " << *(braggPeakNearestVector.begin() + 0.05 * braggPeakNearestVector.size()) << std::endl;
    //std::cout << "Bragg peak 0.95 cut value: " << *(braggPeakNearestVector.begin() + 0.95 * braggPeakNearestVector.size()) << std::endl;

    float braggPeakLowerBound(*(braggPeakNearestVector.begin() + 0.05 * braggPeakNearestVector.size())), braggPeakUpperBound(*(braggPeakNearestVector.begin() + 0.95 * braggPeakNearestVector.size()));

    //std::cout << "Bragg peak +/- N 0.05 cut value: " << *(braggPeakNearestVector.begin() + 0.05 * braggPeakNearestVector.size()) << std::endl;
    //std::cout << "Bragg peak +/- N 0.95 cut value: " << *(braggPeakNearestVector.begin() + 0.95 * braggPeakNearestVector.size()) << std::endl;

    float braggPeakPlusMinusNLowerBound(*(braggPeakPlusMinusNVector.begin() + 0.05 * braggPeakPlusMinusNVector.size())), braggPeakPlusMinusNUpperBound(*(braggPeakPlusMinusNVector.begin() + 0.95 * braggPeakPlusMinusNVector.size()));

    std::cout << "Bragg peak 0.05 Q/w - Q/w_av. cut value: " << *(braggPeakDistanceQWVector.begin() + 0.05 * braggPeakDistanceQWVector.size()) << std::endl;
    //std::cout << "Bragg peak Q/w - Q/w_av. 0.95 cut value: " << *(braggPeakDistanceQWVector.begin() + 0.95 * braggPeakDistanceQWVector.size()) << std::endl;

    float distanceFromQWLowerBound(*(braggPeakDistanceQWVector.begin() + 0.05 * braggPeakDistanceQWVector.size())), distanceFromQWUpperBound(*(braggPeakDistanceQWVector.begin() + 0.95 * braggPeakDistanceQWVector.size()));


    //Find first bin for which there are more bad than good entries
    float nearestSymmetricCut(GetSymmetricCut(braggPeakNearestRatio, nonBraggPeakNearestRatio, true));
    std::cout << "nearestSymmetricCut: " << nearestSymmetricCut << std::endl;
    std::cout << "Bragg peak fraction of events retained: " << braggPeakNearestRatio->Integral(braggPeakNearestRatio->FindBin(1.0 - nearestSymmetricCut), braggPeakNearestRatio->FindBin(1.0 + nearestSymmetricCut)) << std::endl; 
    std::cout << "Non-Bragg peak fraction of events retained: " << nonBraggPeakNearestRatio->Integral(nonBraggPeakNearestRatio->FindBin(1.0 - nearestSymmetricCut), nonBraggPeakNearestRatio->FindBin(1.0 + nearestSymmetricCut)) << std::endl; 
    float plusMinusNSymmetricCut(GetSymmetricCut(braggPeakPlusMinusNRatio, nonBraggPeakPlusMinusNRatio, true));
    std::cout << "plusMinusNSymmetricCut: " << plusMinusNSymmetricCut << std::endl;
    std::cout << "Bragg peak fraction of events retained: " << braggPeakPlusMinusNRatio->Integral(braggPeakPlusMinusNRatio->FindBin(1.0 - plusMinusNSymmetricCut), braggPeakPlusMinusNRatio->FindBin(1.0 + plusMinusNSymmetricCut)) << std::endl; 
    std::cout << "Non-Bragg peak fraction of events retained: " << nonBraggPeakPlusMinusNRatio->Integral(nonBraggPeakPlusMinusNRatio->FindBin(1.0 - plusMinusNSymmetricCut), nonBraggPeakPlusMinusNRatio->FindBin(1.0 + plusMinusNSymmetricCut)) << std::endl; 
    float deviationFromAverageChargeCut(GetSymmetricCut(braggPeakDistanceFromBodyQoverW, nonBraggPeakDistanceFromBodyQoverW, false));
    std::cout << "deviationFromAverageChargeCut: " << deviationFromAverageChargeCut << std::endl;
    std::cout << "Bragg peak fraction of events retained: " << braggPeakDistanceFromBodyQoverW->Integral(braggPeakDistanceFromBodyQoverW->FindBin(0.0), braggPeakDistanceFromBodyQoverW->FindBin(deviationFromAverageChargeCut)) << std::endl; 
    std::cout << "Non-Bragg peak fraction of events retained: " << nonBraggPeakDistanceFromBodyQoverW->Integral(nonBraggPeakDistanceFromBodyQoverW->FindBin(0.0), nonBraggPeakDistanceFromBodyQoverW->FindBin(deviationFromAverageChargeCut)) << std::endl; 


    //--------------------------------------------------------------------------------------------------------------------------------------
    
    TCanvas *c1 = new TCanvas("braggPeakNearestRatio", "braggPeakNearestRatio", 900, 600);
    braggPeakNearestRatio->SetXTitle("max(#hat{Q}_{i}/#hat{Q}_{i#pm1})");
    braggPeakNearestRatio->SetYTitle("Number of entries");
    braggPeakNearestRatio->SetTitle("Bragg Peak and non-Bragg Peak Ratios: max(#hat{Q}_{i}/#hat{Q}_{i#pm1})");
    braggPeakNearestRatio->GetYaxis()->SetTitleOffset(1.5);
    //braggPeakNearestRatio->GetYaxis()->SetRangeUser(0.0, 0.2);
    braggPeakNearestRatio->SetLineColor(kGreen+3);
    braggPeakNearestRatio->SetStats(kFALSE);
    braggPeakNearestRatio->Draw("hist");
    nonBraggPeakNearestRatio->SetLineColor(kRed);
    nonBraggPeakNearestRatio->Draw("histsame");

    TH1F *c_clone1 = (TH1F*)braggPeakNearestRatio->Clone(); 
    int lower1(c_clone1->GetXaxis()->FindBin(1.0 - nearestSymmetricCut)), upper1(c_clone1->GetXaxis()->FindBin(1.0 + nearestSymmetricCut));
    c_clone1->SetFillColor(kGreen+3);
    c_clone1->SetFillStyle(3003);
    c_clone1->GetXaxis()->SetRange(lower1,upper1);
    c_clone1->Draw("histsame");

    auto legend1 = new TLegend(0.5,0.65,0.75,0.85);
    legend1->SetHeader("Colour Legend"); 
    legend1->AddEntry(braggPeakNearestRatio,"Bragg Peak Hits","l");
    legend1->AddEntry(nonBraggPeakNearestRatio,"Non-Bragg Peak Hits","l");
    legend1->AddEntry(c_clone1,"Bragg Peak Entries Preserved by Cut","f");
    legend1->Draw("histsame");

    c1->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/TrackEndFilter/braggPeakNearestRatio.pdf");

    TCanvas *c2 = new TCanvas("braggPeakPlusMinusNRatio", "braggPeakPlusMinusNRatio", 900, 600);
    braggPeakPlusMinusNRatio->SetXTitle("max(#hat{Q}_{i}/#hat{Q}_{i#pm3})");
    braggPeakPlusMinusNRatio->SetYTitle("Number of entries");
    braggPeakPlusMinusNRatio->SetTitle("Bragg Peak and non-Bragg Peak Ratios: max(#hat{Q}_{i}/#hat{Q}_{i#pm3})");
    braggPeakPlusMinusNRatio->GetYaxis()->SetTitleOffset(1.5);
    //braggPeakPlusMinusNRatio->GetYaxis()->SetRangeUser(0.0, 0.12);
    braggPeakPlusMinusNRatio->SetLineColor(kGreen+3);
    braggPeakPlusMinusNRatio->SetStats(kFALSE);
    braggPeakPlusMinusNRatio->Draw("hist");
    nonBraggPeakPlusMinusNRatio->SetLineColor(kRed);
    nonBraggPeakPlusMinusNRatio->Draw("histsame");

    TH1F *c_clone2 = (TH1F*)braggPeakPlusMinusNRatio->Clone(); 
    int lower2(c_clone2->GetXaxis()->FindBin(1.0 - plusMinusNSymmetricCut)), upper2(c_clone2->GetXaxis()->FindBin(1.0 + plusMinusNSymmetricCut));
    c_clone2->SetFillColor(kGreen+3);
    c_clone2->SetFillStyle(3002);
    c_clone2->GetXaxis()->SetRange(lower2,upper2);
    c_clone2->Draw("histsame");

    auto legend2 = new TLegend(0.5,0.65,0.75,0.85);
    legend2->SetHeader("Colour Legend"); 
    legend2->AddEntry(braggPeakPlusMinusNRatio,"Bragg Peak Hits","l");
    legend2->AddEntry(nonBraggPeakPlusMinusNRatio,"Non-Bragg Peak Hits","l");
    legend2->AddEntry(c_clone2,"Bragg Peak Entries Preserved by Cut","f");
    legend2->Draw("histsame");

    c2->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/TrackEndFilter/braggPeakPlusMinusNRatio.pdf");

    TCanvas *c3 = new TCanvas("braggPeakDistanceFromBodyQoverW", "braggPeakDistanceFromBodyQoverW", 900, 600);
    braggPeakDistanceFromBodyQoverW->SetXTitle("|#hat{Q}_{i} - #hat{Q}_{av.}|");
    braggPeakDistanceFromBodyQoverW->SetYTitle("Number of entries");
    braggPeakDistanceFromBodyQoverW->SetTitle("Bragg Peak and non-Bragg Peak Distance from #hat{Q}_{av.}: |#hat{Q}_{i} - #hat{Q}_{av.}|");
    braggPeakDistanceFromBodyQoverW->GetYaxis()->SetTitleOffset(1.5);
    braggPeakDistanceFromBodyQoverW->SetLineColor(kGreen+3);
    braggPeakDistanceFromBodyQoverW->SetStats(kFALSE);
    braggPeakDistanceFromBodyQoverW->Draw("hist");
    nonBraggPeakDistanceFromBodyQoverW->SetLineColor(kRed);
    nonBraggPeakDistanceFromBodyQoverW->Draw("histsame");

    TH1F *c_clone3 = (TH1F*)braggPeakDistanceFromBodyQoverW->Clone(); 
    int lower3(c_clone3->GetXaxis()->FindBin(0.0)), upper3(c_clone3->GetXaxis()->FindBin(deviationFromAverageChargeCut));
    c_clone3->SetFillColor(kGreen+3);
    c_clone3->SetFillStyle(3002);
    c_clone3->GetXaxis()->SetRange(0.f,upper3);
    c_clone3->Draw("histsame");

    auto legend3 = new TLegend(0.5,0.65,0.75,0.85);
    legend3->SetHeader("Colour Legend"); 
    legend3->AddEntry(braggPeakDistanceFromBodyQoverW,"Bragg Peak Hits","l");
    legend3->AddEntry(nonBraggPeakDistanceFromBodyQoverW,"Non-Bragg Peak Hits","l");
    legend3->AddEntry(c_clone3,"Lowest 95% of Bragg Peak Entries","f");
    legend3->Draw("histsame");

    c3->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/TrackEndFilter/braggPeakDistanceFromBodyQoverW.pdf");

}
