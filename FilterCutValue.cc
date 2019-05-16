#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>

//------------------------------------------------------------------------------------------------------------------------------------------

float GetCrossoverCut(TH1F* goodDistribution, TH1F* badDistribution)
{
    float nearestCrossoverCut(0.f);
    for (int bin = 1; bin <= goodDistribution->GetXaxis()->GetNbins(); bin++)
    {   
        float goodBinEntry(goodDistribution->GetBinContent(bin));
        float badBinEntry(badDistribution->GetBinContent(bin));

        if (badBinEntry > goodBinEntry)
        {   
            nearestCrossoverCut = goodDistribution->GetBinCenter(bin);
            break;
        }   
    }   

    return nearestCrossoverCut;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void new_filtertest( void )
{
    TFile *f = TFile::Open("/usera/jjd49/new_LAr/CondorUtilities/saved_results/filtering/filter_02Nov.root");
    
    TTree *t1 = (TTree*)f->Get("Tree");
    int nentries = (Int_t)t1->GetEntries();

    TTree *t3 = (TTree*)f->Get("Tree3");
    int nentries3 = (Int_t)t3->GetEntries();

    TH1F *nHitFiltered = new TH1F("nHitFiltered","", 100, 0.0, 400.0);   
    TH1F *hFractionFilteredBelongingToMuon = new TH1F("hFractionFilteredBelongingToMuon","", 100, 0.0, 1.1);
    TH1F *hFractionFilteredBelongingToNonMuon = new TH1F("hFractionFilteredBelongingToNonMuon","", 75, 0.0, 1.1);
    TH1F *hFractionNonMuonHitsFiltered = new TH1F("hFractionNonMuonHitsFiltered","", 100, 0.0, 1.1);
    TH1F *hFractionMuonHitsFiltered = new TH1F("hFractionMuonHitsFiltered","", 100, 0.0, 1.1);
    TH1F *hFractionClusterHitsFiltered = new TH1F("hFractionClusterHitsFiltered","", 100, 0.0, 1.1);
    TH1F *hFractionImpureHits = new TH1F("hFractionImpureHits","", 100, 0.0, 1.1);

    TH2F *fractionImpureHitsNumberHits = new TH2F("fractionImpureHitsNumberHits","", 50, 0.0, 0.5, 50, 0, 1200);
    TH2F *fractionImpureHitsTrueEnergy = new TH2F("fractionImpureHitsTrueEnergy","", 50, 0.0, 0.5, 50, 0, 2.0);

    int     totalNumberHits;
    float   trueEnergy;
    float     trackLength;
    int    nTotalHitsFiltered;
    int     nMuonHits;
    int     nNonMuonHits;
    int   nMuonHitsFiltered;
    int   nNonMuonHitsFiltered;
    float   fractionImpureHits;

    int fileIdentifier;
    int eventNumber;
    
    t1->SetBranchAddress("totalNumberHits",                                   &totalNumberHits);
    t1->SetBranchAddress("globalTrueTargetParticleEnergy",                                   &trueEnergy);
    t1->SetBranchAddress("nTotalHitsFiltered",                                &nTotalHitsFiltered);
    t1->SetBranchAddress("trackLength",                                       &trackLength);
    t1->SetBranchAddress("nMuonHits",                   &nMuonHits);
    t1->SetBranchAddress("nNonMuonHits",                &nNonMuonHits);
    t1->SetBranchAddress("nMuonHitsFiltered",                       &nMuonHitsFiltered);
    t1->SetBranchAddress("nNonMuonHitsFiltered",                          &nNonMuonHitsFiltered);
    t1->SetBranchAddress("fractionImpureHits",                          &fractionImpureHits);

    t1->SetBranchAddress("fileIdentifier",                          &fileIdentifier);
    t1->SetBranchAddress("eventNumber",                          &eventNumber);


    //*********************************************************************************************************


    float globalTrueTargetParticleEnergy;
    int nNeighboursToConsider;
    int isMuonHit;
    float nearestNeighboursDistanceSum;
    int numberHits;
    //int fileIdentifier;
    //int eventNumber;

    
    t3->SetBranchAddress("globalTrueTargetParticleEnergy",                &globalTrueTargetParticleEnergy);
    t3->SetBranchAddress("nNeighboursToConsider",                       &nNeighboursToConsider);
    t3->SetBranchAddress("isMuonHit",                          &isMuonHit);
    t3->SetBranchAddress("nearestNeighboursDistanceSum",                          &nearestNeighboursDistanceSum);
    t3->SetBranchAddress("numberHits",                          &numberHits);
    //t3->SetBranchAddress("fileIdentifier",                          &fileIdentifier);
    //t3->SetBranchAddress("eventNumber",                          &eventNumber);

    TH1F *nearestDistanceSum = new TH1F("nearestDistanceSum","", 75, 0.0, 1.0);   
    TH1F *nearestDistanceSumImpure = new TH1F("nearestDistanceSumImpure","", 75, 0.0, 1.0);   
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::cout << "There are " << nentries << " entries in the tree" << std::endl;
    
    for (int i = 0; i < nentries; i++)
    {
        t1->GetEntry(i);

        if (nNonMuonHits == 0)
            continue;

        float fractionFilteredBelongingToMuon(((float)nMuonHitsFiltered)/nTotalHitsFiltered), fractionFilteredBelongingToNonMuon(((float)nNonMuonHitsFiltered)/nTotalHitsFiltered),
        fractionMuonHitsFiltered(((float)nMuonHitsFiltered)/nMuonHits), fractionNonMuonHitsFiltered(((float)nNonMuonHitsFiltered)/nNonMuonHits), fractionClusterHitsFiltered(((float)nTotalHitsFiltered)/totalNumberHits);
        
        if (fractionMuonHitsFiltered >= 0.5)
            std::cout << "Filter error: " << fileIdentifier << ": " << eventNumber << std::endl;

        nHitFiltered->Fill(nTotalHitsFiltered);
        hFractionFilteredBelongingToMuon->Fill(fractionFilteredBelongingToMuon);
        hFractionFilteredBelongingToNonMuon->Fill(fractionFilteredBelongingToNonMuon);
        hFractionNonMuonHitsFiltered->Fill(fractionNonMuonHitsFiltered);
        hFractionMuonHitsFiltered->Fill(fractionMuonHitsFiltered);
        hFractionClusterHitsFiltered->Fill(fractionClusterHitsFiltered);
        hFractionImpureHits->Fill(fractionImpureHits);

        fractionImpureHitsNumberHits->Fill(fractionImpureHits, totalNumberHits);
        fractionImpureHitsTrueEnergy->Fill(fractionImpureHits, trueEnergy);
    }

    std::vector<float> nearestDistanceSumVector;
    for (int i = 0; i < nentries3; i++)
    {
        t3->GetEntry(i);

//        if (numberHits <= 300 || numberHits > 400)
 //           continue;

        if (isMuonHit == 1)
        {
            nearestDistanceSum->Fill(nearestNeighboursDistanceSum);
            nearestDistanceSumVector.push_back(nearestNeighboursDistanceSum);
        }
        if (isMuonHit == 0)
            nearestDistanceSumImpure->Fill(nearestNeighboursDistanceSum);
    }

    Double_t scale1 = 1/nearestDistanceSum->Integral();
    nearestDistanceSum->Scale(scale1);
    Double_t scale2 = 1/nearestDistanceSumImpure->Integral();
    nearestDistanceSumImpure->Scale(scale2);

    std::sort(nearestDistanceSumVector.begin(), nearestDistanceSumVector.end());
    //float cutValue(nearestDistanceSumVector.at(0.9 * nearestDistanceSumVector.size()));
    float cutValue = GetCrossoverCut(nearestDistanceSum, nearestDistanceSumImpure);

    std::cout << "Pure hits cut: " << cutValue << std::endl;
    std::cout << "Fraction of pure hits retained by cut: " << nearestDistanceSum->Integral(0.0, nearestDistanceSum->GetXaxis()->FindBin(cutValue)) << std::endl;

    TCanvas *c0 = new TCanvas("nearestDistanceSum", "nearestDistanceSum", 900, 600);
    nearestDistanceSum->SetXTitle("Fractional position p_{f} of hits ranked by M_{5}");
    nearestDistanceSum->SetYTitle("Number of entries");
    nearestDistanceSum->SetTitle("Fractional position p_{f} of hits ranked by M_{5} for pure & impure Hits (k = 5) for all N");
    nearestDistanceSum->GetYaxis()->SetTitleOffset(1.0);
    nearestDistanceSum->GetYaxis()->SetRangeUser(0.0, 0.07);
    nearestDistanceSum->Draw();
    nearestDistanceSumImpure->SetLineColor(kRed);
    nearestDistanceSumImpure->Draw("same");

    TH1F *c_clone1 = (TH1F*)nearestDistanceSum->Clone();
    int lower1(c_clone1->GetXaxis()->FindBin(0.0)), upper1(c_clone1->GetXaxis()->FindBin(cutValue) - 1);
    c_clone1->SetFillColor(kBlue);
    c_clone1->SetFillStyle(3001);
    c_clone1->GetXaxis()->SetRange(lower1,upper1);
    c_clone1->Draw("same");
 
    auto legend1 = new TLegend(0.5,0.65,0.75,0.85);
    legend1->SetHeader("Colour Legend");
    legend1->AddEntry(nearestDistanceSum,"Pure Hits","l");
    legend1->AddEntry(nearestDistanceSumImpure,"Impure Hits","l");
    legend1->AddEntry(c_clone1,"Pure Hits Retained by Cut","f");
    legend1->Draw("same");

    c0->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/nearestDistanceSumImpure_allN.png");
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    TCanvas *c1 = new TCanvas("nHitFiltered", "nHitFiltered", 900, 600);
    nHitFiltered->SetXTitle("Total Number of Hits Removed By Filter");
    nHitFiltered->SetYTitle("Number of entries");
    nHitFiltered->SetTitle("Number of Hits Filtered");
    nHitFiltered->GetYaxis()->SetTitleOffset(1.0);
    nHitFiltered->Draw();
    c1->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/nHitsFiltered.png");
    
    TCanvas *c2 = new TCanvas("hFractionFilteredBelongingToMuon", "hFractionFilteredBelongingToMuon", 900, 600);
    hFractionFilteredBelongingToMuon->SetXTitle("Fraction of Filtered Hits that belong to the Muon");
    hFractionFilteredBelongingToMuon->SetYTitle("Number of entries");
    hFractionFilteredBelongingToMuon->SetTitle("Fraction of Filtered Hits that belong to the Muon");
    hFractionFilteredBelongingToMuon->GetYaxis()->SetTitleOffset(1.0);
    hFractionFilteredBelongingToMuon->Draw();
    c2->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/fractionFilteredBelongingToMuon.png");
    
    TCanvas *c3 = new TCanvas("hFractionFilteredBelongingToNonMuon", "hFractionFilteredBelongingToNonMuon", 900, 600);
    hFractionFilteredBelongingToNonMuon->SetXTitle("Filter Purity P");
    hFractionFilteredBelongingToNonMuon->SetYTitle("Number of entries");
    hFractionFilteredBelongingToNonMuon->SetTitle("Filter Purity P (Fraction of Filtered Hits that are Impure Hits)");
    hFractionFilteredBelongingToNonMuon->GetYaxis()->SetTitleOffset(1.0);
    hFractionFilteredBelongingToNonMuon->Draw();
    c3->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/filterPurity.png");
    
    TCanvas *c4 = new TCanvas("hFractionNonMuonHitsFiltered", "hFractionNonMuonHitsFiltered", 900, 600);
    hFractionNonMuonHitsFiltered->SetXTitle("Filter Completeness C");
    hFractionNonMuonHitsFiltered->SetYTitle("Number of entries");
    hFractionNonMuonHitsFiltered->SetTitle("Filter Completeness C (Fraction of Total Impure Hits Filtered)");
    hFractionNonMuonHitsFiltered->GetYaxis()->SetTitleOffset(1.0);
    hFractionNonMuonHitsFiltered->Draw();
    c4->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/filterCompleteness.png");
    
    TCanvas *c5 = new TCanvas("hFractionMuonHitsFiltered", "hFractionMuonHitsFiltered", 900, 600);
    hFractionMuonHitsFiltered->SetXTitle("Fraction of Total Muon Hits Filtered");
    hFractionMuonHitsFiltered->SetYTitle("Number of entries");
    hFractionMuonHitsFiltered->SetTitle("Fraction of Total Muon Hits Filtered");
    hFractionMuonHitsFiltered->GetYaxis()->SetTitleOffset(1.0);
    hFractionMuonHitsFiltered->Draw();
    c5->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/fractionMuonHitsFiltered.png");

    TCanvas *c55 = new TCanvas("hFractionClusterHitsFiltered", "hFractionClusterHitsFiltered", 900, 600);
    hFractionClusterHitsFiltered->SetXTitle("Fraction of Total Cluster Hits Filtered");
    hFractionClusterHitsFiltered->SetYTitle("Number of entries");
    hFractionClusterHitsFiltered->SetTitle("Fraction of Total Cluster Hits Filtered");
    hFractionClusterHitsFiltered->GetYaxis()->SetTitleOffset(1.0);
    hFractionClusterHitsFiltered->Draw();
    c55->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/fractionClusterHitsFiltered.png");

    TCanvas *c555 = new TCanvas("hFractionImpureHits", "hFractionImpureHits", 900, 600);
    hFractionImpureHits->SetXTitle("Fraction of Impure Hits F");
    hFractionImpureHits->SetYTitle("Number of entries");
    hFractionImpureHits->SetTitle("Fraction of Impure Hits F in Cluster");
    hFractionImpureHits->GetYaxis()->SetTitleOffset(1.0);
    hFractionImpureHits->Draw();
    c555->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/fractionImpureHits.png");

    TCanvas *c6 = new TCanvas("fractionImpureHitsNumberHits", "fractionImpureHitsNumberHits", 900, 600);
    fractionImpureHitsNumberHits->SetXTitle("Fraction of Impure Hits F (2nd MC parent contributes #geq 15% of charge)");
    fractionImpureHitsNumberHits->SetYTitle("Number of Cluster Hits N");
    fractionImpureHitsNumberHits->SetTitle("Fraction of Impure Hits F with Number of Cluster Hits N");
    fractionImpureHitsNumberHits->GetYaxis()->SetTitleOffset(1.2);
    fractionImpureHitsNumberHits->Draw("COLZ");
    c6->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/fractionImpureHitsNumberHits.png");

    TCanvas *c7 = new TCanvas("fractionImpureHitsTrueEnergy", "fractionImpureHitsTrueEnergy", 900, 600);
    fractionImpureHitsTrueEnergy->SetXTitle("Fraction of Impure Hits F (2nd MC parent contributes #geq 15% of charge)");
    fractionImpureHitsTrueEnergy->SetYTitle("True Particle Energy E_{true}");
    fractionImpureHitsTrueEnergy->SetTitle("Fraction of Impure Hits F with True Particle Energy E_{true}");
    fractionImpureHitsTrueEnergy->GetYaxis()->SetTitleOffset(1.2);
    fractionImpureHitsTrueEnergy->Draw("COLZ");
    c7->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/filtering/fractionImpureHitsTrueEnergy.png");
}
