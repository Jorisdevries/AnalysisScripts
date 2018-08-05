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

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>

//------------------------------------------------------------------------------------------------------------------------------------------

void filtertest( void )
{
    TFile *f = TFile::Open("/usera/jjd49/new_LAr/CondorUtilities/saved_results/filtering/filtering_31Aug.root");
    
    TTree *t1 = (TTree*)f->Get("Tree");
    int nentries = (Int_t)t1->GetEntries();

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
        
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    std::cout << "There are " << nentries << " entries in the tree" << std::endl;
    
    for (int i = 1; i < nentries; i++)
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
