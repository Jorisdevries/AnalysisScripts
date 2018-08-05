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
#include <algorithm>

//------------------------------------------------------------------------------------------------------------------------------------------

void normaliseHistogram(TH1F *histogram)
{
    Double_t scale1 = 1/histogram->Integral();
    histogram->Scale(scale1);
}

void fittingtest( void )
{
    TFile *f = TFile::Open("/usera/jjd49/new_LAr/CondorUtilities/saved_results/fitting/fitting.root");

    TTree *t1 = (TTree*)f->Get("Tree");
    int nentries = (Int_t)t1->GetEntries();

    TH1F *hForwardsParameterZero = new TH1F("hForwardsParameterZero","", 100, 0.0, 100.0);
    TH1F *hForwardsParameterOne = new TH1F("hForwardsParameterOne","", 100, 0.0, 5.0);
    TH1F *hForwardsParameterTwo = new TH1F("hForwardsParameterTwo","", 100, 0.9, 1.1);
    TH1F *hForwardsParameterThree = new TH1F("hForwardsParameterThree","", 100, 0.0, 2.0);

    TH1F *hBackwardsParameterZero = new TH1F("hBackwardsParameterZero","", 100, 0.0, 100.0);
    TH1F *hBackwardsParameterOne = new TH1F("hBackwardsParameterOne","", 100, 0.0, 5.0);
    TH1F *hBackwardsParameterTwo = new TH1F("hBackwardsParameterTwo","", 100, 0.9, 1.1);
    TH1F *hBackwardsParameterThree = new TH1F("hBackwardsParameterThree","", 100, 0.0, 2.0);

    TH1F *hFitStatusOne = new TH1F("hFitStatusOne","", 10, 0, 10);
    TH1F *hFitStatusTwo = new TH1F("hFitStatusTwo","", 10, 0, 10);
    TH2F *hFitEnergyWithCharge = new TH2F("hFitEnergyWithCharge","", 40, 0.0, 300.0, 40, 0.0, 300.0);

    TH2F *hForwardsFitEnergyWithTrueEnergy = new TH2F("hForwardsFitEnergyWithTrueEnergy","", 40, 0.0, 1400.0, 40, 0.0, 1400.0);
    TH2F *hBackwardsFitEnergyWithTrueEnergy = new TH2F("hBackwardsFitEnergyWithTrueEnergy","", 40, 0.0, 1.0, 40, 0.0, 1000.0);

    int     trueParticleDirection;
    float   trueParticleEnergy;
    int     numberHits;
    float   totalCharge;
    float   forwardsChiSquared;
    float   backwardsChiSquared;
    int     fitStatus1;
    int     fitStatus2;

    float forwardsParameter0;
    float forwardsParameter1;
    float forwardsParameter2;
    //float forwardsParameter3;
    float backwardsParameter0;
    float backwardsParameter1;
    float backwardsParameter2;
    //float backwardsParameter3;

    int fileIdentifier;
    int eventNumber;

    float forwardsEnergy;
    float backwardsEnergy;

    t1->SetBranchAddress("trueParticleDirection", &trueParticleDirection);
    t1->SetBranchAddress("trueParticleEnergy", &trueParticleEnergy);
    t1->SetBranchAddress("numberHits", &numberHits);
    t1->SetBranchAddress("totalCharge", &totalCharge);
    t1->SetBranchAddress("forwardsChiSquared", &forwardsChiSquared);
    t1->SetBranchAddress("backwardsChiSquared", &backwardsChiSquared);
    t1->SetBranchAddress("fitStatus1", &fitStatus1);
    t1->SetBranchAddress("fitStatus2", &fitStatus2);

    t1->SetBranchAddress("forwardsParameter0", &forwardsParameter0);
    t1->SetBranchAddress("forwardsParameter1", &forwardsParameter1);
    t1->SetBranchAddress("forwardsParameter2", &forwardsParameter2);
    //t1->SetBranchAddress("forwardsParameter3", &forwardsParameter3);
    t1->SetBranchAddress("backwardsParameter0", &backwardsParameter0);
    t1->SetBranchAddress("backwardsParameter1", &backwardsParameter1);
    t1->SetBranchAddress("backwardsParameter2", &backwardsParameter2);
    //t1->SetBranchAddress("backwardsParameter3", &backwardsParameter3);

    t1->SetBranchAddress("fileIdentifier", &fileIdentifier);
    t1->SetBranchAddress("eventNumber", &eventNumber);

    t1->SetBranchAddress("forwardsEnergy", &forwardsEnergy);
    t1->SetBranchAddress("backwardsEnergy", &backwardsEnergy);

    //--------------------------------------------------------------------------------------------------------------------------------------

    std::cout << "There are " << nentries << " entries in the tree" << std::endl;

    float chiSquaredCut(5.0);

    for (int i = 1; i < nentries; i++)
    {
        t1->GetEntry(i);

        if (numberHits == 0)
            continue;
        
        if (trueParticleDirection == 1 && (forwardsChiSquared - backwardsChiSquared)/numberHits <= 0.0)
        {
            hFitEnergyWithCharge->Fill(totalCharge, forwardsParameter1);
            hForwardsFitEnergyWithTrueEnergy->Fill(trueParticleEnergy, forwardsEnergy + 105.7);
            hFitStatusOne->Fill(fitStatus1);
        }
        if (trueParticleDirection == 0 && (forwardsChiSquared - backwardsChiSquared)/numberHits > 0.0)
            hFitStatusTwo->Fill(fitStatus2);

        if (trueParticleDirection == 1 && (forwardsChiSquared - backwardsChiSquared)/numberHits <= 0.0 && (forwardsChiSquared/numberHits) <= chiSquaredCut)
        {
            hForwardsParameterZero->Fill(forwardsParameter0);
            hForwardsParameterOne->Fill(forwardsParameter1);
            hForwardsParameterTwo->Fill(forwardsParameter2);
            //hForwardsParameterThree->Fill(forwardsParameter3);
        }
        if (trueParticleDirection == 0 && (forwardsChiSquared - backwardsChiSquared)/numberHits >= 0.0 && (backwardsChiSquared/numberHits) <= chiSquaredCut)
        {
            //if (backwardsParameter1 <= 0.1)
            //    std::cout << "Very low scaling parameter in event: " << fileIdentifier << ": " << eventNumber << std::endl;

            hBackwardsParameterZero->Fill(backwardsParameter0);
            hBackwardsParameterOne->Fill(backwardsParameter1);
            hBackwardsParameterTwo->Fill(backwardsParameter2);
            //hBackwardsParameterThree->Fill(backwardsParameter3);
            hFitEnergyWithCharge->Fill(totalCharge, backwardsParameter1);
            //hBackwardsFitEnergyWithTrueEnergy->Fill(trueParticleEnergy, backwardsEnergy);
        }
    }

    normaliseHistogram(hForwardsParameterZero);
    normaliseHistogram(hForwardsParameterOne);
    normaliseHistogram(hForwardsParameterTwo);
    normaliseHistogram(hForwardsParameterThree);

    normaliseHistogram(hBackwardsParameterZero);
    normaliseHistogram(hBackwardsParameterOne);
    normaliseHistogram(hBackwardsParameterTwo);
    normaliseHistogram(hBackwardsParameterThree);

    std::cout << "Forwards Muons: fraction correct forwards fits: " << (float)hFitStatusOne->GetBinContent(1)/(hFitStatusOne->GetBinContent(1) + hFitStatusOne->GetBinContent(5)) << std::endl;
    std::cout << "Backwards Muons: fraction correct backwards fits: " << (float)hFitStatusTwo->GetBinContent(1)/(hFitStatusTwo->GetBinContent(1) + hFitStatusTwo->GetBinContent(5)) << std::endl;

    //--------------------------------------------------------------------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas("hForwardsParameterZero", "hForwardsParameterZero", 900, 600);
    hForwardsParameterZero->SetXTitle("Fit Parameter [0] (E_{end})");
    hForwardsParameterZero->SetYTitle("Number of entries");
    hForwardsParameterZero->SetTitle("Forwards & Bsckwards Muons: Fit Parameter [0] (E_{end})");
    hForwardsParameterZero->GetYaxis()->SetTitleOffset(1.4);
    hForwardsParameterZero->GetYaxis()->SetRangeUser(0., 0.1);
    hForwardsParameterZero->SetLineColor(kBlue);
    hForwardsParameterZero->Draw();
    hBackwardsParameterZero->SetLineColor(kRed);
    hBackwardsParameterZero->Draw("same");
    c1->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hParameterZero.png");

    TCanvas *c2 = new TCanvas("hForwardsParameterOne", "hForwardsParameterOne", 900, 600);
    hForwardsParameterOne->SetXTitle("Fit Parameter [1] (2D #rightarrow 3D Scaling Factor)");
    hForwardsParameterOne->SetYTitle("Number of entries");
    hForwardsParameterOne->SetTitle("Forwards & Backwards Muons: Fit Parameter [1] (2D #rightarrow 3D Scaling Factor)");
    hForwardsParameterOne->GetYaxis()->SetTitleOffset(1.4);
    hForwardsParameterOne->Draw();
    hBackwardsParameterOne->SetLineColor(kRed);
    hBackwardsParameterOne->Draw("same");
    c2->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hParameterOne.png");

    TCanvas *c3 = new TCanvas("hForwardsParameterTwo", "hForwardsParameterTwo", 900, 600);
    hForwardsParameterTwo->SetXTitle("Fit Parameter [2] (Constant Term)");
    hForwardsParameterTwo->SetYTitle("Number of entries");
    hForwardsParameterTwo->SetTitle("Forwards & Backwards Muons: Fit Parameter [2] (Constant Term)");
    hForwardsParameterTwo->GetYaxis()->SetTitleOffset(1.0);
    hForwardsParameterTwo->Draw();
    hBackwardsParameterTwo->SetLineColor(kRed);
    hBackwardsParameterTwo->Draw("same");
    c3->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hParameterTwo.png");

    /*TCanvas *c4 = new TCanvas("hForwardsParameterThree", "hForwardsParameterThree", 900, 600);
    hForwardsParameterThree->SetXTitle("Forwards Muons: Fit Parameter [3] (Charge to Energy Scaling Factor)");
    hForwardsParameterThree->SetYTitle("Number of entries");
    hForwardsParameterThree->SetTitle("Forwards Muons: Fit Parameter [3] (Charge to Energy Scaling Factor)");
    hForwardsParameterThree->GetYaxis()->SetTitleOffset(1.0);
    hForwardsParameterThree->Draw();
    c4->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hForwardsParameterThree.png"); */


    TCanvas *c333 = new TCanvas("hFitStatusOne", "hFitStatusOne", 900, 600);
    hFitStatusOne->SetXTitle("Forwards Muons: Forwards Fit Output Status");
    hFitStatusOne->SetYTitle("Number of entries");
    hFitStatusOne->SetTitle("Forwards Muons: Forwards Fit Output Status (0: normal, 4: abnormal termination)");
    hFitStatusOne->GetYaxis()->SetTitleOffset(1.4);
    hFitStatusOne->Draw();
    c333->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hFitStatusOne.png");

    TCanvas *c444 = new TCanvas("hFitStatusTwo", "hFitStatusOne", 900, 600);
    hFitStatusTwo->SetXTitle("Forwards Muons: Backwards Fit Output Status");
    hFitStatusTwo->SetYTitle("Number of entries");
    hFitStatusTwo->SetTitle("Forwards Muons: Backwards Fit Output Status (0: normal, 4: abnormal termination)");
    hFitStatusTwo->GetYaxis()->SetTitleOffset(1.4);
    hFitStatusTwo->Draw();
    c444->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hFitStatusTwo.png");

    /*TCanvas *c5 = new TCanvas("hFitEnergyWithCharge", "hFitEnergyWithCharge", 900, 600);
    hFitEnergyWithCharge->SetXTitle("Total Charge Hit Collection");
    hFitEnergyWithCharge->SetYTitle("Fit Energy (Start Energy  [1] + Integral + Mass [2]) (MeV)");
    hFitEnergyWithCharge->SetTitle("Forwards Muons: Total Fit Energy & Total Charge Hit Collection");
    hFitEnergyWithCharge->GetYaxis()->SetTitleOffset(1.4);
    hFitEnergyWithCharge->Draw("COLZ");
    hFitEnergyWithCharge->SetStats(kFALSE);
    c5->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hFitEnergyWithCharge.png");*/

    TCanvas *c55 = new TCanvas("hForwardsFitEnergyWithTrueEnergy", "hForwardsFitEnergyWithTrueEnergy", 900, 600);
    hForwardsFitEnergyWithTrueEnergy->SetXTitle("E_{true} (MeV)");
    hForwardsFitEnergyWithTrueEnergy->SetYTitle("E_{fit} = Start Energy [0] + #DeltaE (from [1]) + Muon Mass (MeV)");
    hForwardsFitEnergyWithTrueEnergy->SetTitle("Reconstructed Muon Energy E_{fit} & True Muon Energy E_{true}");
    hForwardsFitEnergyWithTrueEnergy->GetYaxis()->SetTitleOffset(1.4);
    hForwardsFitEnergyWithTrueEnergy->Draw("COLZ");
    hForwardsFitEnergyWithTrueEnergy->SetStats(kFALSE);
    c55->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hForwardsFitEnergyWithTrueEnergy.png");

    /*TCanvas *c56 = new TCanvas("hBackwardsFitEnergyWithTrueEnergy", "hFitEnergyWithCharge", 900, 600);
    hBackwardsFitEnergyWithTrueEnergy->SetXTitle("True Muon Energy (GeV)");
    hBackwardsFitEnergyWithTrueEnergy->SetYTitle("Fit Energy (Start Energy  [1] + Integral + Mass [2]) (MeV)");
    hBackwardsFitEnergyWithTrueEnergy->SetTitle("Backwards Muons: Fitted Energy & True Muon Energy");
    hBackwardsFitEnergyWithTrueEnergy->GetYaxis()->SetTitleOffset(1.4);
    hBackwardsFitEnergyWithTrueEnergy->Draw("COLZ");
    hBackwardsFitEnergyWithTrueEnergy->SetStats(kFALSE);
    c56->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hBackwardsFitEnergyWithTrueEnergy.png");
    */

    /*TCanvas *c11 = new TCanvas("hBackwardsParameterZero", "hBackwardsParameterZero", 900, 600);
    hBackwardsParameterZero->SetXTitle("Backwards Muons: Fit Parameter [0] (2D to 3D scaling factor)");
    hBackwardsParameterZero->SetYTitle("Number of entries");
    hBackwardsParameterZero->SetTitle("Backwards Muons: Fit Parameter [0] (2D to 3D scaling factor)");
    hBackwardsParameterZero->GetYaxis()->SetTitleOffset(1.0);
    hBackwardsParameterZero->Draw();
    c11->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hBackwardsParameterZero.png");

    TCanvas *c22 = new TCanvas("hBackwardsParameterOne", "hBackwardsParameterOne", 900, 600);
    hBackwardsParameterOne->SetXTitle("Backwards Muons: Fit Parameter [1] (Particle Kinetic Energy)");
    hBackwardsParameterOne->SetYTitle("Number of entries");
    hBackwardsParameterOne->SetTitle("Backwards Muons: Fit Parameter [1] (Particle Kinetic Energy)");
    hBackwardsParameterOne->GetYaxis()->SetTitleOffset(1.0);
    hBackwardsParameterOne->Draw();
    c22->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hBackwardsParameterOne.png");

    TCanvas *c33 = new TCanvas("hBackwardsParameterTwo", "hBackwardsParameterTwo", 900, 600);
    hBackwardsParameterTwo->SetXTitle("Backwards Muons: Fit Parameter [2] (Density Correction)");
    hBackwardsParameterTwo->SetYTitle("Number of entries");
    hBackwardsParameterTwo->SetTitle("Backwards Muons: Fit Parameter [2] (Density Correction)");
    hBackwardsParameterTwo->GetYaxis()->SetTitleOffset(1.0);
    hBackwardsParameterTwo->Draw();
    c33->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hBackwardsParameterTwo.png");

    TCanvas *c44 = new TCanvas("hBackwardsParameterThree", "hBackwardsParameterThree", 900, 600);
    hBackwardsParameterThree->SetXTitle("Backwards Muons: Fit Parameter [3] (Charge to Energy Scaling Factor)");
    hBackwardsParameterThree->SetYTitle("Number of entries");
    hBackwardsParameterThree->SetTitle("Backwards Muons: Fit Parameter [3] (Charge to Energy Scaling Factor)");
    hBackwardsParameterThree->GetYaxis()->SetTitleOffset(1.0);
    hBackwardsParameterThree->Draw();
    c44->SaveAs("/usera/jjd49/new_LAr/CondorUtilities/figures/fitting/hBackwardsParameterThree.png");*/

}
