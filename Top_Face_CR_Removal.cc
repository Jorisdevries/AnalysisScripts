#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TMath.h"
#include "TColor.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TGaxis.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "Drawing_Helper_Functions.h"

//------------------------------------------------------------------------------------------------------------------------------------------

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/topface_cr_removal/topface_cr_removal.root");
TTree *t1 = (TTree*)f1->Get("TopFaceCosmicRemoval");

//TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/topface_cr_removal/neutrino_id_unmodified.root");
//TTree *t2 = (TTree*)f2->Get("NeutrinoIdCheck");

//TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/topface_cr_removal/neutrino_id_400_upwardsNu_removed.root");
//TTree *t3 = (TTree*)f3->Get("NeutrinoIdCheck");

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateCosmicRemovalEfficiencyCurves(TTree* pTree, int nDataPoints)
{
    float cosmicProbability, deltaChiSquaredUpDownPerHit, polarAngle, pfoHighY, fitEndpointEnergy, deltaChiSquaredUpDown;
    int mcTargetCosmic, targetPfo, intersectsTopFace, fiducialLowY, highTopY, fileIdentifier, eventNumber, mcCosmicRay, trueDownwards, highYZLargerThanLowYZ, sufficientThreeDHits, removedByRegularTagging, pfoPolarAngleIsSomewhatSmall, highYCut, polarAngleCut, pfoChargeCut, neutrinoMomentumZCut, pfoAzimuthalAngleCut, numberConnectedPfosNotNearEndpoint;

    pTree->SetBranchAddress("TrueDownwards", &trueDownwards); 
    pTree->SetBranchAddress("CosmicProbability", &cosmicProbability); 
    pTree->SetBranchAddress("MCTargetCosmic", &mcTargetCosmic); 
    pTree->SetBranchAddress("MCCosmicRay", &mcCosmicRay); 
    pTree->SetBranchAddress("TargetPFO", &targetPfo); 
    pTree->SetBranchAddress("PfoIntersectsTopFace", &intersectsTopFace); 
    pTree->SetBranchAddress("PfoFiducialLowY", &fiducialLowY); 
    pTree->SetBranchAddress("PfoHighTopY", &highTopY); 
    pTree->SetBranchAddress("PfoPolarAngle", &polarAngle); 
    pTree->SetBranchAddress("DeltaChiSquaredUpDownPerHit", &deltaChiSquaredUpDownPerHit); 
    pTree->SetBranchAddress("DeltaChiSquaredUpDown", &deltaChiSquaredUpDown); 
    pTree->SetBranchAddress("HighYZLargerThanLowYZ", &highYZLargerThanLowYZ); 
    pTree->SetBranchAddress("SufficientThreeDHits", &sufficientThreeDHits); 
    pTree->SetBranchAddress("RemovedByRegularTagging", &removedByRegularTagging); 
    pTree->SetBranchAddress("PfoHighY", &pfoHighY); 
    pTree->SetBranchAddress("FitEndpointEnergy", &fitEndpointEnergy); 
    pTree->SetBranchAddress("PfoPolarAngleIsSomewhatSmall", &pfoPolarAngleIsSomewhatSmall); 
    pTree->SetBranchAddress("HighYCut", &highYCut); 
    pTree->SetBranchAddress("PolarAngleCut", &polarAngleCut); 
    pTree->SetBranchAddress("PfoChargeCut", &pfoChargeCut); 
    pTree->SetBranchAddress("NeutrinoMomentumZCut", &neutrinoMomentumZCut); 
    pTree->SetBranchAddress("PfoAzimuthalAngleCut", &pfoAzimuthalAngleCut); 
    pTree->SetBranchAddress("NumberConnectedPfosNotNearEndpoint", &numberConnectedPfosNotNearEndpoint); 
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier); 
    pTree->SetBranchAddress("EventNumber", &eventNumber); 

    TGraph* pFractionTargetCosmicRaysRemovedGraph_CP = new TGraph(nDataPoints);
    pFractionTargetCosmicRaysRemovedGraph_CP->SetMarkerStyle(6);
    pFractionTargetCosmicRaysRemovedGraph_CP->SetMarkerColor(kMagenta);
    pFractionTargetCosmicRaysRemovedGraph_CP->SetLineColor(kMagenta);
    pFractionTargetCosmicRaysRemovedGraph_CP->SetLineStyle(9);

    TGraph* pFractionUpwardsNuPfosRetainedGraph_CP = new TGraph(nDataPoints);
    pFractionUpwardsNuPfosRetainedGraph_CP->SetMarkerStyle(6);
    pFractionUpwardsNuPfosRetainedGraph_CP->SetMarkerColor(kBlue);
    pFractionUpwardsNuPfosRetainedGraph_CP->SetLineColor(kBlue);
    pFractionUpwardsNuPfosRetainedGraph_CP->SetLineStyle(9);

    TGraph* pFractionDownwardsNuPfosRetainedGraph_CP = new TGraph(nDataPoints);
    pFractionDownwardsNuPfosRetainedGraph_CP->SetMarkerStyle(6);
    pFractionDownwardsNuPfosRetainedGraph_CP->SetMarkerColor(kRed);
    pFractionDownwardsNuPfosRetainedGraph_CP->SetLineColor(kRed);
    pFractionDownwardsNuPfosRetainedGraph_CP->SetLineStyle(9);

    TGraph* pFractionTargetCosmicRaysRemovedGraph_DCHS = new TGraph(nDataPoints);
    pFractionTargetCosmicRaysRemovedGraph_DCHS->SetMarkerStyle(6);
    pFractionTargetCosmicRaysRemovedGraph_DCHS->SetMarkerColor(kMagenta);
    pFractionTargetCosmicRaysRemovedGraph_DCHS->SetLineColor(kMagenta);
    //pFractionTargetCosmicRaysRemovedGraph_DCHS->SetLineStyle(9);

    TGraph* pFractionUpwardsNuPfosRetainedGraph_DCHS = new TGraph(nDataPoints);
    pFractionUpwardsNuPfosRetainedGraph_DCHS->SetMarkerStyle(6);
    pFractionUpwardsNuPfosRetainedGraph_DCHS->SetMarkerColor(kBlue);
    pFractionUpwardsNuPfosRetainedGraph_DCHS->SetLineColor(kBlue);

    TGraph* pFractionDownwardsNuPfosRetainedGraph_DCHS = new TGraph(nDataPoints);
    pFractionDownwardsNuPfosRetainedGraph_DCHS->SetMarkerStyle(6);
    pFractionDownwardsNuPfosRetainedGraph_DCHS->SetMarkerColor(kRed);
    pFractionDownwardsNuPfosRetainedGraph_DCHS->SetLineColor(kRed);
    //pFractionDownwardsNuPfosRetainedGraph_DCHS->SetLineStyle(10);


    TGraph* pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN = new TGraph(nDataPoints);
    pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN->SetMarkerStyle(6);
    pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN->SetMarkerColor(kMagenta);
    pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN->SetLineColor(kMagenta);
    //pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN->SetLineStyle(9);

    TGraph* pFractionUpwardsNuPfosRetainedGraph_DCHS_NoN = new TGraph(nDataPoints);
    pFractionUpwardsNuPfosRetainedGraph_DCHS_NoN->SetMarkerStyle(6);
    pFractionUpwardsNuPfosRetainedGraph_DCHS_NoN->SetMarkerColor(kBlue);
    pFractionUpwardsNuPfosRetainedGraph_DCHS_NoN->SetLineColor(kBlue);

    TGraph* pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN = new TGraph(nDataPoints);
    pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN->SetMarkerStyle(6);
    pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN->SetMarkerColor(kRed);
    pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN->SetLineColor(kRed);
    //pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN->SetLineStyle(10);

    for (int i = 0; i <= nDataPoints; ++i)
    {
        float cosmicProbabilityCut(0.5 + i * (0.5/nDataPoints)), deltaChiSquaredCut(i * (-5.0/nDataPoints)), deltaChiSquaredCut_NoN(i * (-500.0/nDataPoints));

        int nTotalCRs(0), nTotalUpwardsNeutrinoPfos(0), nTotalDownwardsNeutrinoPfos(0); 
        int nTotalCRs_HitCut(0), nTotalUpwardsNeutrinoPfos_HitCut(0), nTotalDownwardsNeutrinoPfos_HitCut(0); 
        int nTotalCRs_NotRemovedByNormalTagging(0), nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging(0), nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging(0); 
        int nTotalCRs_NotRemovedByNormalTagging_HighYCut(0), nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut(0), nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut(0); 
        int nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut(0), nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut(0), nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut(0); 
        int nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut(0), nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut(0), nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut(0); 
        int nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut(0), nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut(0), nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut(0); 
        int nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF(0), nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF(0), nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF(0); 

        int nTotalTargetCRs(0), nTotalTargetUpwardsNeutrinoPfos(0), nTotalTargetDownwardsNeutrinoPfos(0);
        int nRemovedPfos_CP(0), nRemovedPfos_DCHS(0), nRemovedPfos_DCHS_NoN(0);
        int nCosmicsRemoved_CP(0), nCosmicsRemoved_DCHS(0), nCosmicsRemoved_DCHS_NoN(0);
        int nUpwardsNeutrinoPfosRemoved_CP(0), nUpwardsNeutrinoPfosRemoved_DCHS(0), nUpwardsNeutrinoPfosRemoved_DCHS_NoN(0);
        int nDownwardsNeutrinoPfosRemoved_CP(0), nDownwardsNeutrinoPfosRemoved_DCHS(0), nDownwardsNeutrinoPfosRemoved_DCHS_NoN(0);

        for (int i = 1; i < t1->GetEntries(); i++)
        {    
            t1->GetEntry(i);

            if (mcCosmicRay == 1)
                ++nTotalCRs;
            if (mcCosmicRay == 0 && trueDownwards == 0) 
                ++nTotalUpwardsNeutrinoPfos;
            if (mcCosmicRay == 0 && trueDownwards == 1) 
                ++nTotalDownwardsNeutrinoPfos;

            if (sufficientThreeDHits == 0)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs_HitCut;
            if (mcCosmicRay == 0 && trueDownwards == 0)
                ++nTotalUpwardsNeutrinoPfos_HitCut;
            if (mcCosmicRay == 0 && trueDownwards == 1)
                ++nTotalDownwardsNeutrinoPfos_HitCut;

            if (removedByRegularTagging == 1)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs_NotRemovedByNormalTagging;
            if (mcCosmicRay == 0 && trueDownwards == 0)
                ++nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging;
            if (mcCosmicRay == 0 && trueDownwards == 1)
                ++nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging;

            //continue if pfo is not of target class or if direction fit fails
            if (highYCut == 1)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs_NotRemovedByNormalTagging_HighYCut;
            if (mcCosmicRay == 0 && trueDownwards == 0)
                ++nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut;
            if (mcCosmicRay == 0 && trueDownwards == 1)
                ++nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut;

            if (polarAngleCut == 1)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut;
            if (mcCosmicRay == 0 && trueDownwards == 0)
                ++nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut;
            if (mcCosmicRay == 0 && trueDownwards == 1)
                ++nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut;

            if (pfoChargeCut == 1)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut;
            if (mcCosmicRay == 0 && trueDownwards == 0)
                ++nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut;
            if (mcCosmicRay == 0 && trueDownwards == 1)
                ++nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut;

            if (neutrinoMomentumZCut == 1)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut;
            if (mcCosmicRay == 0 && trueDownwards == 0)
                ++nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut;
            if (mcCosmicRay == 0 && trueDownwards == 1)
                ++nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut;

            if (numberConnectedPfosNotNearEndpoint == 0)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF;
            if (mcCosmicRay == 0 && trueDownwards == 0)
                ++nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF;
            if (mcCosmicRay == 0 && trueDownwards == 1)
                ++nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF;

            //if (fitEndpointEnergy < 45.0 || highYZLargerThanLowYZ == 0)
            //    continue;

            if (mcCosmicRay == 1)
                ++nTotalTargetCRs;
            if (mcCosmicRay == 0 && trueDownwards == 0) 
                ++nTotalTargetUpwardsNeutrinoPfos;
            if (mcCosmicRay == 0 && trueDownwards == 1) 
                ++nTotalTargetDownwardsNeutrinoPfos;

            if (cosmicProbability >= cosmicProbabilityCut)
            {
                ++nRemovedPfos_CP;
            
                if (mcCosmicRay == 1)
                    ++nCosmicsRemoved_CP;
                if (mcCosmicRay == 0 && trueDownwards == 0) 
                    ++nUpwardsNeutrinoPfosRemoved_CP;
                if (mcCosmicRay == 0 && trueDownwards == 1) 
                    ++nDownwardsNeutrinoPfosRemoved_CP;
            }

            if (deltaChiSquaredUpDownPerHit <= deltaChiSquaredCut)
            {
                ++nRemovedPfos_DCHS;

                if (mcCosmicRay == 1)
                    ++nCosmicsRemoved_DCHS;
                if (mcCosmicRay == 0 && trueDownwards == 0) 
                    ++nUpwardsNeutrinoPfosRemoved_DCHS;
                if (mcCosmicRay == 0 && trueDownwards == 1) 
                    ++nDownwardsNeutrinoPfosRemoved_DCHS;
            }

            if (deltaChiSquaredUpDown <= deltaChiSquaredCut_NoN)
            {
                ++nRemovedPfos_DCHS_NoN;

                if (mcCosmicRay == 1)
                    ++nCosmicsRemoved_DCHS_NoN;
                if (mcCosmicRay == 0 && trueDownwards == 0) 
                    ++nUpwardsNeutrinoPfosRemoved_DCHS_NoN;
                if (mcCosmicRay == 0 && trueDownwards == 1) 
                    ++nDownwardsNeutrinoPfosRemoved_DCHS_NoN;
            }
        }

        if (i == 0)
        {
            std::cout << "****************************************" << std::endl;
            std::cout << "nTotalCRs: " << nTotalCRs << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos: " << nTotalUpwardsNeutrinoPfos << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos: " << nTotalDownwardsNeutrinoPfos << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalCRs_HitCut: " << nTotalCRs_HitCut << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos_HitCut: " << nTotalUpwardsNeutrinoPfos_HitCut << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos_HitCut: " << nTotalDownwardsNeutrinoPfos_HitCut << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalCRs_NotRemovedByNormalTagging: " << nTotalCRs_NotRemovedByNormalTagging << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging: " << nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging: " << nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalCRs_NotRemovedByNormalTagging_HighYCut: " << nTotalCRs_NotRemovedByNormalTagging_HighYCut << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut: " << nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut: " << nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut: " << nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut: " << nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut: " << nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut: " << nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut: " << nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut: " << nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut: " << nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut: " << nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut: " << nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF: " << nTotalCRs_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF: " << nTotalUpwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF: " << nTotalDownwardsNeutrinoPfos_NotRemovedByNormalTagging_HighYCut_PolarAngleCut_ChargeCut_MomentumCut_BPF << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "nTotalTargetCRs: " << nTotalTargetCRs << std::endl;
            std::cout << "nTotalTargetUpwardsNeutrinoPfos: " << nTotalTargetUpwardsNeutrinoPfos << std::endl;
            std::cout << "nTotalTargetDownwardsNeutrinoPfos: " << nTotalTargetDownwardsNeutrinoPfos << std::endl;
            std::cout << "****************************************" << std::endl;
        }

        std::cout << "cosmicProbabilityCut: " << cosmicProbabilityCut << std::endl;
        std::cout << "Additional CR removed: " << nCosmicsRemoved_CP << " (" << static_cast<float>(nCosmicsRemoved_CP)/nTotalTargetCRs << "%)" << std::endl;
        std::cout << "Upwards neutrinos retained: " << " (" << static_cast<float>(nTotalTargetUpwardsNeutrinoPfos - nUpwardsNeutrinoPfosRemoved_CP)/nTotalTargetUpwardsNeutrinoPfos << "%)" << std::endl;
        std::cout << "Downwards neutrinos retained: " << " (" << static_cast<float>(nTotalTargetDownwardsNeutrinoPfos - nDownwardsNeutrinoPfosRemoved_CP)/nTotalTargetDownwardsNeutrinoPfos << "%)" << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        std::cout << "deltaChiSquaredCut: " << deltaChiSquaredCut << std::endl;
        std::cout << "Additional CR removed: " << nCosmicsRemoved_DCHS << " (" << static_cast<float>(nCosmicsRemoved_DCHS)/nTotalTargetCRs << "%)" << std::endl;
        std::cout << "Upwards neutrinos retained: " << (nTotalTargetUpwardsNeutrinoPfos - nUpwardsNeutrinoPfosRemoved_DCHS) << " (" << static_cast<float>(nTotalTargetUpwardsNeutrinoPfos - nUpwardsNeutrinoPfosRemoved_DCHS)/nTotalTargetUpwardsNeutrinoPfos << "%)" << std::endl;
        std::cout << "Downwards neutrinos retained: " << (nTotalTargetDownwardsNeutrinoPfos - nDownwardsNeutrinoPfosRemoved_DCHS) << " (" << static_cast<float>(nTotalTargetDownwardsNeutrinoPfos - nDownwardsNeutrinoPfosRemoved_DCHS)/nTotalTargetDownwardsNeutrinoPfos << "%)" << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        float fractionTargetCRRemoved_CP(static_cast<float>(nCosmicsRemoved_CP)/nTotalTargetCRs);
        float fractionUpwardsNuPfosRetained_CP(1.0 - static_cast<float>(nUpwardsNeutrinoPfosRemoved_CP)/nTotalTargetUpwardsNeutrinoPfos);
        float fractionDownwardsNuPfosRetained_CP(1.0 - static_cast<float>(nDownwardsNeutrinoPfosRemoved_CP)/nTotalTargetDownwardsNeutrinoPfos);

        float fractionTargetCRRemoved_DCHS(static_cast<float>(nCosmicsRemoved_DCHS)/nTotalTargetCRs);
        float fractionUpwardsNuPfosRetained_DCHS(1.0 - static_cast<float>(nUpwardsNeutrinoPfosRemoved_DCHS)/nTotalTargetUpwardsNeutrinoPfos);
        float fractionDownwardsNuPfosRetained_DCHS(1.0 - static_cast<float>(nDownwardsNeutrinoPfosRemoved_DCHS)/nTotalTargetDownwardsNeutrinoPfos);

        float fractionTargetCRRemoved_DCHS_NoN(static_cast<float>(nCosmicsRemoved_DCHS_NoN)/nTotalTargetCRs);
        float fractionUpwardsNuPfosRetained_DCHS_NoN(1.0 - static_cast<float>(nUpwardsNeutrinoPfosRemoved_DCHS_NoN)/nTotalTargetUpwardsNeutrinoPfos);
        float fractionDownwardsNuPfosRetained_DCHS_NoN(1.0 - static_cast<float>(nDownwardsNeutrinoPfosRemoved_DCHS_NoN)/nTotalTargetDownwardsNeutrinoPfos);

        pFractionTargetCosmicRaysRemovedGraph_CP->SetPoint(i, cosmicProbabilityCut, fractionTargetCRRemoved_CP); 
        pFractionUpwardsNuPfosRetainedGraph_CP->SetPoint(i, cosmicProbabilityCut, fractionUpwardsNuPfosRetained_CP); 
        pFractionDownwardsNuPfosRetainedGraph_CP->SetPoint(i, cosmicProbabilityCut, fractionDownwardsNuPfosRetained_CP); 

        pFractionTargetCosmicRaysRemovedGraph_DCHS->SetPoint(i, cosmicProbabilityCut, fractionTargetCRRemoved_DCHS); 
        pFractionUpwardsNuPfosRetainedGraph_DCHS->SetPoint(i, cosmicProbabilityCut, fractionUpwardsNuPfosRetained_DCHS); 
        pFractionDownwardsNuPfosRetainedGraph_DCHS->SetPoint(i, cosmicProbabilityCut, fractionDownwardsNuPfosRetained_DCHS); 

        pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN->SetPoint(i, cosmicProbabilityCut, fractionTargetCRRemoved_DCHS_NoN); 
        pFractionUpwardsNuPfosRetainedGraph_DCHS_NoN->SetPoint(i, cosmicProbabilityCut, fractionUpwardsNuPfosRetained_DCHS_NoN); 
        pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN->SetPoint(i, cosmicProbabilityCut, fractionDownwardsNuPfosRetained_DCHS_NoN); 

    }

    /*
    TCanvas *c1 = new TCanvas("c1","c1",10,10,700,500);
    pCompletenessDeltaChiSquaredGraph->GetXaxis()->SetLimits(-5.0, 0.0);
    pCompletenessDeltaChiSquaredGraph->GetYaxis()->SetLimits(0.0, 1.0);
    pCompletenessDeltaChiSquaredGraph->Draw("ALP");
    pPurityDeltaChiSquaredGraph->Draw("LPsame");
    */

    TCanvas *c2 = new TCanvas("c2","c2",10,10,1200,800);
    TH1 *frame = c2->DrawFrame(0.5,0,1,1);
    frame->GetXaxis()->SetTitle("P_{cosmic} Cut");
    frame->GetYaxis()->SetTitle("Fraction Target CR Pfos Removed");   
    c2->Update();

    pFractionTargetCosmicRaysRemovedGraph_CP->GetXaxis()->SetLimits(0.5, 1.0);
    pFractionTargetCosmicRaysRemovedGraph_CP->GetYaxis()->SetLimits(0.0, 1.0);

    pFractionTargetCosmicRaysRemovedGraph_CP->Draw("LP");
    pFractionUpwardsNuPfosRetainedGraph_CP->Draw("LPsame");
    pFractionDownwardsNuPfosRetainedGraph_CP->Draw("LPsame");

    pFractionTargetCosmicRaysRemovedGraph_DCHS->Draw("LPsame");
    pFractionUpwardsNuPfosRetainedGraph_DCHS->Draw("LPsame");
    pFractionDownwardsNuPfosRetainedGraph_DCHS->Draw("LPsame");

    //pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN->Draw("LPsame");
    //pFractionUpwardsNuPfosRetainedGraph_DCHS_NoN->Draw("LPsame");
    //pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN->Draw("LPsame");

    TGaxis *A0 = new TGaxis(1.0,1,0.5,1,-500.0, 0.0, 510,"+");
    A0->SetTitle("#Delta#chi^{2}_{DU} Cut");
    A0->SetLabelOffset(0.05);
    A0->SetLabelSize(0.04);
    A0->SetTitleSize(0.04);
    //A0->Draw("same");

    TGaxis *A1 = new TGaxis(1.0,1,0.5,1,-5.0, 0.0, 510,"+");
    A1->SetTitle("#Delta#chi^{2}_{DU}/N Cut");
    A1->SetLabelOffset(0.05);
    A1->SetLabelSize(0.04);
    A1->SetTitleSize(0.04);
    A1->Draw("same");

    TGaxis *A2 = new TGaxis(1.0,0,1.0,1, 0.0, 1.0, 510,"+");
    A2->SetTitle("Fraction Target Neutrino Pfos Retained");
    A2->SetLabelOffset(0.05);
    A2->SetLabelSize(0.04);
    A2->SetTitleSize(0.04);
    A2->Draw("same");

    auto legend = new TLegend(0.45,0.4,0.875,0.6);
    legend->SetHeader("Event Type"); 
    legend->AddEntry(pFractionTargetCosmicRaysRemovedGraph_CP, "Fraction of Target CR Removed (P_{c} cut)","l");
    legend->AddEntry(pFractionUpwardsNuPfosRetainedGraph_CP, "Fraction of Target True Upwards #nu Pfos Retained (P_{c} cut)","l");
    legend->AddEntry(pFractionDownwardsNuPfosRetainedGraph_CP, "Fraction of Target True Downwards #nu Pfos Retained (P_{c} cut)","l");
    legend->AddEntry(pFractionTargetCosmicRaysRemovedGraph_DCHS, "Fraction of Target CR Removed (#Delta#chi^{2}/N cut)","l");
    legend->AddEntry(pFractionUpwardsNuPfosRetainedGraph_DCHS, "Fraction of Target True Upwards #nu Pfos Retained (#Delta#chi^{2}/N cut)","l");
    legend->AddEntry(pFractionDownwardsNuPfosRetainedGraph_DCHS, "Fraction of Target True Downwards #nu Pfos Retained (#Delta#chi^{2}/N cut)","l");
    //legend->AddEntry(pFractionTargetCosmicRaysRemovedGraph_DCHS_NoN, "Fraction of Target CR Removed (#Delta#chi^{2} cut)","l");
    //legend->AddEntry(pFractionUpwardsNuPfosRetainedGraph_DCHS_NoN, "Fraction of Target True Upwards #nu Pfos Retained (#Delta#chi^{2} cut)","l");
    //legend->AddEntry(pFractionDownwardsNuPfosRetainedGraph_DCHS_NoN, "Fraction of Target True Downwards #nu Pfos Retained (#Delta#chi^{2} cut)","l");
    legend->Draw("same");

    c2->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/TopFaceCRRemoval/CosmicRemoval_FractionTargetCRRemoved_FractionUpwardsNuPfosRetained_CP_vs_DCHS_small.pdf");
}

//------------------------------------------------------------------------------------------------------------------------------------------

int CreateCosmicRemovalEfficiencyCurves(TTree* pTree, std::vector<std::pair<std::string, int>> filterValues, int nDataPoints, int maximalUpwardsNeutrinoLoss)
{
    float cosmicProbability, deltaChiSquaredUpDownPerHit, fitEndpointEnergy;
    int fileIdentifier, eventNumber, mcCosmicRay, trueDownwards, sufficientThreeDHits, removedByRegularTagging, highYCut, polarAngleCut, pfoChargeCut, neutrinoMomentumZCut, numberConnectedPfosNotNearEndpoint, numberConnectedPfosNearEndpoint;

    pTree->SetBranchAddress("TrueDownwards", &trueDownwards); 
    pTree->SetBranchAddress("CosmicProbability", &cosmicProbability); 
    pTree->SetBranchAddress("MCCosmicRay", &mcCosmicRay); 
    pTree->SetBranchAddress("DeltaChiSquaredUpDownPerHit", &deltaChiSquaredUpDownPerHit); 
    pTree->SetBranchAddress("SufficientThreeDHits", &sufficientThreeDHits); 
    pTree->SetBranchAddress("RemovedByRegularTagging", &removedByRegularTagging); 
    pTree->SetBranchAddress("FitEndpointEnergy", &fitEndpointEnergy); 
    pTree->SetBranchAddress("HighYCut", &highYCut); 
    pTree->SetBranchAddress("PolarAngleCut", &polarAngleCut); 
    pTree->SetBranchAddress("PfoChargeCut", &pfoChargeCut); 
    pTree->SetBranchAddress("NeutrinoMomentumZCut", &neutrinoMomentumZCut); 
    pTree->SetBranchAddress("NumberConnectedPfosNotNearEndpoint", &numberConnectedPfosNotNearEndpoint); 
    pTree->SetBranchAddress("NumberConnectedPfosNearEndpoint", &numberConnectedPfosNearEndpoint); 
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier); 
    pTree->SetBranchAddress("EventNumber", &eventNumber); 

    int filterVariables[99], filterVariableTargetValues[99], filterVariableCounter(0);

    for (int i = 0; i < static_cast<int>(filterValues.size()); ++i)
    {   
        std::cout << "FILTER: " << filterValues.at(i).first << std::endl;
        pTree->SetBranchAddress(filterValues.at(i).first.c_str(), &filterVariables[i]);
        filterVariableTargetValues[i] = filterValues.at(i).second;
        ++filterVariableCounter;
    }   

    TGraph* pFractionCosmicRaysRemovedGraph_DCHS = new TGraph(nDataPoints);
    pFractionCosmicRaysRemovedGraph_DCHS->SetMarkerStyle(6);
    pFractionCosmicRaysRemovedGraph_DCHS->SetMarkerColor(kMagenta);
    pFractionCosmicRaysRemovedGraph_DCHS->SetLineColor(kMagenta);
    //pFractionCosmicRaysRemovedGraph_DCHS->SetLineStyle(9);

    TGraph* pFractionUpwardsNuPfosRetainedGraph_DCHS = new TGraph(nDataPoints);
    pFractionUpwardsNuPfosRetainedGraph_DCHS->SetMarkerStyle(6);
    pFractionUpwardsNuPfosRetainedGraph_DCHS->SetMarkerColor(kBlue);
    pFractionUpwardsNuPfosRetainedGraph_DCHS->SetLineColor(kBlue);

    TGraph* pFractionDownwardsNuPfosRetainedGraph_DCHS = new TGraph(nDataPoints);
    pFractionDownwardsNuPfosRetainedGraph_DCHS->SetMarkerStyle(6);
    pFractionDownwardsNuPfosRetainedGraph_DCHS->SetMarkerColor(kRed);
    pFractionDownwardsNuPfosRetainedGraph_DCHS->SetLineColor(kRed);
    //pFractionDownwardsNuPfosRetainedGraph_DCHS->SetLineStyle(10);

    bool cutNotYetFound(true);

    for (int i = 0; i <= nDataPoints; ++i)
    {
        float deltaChiSquaredCut(0.0 - (i * (5.0/nDataPoints)));
        int nTotalCRs(0), nTotalUpwardsNeutrinoPfos(0), nTotalDownwardsNeutrinoPfos(0);
        int nRemovedPfos_DCHS(0), nCosmicsRemoved_DCHS(0), nUpwardsNeutrinoPfosRemoved_DCHS(0), nDownwardsNeutrinoPfosRemoved_DCHS(0);

        for (int i = 1; i < t1->GetEntries(); i++)
        {    
            t1->GetEntry(i);

            bool skipThisEvent(false);
       
            for (int i = 0; i < static_cast<int>(filterValues.size()); ++i)
            {   
                if (filterVariables[i] != filterVariableTargetValues[i])
                {   
                    skipThisEvent = true;
                    break;
                }   
            }   

            if (skipThisEvent)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs;
            if (mcCosmicRay == 0 && trueDownwards == 0) 
                ++nTotalUpwardsNeutrinoPfos;
            if (mcCosmicRay == 0 && trueDownwards == 1) 
                ++nTotalDownwardsNeutrinoPfos;

            if (deltaChiSquaredUpDownPerHit <= deltaChiSquaredCut)
            {
                ++nRemovedPfos_DCHS;

                if (mcCosmicRay == 1)
                    ++nCosmicsRemoved_DCHS;
                if (mcCosmicRay == 0 && trueDownwards == 0) 
                    ++nUpwardsNeutrinoPfosRemoved_DCHS;
                if (mcCosmicRay == 0 && trueDownwards == 1) 
                    ++nDownwardsNeutrinoPfosRemoved_DCHS;
            }
        }

        if (i == 0)
        {
            std::cout << "****************************************" << std::endl;
            std::cout << "nTotalCRs: " << nTotalCRs << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos: " << nTotalUpwardsNeutrinoPfos << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos: " << nTotalDownwardsNeutrinoPfos << std::endl;
            std::cout << "****************************************" << std::endl;
        }

        /*
        std::cout << "deltaChiSquaredCut: " << deltaChiSquaredCut << std::endl;
        std::cout << "Additional CR removed: " << nCosmicsRemoved_DCHS << " (" << static_cast<float>(nCosmicsRemoved_DCHS)/nTotalCRs << "%)" << std::endl;
        std::cout << "Upwards neutrinos retained: " << (nTotalUpwardsNeutrinoPfos - nUpwardsNeutrinoPfosRemoved_DCHS) << " (" << static_cast<float>(nTotalUpwardsNeutrinoPfos - nUpwardsNeutrinoPfosRemoved_DCHS)/nTotalUpwardsNeutrinoPfos << "%)" << std::endl;
        std::cout << "Downwards neutrinos retained: " << (nTotalDownwardsNeutrinoPfos - nDownwardsNeutrinoPfosRemoved_DCHS) << " (" << static_cast<float>(nTotalDownwardsNeutrinoPfos - nDownwardsNeutrinoPfosRemoved_DCHS)/nTotalDownwardsNeutrinoPfos << "%)" << std::endl;
        */

        if ((cutNotYetFound && nUpwardsNeutrinoPfosRemoved_DCHS <= maximalUpwardsNeutrinoLoss) || deltaChiSquaredCut == -5.0)
        {
            cutNotYetFound = false;

            if (deltaChiSquaredCut == -5.0)
                std::cout << "Failed to find cut. Returning values for deltaChiSquaredCut == -5.0." << std::endl;

            std::cout << "****************************************" << std::endl;
            std::cout << "First cut that retains target upwards neutrinos is at: " << deltaChiSquaredCut << std::endl;
            std::cout << "Number of additional cosmics removed: " << nCosmicsRemoved_DCHS << std::endl;
            std::cout << "Number of upwards neutrinos removed: " << nUpwardsNeutrinoPfosRemoved_DCHS << " (" << static_cast<float>(nUpwardsNeutrinoPfosRemoved_DCHS)/maximalUpwardsNeutrinoLoss << "% of maximum, maximum was " << maximalUpwardsNeutrinoLoss << ")" << std::endl;
            std::cout << "Number of downwards neutrinos removed: " << nDownwardsNeutrinoPfosRemoved_DCHS << std::endl;
            std::cout << "****************************************" << std::endl;

            return nCosmicsRemoved_DCHS;
        }

        //std::cout << "----------------------------------------" << std::endl;

        float fractionCRRemoved_DCHS(static_cast<float>(nCosmicsRemoved_DCHS)/nTotalCRs);
        float fractionUpwardsNuPfosRetained_DCHS(1.0 - static_cast<float>(nUpwardsNeutrinoPfosRemoved_DCHS)/nTotalUpwardsNeutrinoPfos);
        float fractionDownwardsNuPfosRetained_DCHS(1.0 - static_cast<float>(nDownwardsNeutrinoPfosRemoved_DCHS)/nTotalDownwardsNeutrinoPfos);

        pFractionCosmicRaysRemovedGraph_DCHS->SetPoint(i, -deltaChiSquaredCut, fractionCRRemoved_DCHS); 
        pFractionUpwardsNuPfosRetainedGraph_DCHS->SetPoint(i, -deltaChiSquaredCut, fractionUpwardsNuPfosRetained_DCHS); 
        pFractionDownwardsNuPfosRetainedGraph_DCHS->SetPoint(i, -deltaChiSquaredCut, fractionDownwardsNuPfosRetained_DCHS); 

    }

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,800);
    pFractionCosmicRaysRemovedGraph_DCHS->Draw("ALP");
    pFractionCosmicRaysRemovedGraph_DCHS->GetXaxis()->SetLimits(0.0,5.0);
    pFractionCosmicRaysRemovedGraph_DCHS->GetYaxis()->SetRangeUser(0.0,1.0);

    pFractionCosmicRaysRemovedGraph_DCHS->SetTitle("#Delta#chi^{2}_{DU}/N Cut Performance;;Fraction Untagged CR Removed");

    pFractionCosmicRaysRemovedGraph_DCHS->GetXaxis()->SetLabelOffset(999);
    pFractionCosmicRaysRemovedGraph_DCHS->GetXaxis()->SetLabelSize(0);
    pFractionCosmicRaysRemovedGraph_DCHS->GetXaxis()->SetTickLength(0.0);

    TGaxis *A1 = new TGaxis(5.0,0,0.0,0, -5.0, 0.0, 510,"-");
    A1->SetTitle("#Delta#chi^{2}_{DU}/N Cut");
    A1->SetLabelOffset(-0.04);
    A1->SetLabelSize(0.04);
    A1->SetTitleSize(0.04);
    A1->Draw("same");

    TGaxis *A2 = new TGaxis(5.0,0,5.0,1, 0.0, 1.0, 510,"+");
    A2->SetTitle("Fraction Target Neutrino Pfos Retained");
    A2->SetLabelOffset(0.05);
    A2->SetLabelSize(0.04);
    A2->SetTitleSize(0.04);
    A2->Draw("same");

    pFractionUpwardsNuPfosRetainedGraph_DCHS->Draw("LPsame");
    pFractionDownwardsNuPfosRetainedGraph_DCHS->Draw("LPsame");

    auto legend = new TLegend(0.45,0.4,0.875,0.6);
    legend->SetHeader("Event Type"); 
    legend->AddEntry(pFractionCosmicRaysRemovedGraph_DCHS, "Fraction of  CR Removed (#Delta#chi^{2}/N cut)","l");
    legend->AddEntry(pFractionUpwardsNuPfosRetainedGraph_DCHS, "Fraction of  True Upwards #nu Pfos Retained (#Delta#chi^{2}/N cut)","l");
    legend->AddEntry(pFractionDownwardsNuPfosRetainedGraph_DCHS, "Fraction of  True Downwards #nu Pfos Retained (#Delta#chi^{2}/N cut)","l");
    legend->Draw("same");

    std::string plotName("CosmicRemoval_FractionCRRemoved_FractionUpwardsNuPfosRetained_");

    for (const auto &entry : filterValues)
    {
        plotName += entry.first;
        plotName += "_";
    }

    std::string filePath("/usera/jjd49/pandora_direction/Scripts/Figures/TopFaceCRRemoval/" + plotName + "_DCHS.pdf");
    c1->SaveAs(filePath.c_str());

    return 9999;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateCosmicRemovalROCCurve(TTree* pTree, std::vector<std::pair<std::string, int>> filterValues, int nDataPoints)
{
    float cosmicProbability, deltaChiSquaredUpDownPerHit, fitEndpointEnergy;
    int fileIdentifier, eventNumber, mcCosmicRay, trueDownwards, sufficientThreeDHits, removedByRegularTagging, highYCut, polarAngleCut, pfoChargeCut, neutrinoMomentumZCut, numberConnectedPfosNotNearEndpoint, numberConnectedPfosNearEndpoint;

    pTree->SetBranchAddress("TrueDownwards", &trueDownwards); 
    pTree->SetBranchAddress("CosmicProbability", &cosmicProbability); 
    pTree->SetBranchAddress("MCCosmicRay", &mcCosmicRay); 
    pTree->SetBranchAddress("DeltaChiSquaredUpDownPerHit", &deltaChiSquaredUpDownPerHit); 
    pTree->SetBranchAddress("SufficientThreeDHits", &sufficientThreeDHits); 
    pTree->SetBranchAddress("RemovedByRegularTagging", &removedByRegularTagging); 
    pTree->SetBranchAddress("FitEndpointEnergy", &fitEndpointEnergy); 
    pTree->SetBranchAddress("HighYCut", &highYCut); 
    pTree->SetBranchAddress("PolarAngleCut", &polarAngleCut); 
    pTree->SetBranchAddress("PfoChargeCut", &pfoChargeCut); 
    pTree->SetBranchAddress("NeutrinoMomentumZCut", &neutrinoMomentumZCut); 
    pTree->SetBranchAddress("NumberConnectedPfosNotNearEndpoint", &numberConnectedPfosNotNearEndpoint); 
    pTree->SetBranchAddress("NumberConnectedPfosNearEndpoint", &numberConnectedPfosNearEndpoint); 
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier); 
    pTree->SetBranchAddress("EventNumber", &eventNumber); 

    int filterVariables[99], filterVariableTargetValues[99], filterVariableCounter(0);

    for (int i = 0; i < static_cast<int>(filterValues.size()); ++i)
    {   
        std::cout << "FILTER: " << filterValues.at(i).first << std::endl;
        pTree->SetBranchAddress(filterValues.at(i).first.c_str(), &filterVariables[i]);
        filterVariableTargetValues[i] = filterValues.at(i).second;
        ++filterVariableCounter;
    }   

    TGraph* pROCCurve = new TGraph(nDataPoints);
    pROCCurve->SetMarkerStyle(6);
    pROCCurve->SetMarkerColor(kMagenta);
    pROCCurve->SetLineColor(kMagenta);
    //pROCCurve->SetLineStyle(9);

    bool cutNotYetFound(true);

    for (int i = 0; i <= nDataPoints; ++i)
    {
        float deltaChiSquaredCut(-5.0 + (i * (5.0/nDataPoints)));
        int nTotalCRs(0), nTotalUpwardsNeutrinoPfos(0), nTotalDownwardsNeutrinoPfos(0);
        int nRemovedPfos_DCHS(0), nCosmicsRemoved_DCHS(0), nUpwardsNeutrinoPfosRemoved_DCHS(0), nDownwardsNeutrinoPfosRemoved_DCHS(0);

        for (int i = 1; i < t1->GetEntries(); i++)
        {    
            t1->GetEntry(i);

            bool skipThisEvent(false);
       
            for (int i = 0; i < static_cast<int>(filterValues.size()); ++i)
            {   
                if (filterVariables[i] != filterVariableTargetValues[i])
                {   
                    skipThisEvent = true;
                    break;
                }   
            }   

            if (skipThisEvent)
                continue;

            if (mcCosmicRay == 1)
                ++nTotalCRs;
            if (mcCosmicRay == 0 && trueDownwards == 0) 
                ++nTotalUpwardsNeutrinoPfos;
            if (mcCosmicRay == 0 && trueDownwards == 1) 
                ++nTotalDownwardsNeutrinoPfos;

            if (deltaChiSquaredUpDownPerHit <= deltaChiSquaredCut)
            {
                ++nRemovedPfos_DCHS;

                if (mcCosmicRay == 1)
                    ++nCosmicsRemoved_DCHS;
                if (mcCosmicRay == 0 && trueDownwards == 0) 
                    ++nUpwardsNeutrinoPfosRemoved_DCHS;
                if (mcCosmicRay == 0 && trueDownwards == 1) 
                    ++nDownwardsNeutrinoPfosRemoved_DCHS;
            }
        }

        if (i == 0)
        {
            std::cout << "****************************************" << std::endl;
            std::cout << "nTotalCRs: " << nTotalCRs << std::endl;
            std::cout << "nTotalUpwardsNeutrinoPfos: " << nTotalUpwardsNeutrinoPfos << std::endl;
            std::cout << "nTotalDownwardsNeutrinoPfos: " << nTotalDownwardsNeutrinoPfos << std::endl;
            std::cout << "****************************************" << std::endl;
        }

        float fractionCRRemoved_DCHS(static_cast<float>(nCosmicsRemoved_DCHS)/nTotalCRs);
        float fractionUpwardsNuPfosRetained_DCHS(1.0 - static_cast<float>(nUpwardsNeutrinoPfosRemoved_DCHS)/nTotalUpwardsNeutrinoPfos);
        pROCCurve->SetPoint(i, fractionUpwardsNuPfosRetained_DCHS, fractionCRRemoved_DCHS); 

    }

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,800);
    pROCCurve->GetXaxis()->SetLimits(0.9,0.99);
    pROCCurve->GetYaxis()->SetRangeUser(0.0,0.3);

    pROCCurve->SetTitle("#Delta#chi^{2}_{DU}/N Cut ROC Curve;Fraction Upwards #nu Retained (Signal Efficiency); Fraction Untagged CR Removed (Background Rejection)");
    pROCCurve->Draw("ALP");

    std::string plotName("CosmicRemoval_ROC_Curve");

    for (const auto &entry : filterValues)
    {
        plotName += entry.first;
        plotName += "_";
    }

    std::string filePath("/usera/jjd49/pandora_direction/Scripts/Figures/TopFaceCRRemoval/" + plotName + "_NegativeDCHS.pdf");
    c1->SaveAs(filePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CountEventNumbers(TTree* pTree, std::vector<std::pair<std::string, int>> filterValues)
{
    int mcCosmicRay, trueDownwards, sufficientThreeDHits, removedByRegularTagging, highYCut, polarAngleCut, pfoChargeCut, neutrinoMomentumZCut, numberConnectedPfosNotNearEndpoint, numberConnectedPfosNearEndpoint;

    pTree->SetBranchAddress("TrueDownwards", &trueDownwards); 
    pTree->SetBranchAddress("MCCosmicRay", &mcCosmicRay); 
    pTree->SetBranchAddress("RemovedByRegularTagging", &removedByRegularTagging); 
    pTree->SetBranchAddress("HighYCut", &highYCut); 
    pTree->SetBranchAddress("PolarAngleCut", &polarAngleCut); 
    pTree->SetBranchAddress("PfoChargeCut", &pfoChargeCut); 
    pTree->SetBranchAddress("NeutrinoMomentumZCut", &neutrinoMomentumZCut); 
    pTree->SetBranchAddress("NumberConnectedPfosNotNearEndpoint", &numberConnectedPfosNotNearEndpoint); 
    pTree->SetBranchAddress("NumberConnectedPfosNearEndpoint", &numberConnectedPfosNearEndpoint); 

    int filterVariables[99], filterVariableTargetValues[99], filterVariableCounter(0);

    for (int i = 0; i < static_cast<int>(filterValues.size()); ++i)
    {   
        std::cout << "FILTER: " << filterValues.at(i).first << std::endl;
        pTree->SetBranchAddress(filterValues.at(i).first.c_str(), &filterVariables[i]);
        filterVariableTargetValues[i] = filterValues.at(i).second;
        ++filterVariableCounter;
    }   

    int nTotalCRs(0), nTotalUpwardsNeutrinoPfos(0), nTotalDownwardsNeutrinoPfos(0);

    for (int i = 1; i < t1->GetEntries(); i++)
    {    
        t1->GetEntry(i);

        bool skipThisEvent(false);
   
        for (int i = 0; i < static_cast<int>(filterValues.size()); ++i)
        {   
            if (filterVariables[i] != filterVariableTargetValues[i])
            {   
                skipThisEvent = true;
                break;
            }   
        }   

        if (skipThisEvent)
            continue;

        if (mcCosmicRay == 1)
            ++nTotalCRs;
        if (mcCosmicRay == 0 && trueDownwards == 0) 
            ++nTotalUpwardsNeutrinoPfos;
        if (mcCosmicRay == 0 && trueDownwards == 1) 
            ++nTotalDownwardsNeutrinoPfos;
    }

    std::cout << "****************************************" << std::endl;
    std::cout << "nTotalCRs: " << nTotalCRs << std::endl;
    std::cout << "nTotalUpwardsNeutrinoPfos: " << nTotalUpwardsNeutrinoPfos << std::endl;
    std::cout << "nTotalDownwardsNeutrinoPfos: " << nTotalDownwardsNeutrinoPfos << std::endl;
    std::cout << "****************************************" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int CalculateNumberNuEventsLost(TTree* pTree, float cutValue, float polarAngleCut)
{
    //Calculates the number of events that are true nu-induced and flagged as ambiguous (not an obvious CR) at the CR tagging stage
    float cosmicProbability, deltaChiSquaredUpDownPerHit, polarAngle;
    int targetPfo, mcCosmicRay;

    pTree->SetBranchAddress("CosmicProbability", &cosmicProbability); 
    pTree->SetBranchAddress("TargetPFO", &targetPfo); 
    pTree->SetBranchAddress("MCCosmicRay", &mcCosmicRay); 
    pTree->SetBranchAddress("DeltaChiSquaredUpDownPerHit", &deltaChiSquaredUpDownPerHit); 
    pTree->SetBranchAddress("PfoPolarAngle", &polarAngle); 

    int eventCounter(0);

    for (int i = 1; i < t1->GetEntries(); i++)
    {    
        t1->GetEntry(i);

        if (targetPfo == 0 || cosmicProbability < 0.5 || polarAngle >= polarAngleCut)
            continue;

        if (mcCosmicRay == 0)
            ++eventCounter;
    }

    return eventCounter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AssessNeutrinoIdPerformance(TTree* pTreeWithoutCut, TTree* pTreeWithCut)
{
    //Filter events here
    //std::vector<std::pair<std::string, int>> filterValues = {};
    std::vector<std::pair<std::string, int>> filterValues = {}; 

    TH1F* pHitPurityHistogram_NoCut = new TH1F("pHitPurityHistogram_NoCut", "", 100, 0.0, 1.0); 
    CreateFilteredHistogram(pTreeWithoutCut, pHitPurityHistogram_NoCut, "NuRecoHitPurityFraction", filterValues, "float");

    TH1F* pHitPurityHistogram_WithCut = new TH1F("pHitPurityHistogram_WithCut", "", 100, 0.0, 1.0); 
    CreateFilteredHistogram(pTreeWithCut, pHitPurityHistogram_WithCut, "NuRecoHitPurityFraction", filterValues, "float");

    Draw(pHitPurityHistogram_NoCut, pHitPurityHistogram_WithCut, kBlue, kRed, true, "NuId_HitPurity_Comparison", "Hit Purity Fraction", "Fraction of Events", "Normalised Hit Purity Fraction Distributions: Before & After #Delta#chi^{2}_{DU}/N Cut", "No Cut", "With Cut");

    TH1F* pNuIdFailureHistogram_NoCut = new TH1F("pNuIdFailureHistogram_NoCut", "", 2, 0, 2); 
    CreateFilteredHistogram(pTreeWithoutCut, pNuIdFailureHistogram_NoCut, "NuRecoTrueCosmicRay", filterValues, "int");

    TH1F* pNuIdFailureHistogram_WithCut = new TH1F("pNuIdFailureHistogram_WithCut", "", 2, 0, 2); 
    CreateFilteredHistogram(pTreeWithCut, pNuIdFailureHistogram_WithCut, "NuRecoTrueCosmicRay", filterValues, "int");

    //Draw(pNuIdFailureHistogram_NoCut, pNuIdFailureHistogram_NoCut, kBlue, kRed, false, "NuId_Failure_Comparison", "Neutrino ID Failure? (bool)", "Number of Events", "Neutrino ID Failure Distributions: Before & After #Delta#chi^{2}_{DU}/N Cut", "No Cut", "With Cut");

    std::cout << "Number of Events without cut: " << pNuIdFailureHistogram_NoCut->GetEntries() << std::endl; 
    std::cout << "Number of Events with cut: " << pNuIdFailureHistogram_WithCut->GetEntries() << std::endl; 
    std::cout << "-----------------------------------------" << std::endl; 
    std::cout << "Number of Successes without cut: " << pNuIdFailureHistogram_NoCut->GetBinContent(1) << std::endl; 
    std::cout << "Number of Successes with cut: " << pNuIdFailureHistogram_WithCut->GetBinContent(1) << std::endl; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FindBadEvents(TTree* pTree)
{
    float deltaChiSquaredUpDownPerHit;
    int targetPfo, trueDownwards, mcCosmicRay, fileIdentifier, eventNumber;

    pTree->SetBranchAddress("DeltaChiSquaredUpDownPerHit", &deltaChiSquaredUpDownPerHit); 
    pTree->SetBranchAddress("TargetPFO", &targetPfo); 
    pTree->SetBranchAddress("TrueDownwards", &trueDownwards); 
    pTree->SetBranchAddress("MCCosmicRay", &mcCosmicRay); 
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier); 
    pTree->SetBranchAddress("EventNumber", &eventNumber); 

    for (int i = 0; i < pTree->GetEntries(); i++)
    {    
        t1->GetEntry(i);

        if (targetPfo == 1 && mcCosmicRay == 0 && trueDownwards == 0 && deltaChiSquaredUpDownPerHit < 0.f && eventNumber < 5)
            std::cout << "Bad event: " << fileIdentifier << ":" << eventNumber << " with delta chi squared: " << deltaChiSquaredUpDownPerHit << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GetNinetyPercentBoundaries(TH1F* pCosmicHistogram, TH1F* pNeutrinoHistogram)
{
    float cr_integral(0.f), nu_integral(0.f);
    float cr_reverse_integral(0.f), nu_reverse_integral(0.f);

    for (int i = 0; i < pCosmicHistogram->GetNbinsX(); ++i)
    {
        if (cr_integral >= 0.9)
        {
            std::cout << "Bin containing " << cr_integral << "% of CR entries: " << pCosmicHistogram->GetBinCenter(i) << std::endl;
            std::cout << "This cut contains " << nu_integral * 100.0 << "% of neutrino entries" << std::endl;
            break;
        }

        cr_integral += (static_cast<float>(pCosmicHistogram->GetBinContent(i))/static_cast<float>(pCosmicHistogram->GetEntries()));
        nu_integral += (static_cast<float>(pNeutrinoHistogram->GetBinContent(i))/static_cast<float>(pNeutrinoHistogram->GetEntries()));
    }

    for (int i = pCosmicHistogram->GetNbinsX(); i > 0; --i)
    {
        if (cr_reverse_integral >= 0.9)
        {
            std::cout << "Reverse direction bin containing " << cr_reverse_integral << "% of CR entries: " << pCosmicHistogram->GetBinCenter(i) << std::endl;
            std::cout << "This cut contains " << nu_reverse_integral * 100.0 << "% of neutrino entries" << std::endl;
            break;
        }

        cr_reverse_integral += (static_cast<float>(pCosmicHistogram->GetBinContent(i))/static_cast<float>(pCosmicHistogram->GetEntries()));
        nu_reverse_integral += (static_cast<float>(pNeutrinoHistogram->GetBinContent(i))/static_cast<float>(pNeutrinoHistogram->GetEntries()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Neutrino_CR_Comparison(TTree* pTree, std::vector<std::pair<std::string, int>> additionalFilterValues, std::string variableName, int numberBins, float lowerBound, float upperBound, std::string typeName)
{
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("MCCosmicRay", 1)};
    std::vector<std::pair<std::string, int>> filterValues_Neutrino = {std::make_pair("MCCosmicRay", 0)};

    for (const auto &pair : additionalFilterValues)
    {
        filterValues.push_back(pair);
        filterValues_Neutrino.push_back(pair);
    }

    TH1F* pHistogram = new TH1F("pHistogram", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pHistogram, variableName, filterValues, typeName);

    TH1F* pHistogram_Neutrino = new TH1F("pHistogram_Neutrino", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pHistogram_Neutrino, variableName, filterValues_Neutrino, typeName);

    //Counting for daughter-like and delta-ray-like PFOs
    if (variableName == "NumberConnectedPfosNearEndpoint")
    {
        //std::cout << 1.0 - static_cast<float>(pHistogram->GetBinContent(1))/pHistogram->GetEntries() << std::endl;
        //std::cout << 1.0 - static_cast<float>(pHistogram_Neutrino->GetBinContent(1))/pHistogram_Neutrino->GetEntries() << std::endl;

        std::cout << "Fraction of cosmic rays in bins 3-6: " << 1.0 - (static_cast<float>(pHistogram->GetBinContent(4))/pHistogram->GetEntries() + (pHistogram->GetBinContent(5))/pHistogram->GetEntries() + (pHistogram->GetBinContent(6))/pHistogram->GetEntries()) << std::endl;
        std::cout << "Fraction of neutrino PFOS in bins 3-6: " << 1.0 - (static_cast<float>(pHistogram_Neutrino->GetBinContent(4))/pHistogram_Neutrino->GetEntries() + (pHistogram_Neutrino->GetBinContent(5))/pHistogram_Neutrino->GetEntries() + (pHistogram_Neutrino->GetBinContent(6))/pHistogram_Neutrino->GetEntries()) << std::endl;
    }
    
    GetNinetyPercentBoundaries(pHistogram, pHistogram_Neutrino);

    std::string plotName(variableName+"_Comparison_Nu_CR_");

    for (const auto &entry : additionalFilterValues)
    {
        plotName += entry.first;
        plotName += "_";
    }

    Draw(pHistogram, pHistogram_Neutrino, kBlue, kRed, true, plotName, variableName, "Number of Events", variableName + " Distribution Comparing #nu-Induced/CR PFOs", "True Cosmic", "True Neutrino-Induced");

    delete pHistogram;
    delete pHistogram_Neutrino;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Neutrino_CR_Comparison(TTree* pTree, std::vector<std::pair<std::string, int>> additionalFilterValues, std::string variableName, int numberBins, float lowerBound, float upperBound, std::string typeName, std::string plotTitle, std::string xAxisTitle, std::string yAxisTitle)
{
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("MCCosmicRay", 1)};
    std::vector<std::pair<std::string, int>> filterValues_Neutrino = {std::make_pair("MCCosmicRay", 0)};

    for (const auto &pair : additionalFilterValues)
    {
        filterValues.push_back(pair);
        filterValues_Neutrino.push_back(pair);
    }

    TH1F* pHistogram = new TH1F("pHistogram", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pHistogram, variableName, filterValues, typeName);

    TH1F* pHistogram_Neutrino = new TH1F("pHistogram_Neutrino", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pHistogram_Neutrino, variableName, filterValues_Neutrino, typeName);

    //Counting for daughter-like and delta-ray-like PFOs
    if (variableName == "NumberConnectedPfosNearEndpoint")
    {
        //std::cout << 1.0 - static_cast<float>(pHistogram->GetBinContent(1))/pHistogram->GetEntries() << std::endl;
        //std::cout << 1.0 - static_cast<float>(pHistogram_Neutrino->GetBinContent(1))/pHistogram_Neutrino->GetEntries() << std::endl;

        std::cout << "Fraction of cosmic rays in bins 3-6: " << 1.0 - (static_cast<float>(pHistogram->GetBinContent(4))/pHistogram->GetEntries() + (pHistogram->GetBinContent(5))/pHistogram->GetEntries() + (pHistogram->GetBinContent(6))/pHistogram->GetEntries()) << std::endl;
        std::cout << "Fraction of neutrino PFOS in bins 3-6: " << 1.0 - (static_cast<float>(pHistogram_Neutrino->GetBinContent(4))/pHistogram_Neutrino->GetEntries() + (pHistogram_Neutrino->GetBinContent(5))/pHistogram_Neutrino->GetEntries() + (pHistogram_Neutrino->GetBinContent(6))/pHistogram_Neutrino->GetEntries()) << std::endl;
    }
    
    GetNinetyPercentBoundaries(pHistogram, pHistogram_Neutrino);

    std::string plotName(variableName+"_Comparison_Nu_CR_");

    for (const auto &entry : additionalFilterValues)
    {
        plotName += entry.first;
        plotName += "_";
    }

    Draw(pHistogram, pHistogram_Neutrino, kBlue, kRed, true, plotName, xAxisTitle, yAxisTitle, plotTitle, "True Cosmic", "True Neutrino-Induced");

    delete pHistogram;
    delete pHistogram_Neutrino;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Upwards_Downwards_Neutrino_Comparison(TTree* pTree, std::string variableName, int numberBins, float lowerBound, float upperBound, std::string typeName)
{
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("MCCosmicRay", 0), std::make_pair("TrueDownwards", 0)};
    std::vector<std::pair<std::string, int>> filterValues_Downwards = {std::make_pair("MCCosmicRay", 0), std::make_pair("TrueDownwards", 1)};
    //std::vector<std::pair<std::string, int>> filterValues = {};

    TH1F* pHistogram = new TH1F("pHistogram", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pHistogram, variableName, filterValues, typeName);

    TH1F* pHistogram_Downwards = new TH1F("pHistogram_Downwards", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pHistogram_Downwards, variableName, filterValues_Downwards, typeName);

    GetNinetyPercentBoundaries(pHistogram, pHistogram_Downwards);

    Draw(pHistogram, pHistogram_Downwards, kBlue, kRed, true, variableName+"_Comparison_Nu_Up_Down", variableName, "Number of Events", variableName + " Distribution Comparing True Upwards/Downwards #nu-Induced Particles", "True Upwards", "True Downwards");

    delete pHistogram;
    delete pHistogram_Downwards;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateDistributions(TTree* pTree)
{
    /*
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("MCCosmicRay", 1), std::make_pair("MCPDG", 13)};
    std::vector<std::pair<std::string, int>> filterValues_Neutrino = {std::make_pair("MCCosmicRay", 0), std::make_pair("MCPDG", 13)};

    TH1F* pHighYHistogram = new TH1F("pHighYHistogram", "", 100, 80.0, 150.0); 
    CreateFilteredHistogram(t1, pHighYHistogram, "PfoHighY", filterValues, "float");

    TH1F* pHighYHistogram_Neutrino = new TH1F("pHighYHistogram_Neutrino", "", 100, 80.0, 150.0); 
    CreateFilteredHistogram(t1, pHighYHistogram_Neutrino, "PfoHighY", filterValues_Neutrino, "float");

    Draw(pHighYHistogram, pHighYHistogram_Neutrino, kBlue, kRed, true, "pHighYHistogram_Comparison", "Reconstucted High Y Coordinate (cm)", "Number of Events", "Y_{high} Distribution (With SCE) For True Primary Muons", "True Cosmic", "True Neutrino-Induced");

    std::vector<std::pair<std::string, int>> filterValues_Cosmic_Neutrino = {std::make_pair("MCCosmicRay", 0), std::make_pair("MCPDG", 13)};

    TH1F* pNumberThreeDHitsHistogram = new TH1F("pNumberThreeDHitsHistogram", "", 100, 0.0, 2000.0); 
    CreateFilteredHistogram(t1, pNumberThreeDHitsHistogram, "NumberThreeDHits", filterValues_Cosmic_Muon, "int");

    TH1F* pNumberThreeDHitsHistogram_Neutrino = new TH1F("pNumberThreeDHitsHistogram_Neutrino", "", 100, 0.0, 2000.0); 
    CreateFilteredHistogram(t1, pNumberThreeDHitsHistogram_Neutrino, "NumberThreeDHits", filterValues_Cosmic_Neutrino, "int");

    Draw(pNumberThreeDHitsHistogram, pNumberThreeDHitsHistogram_Neutrino, kBlue, kRed, true, "pNumberThreeDHitsHistogram_Comparison", "Number of 3D Hits", "Number of Events", "Number of 3D Hits Distribution (Cosmic Muons & Neutrinos)", "True Cosmic Muon", "True Neutrino-induced Muon");
    */

    //TH1F* pLowYHistogram = new TH1F("pLowYHistogram", "", 100, -120.0, 120.0); 
    //CreateFilteredHistogram(t1, pLowYHistogram, "PfoLowY", filterValues, "float");
    //Draw(pLowYHistogram, kBlue, "pLowYHistogram", "Reconstucted Low Y Coordinate (cm)", "Number of Events", "True CR && Target PFO Reconstructed Low Y Coordinate Distribution");

    //TH1F* pPolarAngleHistogram = new TH1F("pPolarAngleHistogram", "", 100, 0.0, 3.1415); 
    //CreateFilteredHistogram(t1, pPolarAngleHistogram, "PfoPolarAngle", filterValues, "float");
    //Draw(pPolarAngleHistogram, kBlue, "pPolarAngleHistogram", "Polar Angle #theta (rad)", "Number of Events", "Polar Angle #theta Distribution");

    /*
    std::vector<std::pair<std::string, int>> filterValuesBraggFinder = {std::make_pair("MCFiducialLowY", 1), std::make_pair("MCCosmicRay", 1)};
    std::vector<std::pair<std::string, int>> filterValuesBraggFinder_notContained = {std::make_pair("MCFiducialLowY", 0), std::make_pair("MCCosmicRay", 1)};

    TH1F* pExtentHistogram = new TH1F("pExtentHistogram", "", 100, 0.0, 300.0); 
    CreateFilteredHistogram(t1, pExtentHistogram, "TrueExtent", filterValues, "float");
    Draw(pExtentHistogram, kBlue, "pExtentHistogram", "Length L (cm)", "Number of Events", "Length L Distribution");

    TH1F* pFitEndpointEnergyHistogram = new TH1F("pFitEndpointEnergyHistogram", "", 100, 0.0, 100.0); 
    CreateFilteredHistogram(t1, pFitEndpointEnergyHistogram, "FitEndpointEnergy", filterValuesBraggFinder, "float");

    TH1F* pFitEndpointEnergyHistogram_notContained = new TH1F("pFitEndpointEnergyHistogram_notContained", "", 100, 0.0, 1000.0); 
    CreateFilteredHistogram(t1, pFitEndpointEnergyHistogram_notContained, "FitEndpointEnergy", filterValuesBraggFinder_notContained, "float");

    Draw(pFitEndpointEnergyHistogram, pFitEndpointEnergyHistogram_notContained, kBlue, kRed, true, "FitEndpointEnergy_Comparison", "Fit End Energy (MeV)", "Number of Events", "Fit End Energy Distribution", "True Contained Endpoint", "True Uncontained Endpoint");

    TH1F* pEndpointFitChargeRatioHistogram = new TH1F("pEndpointFitChargeRatioHistogram", "", 100, 0.0, 2.0); 
    CreateFilteredHistogram(t1, pEndpointFitChargeRatioHistogram, "EndpointFitChargeRatio", filterValuesBraggFinder, "float");

    TH1F* pEndpointFitChargeRatioHistogram_notContained = new TH1F("pEndpointFitChargeRatioHistogram_notContained", "", 100, 0.0, 2.0); 
    CreateFilteredHistogram(t1, pEndpointFitChargeRatioHistogram_notContained, "EndpointFitChargeRatio", filterValuesBraggFinder_notContained, "float");

    Draw(pEndpointFitChargeRatioHistogram, kBlue, "pEndpointFitChargeRatioHistogram", "First/Last 5 Fit Charge Ratio", "Number of Events", "Fit Charge Ratio Distribution");
    Draw(pEndpointFitChargeRatioHistogram_notContained, kRed, "pEndpointFitChargeRatioHistogram_notContained", "First/Last 5 Fit Charge Ratio", "Number of Events", "Fit Charge Ratio Distribution");

    Draw(pEndpointFitChargeRatioHistogram, pEndpointFitChargeRatioHistogram_notContained, kBlue, kRed, false, "EndpointFitChargeRatio_Comparison", "First/Last 5 Fit Charge Ratio", "Number of Events", "Fit Charge Ratio Distribution", "True Contained Endpoint", "True Uncontained Endpoint");
    */

    TH1F* pCosmicDCHSHistogram = new TH1F("pCosmicDCHSHistogram", "", 100, -5.0, 5.0); 
    std::vector<std::pair<std::string, int>> filterValuesCosmic = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("PfoHighTopY", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("FitEndpointEnergyIsSmall", 1), std::make_pair("MCCosmicRay", 1)};
    CreateFilteredHistogram(t1, pCosmicDCHSHistogram, "DeltaChiSquaredUpDownPerHit", filterValuesCosmic, "float");

    TH1F* pNeutrinoDCHSHistogram = new TH1F("pNeutrinoDCHSHistogram", "", 100, -5.0, 5.0); 
    std::vector<std::pair<std::string, int>> filterValuesNu = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("PfoHighTopY", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("FitEndpointEnergyIsSmall", 1), std::make_pair("MCCosmicRay", 0)};
    CreateFilteredHistogram(t1, pNeutrinoDCHSHistogram, "DeltaChiSquaredUpDownPerHit", filterValuesNu, "float");

    std::string plotName("DCHS_Distribution_Comparison_");

    for (const auto &pair : filterValuesCosmic)
    {
        plotName += pair.first;
        plotName += "_";
    }

    //Draw(pNeutrinoDCHSHistogram, pCosmicDCHSHistogram, kBlue, kRed, false, plotName.c_str(), "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "#Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");
    Draw(pNeutrinoDCHSHistogram, pCosmicDCHSHistogram, kBlue, kRed, true, plotName.c_str(), "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised #Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateDirectionalDistributions(TTree* pTree, std::string variableName, std::string titleString, int numberBins, float lowerBound, float upperBound)
{
    std::vector<std::pair<std::string, int>> filterValuesCosmic = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("PfoHighTopY", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("FitEndpointEnergyIsSmall", 1), std::make_pair("MCCosmicRay", 1), std::make_pair("TrueDownwards", 1)};
    std::vector<std::pair<std::string, int>> filterValuesUpwardsNeutrino = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("PfoHighTopY", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("FitEndpointEnergyIsSmall", 1), std::make_pair("MCCosmicRay", 0), std::make_pair("TrueDownwards", 0)};
    std::vector<std::pair<std::string, int>> filterValuesDownwardsNeutrino = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("PfoHighTopY", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("FitEndpointEnergyIsSmall", 1), std::make_pair("MCCosmicRay", 0), std::make_pair("TrueDownwards", 1)};

    TH1F* pCosmicDirectionHistogram = new TH1F("pCosmicDirectionHistogram", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pCosmicDirectionHistogram, variableName, filterValuesCosmic, "float");

    TH1F* pUpwardsNeutrinoDirectionHistogram = new TH1F("pUpwardsNeutrinoDirectionHistogram", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pUpwardsNeutrinoDirectionHistogram, variableName, filterValuesUpwardsNeutrino, "float");

    TH1F* pDownwardsNeutrinoDirectionHistogram = new TH1F("pDownwardsNeutrinoDirectionHistogram", "", numberBins, lowerBound, upperBound); 
    CreateFilteredHistogram(t1, pDownwardsNeutrinoDirectionHistogram, variableName, filterValuesDownwardsNeutrino, "float");

    std::string plotName("Direction_Distribution_Comparison_");
    plotName += variableName;
    plotName += "_";

    for (const auto &pair : filterValuesCosmic)
    {
        plotName += pair.first;
        plotName += "_";
    }

    std::vector<TH1F*> histogramVector = {pCosmicDirectionHistogram, pUpwardsNeutrinoDirectionHistogram, pDownwardsNeutrinoDirectionHistogram};
    std::vector<EColor> colourVector = {kMagenta, kBlue, kRed};
    std::vector<std::string> legendVector = {"True Cosmic", "True Upwards #nu", "True Downwards #nu"};

    std::string plotTitle(titleString + " Distributions: Cosmic And Neutrino-Induced Muons");
    Draw(histogramVector, colourVector, legendVector, true, true, plotName, titleString, "Fraction of Events", plotTitle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SequentiallyApplyCutsToDCHSDistribution(TTree* pTree)
{
    std::vector<std::pair<std::string, int>> filterValuesCosmic = {std::make_pair("MCCosmicRay", 1), std::make_pair("TrueDownwards", 1), std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0), std::make_pair("FitEndpointEnergyIsSmall", 1)};
    std::vector<std::pair<std::string, int>> filterValuesNu_Up = {std::make_pair("MCCosmicRay", 0), std::make_pair("TrueDownwards", 0), std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0), std::make_pair("FitEndpointEnergyIsSmall", 1)};
    std::vector<std::pair<std::string, int>> filterValuesNu_Down = {std::make_pair("MCCosmicRay", 0), std::make_pair("TrueDownwards", 1), std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0), std::make_pair("FitEndpointEnergyIsSmall", 1)};

    while (filterValuesCosmic.size() > 2)
    {
        TH1F* pCosmicDCHSHistogram = new TH1F("pCosmicDCHSHistogram", "", 100, -5.0, 5.0); 
        CreateFilteredHistogram(t1, pCosmicDCHSHistogram, "DeltaChiSquaredUpDownPerHit", filterValuesCosmic, "float");

        TH1F* pUpwardsNeutrinoDCHSHistogram = new TH1F("pUpwardsNeutrinoDCHSHistogram", "", 100, -5.0, 5.0); 
        CreateFilteredHistogram(t1, pUpwardsNeutrinoDCHSHistogram, "DeltaChiSquaredUpDownPerHit", filterValuesNu_Up, "float");

        TH1F* pDownwardsNeutrinoDCHSHistogram = new TH1F("pDownwardsNeutrinoDCHSHistogram", "", 100, -5.0, 5.0); 
        CreateFilteredHistogram(t1, pDownwardsNeutrinoDCHSHistogram, "DeltaChiSquaredUpDownPerHit", filterValuesNu_Down, "float");

        std::string nameString("DCHS_Distribution_Comparison_");

        for (const auto &entry : filterValuesCosmic)
        {
            nameString += entry.first;
            nameString += "_";
        }

        //Draw(pNeutrinoDCHSHistogram, pCosmicDCHSHistogram, kBlue, kRed, false, "UpDownDeltaChi_Distribution_Comparison_Neutrino_Cosmic_PolarAngleCut_NotNormalised", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "#Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");
        //Draw(pNeutrinoDCHSHistogram, pCosmicDCHSHistogram, kBlue, kRed, true, nameString, "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised #Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");

        std::vector<TH1F*> histogramVector = {pCosmicDCHSHistogram, pUpwardsNeutrinoDCHSHistogram, pDownwardsNeutrinoDCHSHistogram};
        std::vector<EColor> colourVector = {kMagenta, kBlue, kRed};
        std::vector<std::string> legendVector = {"True Cosmic", "True Upwards #nu", "True Downwards #nu"};

        Draw(histogramVector, colourVector, legendVector, true, true, nameString, "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "#Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons");

        filterValuesCosmic.pop_back();
        filterValuesNu_Up.pop_back();
        filterValuesNu_Down.pop_back();

        delete pCosmicDCHSHistogram;
        delete pUpwardsNeutrinoDCHSHistogram;
        delete pDownwardsNeutrinoDCHSHistogram;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BraggPeakFinding(TTree* pTree)
{
    std::vector<std::pair<std::string, int>> filterValuesBraggFinder = {std::make_pair("MCFiducialLowY", 1), std::make_pair("MCCosmicRay", 1), std::make_pair("SufficientThreeDHits", 1)};
    std::vector<std::pair<std::string, int>> filterValuesBraggFinder_notContained = {std::make_pair("MCFiducialLowY", 0), std::make_pair("MCCosmicRay", 1), std::make_pair("SufficientThreeDHits", 1)};

    TH1F* pFitEndpointEnergyHistogram = new TH1F("pFitEndpointEnergyHistogram", "", 100, 0.0, 300.0); 
    CreateFilteredHistogram(pTree, pFitEndpointEnergyHistogram, "FitEndpointEnergy", filterValuesBraggFinder, "float");

    TH1F* pFitEndpointEnergyHistogram_notContained = new TH1F("pFitEndpointEnergyHistogram_notContained", "", 100, 0.0, 300.0); 
    CreateFilteredHistogram(pTree, pFitEndpointEnergyHistogram_notContained, "FitEndpointEnergy", filterValuesBraggFinder_notContained, "float");

    Draw(pFitEndpointEnergyHistogram, pFitEndpointEnergyHistogram_notContained, kBlue, kRed, true, "FitEndpointEnergy_Comparison", "Fit End Energy (MeV)", "Number of Events", "Fit End Energy Distribution", "True Contained Endpoint", "True Uncontained Endpoint");

    std::vector<std::pair<std::string, int>> filterValuesBraggFinder_smallEndCharge = {std::make_pair("MCFiducialLowY", 1), std::make_pair("MCCosmicRay", 1), std::make_pair("SufficientThreeDHits", 1), std::make_pair("FitEndpointEnergyIsSmall", 1)};
    std::vector<std::pair<std::string, int>> filterValuesBraggFinder_notContained_smallEndCharge = {std::make_pair("MCFiducialLowY", 0), std::make_pair("MCCosmicRay", 1), std::make_pair("SufficientThreeDHits", 1), std::make_pair("FitEndpointEnergyIsSmall", 1)};

    TH1F* pEndpointFitChargeRatioHistogram = new TH1F("pEndpointFitChargeRatioHistogram", "", 100, 0.0, 3.0); 
    CreateFilteredHistogram(pTree, pEndpointFitChargeRatioHistogram, "EndpointFitChargeRatio", filterValuesBraggFinder_smallEndCharge, "float");

    TH1F* pEndpointFitChargeRatioHistogram_notContained = new TH1F("pEndpointFitChargeRatioHistogram_notContained", "", 100, 0.0, 3.0); 
    CreateFilteredHistogram(pTree, pEndpointFitChargeRatioHistogram_notContained, "EndpointFitChargeRatio", filterValuesBraggFinder_notContained_smallEndCharge, "float");

    //Draw(pEndpointFitChargeRatioHistogram, kBlue, "pEndpointFitChargeRatioHistogram_smallEndCharge", "First/Last 5 Fit Charge Ratio", "Number of Events", "Fit Charge Ratio Distribution");
    //Draw(pEndpointFitChargeRatioHistogram_notContained, kRed, "pEndpointFitChargeRatioHistogram_notContained_smallEndCharge", "First/Last 5 Fit Charge Ratio", "Number of Events", "Fit Charge Ratio Distribution");

    Draw(pEndpointFitChargeRatioHistogram, pEndpointFitChargeRatioHistogram_notContained, kBlue, kRed, true, "EndpointFitChargeRatio_Comparison_smallEndCharge", "First/Last 5 Fit Charge Ratio", "Number of Events", "Fit Charge Ratio Distribution", "True Contained Endpoint", "True Uncontained Endpoint");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RankCuts(TTree* pTree)
{
    std::vector<std::pair<std::string, int>> additionalFilterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0)};
    //std::vector<std::pair<std::string, int>> additionalFilterValues = {};

    //Neutrino_CR_Comparison(t1, additionalFilterValues, "PfoHighY", 100, -120.0, 120.0, "float");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "PfoPolarAngle", 100, 0.0, 3.14/2.0, "float");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "PfoCharge", 100, 0.0, 3e5, "float");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "ApproximateNeutrinoMomentumZ", 100, 0.0, 1.0, "float");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "NumberConnectedPfosNearEndpoint", 8, 0, 8, "int");
    Neutrino_CR_Comparison(t1, additionalFilterValues, "NumberConnectedPfosNotNearEndpoint", 20, 0, 20, "int");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ApplyCutsSequentially(TTree* pTree)
{
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0)};

    std::cout << "VARIABLE: NumberConnectedPfosNearEndpoint" << std::endl;
    Neutrino_CR_Comparison(t1, filterValues, "NumberConnectedPfosNearEndpoint", 8, 0, 8, "int");
    filterValues.pop_back();

    std::cout << "VARIABLE: ApproximateNeutrinoMomentumZ" << std::endl;
    Neutrino_CR_Comparison(t1, filterValues, "ApproximateNeutrinoMomentumZ", 100, 0.0, 1.0, "float");
    filterValues.pop_back();

    std::cout << "VARIABLE: PfoCharge" << std::endl;
    Neutrino_CR_Comparison(t1, filterValues, "PfoCharge", 100, 0.0, 3e5, "float");
    filterValues.pop_back();

    std::cout << "VARIABLE: PfoPolarAngle" << std::endl;
    Neutrino_CR_Comparison(t1, filterValues, "PfoPolarAngle", 100, 0.0, 3.14/2.0, "float");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateNumberCosmicsWithCutNumberPlot(TTree* t1, int numberDirectionCuts, int maximumUpwardsNuLoss)
{
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0), std::make_pair("NumberConnectedPfosNearEndpointCut", 0)};

    int nCuts(static_cast<int>(filterValues.size()) - 2);
    TGraph* pNumberCosmicsRemovedWithCutNumberGraph = new TGraph(nCuts);
    pNumberCosmicsRemovedWithCutNumberGraph->GetXaxis()->SetNdivisions(nCuts + 2);
    
    for (int i = 1; i <= nCuts; ++i)
    {
        int cutValue = CreateCosmicRemovalEfficiencyCurves(t1, filterValues, numberDirectionCuts, maximumUpwardsNuLoss);
        pNumberCosmicsRemovedWithCutNumberGraph->SetPoint(nCuts - i, nCuts + 1 - i, cutValue);
        filterValues.pop_back();
    }

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,800);
    pNumberCosmicsRemovedWithCutNumberGraph->SetTitle("Number of Additional CR Removed With Cut Number; Cut Number; Number of Additional CR Removed");
    pNumberCosmicsRemovedWithCutNumberGraph->Draw("ALP");
    c1->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/TopFaceCRRemoval/AdditionalCosmics_CutNumber_Graph_10UpwardsNuLost.pdf");
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GetCutValue(TTree* t1, int numberDirectionCuts, int maximumUpwardsNuLoss)
{
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0), std::make_pair("NumberConnectedPfosNearEndpointCut", 0)};

    int cutValue = CreateCosmicRemovalEfficiencyCurves(t1, filterValues, numberDirectionCuts, maximumUpwardsNuLoss);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateEventNumbersTable(TTree* t1)
{
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0), std::make_pair("NumberConnectedPfosNearEndpointCut", 0)};

    while (filterValues.size() > 0)
    {
        CountEventNumbers(t1, filterValues);
        filterValues.pop_back();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Top_Face_CR_Removal(void)
{
    //19/02/19 is a good reference presentation

    //AssessNeutrinoIdPerformance(t2, t3);

    //use this to create slide 9
    //CreateCosmicRemovalEfficiencyCurves(t1, 20);

    //std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0)};
    //CreateCosmicRemovalEfficiencyCurves(t1, filterValues, 50, 2625);

    //use this to create a table of CR, nu up and nu down per cut
    //CreateEventNumbersTable(t1);

    //use this to create slide 11
    CreateNumberCosmicsWithCutNumberPlot(t1, 200, 400);

    //GetCutValue(t1, 1000, 1000);

    //std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0), std::make_pair("PolarAngleCut", 0), std::make_pair("PfoChargeCut", 0), std::make_pair("NeutrinoMomentumZCut", 0), std::make_pair("NumberConnectedPfosNearEndpointCut", 0)};
    //CreateCosmicRemovalROCCurve(t1, filterValues, 100);

    //RankCuts(t1);

    //use this to create slide 8
    //ApplyCutsSequentially(t1);

    //CreateDistributions(t1);

    //use this to create slide 10
    //CreateDirectionalDistributions(t1, "DeltaChiSquaredUpDownPerHit", "#Delta#chi^{2}_{DU}/N", 100, -5.0, 5.0);
    //CreateDirectionalDistributions(t1, "DeltaChiSquaredUpDown", "#Delta#chi^{2}_{DU}", 100, -500.0, 500.0);

    //std::vector<std::pair<std::string, int>> additionalFilterValues = {std::make_pair("SufficientThreeDHits", 1), std::make_pair("RemovedByRegularTagging", 0), std::make_pair("HighYCut", 0)};
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "PfoHighY", 100, 80.0, 150.0, "float", "Distribution of y_{high}", "y_{high} (cm)", "Number of Entries");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "PfoPolarAngle", 100, 0.0, 3.14/2.0, "float", "Distribution of Polar Angle #theta", "Polar Angle #theta (rad)", "Number of Entries");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "PfoCharge", 100, 0.0, 3e5, "float", "Distribution of Q_{W}", "Q_{W} (Integrated ADC Counts)", "Number of Entries");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "ApproximateNeutrinoMomentumZ", 100, 0.0, 1.0, "float", "Distribution of Approximate #nu Momentum Z-Component", "#nu Momentum Z-Component |#hat{p}_{#nu} #bullet #hat{z}|", "Number of Entries");
    //Neutrino_CR_Comparison(t1, additionalFilterValues, "NumberConnectedPfosNearEndpoint", 8, 0, 8, "int", "Distribution of The Number of Daughter-Like Particles N_{D}", "Number of Daughter-Like Particles N_{D}", "Number of Entries");

    //Neutrino_CR_Comparison(t1, additionalFilterValues, "NumberConnectedPfosNotNearEndpoint", 20, 0, 20, "int");

    //Neutrino_CR_Comparison(t1, "PfoAzimuthalAngle", 100, 0.0, 3.14, "float");
    //Upwards_Downwards_Neutrino_Comparison(t1, "PfoAzimuthalAngle", 100, 0.0, 3.14, "float");
    //Upwards_Downwards_Neutrino_Comparison(t1, "HighYZLargerThanLowYZ", 2, 0, 2, "int");

    //SequentiallyApplyCutsToDCHSDistribution(t1);
    
    //CreateDistributions(t1);


    //FindBadEvents(t1);

    //BraggPeakFinding(t1);
}

//------------------------------------------------------------------------------------------------------------------------------------------
