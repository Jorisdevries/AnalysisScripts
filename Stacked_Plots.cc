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

#include "InteractionType_Definitions.h"
#include "Drawing_Helper_Functions.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

//------------------------------------------------------------------------------------------------------------------------------------------

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/0.root");
TTree *t1 = (TTree*)f1->Get("EventSelection");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/1.root");
TTree *t2 = (TTree*)f2->Get("EventSelection");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/2.root");
TTree *t3 = (TTree*)f3->Get("EventSelection");

TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/3.root");
TTree *t4 = (TTree*)f4->Get("EventSelection");

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateIntStackHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, int targetMultiplicity, int targetInteractionType, bool onlyContained)
{
    int interactionType;
    int variable;
    int particleMultiplicity;
    int nCosmics;
    int nuRecoContained;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress("ModifiedInteractionType", &interactionType);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmics);
    pTree->SetBranchAddress("NuRecoContained", &nuRecoContained);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {   
        pTree->GetEntry(i);

        if (onlyContained && nuRecoContained == 0)
            continue;

        if (particleMultiplicity == 2 && interactionType == 0 && nCosmics == 1)
            interactionType = 169;

        if (particleMultiplicity == 2 && interactionType == 163 && nCosmics == 1)
            interactionType = 170;

        if (interactionType != targetInteractionType || particleMultiplicity != targetMultiplicity)
            continue;

        pHistogram->Fill(variable);
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateFloatStackHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, int targetMultiplicity, int targetInteractionType, bool onlyContained)
{
    float variable;
    int interactionType;
    int particleMultiplicity;
    int nCosmics;
    int nuRecoContained;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress("ModifiedInteractionType", &interactionType);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmics);
    pTree->SetBranchAddress("NuRecoContained", &nuRecoContained);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {   
        pTree->GetEntry(i);

        if (onlyContained && nuRecoContained == 0)
            continue;

        if (particleMultiplicity == 2 && interactionType == 0 && nCosmics == 1)
            interactionType = 169;

        if (particleMultiplicity == 2 && interactionType == 163 && nCosmics == 1)
            interactionType = 170;

        if (interactionType != targetInteractionType || particleMultiplicity != targetMultiplicity)
            continue;

        pHistogram->Fill(variable);
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStack(TTree* pTree, std::string variableName, std::vector<int> &targetInteractionTypes, int targetNumberParticles, bool onlyContained, int nBins, float minRange, float maxRange, std::string xTitle, std::string yTitle, std::string plotTitle, std::string legendPosition, std::string typeName)
{
    THStack *pStack = new THStack("pStack","");
    std::vector<TH1F*> histogramVector;

    //TLegend* pLegend = new TLegend(0.68,0.72,0.98,0.92);
    TLegend* pLegend = new TLegend(0.6,0.55,0.85,0.85);
    TLegend* pLeftLegend = new TLegend(0.15,0.55,0.4,0.85);

    if (legendPosition == "L")
        pLegend = pLeftLegend;

    for (const auto &interactionType : targetInteractionTypes)
    {
        std::string histogramName("pHistogram" + std::to_string(interactionType));
        TH1F *pHistogram = new TH1F(histogramName.c_str(),"", nBins, minRange, maxRange);
        pHistogram->SetFillStyle(3004);

        if (typeName == "int")
            CreateIntStackHistogram(pTree, pHistogram, variableName, targetNumberParticles, interactionType, onlyContained);
        else if (typeName == "float")
            CreateFloatStackHistogram(pTree, pHistogram, variableName, targetNumberParticles, interactionType, onlyContained);
        else
            std::cout << "Unrecognised typename" << std::endl;

        if (interactionType == 169)
            pLegend->AddEntry(pHistogram, "CCQEL_MU_1_CR");
        else if (interactionType == 170)
            pLegend->AddEntry(pHistogram, "OTHER_INTERACTION_1_CR");
        else
            pLegend->AddEntry(pHistogram, ToString(static_cast<InteractionType>(interactionType)).c_str());

        histogramVector.push_back(pHistogram);
    }

    for (const auto pHistogram : histogramVector)
        pStack->Add(pHistogram);

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    std::string fullTitle(plotTitle + ";" + xTitle + ";" + yTitle);
    pStack->SetTitle(fullTitle.c_str());
    pStack->Draw("pfc plc");
    pLegend->Draw("same");

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Stacked_Plots/N" + std::to_string(targetNumberParticles) + "_" + variableName + "_Stack.pdf");
    canvas->SaveAs(savePath.c_str());
    delete canvas;

    for (const auto pHistogram : histogramVector)
       delete pHistogram; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStackedPlots(TTree* pTree, std::vector<int> &interactionTypes, std::vector<std::string> &variableNames, std::vector<std::string> &titleNames, std::vector<float> lowerBounds, std::vector<float> upperBounds, std::vector<std::string> legendPositions, int particleMultiplicity, bool onlyContained)
{
    gStyle->SetPalette(kPastel);

    for (int i = 0; i < static_cast<int>(variableNames.size()); ++i)
        CreateStack(pTree, variableNames.at(i), interactionTypes, particleMultiplicity, onlyContained, 100, lowerBounds.at(i), upperBounds.at(i), titleNames.at(i), "Number of Events", "Distribution of " + titleNames.at(i) + " for N=" + std::to_string(particleMultiplicity), legendPositions.at(i), "float");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Stacked_Plots(void)
{
    /*
    std::vector<int> interactionTypes = {167, 165, 1, 163, 12, 0};
    std::vector<std::string> variableNames = {"LongestPfoCharge", "LongestPfoLength", "LongestPfoPolarAngle", "LongestPfoAzimuthalAngle"};
    std::vector<std::string> titleNames = {"Total PFO Charge Q_{PFO}", "Total PFO Length L_{PFO}", "PFO Polar Angle #theta", "PFO Azimuthal Angle #phi"};
    std::vector<float> lowerBounds = {0.0, 0.0, 0.0, 0.0};
    std::vector<float> upperBounds = {3e5, 3e2, 3.14, 3.14};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, 1);
    */

    /*
    std::vector<int> interactionTypes = {167, 165, 23, 18, 100, 0, 163, 1}; 
    std::vector<std::string> variableNames = {"ShortestPfoCharge", "ShortestPfoLength"};
    std::vector<std::string> titleNames = {"Total PFO Charge Q_{PFO}", "Total PFO Length L_{PFO}"};
    std::vector<float> lowerBounds = {0.0, 0.0};
    std::vector<float> upperBounds = {5e4, 5e1};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, 2);
    */

    /*
    std::vector<int> interactionTypes = {167, 165, 169, 170, 23, 18, 100, 0, 163, 1}; 
    std::vector<std::string> variableNames = {"OpeningAngle", "NeutrinoMomentumX", "NeutrinoMomentumY", "NeutrinoMomentumZ"};
    std::vector<std::string> titleNames = {"Opening Angle #zeta (rad)", "Neutrino Momentum p_{#nu} X Component", "Neutrino Momentum p_{#nu} Y Component", "Neutrino Momentum p_{#nu} Z Component"};
    std::vector<float> lowerBounds = {0.0, 0.0, 0.0, 0.0};
    std::vector<float> upperBounds = {3.14, 1.0, 1.0, 1.0};
    std::vector<std::string> legendPositions = {"R", "R", "R", "L"};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, legendPositions, 2);
    */

    /*
    std::vector<int> interactionTypes = {167, 165, 1, 163, 12, 0};
    std::vector<std::string> variableNames = {"LongestPfoCorrectedMaxY", "LongestPfoCorrectedMaxZ"};
    std::vector<std::string> titleNames = {"PFO Endpoint Maximum Y Coordinate (cm)", "PFO Endpoint Maximum Z Coordinate (cm)"};
    std::vector<float> lowerBounds = {-116.5, 0.0};
    std::vector<float> upperBounds = {116.5, 256.0}; 
    std::vector<std::string> legendPositions = {"L", "R"};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, legendPositions, 1);
    */

    /*
    std::vector<int> interactionTypes = {167, 165, 1, 163, 12, 0};
    std::vector<std::string> variableNames = {"LongestPfoDeltaChiSquaredPerHit", "LongestPfoUpDownDeltaChiSquaredPerHit"};
    std::vector<std::string> titleNames = {"PFO #Delta#chi^{2}_{FB}/N", "PFO PFO #Delta#chi^{2}_{DU}/N"};
    std::vector<float> lowerBounds = {-5.0, -5.0};
    std::vector<float> upperBounds = {5.0, 5.0}; 
    std::vector<std::string> legendPositions = {"R", "R"};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, legendPositions, 1);
    */

    /*
    std::vector<int> interactionTypes = {167, 165, 3, 12, 0};
    std::vector<std::string> variableNames = {"LongestPfoUpDownDeltaChiSquaredPerHit"};
    std::vector<std::string> titleNames = {"PFO #Delta#chi^{2}_{DU}/N"};
    std::vector<float> lowerBounds = {-5.0};
    std::vector<float> upperBounds = {5.0}; 
    std::vector<std::string> legendPositions = {"R"};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, legendPositions, 1, true);
    */

    /*
    std::vector<int> interactionTypes = {23, 100, 1, 18};
    std::vector<std::string> variableNames = {"ShortestPfoFitChargeRangeOverLength", "ShortestPfoMinSumQW", "ShortestPfoMaxSumQW", "ShortestPfoSumQWRatio", "ShortestPfoMeanQW"};
    std::vector<std::string> titleNames = {"PFO Q_{range}/L_{2D}", "PFO Q_{min}", "PFO Q_{max}", "PFO Q_{ratio}", "PFO <Q_{fit}>"};
    std::vector<float> lowerBounds = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<float> upperBounds = {0.2, 20.0, 30.0, 3.0, 10.0}; 
    std::vector<std::string> legendPositions = {"R", "R", "R", "R", "R"};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, legendPositions, 2, false);
    */

    std::vector<int> interactionTypes = {167, 165, 1, 163, 12, 0};
    std::vector<std::string> variableNames = {"LongestPfoPolarAngle", "NeutrinoMomentumY", "LongestPfoMinSumQW", "NeutrinoMomentumZ", "LongestPfoDeltaChiSquaredPerHit"};
    std::vector<std::string> titleNames = {"Polar Angle #theta", "p_{#nu,y}", "PFO Q_{min}", "p_{#nu,z}", "PFO #Delta#chi^{2}_{FB}/N"};
    std::vector<float> lowerBounds = {0.0, 0.0, 0.0, 0.0, -5.0};
    std::vector<float> upperBounds = {3.14, 1.0, 30.0, 1.0, 5.0}; 
    std::vector<std::string> legendPositions = {"R", "L", "R", "L", "R"};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, legendPositions, 1, true);
}

//------------------------------------------------------------------------------------------------------------------------------------------
