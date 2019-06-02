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

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/eventselection_III/0.root");
TTree *t1 = (TTree*)f1->Get("EventSelection");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/fullreco_bdtresponse.root");
TTree *t2 = (TTree*)f2->Get("EventSelection");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/fullreco_bdtresponse.root");
TTree *t3 = (TTree*)f3->Get("EventSelection");

TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/fullreco_bdtresponse.root");
TTree *t4 = (TTree*)f4->Get("EventSelection");

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStackHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, int targetMultiplicity, int targetInteractionType, bool onlyContained, std::vector<std::pair<std::string, int>> &filterValues, std::vector<std::pair<std::string, float>> &bdtVector, std::string typeName)
{
    int intVariable;
    float floatVariable;
    int interactionType;
    int particleMultiplicity;
    int nCosmics;
    int nuRecoContained;
    float bdtResponse[5];

    if (typeName == "int")
        pTree->SetBranchAddress(variableName.c_str(), &intVariable);
    else if (typeName == "float")
        pTree->SetBranchAddress(variableName.c_str(), &floatVariable);
    else
        std::cout << "UNKNOWN TYPENAME" << std::endl;

    pTree->SetBranchAddress("ModifiedInteractionType", &interactionType);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmics);
    pTree->SetBranchAddress("NuRecoContained", &nuRecoContained);

    //filter values
    int filterVariables[99], filterVariableTargetValues[99], filterVariableCounter(0);

    for (int i = 0; i < static_cast<int>(filterValues.size()); ++i)
    {   
        pTree->SetBranchAddress(filterValues.at(i).first.c_str(), &filterVariables[i]);
        filterVariableTargetValues[i] = filterValues.at(i).second;
        ++filterVariableCounter;
    }   

    //Create a vector of BDTs for which the response branch can actually be found
    std::vector<std::pair<std::string, float>> validBdtVector;    

    int i(0);
    for (const auto &pair : bdtVector)
    {   
        std::string bdtBranchName(pair.first);
        static TString invalidBranch(bdtBranchName.c_str());
        TBranch* reponseBranch = static_cast<TBranch*>(pTree->GetListOfBranches()->FindObject(invalidBranch));

        if (reponseBranch) 
        {   
            pTree->SetBranchAddress(bdtBranchName.c_str(), &bdtResponse[i]);
            validBdtVector.emplace_back(pair.first, pair.second);
            ++i;
        }   
    }   

    //Report which BDTs are used
    std::cout << "Using BDT cuts: ";

    for (const auto &pair : validBdtVector)
        std::cout << pair.first << " ";

    std::cout << endl;

    for (int i = 0; i < pTree->GetEntries(); i++)
    {   
        pTree->GetEntry(i);

        if (onlyContained && nuRecoContained == 0)
            continue;

        if (interactionType != targetInteractionType || particleMultiplicity != targetMultiplicity)
            continue;

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

        bool passedBdtCuts(true);

        for (unsigned int i = 0; i < validBdtVector.size(); ++i)
        {   
            if (bdtResponse[i] < validBdtVector.at(i).second)
                passedBdtCuts = false;
        }   

        if ((!passedBdtCuts) || skipThisEvent)
            continue;

        if (typeName == "int")
            pHistogram->Fill(intVariable);
        else if (typeName == "float")
            pHistogram->Fill(floatVariable);       
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStack(TTree* pTree, std::string variableName, std::vector<int> &targetInteractionTypes, int targetNumberParticles, bool onlyContained, std::vector<std::pair<std::string, int>> &filterValues, std::vector<std::pair<std::string, float>> &bdtVector, int nBins, float minRange, float maxRange, std::string xTitle, std::string yTitle, std::string plotTitle, std::string legendPosition, std::string typeName)
{
    THStack *pStack = new THStack("pStack","");
    std::vector<TH1F*> histogramVector;

    //TLegend* pLegend = new TLegend(0.68,0.72,0.98,0.92);
    TLegend* pLegend = new TLegend(0.63,0.68,0.88,0.88);
    TLegend* pLeftLegend = new TLegend(0.12,0.68,0.37,0.88);

    if (legendPosition == "L")
        pLegend = pLeftLegend;

    for (const auto &interactionType : targetInteractionTypes)
    {
        std::string histogramName("pHistogram" + std::to_string(interactionType));
        TH1F *pHistogram = new TH1F(histogramName.c_str(),"", nBins, minRange, maxRange);
        pHistogram->SetFillStyle(3004);

        CreateStackHistogram(pTree, pHistogram, variableName, targetNumberParticles, interactionType, onlyContained, filterValues, bdtVector, typeName);

        //custom legend entries
        if (interactionType == 998)
            pLegend->AddEntry(pHistogram, "CCQEL_MU_1_CR");
        else if (interactionType == 999)
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

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Stacked_Plots/N" + std::to_string(targetNumberParticles) + "_" + variableName + "_Stack_BDT.pdf");
    canvas->SaveAs(savePath.c_str());
    delete canvas;

    for (const auto pHistogram : histogramVector)
       delete pHistogram; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStackedPlots(TTree* pTree, std::vector<int> &interactionTypes, std::vector<std::string> &variableNames, std::vector<std::string> &titleNames, std::vector<float> lowerBounds, std::vector<float> upperBounds, std::vector<std::string> legendPositions, int particleMultiplicity, bool onlyContained, std::vector<std::pair<std::string, int>> &filterValues, std::vector<std::pair<std::string, float>> &bdtVector, std::vector<std::string> &typeNames)
{
    gStyle->SetPalette(kPastel);

    for (int i = 0; i < static_cast<int>(variableNames.size()); ++i)
        CreateStack(pTree, variableNames.at(i), interactionTypes, particleMultiplicity, onlyContained, filterValues, bdtVector, 100, lowerBounds.at(i), upperBounds.at(i), titleNames.at(i), "Number of Events", "Distribution of " + titleNames.at(i) + " for N=" + std::to_string(particleMultiplicity), legendPositions.at(i), typeNames.at(i));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Stacked_Plots_BDT(void)
{
    std::vector<int> interactionTypes = {182, 181, 0, 171, 1, 12}; 
    std::vector<std::string> variableNames = {"LongestPfoMinChiSquaredPerHit3D", "LongestPfoPolarAngle", "NeutrinoMomentumY", "LongestPfoMinSumQW", "NeutrinoMomentumZ", "LongestPfoDeltaChiSquaredPerHit", "LongestPfoRecoLowestTenCmTotalCharge"};
    std::vector<std::string> titleNames = {"#chi^{2}_{p}", "Polar Angle #theta", "p_{#nu,y}", "PFO Q_{min}", "p_{#nu,z}", "PFO #Delta#chi^{2}_{FB}/N", "Q_{low,10cm} (MeV)"};
    std::vector<float> lowerBounds = {0.0, 0.0, 0.0, 0.0, 0.0, -5.0, 0.0};
    std::vector<float> upperBounds = {100.0, 3.14, 1.0, 30.0, 1.0, 5.0, 150.0}; 
    std::vector<std::string> legendPositions = {"R", "R", "L", "R", "L", "R", "R"};
    std::vector<std::pair<std::string, int>> filterValues = {};
    std::vector<std::pair<std::string, float>> bdtVector = {std::make_pair("TMVA_BDT_AllVariables_N2_Contained", -0.05)};
    std::vector<std::string> typeNames = {"float", "float", "float", "float", "float", "float", "float"};

    CreateStackedPlots(t1, interactionTypes, variableNames, titleNames, lowerBounds, upperBounds, legendPositions, 1, true, filterValues, bdtVector, typeNames);
}

//------------------------------------------------------------------------------------------------------------------------------------------
