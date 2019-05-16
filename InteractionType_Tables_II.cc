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
#include <utility>

//------------------------------------------------------------------------------------------------------------------------------------------

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/0.root");
TTree *t1 = (TTree*)f1->Get("EventSelection");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/1.root");
TTree *t2 = (TTree*)f2->Get("EventSelection");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/2.root");
TTree *t3 = (TTree*)f3->Get("EventSelection");

TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/3.root");
TTree *t4 = (TTree*)f4->Get("EventSelection");

TFile *f5 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/fullreco_bdtresponse.root");
TTree *t5 = (TTree*)f5->Get("EventSelection");

//------------------------------------------------------------------------------------------------------------------------------------------

int CalculateNumberMCSignalEvents(TTree* pTree, int targetMultiplicity, bool containedOnly, int variables[5])
{
    int &nuTrueContained = variables[1];
    int &nPrimaries = variables[2];
    int &nMuons = variables[3];
    int &nProtons = variables[4];

    pTree->SetBranchAddress("NuTrueContained", &nuTrueContained);
    pTree->SetBranchAddress("NumberTruePrimaries", &nPrimaries);
    pTree->SetBranchAddress("NumberTrueMuons", &nMuons);
    pTree->SetBranchAddress("NumberTrueProtons", &nProtons);

    int nTrueSignalEvents(0);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if ((containedOnly && nuTrueContained == 0) || (nPrimaries != targetMultiplicity))
            continue;

        if ((targetMultiplicity == 1 && nMuons == 1) || (targetMultiplicity == 2 && nMuons == 1 && nProtons == 1))
            ++nTrueSignalEvents;
    }

    return nTrueSignalEvents;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PassesPidCut(int targetMultiplicity, float shortestPfoProtonChiSquared, float longestPfoProtonChiSquared, float cutValue)
{
    if (targetMultiplicity == 2) 
    {
        if (shortestPfoProtonChiSquared == -1.f)
            return true;

        bool shortestPfoProtonLike(shortestPfoProtonChiSquared <= cutValue);

        if (!shortestPfoProtonLike) 
            return false;
    }
    else if (targetMultiplicity == 1)
    {
        if (longestPfoProtonChiSquared == -1.f)
            return true;

        bool longestPfoProtonLike(longestPfoProtonChiSquared <= cutValue);

        if (longestPfoProtonLike)
            return false;
    }
    else
        return true;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateInteractionTypeTable(TTree* pTree, int targetMultiplicity, bool containedEventsOnly, std::string containmentDefinition, bool applyPidCut, std::vector<std::pair<std::string, float>> &bdtVector)
{
    //containedEventsOnly: whether to only use events contained within the (bipartite) containment volume
    //bdtVector consists of std::pair objects of std::string bdtBranchName and float bdtCut
    //bdtBranchName defines the name of the branch containing the BDT response. Code checks whether branch exists. Leave blank to not use a BDT.
    //bdtCut The cut on the BDT response to use: >= bdtCut is signal

    std::map<std::string, int> eventTypeCount;

    int nuRecoContained;
    int interactionType;
    int particleMultiplicity;
    int nCosmicRays;
    int signal;
    float longestPfoProtonChiSquared; 
    float shortestPfoProtonChiSquared; 
    int fileIdentifier;
    int eventNumber;
    float bdtResponse[5];

    pTree->SetBranchAddress(containmentDefinition.c_str(), &nuRecoContained);
    pTree->SetBranchAddress("ModifiedInteractionType", &interactionType);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmicRays);
    pTree->SetBranchAddress("LongestPfoMinChiSquaredPerHit3D", &longestPfoProtonChiSquared);

    if (targetMultiplicity > 1)
        pTree->SetBranchAddress("ShortestPfoMinChiSquaredPerHit3D", &shortestPfoProtonChiSquared);

    pTree->SetBranchAddress("Signal", &signal);
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier);
    pTree->SetBranchAddress("EventNumber", &eventNumber);

    float percentageCutoff(2.0);

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

    std::cout << "PID: " << applyPidCut << std::endl;

    //Fill eventTypeCount map 
    int signalCount(0); 
    int variables[5];
    int trueSignalCount(CalculateNumberMCSignalEvents(pTree, targetMultiplicity, containedEventsOnly, variables));

    std::cout << "trueSignalCount: " << trueSignalCount << std::endl;

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (particleMultiplicity != targetMultiplicity)
            continue;

        if (containedEventsOnly && nuRecoContained == 0)
            continue;

        if (applyPidCut && shortestPfoProtonChiSquared >= 20.0)
            continue;

        bool passedBdtCuts(true);

        for (unsigned int i = 0; i < validBdtVector.size(); ++i)
        {
            if (bdtResponse[i] < validBdtVector.at(i).second)
                passedBdtCuts = false;
        }

        if (!passedBdtCuts)
            continue;

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));
        std::string nCosmicRaysString(std::to_string(nCosmicRays)+"_CR");
    
        if (interactionType != 165 && interactionType != 167 && nCosmicRays != 0)
            interactionTypeString += "_" + nCosmicRaysString;

        if (signal == 1)
            ++signalCount;

        eventTypeCount[interactionTypeString]++;
    }

    int totalCount(0), otherCount(0);

    for (const auto &pair : eventTypeCount)
        totalCount += pair.second;

    std::cout << "---------------------------------------------" << std::endl;

    typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;
    Comparator compFunctor = [](std::pair<std::string, int> elem1 ,std::pair<std::string, int> elem2) { return elem1.second > elem2.second; };
    std::set<std::pair<std::string, int>, Comparator> interactionTypeSet(eventTypeCount.begin(), eventTypeCount.end(), compFunctor);
 
    for (std::pair<std::string, int> pair : interactionTypeSet)
    {
        if (100.0 * (float)pair.second/totalCount >= percentageCutoff)
            std::cout << pair.first << " (" << FromString(pair.first) << ") : " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;
        else
            otherCount += pair.second; 
    }

    std::cout << "............................................." << std::endl;
    std::cout << "PURITY: " << signalCount << "/" << totalCount << " (" << 100.0 * (float)signalCount/totalCount << "%)" << std::endl;
    std::cout << "SELECTION EFFICIENCY: " << signalCount << "/" << trueSignalCount << " (" << 100.0 * (float)signalCount/trueSignalCount<< "%)" << std::endl;
    std::cout << "OTHER (< " << percentageCutoff << "%): " << otherCount << " (" << 100.0 * (float)otherCount/totalCount << "%)" << std::endl;
    std::cout << "TOTAL: " << totalCount << " (100%)" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << std::setprecision(6);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SubTableAnalysis(TTree* pTree, int targetMultiplicity, int targetModifiedInteractionType, bool onlyContained)
{
    std::map<std::string, int> eventTypeCount;

    int nuRecoContained;
    int interactionType;
    int modifiedInteractionType;
    int particleMultiplicity; 
    int nCosmicRays;
    int taggingFailure;

    pTree->SetBranchAddress("NuRecoContained", &nuRecoContained);
    pTree->SetBranchAddress("TrueInteractionType", &interactionType);
    pTree->SetBranchAddress("ModifiedInteractionType", &modifiedInteractionType);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmicRays);
    pTree->SetBranchAddress("TaggingFailure", &taggingFailure);

    int nCosmicRayCategory(0), nTaggingFailure(0);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if ((particleMultiplicity != targetMultiplicity) || (onlyContained && nuRecoContained != 1) || (modifiedInteractionType != targetModifiedInteractionType))
            continue;

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));

        if (interactionType != 167 && nCosmicRays == targetMultiplicity)
        {
            ++nCosmicRayCategory;
            if (taggingFailure == 1) ++nTaggingFailure;
            eventTypeCount[interactionTypeString]++;
        }
    }

    std::cout << "nCosmicRayCategory: " << nCosmicRayCategory << std::endl;
    std::cout << "nTaggingFailure: " << nTaggingFailure << std::endl;

    int totalCount(0), otherCount(0);

    for (const auto &pair : eventTypeCount)
        totalCount += pair.second;

    std::cout << "---------------------------------------------" << std::endl;

    typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;
    Comparator compFunctor = [](std::pair<std::string, int> elem1 ,std::pair<std::string, int> elem2) { return elem1.second > elem2.second; };
    std::set<std::pair<std::string, int>, Comparator> interactionTypeSet(eventTypeCount.begin(), eventTypeCount.end(), compFunctor);
 
    for (std::pair<std::string, int> pair : interactionTypeSet)
    {
        if (100.0 * (float)pair.second/totalCount > 2)
            std::cout << pair.first << " (" << FromString(pair.first) << ") : " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;
        else
            otherCount += pair.second; 
    }

    std::cout << "............................................." << std::endl;
    std::cout << "OTHER: " << otherCount << " (" << 100.0 * (float)otherCount/totalCount << "%)" << std::endl;
    std::cout << "TOTAL: " << totalCount << " (100%)" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << std::setprecision(6);


    //Number of nu true hits in cosmic ray category
    TH1F* pTrueNuNumberHitsHistogram = new TH1F("pTrueNuNumberHitsHistogram", "", 100, 0, 100); 
    TH1F* pTrueNuNumberHitsInPFOsHistogram = new TH1F("pTrueNuNumberHitsInPFOsHistogram", "", 100, 0, 100); 
    TH1F* pTrueNuNumberHitsInRecoNuHistogram = new TH1F("pTrueNuNumberHitsInRecoNuHistogram", "", 100, 0, 100); 

    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("RecoNeutrinoNumberAssociatedParticles", targetMultiplicity), std::make_pair("RecoNeutrinoNumberCosmicRays", targetMultiplicity), std::make_pair("TaggingFailure", 0)};

    CreateFilteredHistogram(pTree, pTrueNuNumberHitsHistogram, "TrueNeutrinoNumberInducedHits", filterValues);
    CreateFilteredHistogram(pTree, pTrueNuNumberHitsInPFOsHistogram, "PfoAllTrueNeutrinoHits", filterValues);
    CreateFilteredHistogram(pTree, pTrueNuNumberHitsInRecoNuHistogram, "RecoNuAllTrueNeutrinoHits", filterValues);

    //Draw histograms
    std::vector<TH1F*> histogramVector = {pTrueNuNumberHitsHistogram, pTrueNuNumberHitsInPFOsHistogram, pTrueNuNumberHitsInRecoNuHistogram}; 
    std::vector<EColor> colourVector = {kRed, kBlue, kMagenta};
    std::vector<std::string> legendVector = {"Number of #nu_{true} hits in Event", "Number of #nu_{true} hits in PFOs", "Number of #nu_{true} hits in #nu_{reco}"};
    Draw(histogramVector, colourVector, legendVector, false, true, "true_nu_hits_distributions_CR_N_" + std::to_string(targetMultiplicity), "Number of #nu_{true} hits", "Number of Entries", "Distributions of number of #nu_{true} hits (N=" + std::to_string(targetMultiplicity) + " COSMIC_RAY Category)");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OtherInteractionAnalysis(TTree* pTree, int targetMultiplicity, bool onlyContained)
{
    std::cout << "---------------------------------------------" << std::endl;

    int nuRecoContained;
    int nuanceCode;
    int particleMultiplicity;
    int interactionType;
    int nCosmicRays; 
    int trueNumberMuons;
    int trueNumberProtons;

    pTree->SetBranchAddress("NuRecoContained", &nuRecoContained);
    pTree->SetBranchAddress("NeutrinoNuanceCode", &nuanceCode);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("ModifiedInteractionType", &interactionType);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmicRays);
    pTree->SetBranchAddress("NumberTrueProtons", &trueNumberProtons);
    pTree->SetBranchAddress("NumberTrueMuons", &trueNumberMuons);

    std::map<int, int> nuanceCodeCount, trueMuonCount, trueProtonCount;
    int totalCount(0), nSingleMuonNoProton(0), nSingleMuonSingleProton(0);

    for (int i = 0; i < pTree->GetEntries(); i++) 
    {    
        pTree->GetEntry(i);

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));
        std::string nCosmicRaysString(std::to_string(nCosmicRays)+"_CR");
    
        if (interactionType != 165 && interactionType != 167 && nCosmicRays != 0)
            interactionTypeString += "_" + nCosmicRaysString;

        if (particleMultiplicity != targetMultiplicity || interactionTypeString != "OTHER_INTERACTION" || (onlyContained && nuRecoContained != 1)) 
            continue;

        nuanceCodeCount[nuanceCode]++;

        if (nuanceCode == 1000)
        {
            trueMuonCount[trueNumberMuons]++;
            trueProtonCount[trueNumberProtons]++;

            if (trueNumberMuons == 1 && trueNumberProtons == 0) ++ nSingleMuonNoProton;
            if (trueNumberMuons == 1 && trueNumberProtons == 1) ++ nSingleMuonSingleProton;
        }

        ++totalCount;
    }    

    typedef std::function<bool(std::pair<int, int>, std::pair<int, int>)> Comparator;
    Comparator compFunctor = [](std::pair<int, int> elem1 ,std::pair<int, int> elem2) { return elem1.second > elem2.second; };
    std::set<std::pair<int, int>, Comparator> interactionTypeSet(nuanceCodeCount.begin(), nuanceCodeCount.end(), compFunctor);
    std::set<std::pair<int, int>, Comparator> muonSet(trueMuonCount.begin(), trueMuonCount.end(), compFunctor);
    std::set<std::pair<int, int>, Comparator> protonSet(trueProtonCount.begin(), trueProtonCount.end(), compFunctor);

    for (std::pair<int, int> pair : interactionTypeSet)
        std::cout << "Nuance " << pair.first << ": " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;

    std::cout << "Total: " << totalCount << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    std::cout << "Single muon no proton count: " << nSingleMuonNoProton << std::endl;
    std::cout << "Single muon single proton count: " << nSingleMuonSingleProton << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    /*
    for (std::pair<int, int> pair : muonSet)
        std::cout << "Muon count " << pair.first << ": " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;

    std::cout << "---------------------------------------------" << std::endl;

    for (std::pair<int, int> pair : protonSet)
        std::cout << "Proton count " << pair.first << ": " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;

    std::cout << "---------------------------------------------" << std::endl;
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareCheatingConfigurations(int targetMultiplicity, bool containedOnly)
{
    std::vector<std::pair<std::string, float>> emptyBdtVector = {};

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t1, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector);

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t2, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector);

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t3, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector);

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Slicing, Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t4, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareContainmentDefinitions(TTree* pTree, int targetMultiplicity, std::vector<std::string> &containmentDefinitions)
{
    std::vector<std::pair<std::string, float>> emptyBdtVector = {};

    //Table without containment cut for reference
    CreateInteractionTypeTable(pTree, targetMultiplicity, false, "NuRecoContained", false, emptyBdtVector);

    for (std::string &containmentDefinition : containmentDefinitions)
    {
        std::cout << "DEFINITION: " << containmentDefinition << std::endl;
        CreateInteractionTypeTable(pTree, targetMultiplicity, true, containmentDefinition, false, emptyBdtVector);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void InteractionType_Tables(void)
{
    int targetMultiplicity(1);

    //std::vector<std::string> containmentDefinitions = {"NuRecoContained", "NuRecoEndpointsContained", "NuRecoVertexContained"};
    //CompareContainmentDefinitions(t1, targetMultiplicity, containmentDefinitions);
    
    //OtherInteractionAnalysis(t1, targetMultiplicity, true);

    //CompareCheatingConfigurations(targetMultiplicity, true);

    SubTableAnalysis(t2, targetMultiplicity, 165, true);

    //std::vector<std::pair<std::string, float>> emptyBdtVector = {};
    //CreateInteractionTypeTable(t1, targetMultiplicity, true, "NuRecoContained", false, emptyBdtVector);

    /*
    std::vector<std::pair<std::string, float>> bdtVector0 = {};
    std::vector<std::pair<std::string, float>> bdtVector1 = {std::make_pair("BDT_PID_precut_N2", -0.1)};
    //std::vector<std::pair<std::string, float>> bdtVector1 = {std::make_pair("BDT_N2_Contained_CleanedVariables", -0.1)};
    std::vector<std::pair<std::string, float>> bdtVector2 = {std::make_pair("BDT_N2_Contained_CleanedVariables", -0.01)};

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t5, targetMultiplicity, true, true, bdtVector0);
    CreateInteractionTypeTable(t5, targetMultiplicity, true, true, bdtVector1);
    CreateInteractionTypeTable(t5, targetMultiplicity, true, false, bdtVector2);
    */

    /*
    std::vector<std::pair<std::string, float>> emptyBdtVector = {};
    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t5, targetMultiplicity, true, false, emptyBdtVector);

    std::vector<std::pair<std::string, float>> bdtVector0 = {std::make_pair("TMVA_BDT_AllVariables_N2_Contained", -0.05)};
    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t5, targetMultiplicity, true, false, bdtVector0);

    std::vector<std::pair<std::string, float>> bdtVector1 = {std::make_pair("TMVA_BDT_AllVariables_N2_Contained", -0.05)};
    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t5, targetMultiplicity, true, true, bdtVector1);
    */

    /*
    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    std::vector<std::pair<std::string, float>> bdtVector0 = {std::make_pair("TMVA_BDT_AllVariables_N1_Contained", -0.05)};
    CreateInteractionTypeTable(t5, targetMultiplicity, true, bdtVector0);

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    std::vector<std::pair<std::string, float>> bdtVector1 = {std::make_pair("TMVA_BDT_NoDirection_N1_Contained", -0.05)};
    CreateInteractionTypeTable(t5, targetMultiplicity, true, bdtVector1);
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------
