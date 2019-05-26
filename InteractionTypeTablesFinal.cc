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

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/eventselection/0.root");
TTree *t1 = (TTree*)f1->Get("EventSelection");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/eventselection/1.root");
TTree *t2 = (TTree*)f2->Get("EventSelection");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/eventselection/2.root");
TTree *t3 = (TTree*)f3->Get("EventSelection");

TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/eventselection/3.root");
TTree *t4 = (TTree*)f4->Get("EventSelection");

TFile *f5 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/fullreco_bdtresponse.root");
TTree *t5 = (TTree*)f5->Get("EventSelection");

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

void FormattedPrintout(int printoutWidth, std::string title, int denominator1, int selectedCount, int selectedCorrectCount, int denominator2, int totalCount, int totalCorrectCount)
{
    title += ": ";
    std::stringstream selectedStream, totalStream;
    selectedStream << std::fixed << setprecision(2) << selectedCount << " [" << selectedCorrectCount << "] (" << 100.0 * static_cast<float>(selectedCount)/denominator1 << "%)";
    totalStream << std::fixed << setprecision(2) << totalCount << " [" << totalCorrectCount << "] (" << 100.0 * static_cast<float>(totalCount)/denominator2 << "%)";

    std::cout << setw(printoutWidth) << left 
              << title 
              << setw(printoutWidth) << left
              << selectedStream.str() 
              << setw(printoutWidth) << left
              << totalStream.str() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateInteractionTypeTable(TTree* pTree, int targetMultiplicity, bool containedEventsOnly, std::string containmentDefinition, bool applyPidCut, std::vector<std::pair<std::string, float>> &bdtVector, std::pair<bool, int> filterEventClass)
{
    //containedEventsOnly: whether to only use events contained within the (bipartite) containment volume
    //bdtVector consists of std::pair objects of std::string bdtBranchName and float bdtCut
    //bdtBranchName defines the name of the branch containing the BDT response. Code checks whether branch exists. Leave blank to not use a BDT.
    //bdtCut The cut on the BDT response to use: >= bdtCut is signal

    std::map<std::string, int> eventTypeCount, correctEventTypeCount;
    std::map<std::string, int> totalInteractionTypeCountMap, totalCorrectInteractionTypeCountMap;

    int nuRecoContained;
    int interactionType;
    int trueInteractionType;
    int modifiedInteractionType;
    int nuanceCode;
    int particleMultiplicity;
    int nCosmicRays;
    int signal;
    int correctlyReconstructed;
    float longestPfoProtonChiSquared; 
    float shortestPfoProtonChiSquared; 
    int fileIdentifier;
    int eventNumber;
    float bdtResponse[5];

    pTree->SetBranchAddress(containmentDefinition.c_str(), &nuRecoContained);
    pTree->SetBranchAddress("TrueInteractionType", &trueInteractionType);
    pTree->SetBranchAddress("ModifiedInteractionType", &modifiedInteractionType);
    pTree->SetBranchAddress("NeutrinoNuanceCode", &nuanceCode);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmicRays);
    pTree->SetBranchAddress("LongestPfoMinChiSquaredPerHit3D", &longestPfoProtonChiSquared);

    if (targetMultiplicity > 1)
        pTree->SetBranchAddress("ShortestPfoMinChiSquaredPerHit3D", &shortestPfoProtonChiSquared);

    pTree->SetBranchAddress("Signal", &signal);
    pTree->SetBranchAddress("CorrectlyReconstructed", &correctlyReconstructed);
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier);
    pTree->SetBranchAddress("EventNumber", &eventNumber);


    std::cout << "******************************************************************************************" << std::endl;

    if (filterEventClass.first == true)
        std::cout << ">>> Filtering to modified interaction type " << ToString(static_cast<InteractionType>(filterEventClass.second)) << std::endl;

    int totalCount(0), totalCorrectCount(0), totalSignalCount(0), totalCorrectSignalCount(0);

    //Fill total count maps
    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterEventClass.first == true)
            interactionType = trueInteractionType;
        else 
            interactionType = modifiedInteractionType;

        if (containedEventsOnly && nuRecoContained == 0)
            continue;

        if (filterEventClass.first == true && modifiedInteractionType != filterEventClass.second)
            continue;

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));
        std::string nCosmicRaysString(std::to_string(nCosmicRays)+"_CR");

        //exclude _CR suffix if  COSMIC_RAY and only add it if nu reco is contaminated 
        if (interactionType != 181 && nCosmicRays != 0 && nCosmicRays != particleMultiplicity)
            interactionTypeString += "_" + nCosmicRaysString;

        //fill counts
        ++totalCount;
        if (correctlyReconstructed == 1) ++totalCorrectCount;

        if (signal == 1) 
        {
            ++totalSignalCount;
            if (correctlyReconstructed == 1) ++totalCorrectSignalCount;
        }

        //fill maps
        totalInteractionTypeCountMap[interactionTypeString]++;
        if (correctlyReconstructed == 1) totalCorrectInteractionTypeCountMap[interactionTypeString]++;
    }

    float percentageCutoff(1.0);
    std::stringstream percentageStream;
    percentageStream << std::fixed << setprecision(1) << percentageCutoff; 

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
    std::cout << "******************************************************************************************" << std::endl;

    //Fill eventTypeCount map 
    int selectedCount(0), correctSelectedCount(0), selectedSignalCount(0), correctSelectedSignalCount(0);
    int variables[5];

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterEventClass.first == true)
            interactionType = trueInteractionType;
        else 
            interactionType = modifiedInteractionType;

        if (particleMultiplicity != targetMultiplicity)
            continue;

        if (containedEventsOnly && nuRecoContained == 0)
            continue;

        if (filterEventClass.first == true && modifiedInteractionType != filterEventClass.second)
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
    
        //exclude _CR suffix if  COSMIC_RAY and only add it if nu reco is contaminated 
        if (interactionType != 181 && nCosmicRays != 0 && nCosmicRays != particleMultiplicity)
            interactionTypeString += "_" + nCosmicRaysString;

        //fill counts
        ++selectedCount;
        if (correctlyReconstructed == 1) ++correctSelectedCount;

        if (signal == 1) 
        {
            ++selectedSignalCount;
            if (correctlyReconstructed == 1) ++correctSelectedSignalCount;
        }

        //fill maps
        eventTypeCount[interactionTypeString]++;
        
        if (correctlyReconstructed == 1)
            correctEventTypeCount[interactionTypeString]++;
    }

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;

    //sort map entries by count by filling sorting vector
    std::vector<std::pair<std::string, int>> pairs;

    for (const auto &pair : eventTypeCount)
        pairs.push_back(pair);

    sort(pairs.begin(), pairs.end(), [=](std::pair<std::string, int>& a, std::pair<std::string, int>& b) {return a.second > b.second;});

    int printoutWidth(30);

    std::cout << setw(printoutWidth) << left
              << "Interaction Type"
              << setw(printoutWidth) << left 
              << "Selected Events"
              << setw(printoutWidth) << left
              << "Total Events" << std::endl;

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;

    int otherCount(0), otherCorrectCount(0), otherTotalCount(0), otherTotalCorrectCount(0);
 
    for (std::pair<std::string, int> pair : pairs)
    {
        if (100.0 * (float)pair.second/selectedCount >= percentageCutoff)
        {
            std::string title(pair.first + " (" + std::to_string(FromString(pair.first)) + ")");
            int correctChannelCount(correctEventTypeCount.count(pair.first) != 0 ? correctEventTypeCount.at(pair.first) : 0);
            int totalChannelCount(totalInteractionTypeCountMap.at(pair.first) != 0 ? totalInteractionTypeCountMap.at(pair.first) : 0);
            int totalCorrectChannelCount(totalCorrectInteractionTypeCountMap.count(pair.first) != 0 ? totalCorrectInteractionTypeCountMap.at(pair.first) : 0);
            FormattedPrintout(printoutWidth, title, selectedCount, pair.second, correctChannelCount, totalCount, totalChannelCount, totalCorrectChannelCount);

            /*
            std::string interactionTypeCounts(pair.first + " (" + std::to_string(FromString(pair.first)) + ")");
            std::string selectedEvents(std::to_string(pair.second) + " [" + std::to_string(correctEventTypeCount.count(pair.first) != 0 ? correctEventTypeCount.at(pair.first) : 0) + "] (" + std::to_string(100.0 * (float)pair.second/selectedCount) + "%)");
            std::string totalEvents(std::to_string(totalInteractionTypeCountMap.at(pair.first) != 0 ? totalInteractionTypeCountMap.at(pair.first) : 0) + " [" + std::to_string(totalCorrectInteractionTypeCountMap.count(pair.first) != 0 ? totalCorrectInteractionTypeCountMap.at(pair.first) : 0) + "]");

            std::stringstream selectedEventsStream;
            selectedEventsStream << std::fixed << setprecision(2) << pair.second << " [" << (correctEventTypeCount.count(pair.first) != 0 ? correctEventTypeCount.at(pair.first) : 0) << "] (" << 100.0 * (float)pair.second/selectedCount << "%)";

            std::cout << setw(printoutWidth) << left 
                      << interactionTypeCounts
                      << setw(printoutWidth) << left
                      << selectedEventsStream.str()
                      << setw(printoutWidth) << left
                      << totalEvents << std::endl;
            */
        }
        else
        {
            otherCount += pair.second; 
            otherCorrectCount += (correctEventTypeCount.count(pair.first) != 0 ? correctEventTypeCount.at(pair.first) : 0);
            otherTotalCount += (totalInteractionTypeCountMap.count(pair.first) != 0 ? totalInteractionTypeCountMap.at(pair.first) : 0);
            otherTotalCorrectCount += (totalCorrectInteractionTypeCountMap.count(pair.first) != 0 ? totalCorrectInteractionTypeCountMap.at(pair.first) : 0);
        }
    }

    FormattedPrintout(printoutWidth, "OTHER (< " + percentageStream.str() + "%)", selectedCount, otherCount, otherCorrectCount, totalCount, otherTotalCount, otherTotalCorrectCount);
    std::cout << ".........................................................................................." << std::endl;
    FormattedPrintout(printoutWidth, "SIGNAL", selectedCount, selectedSignalCount, correctSelectedSignalCount, totalCount, totalSignalCount, totalCorrectSignalCount);
    FormattedPrintout(printoutWidth, "TOTAL", selectedCount,  selectedCount, correctSelectedCount, totalCount, totalCount, totalCorrectCount);

    /*
    std::string otherInteractionTypeCounts("OTHER (<" + std::to_string(percentageCutoff) + "%)"); 
    std::string otherSelectedEvents(std::to_string(otherCount) + " [" + std::to_string(otherCorrectCount) + "] (" + std::to_string(100.0 * (float)otherCount/selectedCount) + "%)");
    std::string otherTotalEvents(std::to_string(otherTotalCount) + " [" + std::to_string(otherTotalCorrectCount) + "]");

    std::stringstream otherSelectedEventsStream;
    otherSelectedEventsStream << std::fixed << setprecision(2) << otherCount << " [" << otherCorrectCount << "] (" << 100.0 * (float)otherCount/selectedCount << "%)";

    std::cout << setw(printoutWidth) << left 
              << otherInteractionTypeCounts
              << setw(printoutWidth) << left
              << otherSelectedEventsStream.str() 
              << setw(printoutWidth) << left
              << otherTotalEvents << std::endl;

    std::cout << ".........................................................................................." << std::endl;
    //std::cout << "SIGNAL: " << selectedSignalCount << " [" << correctSelectedSignalCount << "]" << std::endl;

    std::string signalInteractionTypeCounts("SIGNAL: "); 
    std::stringstream signalSelectedEventsStream;
    signalSelectedEventsStream << std::fixed << setprecision(2) << selectedSignalCount << " (" << 100 * static_cast<float>(selectedSignalCount)/selectedCount << "%) [" << correctSelectedSignalCount << "]"; 
    std::stringstream signalTotalEventsStream;
    signalTotalEventsStream << std::fixed << setprecision(2) << totalSignalCount << " (" << 100 * static_cast<float>(totalSignalCount)/totalCount << "%) [" << totalCorrectSignalCount << "]"; 

    std::cout << setw(printoutWidth) << left 
              << signalInteractionTypeCounts
              << setw(printoutWidth) << left
              << signalSelectedEventsStream.str() 
              << setw(printoutWidth) << left
              << signalTotalEventsStream.str() << std::endl;

    std::string totalInteractionTypeCounts("TOTAL: "); 
    std::string totalSelectedEvents(std::to_string(selectedCount) + " (100%) [" + std::to_string(correctSelectedCount) + "]"); 
    std::string totalTotalEvents(std::to_string(totalCount) + " (100%) [" + std::to_string(totalCorrectCount) + "]"); 

    std::cout << setw(printoutWidth) << left 
              << totalInteractionTypeCounts
              << setw(printoutWidth) << left
              << totalSelectedEvents 
              << setw(printoutWidth) << left
              << totalTotalEvents << std::endl;
    */

    std::cout << ".........................................................................................." << std::endl;
    std::cout << "PURITY: " << correctSelectedSignalCount << "/" << selectedCount << " (" << 100.0 * (float)correctSelectedSignalCount/selectedCount<< "%)" << std::endl;
    std::cout << "SELECTION CUTS EFFICIENCY: " << correctSelectedSignalCount << "/" << totalCorrectSignalCount << " (" << 100.0 * (float)correctSelectedSignalCount/totalCorrectSignalCount<< "%)" << std::endl;
    std::cout << "OVERALL SELECTION EFFICIENCY: " << correctSelectedSignalCount << "/" << totalSignalCount << " (" << 100.0 * (float)correctSelectedSignalCount/totalSignalCount << "%)" << std::endl;

    //std::cout << "OTHER (< " << percentageCutoff << "%): " << otherCount << " (" << 100.0 * (float)otherCount/selectedCount << "%)" << std::endl;
    //std::cout << "TOTAL: " << selectedCount << " (100%)" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OtherInteractionAnalysis(TTree* pTree, int targetMultiplicity, bool onlyContained)
{
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "               OTHER_INTERACTION             " << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    int nuRecoContained;
    int nuanceCode;
    int particleMultiplicity;
    int interactionType;
    int nCosmicRays; 
    int trueNumberMuons;
    int trueNumberProtons;
    int trueNumberPhotons;
    int signal;
    int correctlyReconstructed;

    pTree->SetBranchAddress("NuRecoContained", &nuRecoContained);
    pTree->SetBranchAddress("NeutrinoNuanceCode", &nuanceCode);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("ModifiedInteractionType", &interactionType);
    pTree->SetBranchAddress("RecoNeutrinoNumberCosmicRays", &nCosmicRays);
    pTree->SetBranchAddress("NumberTrueProtons", &trueNumberProtons);
    pTree->SetBranchAddress("NumberTrueMuons", &trueNumberMuons);
    pTree->SetBranchAddress("NumberTruePhotons", &trueNumberPhotons);
    pTree->SetBranchAddress("Signal", &signal);
    pTree->SetBranchAddress("CorrectlyReconstructed", &correctlyReconstructed);

    std::map<int, int> nuanceCodeCount, trueMuonCount, trueProtonCount;
    std::map<int, std::string> nuanceCodeNameMapping = {{1091, "CCDIS"}, {1092, "NCDIS"}, {5001, "NCRES"}, {5101, "CCRES"}, {1002, "NCQEL"}, {1004, "CCRES_N_PIZERO"}, {1005, "CCRES_N_PIPLUS"}, {1098, "NEUTRINO_E_ELASTIC"}};
    int totalCount(0), nSingleMuonNoProton(0), nSingleMuonSingleProton(0);
    int signalCount(0), correctSignalCount(0);

    for (int i = 0; i < pTree->GetEntries(); i++) 
    {    
        pTree->GetEntry(i);

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));
        std::string nCosmicRaysString(std::to_string(nCosmicRays)+"_CR");
    
        if (interactionType != 181 && interactionType != 183 && nCosmicRays != 0)
            interactionTypeString += "_" + nCosmicRaysString;

        if (particleMultiplicity != targetMultiplicity || interactionTypeString != "OTHER_INTERACTION" || (onlyContained && nuRecoContained != 1)) 
            continue;

        /*
        if (nuanceCode == 5001 || nuanceCode == 5101)
        {
            std::cout << "nuanceCode: " << nuanceCode << std::endl;
            std::cout << "Muons: " << trueNumberMuons << std::endl;
            std::cout << "Protons: " << trueNumberProtons << std::endl;
            std::cout << "Photons: " << trueNumberPhotons << std::endl;
        }
        */

        if (signal == 1) ++signalCount;
        if (signal == 1 && correctlyReconstructed == 1) ++correctSignalCount;

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

    std::vector<std::pair<int, int>> pairs;

    for (const auto &pair : nuanceCodeCount)
        pairs.push_back(pair);

    sort(pairs.begin(), pairs.end(), [=](std::pair<int, int>& a, std::pair<int, int>& b) {return a.second > b.second;});

    int total(0);

    for (std::pair<int, int> pair : pairs)
    {
        std::string nuanceCodeName(nuanceCodeNameMapping.find(pair.first) != nuanceCodeNameMapping.end() ? nuanceCodeNameMapping.at(pair.first) : "UNKNOWN");
        std::cout << "Nuance " << pair.first << " " << nuanceCodeName  << ": " << pair.second << " (" << std::fixed << setprecision(2) << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;
        ++total;
    }

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "TOTAL: " << total << " (100%)" << std::endl;
    std::cout << "SIGNAL: " << signalCount << " (" << 100.0 * (float)signalCount/totalCount << "%)" << std::endl;
    std::cout << "CORRECT SIGNAL: " << correctSignalCount << " (" << 100.0 * (float)correctSignalCount/totalCount << "%)" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NoReconstructableAnalysis(TTree* pTree, int targetMultiplicity, bool onlyContained)
{
    //Energy distribution no_reco
    std::vector<std::pair<std::string, int>> filterValues = {std::make_pair("ModifiedInteractionType", 182), std::make_pair("RecoNeutrinoNumberAssociatedParticles", targetMultiplicity), std::make_pair("NuRecoConatained", onlyContained)};

    TH1F *pNoReconstructableEnergyDistribution = new TH1F("pNoReconstructableEnergyDistribution","", 100, 0.0, 500.0);
    CreateFilteredHistogram(t1, pNoReconstructableEnergyDistribution, "TrueNeutrinoEnergy", filterValues, "float")

    //number of neutrons
    TH1F *pNoReconstructableNumberNeutrons = new TH1F("pNoReconstructableNumberNeutrons","", 6, 0, 5);
    CreateFilteredHistogram(t1, pNoReconstructableNumberNeutrons, "TrueNumberNeutrons", filterValues, "int")
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareCheatingConfigurations(int targetMultiplicity, bool containedOnly)
{
    std::vector<std::pair<std::string, float>> emptyBdtVector = {};

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t1, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector, std::make_pair(false, 181));

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t2, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector, std::make_pair(false, 181));

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t3, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector, std::make_pair(false, 181));

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Slicing, Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t4, targetMultiplicity, containedOnly, "NuRecoContained", false, emptyBdtVector, std::make_pair(false, 181));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareContainmentDefinitions(TTree* pTree, int targetMultiplicity, std::vector<std::string> &containmentDefinitions)
{
    std::vector<std::pair<std::string, float>> emptyBdtVector = {};

    //Table without containment cut for reference
    CreateInteractionTypeTable(pTree, targetMultiplicity, false, "NuRecoContained", false, emptyBdtVector, std::make_pair(false, 181));

    for (std::string &containmentDefinition : containmentDefinitions)
    {
        std::cout << "DEFINITION: " << containmentDefinition << std::endl;
        CreateInteractionTypeTable(pTree, targetMultiplicity, true, containmentDefinition, false, emptyBdtVector, std::make_pair(false, 181));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void InteractionTypeTablesFinal(void)
{
    int targetMultiplicity(1);

    //std::vector<std::string> containmentDefinitions = {"NuRecoContained", "NuRecoEndpointsContained", "NuRecoVertexContained"};
    //CompareContainmentDefinitions(t1, targetMultiplicity, containmentDefinitions);
    
    OtherInteractionAnalysis(t1, targetMultiplicity, true);
    //NoReconstructableAnalysis(t1, targetMultiplicity, true);

    CompareCheatingConfigurations(targetMultiplicity, true);

    std::vector<std::pair<std::string, float>> emptyBdtVector = {};
    CreateInteractionTypeTable(t1, targetMultiplicity, true, "NuRecoContained", false, emptyBdtVector, std::make_pair(true, 181));
}

//------------------------------------------------------------------------------------------------------------------------------------------
