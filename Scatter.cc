#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <sstream>

//------------------------------------------------------------------------------------------------------------------------------------------

float GetCompletenessPurityIntersection(TTree* tree, int NLowerBound, int NUpperBound, std::string variable1, std::string variable2, std::string numberhits)
{
    std::cout << "*" << std::endl;
    int     MCForwards;
    int     MCIsClusterTwoParticles;
    int     NumberHits;
    float   beforeRemovalMinchiSquared;
    float   afterRemovalMinchiSquared;
    
    tree->SetBranchAddress("MCForwards",                                             &MCForwards);
    tree->SetBranchAddress("MCIsClusterTwoParticles",                                             &MCIsClusterTwoParticles);
    tree->SetBranchAddress(numberhits.c_str(),                                             &NumberHits);
    tree->SetBranchAddress(variable1.c_str(),                                        &beforeRemovalMinchiSquared);
    tree->SetBranchAddress(variable2.c_str(),                                         &afterRemovalMinchiSquared);

    int nGoodCuts(0), nBadCuts(0);

    for (int i = 0; i <= tree->GetEntries(); i++)
    {   
        tree->GetEntry(i);

        if (beforeRemovalMinchiSquared == 0)
            continue;

        if (!(NumberHits > NLowerBound && NumberHits <= NUpperBound))
            continue;

        if (MCIsClusterTwoParticles == 0)
            nBadCuts++;

        if (MCIsClusterTwoParticles == 1)
            nGoodCuts++;
    }

    float maxCut(10.0), cutStep(0.1);
    int n(0);

    TGraph *cutValueWithPurity = new TGraph(maxCut/cutStep + 1); 
    TGraph *cutValueWithCompleteness = new TGraph(maxCut/cutStep + 1); 

    for (float cut = 0.0; cut < maxCut; cut += cutStep)
    {   
        int nIncorrect(0), nCorrect(0);    

        for (int i = 0; i <= tree->GetEntries(); i++)
        {   
            tree->GetEntry(i);

            if (beforeRemovalMinchiSquared == 0)
                continue;

            if (!(NumberHits > NLowerBound && NumberHits <= NUpperBound))
                continue;

            if ((beforeRemovalMinchiSquared - afterRemovalMinchiSquared) >= cut)
            {   
                if (MCIsClusterTwoParticles == 0)
                    nIncorrect++;

                if (MCIsClusterTwoParticles == 1)
                    nCorrect++;
            }   
        }   

        float cutPurity((float)nCorrect/(nCorrect + nIncorrect));
        float cutCompleteness((float)nCorrect/nGoodCuts);

        cutValueWithPurity->SetPoint(n, cut, cutPurity);
        cutValueWithCompleteness->SetPoint(n, cut, cutCompleteness);
        n++;
    } 

    float bestCutValue(0.f), bestPCDistance(100000.f);

    for (int j = 0; j < n; j++)
    {
        Double_t cutValue(0.f), completeness(0.f), purity(0.f);       
        cutValueWithPurity->GetPoint(j, cutValue, purity);
        cutValueWithCompleteness->GetPoint(j, cutValue, completeness);

        if (std::abs(completeness - purity) < bestPCDistance)
        {
            bestPCDistance = std::abs(completeness - purity);
            bestCutValue = cutValue;
        }
    }

    delete cutValueWithPurity;
    delete cutValueWithCompleteness;

    return bestCutValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<int, float> CreateCutsArray(TTree* t1, int nLowerBound, int nUpperBound, int nSteps, std::string variable1, std::string variable2, std::string numberhits)
{
    int nBinSize((nUpperBound - nLowerBound)/nSteps);
    std::map<int, float> cutValueMap;

    for (int i = 1; i <= nSteps; ++i)
        cutValueMap[i * nBinSize] = GetCompletenessPurityIntersection(t1, (i - 1) * nBinSize, i * nBinSize, variable1, variable2, numberhits);

    cutValueMap[(nSteps + 1) * nBinSize] = GetCompletenessPurityIntersection(t1, nUpperBound, 1e6, variable1, variable2, numberhits);
    return cutValueMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<int, float> CreateFailuresWithN(TTree* t1, int nLowerBound, int nUpperBound, int nSteps, std::string variable1, std::string variable2, std::string numberhits)
{
    int     MCForwards;
    int     MCIsClusterTwoParticles;
    int     NumberHits;
    float   beforeRemovalMinchiSquared;
    float   afterRemovalMinchiSquared;
    
    t1->SetBranchAddress("MCForwards",                                             &MCForwards);
    t1->SetBranchAddress("MCIsClusterTwoParticles",                                             &MCIsClusterTwoParticles);
    t1->SetBranchAddress(numberhits.c_str(),                                             &NumberHits);
    t1->SetBranchAddress(variable1.c_str(),                                        &beforeRemovalMinchiSquared);
    t1->SetBranchAddress(variable2.c_str(),                                         &afterRemovalMinchiSquared);

    std::map<int, float> badCutCounter, totalCutCounter, failureCounterWithN;
    int nBinSize((nUpperBound - nLowerBound)/nSteps);

    for (int i = 1; i <= nSteps + 1; ++i)
    {
        totalCutCounter[i] = 0.f;
        badCutCounter[i] = 0.f;
        failureCounterWithN[i * nBinSize] = 0.f;
    }

    for (int j = 0; j <= t1->GetEntries(); j++)
    {
        t1->GetEntry(j);

        if (afterRemovalMinchiSquared == 0)
            continue;

        for (int i = 1; i <= nSteps; ++i)
        {
            if (NumberHits > (i - 1) * nBinSize && NumberHits <= i * nBinSize)
            {
                totalCutCounter[i] += 1.0;
            
                if (!MCIsClusterTwoParticles)
                    badCutCounter[i] += 1.0;
            }
        }

        if (NumberHits > nSteps * nBinSize)
        {
            totalCutCounter[nSteps + 1] += 1.0;
        
            if (!MCIsClusterTwoParticles)
                badCutCounter[nSteps + 1] += 1.0;
        }
    }

    for (int i = 1; i <= nSteps + 1; i++)
    {
        if (totalCutCounter[i] != 0)
            failureCounterWithN[i * nBinSize] = badCutCounter[i]/totalCutCounter[i];
    }

    return failureCounterWithN; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

int GetCutMapKey(std::map<int, float> cutValueMap, int numberHits)
{
    int binSize((*(cutValueMap.begin())).first);

    if (numberHits > ((*(cutValueMap.rbegin())).first))
        return cutValueMap.size() * binSize;

    for (int i = 1; i <= (int) cutValueMap.size(); ++i)
    {   
        if (numberHits > (i - 1) * binSize && numberHits <= i * binSize)
            return i * binSize;
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

TGraph* CreateScatterGraph(TTree* t1, TGraph* correctnessDirectionGraph, std::map<int, float> cutValueMap, std::string hitsDefinition, std::string directionDefinition, int alongDirection, std::string correctnessDefinition, int isCorrect, std::string variable1, std::string variable2, EColor colour, int markerStyle, Size_t markerSize)
{
    int     trueParticleDirection;
    int     isTwoParticles;
    int     numberHits;
    float   beforeRemovalMinchiSquared;
    float   afterRemovalMinchiSquared;

    t1->SetBranchAddress(directionDefinition.c_str(),                              &trueParticleDirection);
    t1->SetBranchAddress(correctnessDefinition.c_str(),                            &isTwoParticles);
    t1->SetBranchAddress(hitsDefinition.c_str(),                                   &numberHits);
    t1->SetBranchAddress(variable1.c_str(),                                        &beforeRemovalMinchiSquared);
    t1->SetBranchAddress(variable2.c_str(),                                        &afterRemovalMinchiSquared);

    int nPoints(0);
    float directionCorrectnessBefore[100000] = {}; 
    float directionCorrectnessAfter[100000] = {}; 

    for (int i = 0; i <= t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if (beforeRemovalMinchiSquared == 0)
            continue;

        if (isTwoParticles == isCorrect)
        {
            if (trueParticleDirection == alongDirection)
            {
                if (cutValueMap.size() != 0)
                {
                    if ((beforeRemovalMinchiSquared - afterRemovalMinchiSquared) < cutValueMap.at(GetCutMapKey(cutValueMap, numberHits)))
                        continue;
                }

                directionCorrectnessBefore[nPoints] = beforeRemovalMinchiSquared;                
                directionCorrectnessAfter[nPoints] = afterRemovalMinchiSquared;
                nPoints++;
            }
        }
    }

    for (int i = 1; i <= nPoints; ++i)
        correctnessDirectionGraph->SetPoint(i, directionCorrectnessBefore[i], directionCorrectnessAfter[i]);

    if (colour == kGreen)
        correctnessDirectionGraph->SetMarkerColor(colour+2);
    else    
        correctnessDirectionGraph->SetMarkerColor(colour);

    correctnessDirectionGraph->SetMarkerStyle(markerStyle);
    correctnessDirectionGraph->SetMarkerSize(markerSize);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DrawMap(std::map<int, float> map, int nLowerBound, int nUpperBound, int nSteps, bool fitPlot, std::string plotName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    Int_t nPoints(nSteps+1);
    TGraph* outputGraph = new TGraph(nPoints);

    int nBinSize((nUpperBound - nLowerBound)/nSteps);

    for (int i = 1; i <= nSteps + 1; ++i)
        outputGraph->SetPoint((i - 1), (i - 1), map.at(i * nBinSize));


    TCanvas *canvas = new TCanvas(plotName.c_str(), plotName.c_str(), 900, 600);
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/scatter/" + plotName + ".png");
    outputGraph->SetMarkerStyle(7);
    outputGraph->SetMarkerSize(1.5);
    outputGraph->SetTitle(plotTitle.c_str());
    outputGraph->GetXaxis()->SetTitle(xTitle.c_str());
    outputGraph->GetYaxis()->SetTitle(yTitle.c_str());
    outputGraph->GetYaxis()->SetRangeUser(0, std::ceil(map.at(nBinSize)));

    TAxis *ax = outputGraph->GetHistogram()->GetXaxis();
    Double_t x1 = ax->GetBinLowEdge(1);
    Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
    outputGraph->GetHistogram()->GetXaxis()->Set(nSteps+1,x1,x2);

    for(int k = 1; k <= nSteps; ++k)
    {
        std::string binName = std::to_string((k - 1) * nBinSize) + "-" + std::to_string(k * nBinSize);
        outputGraph->GetHistogram()->GetXaxis()->SetBinLabel(k, binName.c_str());
    }

    std::string lastBinName = std::to_string(nSteps * nBinSize) + "+";
    outputGraph->GetHistogram()->GetXaxis()->SetBinLabel(nSteps + 1, lastBinName.c_str());

    if (fitPlot)
        outputGraph->Fit("pol1");

    outputGraph->Draw("ALP");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
    delete outputGraph;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DrawScatterPlot(TTree* t1, std::string hitsDefinition, std::string directionDefinition, std::map<int, float> &cutValueMap, std::string correctnessDefinition, std::string variable1, std::string variable2, std::string plotName, std::string xTitle, std::string yTitle, std::string plotTitle)
{

    TGraph* forwardsCorrect = new TGraph();
    TGraph* forwardsIncorrect = new TGraph();
    TGraph* backwardsCorrect = new TGraph();
    TGraph* backwardsIncorrect = new TGraph();

    CreateScatterGraph(t1, forwardsCorrect, cutValueMap, hitsDefinition, directionDefinition, 1, correctnessDefinition, 1, variable1, variable2, kGreen, 22, 0.9);
    CreateScatterGraph(t1, forwardsIncorrect, cutValueMap, hitsDefinition, directionDefinition, 1, correctnessDefinition, 0, variable1, variable2, kRed, 22, 0.9);
    CreateScatterGraph(t1, backwardsCorrect, cutValueMap, hitsDefinition, directionDefinition, 0, correctnessDefinition, 1, variable1, variable2, kGreen, 21, 0.7);
    CreateScatterGraph(t1, backwardsIncorrect, cutValueMap, hitsDefinition, directionDefinition, 0, correctnessDefinition, 0, variable1, variable2, kRed, 21, 0.7);

    TCanvas *canvas1 = new TCanvas(plotName.c_str(), plotName.c_str(), 900, 600);
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/scatter/"+plotName+".png");

    forwardsCorrect->GetXaxis()->SetLimits(0.f, 20.f);
    forwardsCorrect->GetYaxis()->SetRangeUser(0.f, 10.f);
    forwardsCorrect->SetTitle(plotTitle.c_str());
    forwardsCorrect->GetXaxis()->SetTitle(xTitle.c_str());
    forwardsCorrect->GetYaxis()->SetTitle(yTitle.c_str());

    forwardsCorrect->Draw("AP");
    forwardsIncorrect->Draw("PSame");
    backwardsCorrect->Draw("PSame");
    backwardsIncorrect->Draw("PSame");

    canvas1->SaveAs(savePath.c_str());
    delete canvas1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Scatter(void)
{
    //--------------Change settings here--------------------:

    bool multiNCut(true);
    int NLowerBound(0), NUpperBound(1e6); //use NLowerBound == 0 && NUpperBound == 1e6 for all N (sets plot title automatically)
    int removalRoutine(2); //0: track end filter, 1:splitting 2: fragment removal

    //-------------------------------------------------------------------------------------------------------------------------------

    TFile *f = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/saved_results/splitting.root");
    TTree *t1 = (TTree*)f->Get("PFO");
    int nentries = (Int_t)t1->GetEntries();
    std::cout << "There are " << nentries << " entries in the tree" << std::endl;

    std::map<int, float> cutValueMap = CreateCutsArray(t1, 0, 400, 8, "SplittingBeforeChiSquaredPerHit", "SplittingAfterChiSquaredPerHit", "SplittingBeforeNumberHits");
    DrawMap(cutValueMap, 0, 400, 8, true, "cutValueMap", "Bins in N", "Best #chi_{B}^{2}/N - #chi_{A}^{2}/N cut value", "Purity-Completeness Intersection Cut Values for Bins in N");

    std::map<int, float> failureCounterWithN = CreateFailuresWithN(t1, 0, 400, 8, "SplittingBeforeChiSquaredPerHit", "SplittingAfterChiSquaredPerHit", "SplittingBeforeNumberHits");
    DrawMap(failureCounterWithN, 0, 400, 8, false, "failureCounterWithN", "Bins in N", "Fraction of Incorrect Splits", "Fraction of Incorrect Splits for Bins in N");

    std::map<int, float> emptyMap;
    DrawScatterPlot(t1, "SplittingBeforeNumberHits", "MCForwards", emptyMap, "MCIsClusterTwoParticles", "SplittingBeforeChiSquaredPerHit", "SplittingAfterChiSquaredPerHit", "splitting_scatter", "Best #chi^{2}_{B}/N", "Best #chi^{2}_{A}/N", "Splitting scatter plot: best #chi^{2}/N before (B) and after (A) splitting");

    DrawScatterPlot(t1, "SplittingBeforeNumberHits", "MCForwards", cutValueMap, "MCIsClusterTwoParticles", "SplittingBeforeChiSquaredPerHit", "SplittingAfterChiSquaredPerHit", "splitting_scatter_cutapplied", "Best #chi^{2}_{B}/N", "Best #chi^{2}_{A}/N", "Splitting scatter plot after cuts: best #chi^{2}/N before (B) and after (A) splitting");

}
