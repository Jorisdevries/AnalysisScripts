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
#include "TF1.h"
#include "TLine.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <sstream>

#include "Drawing_Helper_Functions.h"

TFile *f = TFile::Open("/usera/jjd49/pandora_direction/Scripts/roots/alldata_final.root");
TTree *t1 = (TTree*)f->Get("DirectionFitTree");

//--------------Change settings here--------------------:

std::string prefix("TEF");
float jumpThreshold(-0.6);
float scatterRangeBefore(50.0), scatterRangeAfter(20.0);

//------------------------------------------------------:

//------------------------------------------------------------------------------------------------------------------------------------------

float GetCompletenessPurityIntersection(TTree* tree, int NLowerBound, int NUpperBound, std::string variable1, std::string variable2, std::string numberhits, std::string correctnessDefinition, std::string jumpDefinition)
{
    std::cout << "*" << std::endl;
    int     TrueForwards = 0;
    int     SplittingCorrect = 0;
    int     NumberHits = 0;
    float   beforeRemovalMinchiSquared = 0;
    float   afterRemovalMinchiSquared = 0;
    float   jumpSize = 0;
    int     fileIdentifier = 0;
    int     eventNumber = 0;

    tree->SetBranchAddress("TrueForwards",                                             &TrueForwards);
    tree->SetBranchAddress(correctnessDefinition.c_str(),                                             &SplittingCorrect);
    tree->SetBranchAddress(numberhits.c_str(),                                             &NumberHits);
    tree->SetBranchAddress(variable1.c_str(),                                        &beforeRemovalMinchiSquared);
    tree->SetBranchAddress(variable2.c_str(),                                         &afterRemovalMinchiSquared);
    tree->SetBranchAddress(jumpDefinition.c_str(),                                         &jumpSize);
    tree->SetBranchAddress("FileIdentifier",                                         &fileIdentifier);
    tree->SetBranchAddress("EventNumber",                                         &eventNumber);

    int nGoodCuts(0), nBadCuts(0);

    for (int i = 0; i <= tree->GetEntries(); i++)
    {   
        tree->GetEntry(i);

        if (beforeRemovalMinchiSquared == 0 || jumpSize < jumpThreshold)
            continue;

        //if (afterRemovalMinchiSquared - beforeRemovalMinchiSquared < -2.0 && SplittingCorrect == 0)
        //    std::cout << "Bad split: " << fileIdentifier << ":" << eventNumber << " with change " << afterRemovalMinchiSquared - beforeRemovalMinchiSquared << std::endl;

        if (!(NumberHits > NLowerBound && NumberHits <= NUpperBound))
            continue;

        if (SplittingCorrect == 0)
            nBadCuts++;

        if (SplittingCorrect == 1)
            nGoodCuts++;
    }

    std::cout << "nGoodCuts: " << nGoodCuts << std::endl;
    std::cout << "nBadCuts: " << nBadCuts << std::endl;

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

            if (beforeRemovalMinchiSquared == 0 || jumpSize < jumpThreshold)
                continue;

            if (!(NumberHits > NLowerBound && NumberHits <= NUpperBound))
                continue;

            if ((beforeRemovalMinchiSquared - afterRemovalMinchiSquared) >= cut)
            {   
                if (SplittingCorrect == 0)
                    nIncorrect++;

                if (SplittingCorrect == 1)
                    nCorrect++;
            }   
        }   

        float cutPurity((float)nCorrect/(nCorrect + nIncorrect));
        float cutCompleteness((float)nCorrect/nGoodCuts);

        cutValueWithPurity->SetPoint(n, cut, cutPurity);
        cutValueWithCompleteness->SetPoint(n, cut, cutCompleteness);
        n++;
    } 

    TCanvas *canvas = new TCanvas("c1", "c1", 900, 600);
    std::string substring(";Cut Value T;Purity/Completeness value");
    std::string title("Purity/Completeness Curves for " + std::to_string(NLowerBound) + " #leq N < " + std::to_string(NUpperBound) + substring);
    cutValueWithPurity->SetTitle(title.c_str());

    cutValueWithPurity->SetMarkerColor(kBlue);
    cutValueWithCompleteness->SetMarkerColor(kMagenta);

    cutValueWithPurity->SetMarkerStyle(9);
    cutValueWithCompleteness->SetMarkerStyle(9);

    cutValueWithPurity->SetMarkerSize(0.5);
    cutValueWithCompleteness->SetMarkerSize(0.5);

    cutValueWithPurity->SetLineColor(kBlue);
    cutValueWithCompleteness->SetLineColor(kMagenta);

    cutValueWithPurity->GetXaxis()->SetLimits(0.0, 10.0);
    cutValueWithPurity->GetYaxis()->SetRangeUser(0.0, 1.0);

    auto legend = new TLegend(0.15,0.15,0.35,0.25);
    legend->SetHeader("Type of Curve"); 
    legend->AddEntry(cutValueWithPurity, "Purity","l");
    legend->AddEntry(cutValueWithCompleteness, "Completeness","l");

    cutValueWithPurity->Draw("AP");
    cutValueWithCompleteness->Draw("Psame");
    legend->Draw("same");
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Scatter/" + prefix + "/" + prefix + "_Purity_Completeness_Intersection_" + std::to_string(NLowerBound) + "_" + std::to_string(NUpperBound)  + ".pdf");
    canvas->SaveAs(savePath.c_str());

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
    delete canvas;

    return bestCutValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<int, float> CreateCutsArray(TTree* t1, int nLowerBound, int nUpperBound, int nSteps, std::string variable1, std::string variable2, std::string numberhits, std::string correctnessDefinition, std::string jumpDefinition)
{
    int nBinSize((nUpperBound - nLowerBound)/nSteps);
    std::map<int, float> cutValueMap;

    for (int i = 1; i <= nSteps; ++i)
        cutValueMap[i * nBinSize] = GetCompletenessPurityIntersection(t1, (i - 1) * nBinSize, i * nBinSize, variable1, variable2, numberhits, correctnessDefinition, jumpDefinition);

    cutValueMap[(nSteps + 1) * nBinSize] = GetCompletenessPurityIntersection(t1, nUpperBound, 1e6, variable1, variable2, numberhits, correctnessDefinition, jumpDefinition);
    return cutValueMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<int, float> CreateFailuresWithN(TTree* t1, int nLowerBound, int nUpperBound, int nSteps, std::string variable1, std::string variable2, std::string numberhits, std::string correctnessDefinition, std::string jumpDefinition)
{
    int     TrueForwards;
    int     SplittingCorrect;
    int     NumberHits;
    float   beforeRemovalMinchiSquared;
    float   afterRemovalMinchiSquared;
    float   jumpSize;
    
    t1->SetBranchAddress("TrueForwards",                                             &TrueForwards);
    t1->SetBranchAddress(correctnessDefinition.c_str(),                                             &SplittingCorrect);
    t1->SetBranchAddress(numberhits.c_str(),                                             &NumberHits);
    t1->SetBranchAddress(variable1.c_str(),                                        &beforeRemovalMinchiSquared);
    t1->SetBranchAddress(variable2.c_str(),                                         &afterRemovalMinchiSquared);
    t1->SetBranchAddress(jumpDefinition.c_str(),                                         &jumpSize);
    std::cout << "-" << std::endl;

    std::map<int, float> badCutCounter, totalCutCounter, failureCounterWithN;
    int nBinSize((nUpperBound - nLowerBound)/nSteps);

    for (int i = 1; i <= nSteps + 1; ++i)
    {
        totalCutCounter[i] = 0.f;
        badCutCounter[i] = 0.f;
        failureCounterWithN[i * nBinSize] = 0.f;
    }

    std::cout << "--" << std::endl;

    for (int j = 1; j < t1->GetEntries(); j++)
    {
        t1->GetEntry(j);

        if (afterRemovalMinchiSquared == 0 || jumpSize < jumpThreshold)
            continue;


        for (int i = 1; i <= nSteps; ++i)
        {
            if (NumberHits > (i - 1) * nBinSize && NumberHits <= i * nBinSize)
            {
                totalCutCounter[i] += 1.0;
            
                if (SplittingCorrect == 0)
                    badCutCounter[i] += 1.0;
            }
        }

        if (NumberHits > nSteps * nBinSize)
        {
            totalCutCounter[nSteps + 1] += 1.0;
        
            if (SplittingCorrect == 0)
                badCutCounter[nSteps + 1] += 1.0;
        }
    }

    std::cout << "---" << std::endl;

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

TGraph* CreateScatterGraph(TTree* t1, TGraph* correctnessDirectionGraph, std::map<int, float> cutValueMap, std::string hitsDefinition, int minN, int maxN, std::string directionDefinition, int alongDirection, std::string correctnessDefinition, int isCorrect, std::string variable1, std::string variable2, std::string jumpDefinition, EColor colour, int markerStyle, Size_t markerSize)
{
    int     trueParticleDirection;
    int     isTwoParticles;
    int     numberHits;
    float   beforeRemovalMinchiSquared;
    float   afterRemovalMinchiSquared;
    float   jumpSize;
    int     fileIdentifier;
    int     eventNumber;

    t1->SetBranchAddress(directionDefinition.c_str(),                              &trueParticleDirection);
    t1->SetBranchAddress(correctnessDefinition.c_str(),                            &isTwoParticles);
    t1->SetBranchAddress(hitsDefinition.c_str(),                                   &numberHits);
    t1->SetBranchAddress(variable1.c_str(),                                        &beforeRemovalMinchiSquared);
    t1->SetBranchAddress(variable2.c_str(),                                        &afterRemovalMinchiSquared);
    t1->SetBranchAddress(jumpDefinition.c_str(),                                   &jumpSize);
    t1->SetBranchAddress("FileIdentifier",                                         &fileIdentifier);
    t1->SetBranchAddress("EventNumber",                                            &eventNumber);

    int nPoints(0);
    float directionCorrectnessBefore[100000] = {}; 
    float directionCorrectnessAfter[100000] = {}; 

    int nGoodEntriesAfterFilter(0), nBadEntriesAfterFilter(0);

    for (int i = 0; i <= t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if (beforeRemovalMinchiSquared == 0 || jumpSize < jumpThreshold || (afterRemovalMinchiSquared - beforeRemovalMinchiSquared) > 0)
            continue;

        if (!(numberHits >= minN && numberHits < maxN))
            continue;

        if (cutValueMap.size() != 0)
        {
            if ((beforeRemovalMinchiSquared - afterRemovalMinchiSquared) < cutValueMap.at(GetCutMapKey(cutValueMap, numberHits)))
                continue;
        }

        if (isTwoParticles == 1)
            ++nGoodEntriesAfterFilter;
        else
            ++nBadEntriesAfterFilter;

        //if (beforeRemovalMinchiSquared > 20.0 && afterRemovalMinchiSquared < 5.0)
        //    std::cout << fileIdentifier << ":" << eventNumber << "(" << beforeRemovalMinchiSquared << " -> " << afterRemovalMinchiSquared << ", N = " << numberHits << ")" << std::endl;

        if (isTwoParticles == isCorrect)
        {
            if (trueParticleDirection == alongDirection)
            {
                directionCorrectnessBefore[nPoints] = beforeRemovalMinchiSquared;                
                directionCorrectnessAfter[nPoints] = afterRemovalMinchiSquared;
                nPoints++;
            }
        }
    }

    std::cout << "Correctness: " << isCorrect << std::endl;
    std::cout << "Direction: " << alongDirection << std::endl;
    std::cout << "Filter: " << (cutValueMap.size() != 0) << std::endl;
    std::cout << "Good entries: " << nGoodEntriesAfterFilter << std::endl;
    std::cout << "Bad entries: " << nBadEntriesAfterFilter << std::endl;

    for (int i = 0; i < nPoints; ++i)
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
    int lastAboveThresholdBin(1);

    for (int i = 1; i <= nSteps + 1; ++i)
    {
        outputGraph->SetPoint((i - 1), (i - 1), map.at(i * nBinSize));

        if (map.at(i * nBinSize) > 0.5)
            lastAboveThresholdBin = i;
    }

    TCanvas *canvas = new TCanvas(plotName.c_str(), plotName.c_str(), 900, 600);
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

    std::cout << "lastAboveThresholdBin: " << lastAboveThresholdBin << std::endl;

    outputGraph->Draw("ALP");

    if (fitPlot)
    {
        /*
        Double_t x1(0.f), y1(0.f);
        outputGraph->GetPoint(0, x1, y1);
        TLine* line = new TLine(0, y1, lastAboveThresholdBin - 1, 0.5);
        TLine* line2 = new TLine(lastAboveThresholdBin - 1, 0.5, nSteps, 0.5);
        line->SetLineColor(kRed);
        line2->SetLineColor(kRed);
        line->Draw("same");
        line2->Draw("same");
        */

        TF1 *f1 = new TF1("f1","[0] + exp([1] + [2]*x)", 0, nSteps);
        f1->SetParameters(1.0, 2, -0.5);

        outputGraph->Fit("f1");
    }
    //    outputGraph->Fit("f1");


    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Scatter/" + plotName + ".pdf");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
    delete outputGraph;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DrawScatterPlot(TTree* t1, std::string hitsDefinition, int minN, int maxN, std::string directionDefinition, std::map<int, float> &cutValueMap, std::string correctnessDefinition, std::string variable1, std::string variable2, std::string plotName, std::string xTitle, std::string yTitle, std::string plotTitle, std::string jumpDefinition)
{

    TGraph* forwardsCorrect = new TGraph();
    TGraph* forwardsIncorrect = new TGraph();
    TGraph* backwardsCorrect = new TGraph();
    TGraph* backwardsIncorrect = new TGraph();

    CreateScatterGraph(t1, forwardsCorrect, cutValueMap, hitsDefinition, minN, maxN, directionDefinition, 1, correctnessDefinition, 1, variable1, variable2, jumpDefinition, kGreen, 22, 0.9);
    CreateScatterGraph(t1, forwardsIncorrect, cutValueMap, hitsDefinition, minN, maxN, directionDefinition, 1, correctnessDefinition, 0, variable1, variable2, jumpDefinition, kRed, 22, 0.9);
    CreateScatterGraph(t1, backwardsCorrect, cutValueMap, hitsDefinition, minN, maxN, directionDefinition, 0, correctnessDefinition, 1, variable1, variable2, jumpDefinition, kGreen, 21, 0.7);
    CreateScatterGraph(t1, backwardsIncorrect, cutValueMap, hitsDefinition, minN, maxN, directionDefinition, 0, correctnessDefinition, 0, variable1, variable2, jumpDefinition, kRed, 21, 0.7);

    TCanvas *canvas1 = new TCanvas(plotName.c_str(), plotName.c_str(), 900, 600);
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Scatter/"+plotName+".pdf");

    forwardsCorrect->GetXaxis()->SetLimits(0.f, scatterRangeBefore);
    forwardsCorrect->GetYaxis()->SetRangeUser(0.f, scatterRangeAfter);
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
    int nentries = (Int_t)t1->GetEntries();
    std::cout << "There are " << nentries << " entries in the tree" << std::endl;

    TH1F *pBadSplitJumpSizeHistogram = new TH1F("pBadSplitJumpSizeHistogram", "", 75, 0.0, 2.0);
    TH1F *pGoodSplitJumpSizeHistogram = new TH1F("pGoodSplitJumpSizeHistogram", "", 75, 0.0, 2.0);

    std::vector<std::pair<std::string, int>> badFilterValues = {std::make_pair(prefix + "Correct", 0)};
    std::vector<std::pair<std::string, int>> goodFilterValues = {std::make_pair(prefix + "Correct", 1)};

    CreateFilteredHistogram(t1, pBadSplitJumpSizeHistogram, prefix + "JumpSize", badFilterValues, "float");
    CreateFilteredHistogram(t1, pGoodSplitJumpSizeHistogram, prefix + "JumpSize", goodFilterValues, "float");

    std::string titleString(prefix);
    if (prefix == "TEF") {titleString = "Bragg Peak Filter";}

    Draw(pBadSplitJumpSizeHistogram, pGoodSplitJumpSizeHistogram, kRed, kGreen, true, "JumpSizeComparison", "Jump Size (MeV/cm)", "Fraction of Entries", "Comparison of Jump Size for Good and Bad Splits", "Bad Splits", "Good Splits");

    std::map<int, float> cutValueMap = CreateCutsArray(t1, 0, 400, 8, prefix + "BeforeMinChiSquaredPerHit", prefix + "AfterMinChiSquaredPerHit", prefix + "BeforeNumberHits", prefix + "Correct", prefix + "JumpSize");
    DrawMap(cutValueMap, 0, 400, 8, true, prefix+"_cut_values", "Bins in N", "Best #chi_{B}^{2}/N - #chi_{A}^{2}/N cut value", titleString + ": Purity-Completeness Intersection Cut Values for Bins in N");

    std::map<int, float> failureCounterWithN = CreateFailuresWithN(t1, 0, 400, 8, prefix + "BeforeMinChiSquaredPerHit", prefix + "AfterMinChiSquaredPerHit", prefix + "BeforeNumberHits", prefix + "Correct", prefix + "JumpSize");
    DrawMap(failureCounterWithN, 0, 400, 8, false, prefix+"_failureCounterWithN", "Bins in N", titleString + ": Fraction of Incorrect Splits", titleString + ": Fraction of Incorrect Splits for Bins in N");

    std::map<int, float> emptyMap;
    DrawScatterPlot(t1, prefix + "BeforeNumberHits", 0, 1000, "TrueForwards", emptyMap, prefix + "Correct", prefix + "BeforeMinChiSquaredPerHit", prefix + "AfterMinChiSquaredPerHit", prefix+"_splitting_scatter_150_200", "#chi^{2}_{min,B}/N", "#chi^{2}_{min,A}/N", titleString + " scatter plot: #chi^{2}_{min,B}/N and #chi^{2}_{min,A}/N", prefix + "JumpSize");

    DrawScatterPlot(t1, prefix + "BeforeNumberHits", 0, 1000, "TrueForwards", cutValueMap, prefix + "Correct", prefix + "BeforeMinChiSquaredPerHit", prefix + "AfterMinChiSquaredPerHit", prefix+"_splitting_scatter_cutapplied", "#chi^{2}_{min,B}/N", "#chi^{2}_{min,A}/N", titleString + " scatter plot after cuts: #chi_{min,B}^{2}/N and #chi^{2}_{min,A}/N", prefix + "JumpSize");
}
