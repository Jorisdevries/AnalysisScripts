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
#include "TLegend.h"
#include "TLine.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iomanip> // setprecision
#include <sstream> // stringstream

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH1F* pHistogram, EColor colour, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram->SetXTitle(xTitle.c_str());
    pHistogram->SetYTitle(yTitle.c_str());
    pHistogram->SetTitle(plotTitle.c_str());
    pHistogram->GetYaxis()->SetTitleOffset(1.5);
    //pHistogram->GetYaxis()->SetRangeUser(0.0, 1.0);
    pHistogram->SetLineColor(colour);
    pHistogram->SetStats(kFALSE);
    pHistogram->Draw();
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/probability/" + histogramName + "_TwoDirections.png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DrawNormalisedHistograms(TTree* t1, TCanvas* canvas, TH1F* pCosmicHistogram, TH1F* pNeutrinoHistogram, std::string variableName, float upper_bound, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    pCosmicHistogram->SetBins(100, 0.0, upper_bound);
    pNeutrinoHistogram->SetBins(100, 0.0, upper_bound);

    float variable;
    float recoLength;
    int   intVariable;
    int NeutrinoInduced;

    if (variableName == "NumberHits")
        t1->SetBranchAddress(variableName.c_str(), &intVariable);
    else
        t1->SetBranchAddress(variableName.c_str(), &variable);
        
    t1->SetBranchAddress("NeutrinoInduced", &NeutrinoInduced);
    t1->SetBranchAddress("recoLength", &recoLength);

    for (int i = 0; i < t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if (variableName == "recoLength")
            variable = recoLength;

        if (NeutrinoInduced == 0)
        {
            if (variableName == "NumberHits")
                pCosmicHistogram->Fill(intVariable);
            else
                pCosmicHistogram->Fill(variable);
        }

        if (NeutrinoInduced == 1)
        {
            if (variableName == "NumberHits")
                pNeutrinoHistogram->Fill(intVariable);
            else
                pNeutrinoHistogram->Fill(variable);
        }
    }

    pCosmicHistogram->Scale(1/(pCosmicHistogram->GetEntries()));
    pNeutrinoHistogram->Scale(1/(pNeutrinoHistogram->GetEntries()));

    float largestBinEntry(0.f);

    for (int i = 0; i < pCosmicHistogram->GetNbinsX(); i++)
    {
        if (pCosmicHistogram->GetBinContent(i) > largestBinEntry) 
            largestBinEntry = pCosmicHistogram->GetBinContent(i);
    }

    for (int i = 0; i < pNeutrinoHistogram->GetNbinsX(); i++)
    {
        if (pNeutrinoHistogram->GetBinContent(i) > largestBinEntry) 
            largestBinEntry = pNeutrinoHistogram->GetBinContent(i);
    }

    canvas->Clear();
    pCosmicHistogram->SetXTitle(xTitle.c_str());
    pCosmicHistogram->SetYTitle(yTitle.c_str());
    pCosmicHistogram->SetTitle(plotTitle.c_str());

    pCosmicHistogram->GetYaxis()->SetTitleOffset(1.3);
    pCosmicHistogram->GetYaxis()->SetRangeUser(0.0, 1.25 * largestBinEntry);

    pCosmicHistogram->SetLineColor(kBlue);
    pCosmicHistogram->Draw();
    pNeutrinoHistogram->SetLineColor(kRed);
    pNeutrinoHistogram->Draw("same");

    auto legend = new TLegend(0.15,0.75,0.35,0.85);
    legend->SetHeader("Event Type"); 
    legend->AddEntry(pCosmicHistogram, "Cosmic Events","l");
    legend->AddEntry(pNeutrinoHistogram, "Neutrino Events","l");
    legend->Draw("same");

    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/probability/" + histogramName + "_TwoDirections.png");
    canvas->SaveAs(savePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH2F* pHistogram, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram->SetXTitle(xTitle.c_str());
    pHistogram->SetYTitle(yTitle.c_str());
    pHistogram->SetTitle(plotTitle.c_str());
    pHistogram->GetYaxis()->SetTitleOffset(1.5);
    pHistogram->Draw("COLZ");
    std::string savePath("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/probability/" + histogramName + "_TwoDirections.png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TH1F* ExtendTails(TH1F* pHistogram, float tailDefinition)
{
    float tailStart(pHistogram->FindBin(tailDefinition));
    float averageBinContent(0.f);
    int NBins(0);

    for (int i = 1; i <= pHistogram->GetNbinsX(); i++)
    {
        if (tailDefinition >= 0.f)
        {
            if (i >= tailStart)
            {
                averageBinContent += pHistogram->GetBinContent(i);
                NBins++;
            }
        }
        else
        {
            if (i <= tailStart)
            {
                averageBinContent += pHistogram->GetBinContent(i);
                NBins++;
            }
        }
    }

    averageBinContent /= NBins;

    for (int i = 1; i <= pHistogram->GetNbinsX(); i++)
    {
        if (tailDefinition >= 0.f)
        {
            if (i >= tailStart)
                pHistogram->SetBinContent(i, averageBinContent);
        }
        else
        {
            if (i <= tailStart)
                pHistogram->SetBinContent(i, averageBinContent);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateProbabilityPlot(TTree* tree, std::string name, int neutrinoMode, int yFaceOnly, std::string direction, std::string variableOne, std::string variableTwo, TGraph* ProbabilityPlot, float histogramRange, int &nBins, float tailDefinition, bool enableFiltering, std::string filterVariable, float lowerFilterBound, float upperFilterBound, std::string xTitle1, std::string title1, std::string xTitle2, std::string title2, bool drawUnderlying)
{
    //TFile f("probability.root","RECREATE");

    TH1F *pHistogramOne = new TH1F("pHistogramOne","", nBins, -histogramRange, histogramRange);
    TH1F *pHistogramTwo = new TH1F("pHistogramTwo","", nBins, -histogramRange, histogramRange);

    int   Direction;
    float VariableOne;
    float VariableTwo;
    float FilterVariable;
    int   NeutrinoInduced;
    int   MCIntersectsYFace;
    int   NumberHits;
    float recoLength;
    int fileIdentifier;
    int eventNumber;

    tree->SetBranchAddress(direction.c_str(), &Direction);
    tree->SetBranchAddress(variableOne.c_str(), &VariableOne);
    tree->SetBranchAddress(variableTwo.c_str(), &VariableTwo);
    if (filterVariable != "NumberHits")
        tree->SetBranchAddress(filterVariable.c_str(), &FilterVariable);
    tree->SetBranchAddress("NeutrinoInduced", &NeutrinoInduced);
    tree->SetBranchAddress("MCIntersectsYFace", &MCIntersectsYFace);
    tree->SetBranchAddress("NumberHits", &NumberHits);
    tree->SetBranchAddress("recoLength", &recoLength);
    tree->SetBranchAddress("FileIdentifier", &fileIdentifier);
    tree->SetBranchAddress("EventNumber", &eventNumber);

    int nentries = (Int_t)tree->GetEntries();
    //std::cout << "There are " << nentries << " entries in the tree" << std::endl;

    float notAlongDirectionPositive(0.f), notAlongDirectionNegative(0.f), totalNotAlongDirection(0.f), alongDirectionPositive(0.f), alongDirectionNegative(0.f), totalAlongDirection(0.f);
    float largerThanNinetyFive(0.f), totalEntries(0.f);
    
    for (int i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        if (enableFiltering && filterVariable == "recoLength")
            FilterVariable = recoLength;

        if (enableFiltering && filterVariable == "NumberHits")
            FilterVariable = NumberHits;

        if (enableFiltering && !(FilterVariable >= lowerFilterBound && FilterVariable <= upperFilterBound))
            continue;

        if ((neutrinoMode == 0 && NeutrinoInduced != 1) || (neutrinoMode == 1 && NeutrinoInduced != 0))
            continue;

        if (yFaceOnly == 1 && MCIntersectsYFace != 1)
            continue;

        if (Direction == 0)
        {
            pHistogramTwo->Fill((VariableOne - VariableTwo));
            totalNotAlongDirection += 1.f;

            if ((VariableOne - VariableTwo) <= 0.f)
                notAlongDirectionNegative += 1.f;
            else
                notAlongDirectionPositive += 1.f;
        }
        
        if (Direction == 1)
        {
            pHistogramOne->Fill((VariableOne - VariableTwo));
            totalAlongDirection += 1.f;

            if ((VariableOne - VariableTwo) >= -0.5 && (VariableOne - VariableTwo) <= 0.5 && eventNumber <= 5)
                std::cout << "Bad event: " << fileIdentifier << ":" << eventNumber << std::endl;

            if ((VariableOne - VariableTwo) <= 0.f)
                alongDirectionNegative += 1.f;
            else
                alongDirectionPositive += 1.f;
        }
    }
    
    /*
    std::cout << "Number of entries along direction: " << totalAlongDirection << std::endl;
    std::cout << "Number of entries not along direction: " << totalNotAlongDirection << std::endl;
    std::cout << "Fraction of entries along direction that are positive: " << alongDirectionPositive/totalAlongDirection << std::endl;
    std::cout << "Fraction of entries along direction that are negative: " << alongDirectionNegative/totalAlongDirection << std::endl;
    std::cout << "Fraction of entries not along direction that are positive: " << notAlongDirectionPositive/totalNotAlongDirection << std::endl;
    std::cout << "Fraction of entries not along direction that are negative: " << notAlongDirectionNegative/totalNotAlongDirection << std::endl;
    */

    if (drawUnderlying)
    {
        Draw(pHistogramOne, kBlue, name+variableOne, xTitle1.c_str(), "Number of Entries", title1.c_str());
        Draw(pHistogramTwo, kRed, name+variableTwo, xTitle2.c_str(), "Number of Entries", title2.c_str());
    }

    ExtendTails(pHistogramOne, tailDefinition);
    ExtendTails(pHistogramTwo, -tailDefinition);

    if (drawUnderlying)
    {
        Draw(pHistogramOne, kBlue, name+variableOne+"_extended", (xTitle1 + " with Extended Tails").c_str(), "Number of Entries", (title1 + " with Extended Tails").c_str());
        Draw(pHistogramTwo, kRed, name+variableTwo+"_extended", (xTitle2 + " with Extended Tails").c_str(), "Number of Entries", (title2 + " with Extended Tails").c_str());
    }

    for (int i = 1; i <= nBins; i++)
    {
        //Normalised distributions
        float forwardsBinEntry(pHistogramOne->GetBinContent(i)/pHistogramOne->GetEntries());
        float backwardsBinEntry(pHistogramTwo->GetBinContent(i)/pHistogramTwo->GetEntries());

        float binCenter(pHistogramOne->GetBinCenter(i));
        float forwardsProbability(forwardsBinEntry/(forwardsBinEntry + backwardsBinEntry));

        totalEntries += 1.f;

        if (forwardsProbability >= 0.95 || forwardsProbability <= 0.5)
            largerThanNinetyFive += 1.f;

        ProbabilityPlot->SetPoint(i, binCenter, forwardsProbability);
    }

    //std::cout << "Fraction where probability is larger than 95% or smaller than 5%: " << largerThanNinetyFive/totalEntries << std::endl;
    //std::cout << "--------------------------------------------------------------------------" << std::endl;

    delete pHistogramOne;
    delete pHistogramTwo;

    //std::cout << "----------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProbabilityComparison(std::vector<std::pair<float, float>> cutValues, std::string tokenValue, TTree* tree, std::string name, int neutrinoMode, int yFaceOnly, std::string direction, std::string variableOne, std::string variableTwo, float histogramRange, int &nBins, float tailDefinition, std::string filterVariable, std::string plotTitle)
{
    std::vector<EColor> colourVector = {kRed, kBlue, kMagenta, kBlack, kOrange, kGreen, kYellow};

    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    auto legend = new TLegend(0.55,0.675,0.875,0.875);
    legend->SetHeader("Colour Legend"); 

    TGraph *Probability = new TGraph(nBins);
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(colourVector.at(0));
    Probability->SetMinimum(0.0);
    Probability->SetMaximum(1.0);
    Probability->SetTitle(plotTitle.c_str());

    CreateProbabilityPlot(tree, "unmodified_comparison_", neutrinoMode, yFaceOnly, direction, variableOne, variableTwo, Probability, histogramRange, nBins, tailDefinition, false, filterVariable, 0.0, 1.0, "xTitle1", "title1", "xTitle2", "title2", false);
    legend->AddEntry(Probability,"BNB + Cosmic: All Entries","p");
    Probability->Draw("AP");

    int counter(1);

    for (auto &pair : cutValues)
    {
        TGraph *Probability_Filtered = new TGraph(nBins);
        Probability_Filtered->SetMarkerStyle(6);
        Probability_Filtered->SetMarkerColor(colourVector.at(counter));

        CreateProbabilityPlot(tree, name+to_string(counter)+"_comparison_", neutrinoMode, yFaceOnly, direction, variableOne, variableTwo, Probability_Filtered, histogramRange, nBins, tailDefinition, true, filterVariable, pair.first, pair.second, "xTitle1", "title1", "xTitle2", "title2", false);

        Probability_Filtered->Draw("PSame");

        std::stringstream stream;
        if (filterVariable == "NumberHits" || filterVariable == "recoLength")
            stream << fixed << setprecision(0) << pair.first;
        else
            stream << fixed << setprecision(2) << pair.first;
        std::string lower_bound = stream.str();

        std::stringstream stream2;
        if (filterVariable == "NumberHits" || filterVariable == "recoLength")
            stream2 << fixed << setprecision(0) << pair.second;
        else
            stream2 << fixed << setprecision(2) << pair.second;
        std::string upper_bound = stream2.str();
        if (pair.second == 1e6)
            upper_bound = lower_bound + "+";

        std::string legendEntry("BNB + Cosmic: "+ lower_bound +" #leq "+tokenValue+" #leq " + upper_bound);
        legend->AddEntry(Probability_Filtered, legendEntry.c_str(),"p");

        counter++;
    }

    legend->Draw("same");
    std::string saveName("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/probability/Probability_comparison_" + name + "_TwoDirections.png");
    c0->SaveAs(saveName.c_str());

    delete Probability;
    delete c0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProbabilityYFaceComparison(TTree* tree, std::string name, int neutrinoMode, std::string direction, std::string variableOne, std::string variableTwo, float histogramRange, int &nBins, float tailDefinition, std::string plotTitle)
{
    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    auto legend = new TLegend(0.55,0.675,0.875,0.875);
    legend->SetHeader("Colour Legend"); 

    TGraph *Probability = new TGraph(nBins);
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(kRed);
    Probability->SetMinimum(0.0);
    Probability->SetMaximum(1.0);
    Probability->SetTitle(plotTitle.c_str());

    CreateProbabilityPlot(tree, "_all_", neutrinoMode, 0, direction, variableOne, variableTwo, Probability, histogramRange, nBins, tailDefinition, false, "MCPhi", 0.0, 1.0, "xTitle1", "title1", "xTitle2", "title2", false);
    legend->AddEntry(Probability,"BNB + Cosmic: All Entries","p");
    Probability->Draw("AP");

    TGraph *Probability_YFace = new TGraph(nBins);
    Probability_YFace->SetMarkerStyle(6);
    Probability_YFace->SetMarkerColor(kBlue);

    CreateProbabilityPlot(tree, "_yface_", neutrinoMode, 1, direction, variableOne, variableTwo, Probability_YFace, histogramRange, nBins, tailDefinition, false, "MCPhi", 0.0, 1.0, "xTitle1", "title1", "xTitle2", "title2", false);

    legend->AddEntry(Probability_YFace,"BNB + Cosmic: Tracks Intersect Top Face","p");
    Probability_YFace->Draw("PSame");

    legend->Draw("same");
    std::string saveName("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/probability/Probability_comparison_" + name + "_TwoDirections.png");
    c0->SaveAs(saveName.c_str());

    delete Probability;
    delete c0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SingleProbabilityCurve(TTree* tree, std::string name, int neutrinoMode, std::string direction, std::string variableOne, std::string variableTwo, float histogramRange, int &nBins, float tailDefinition, std::string plotTitle)
{
    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    auto legend = new TLegend(0.65,0.75,0.85,0.85);
    legend->SetHeader("Colour Legend"); 

    TGraph *Probability = new TGraph(nBins);
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(kBlack);
    Probability->SetMinimum(0.0);
    Probability->SetMaximum(1.0);
    Probability->SetTitle(plotTitle.c_str());

    CreateProbabilityPlot(tree, "_all_", neutrinoMode, 0, direction, variableOne, variableTwo, Probability, histogramRange, nBins, tailDefinition, false, "MCPhi", 0.0, 1.0, "xTitle1", "title1", "xTitle2", "title2", true);
    legend->AddEntry(Probability,"BNB + Cosmic: All Entries","p");
    Probability->Draw("AP");
    legend->Draw("same");

    std::string saveName("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/probability/Probability_comparison_" + name + "_TwoDirections.png");
    c0->SaveAs(saveName.c_str());

    delete Probability;
    delete c0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreatePurityCompletenessCurve(TTree* tree, std::string direction, std::string deltaChiSquaredName, std::string variableName, float maxCut, float cutStep, std::string comparison, std::string tokenValue)
{
    int   Direction;
    float Variable;
    float DeltaChiSquaredPerHit;

    tree->SetBranchAddress(direction.c_str(), &Direction);
    tree->SetBranchAddress(variableName.c_str(), &Variable);
    tree->SetBranchAddress(deltaChiSquaredName.c_str(), &DeltaChiSquaredPerHit);

    TGraph *cutValueWithPurity = new TGraph(maxCut/cutStep + 1);
    TGraph *cutValueWithCompleteness = new TGraph(maxCut/cutStep + 1);

    cutValueWithPurity->SetMarkerStyle(6);
    cutValueWithPurity->SetMarkerColor(kBlue);
    std::string title("Purity and Completeness as a function of a cut on " + tokenValue +";Cut T on " + tokenValue + ": " + tokenValue + "#geq T" + ";Cut Purity/Efficiency");
    if (comparison == "lesser")
        title = ("Purity and Completeness as a function of a cut on " + tokenValue +";Cut T on " + tokenValue + ": " + tokenValue + "#leq T" + ";Cut Purity/Efficiency");
    cutValueWithPurity->SetTitle(title.c_str());


    cutValueWithCompleteness->SetMarkerStyle(6);
    cutValueWithCompleteness->SetMarkerColor(kRed);

    int n(0);

    for (float cut = 0.0; cut < maxCut; cut += cutStep)
    {
        int nIncorrect(0), nCorrect(0), nCosmics(0);

        for (int i = 0; i <= tree->GetEntries(); i++)
        {
            tree->GetEntry(i);

            if (Direction == 1)
                nCosmics++;

            if (comparison == "lesser" && Variable <= cut)
            {   
                if ((Direction == 0 && DeltaChiSquaredPerHit < 0) || (Direction == 1 && DeltaChiSquaredPerHit > 0))
                    nIncorrect++;

                if ((Direction == 1 && DeltaChiSquaredPerHit < 0) || (Direction == 0 && DeltaChiSquaredPerHit > 0))
                    nCorrect++;
            }   

            if (comparison == "greater" && Variable >= cut)
            {   
                if ((Direction == 0 && DeltaChiSquaredPerHit < 0) || (Direction == 1 && DeltaChiSquaredPerHit > 0))
                    nIncorrect++;

                if ((Direction == 1 && DeltaChiSquaredPerHit < 0) || (Direction == 0 && DeltaChiSquaredPerHit > 0))
                    nCorrect++;
            }   
        }   

        if (nCorrect == 0)
            continue;

        float cutPurity((float)nCorrect/(nCorrect + nIncorrect));
        float cutCompleteness((float)nCorrect/nCosmics);

        cutValueWithPurity->SetPoint(n, cut, cutPurity);
        cutValueWithCompleteness->SetPoint(n, cut, cutCompleteness);
        n++;
    }  

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    auto legend = new TLegend(0.55,0.375,0.75,0.575);
    legend->SetHeader("Colour Legend"); 
    legend->AddEntry(cutValueWithPurity,"Purity","p");
    legend->AddEntry(cutValueWithCompleteness,"Efficiency","p");
    cutValueWithPurity->Draw("AP");
    cutValueWithCompleteness->Draw("Psame");
    legend->Draw("same");

    std::string saveName("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/probability/Purity_completeness_" + variableName + "_TwoDirections.png");
    canvas->SaveAs(saveName.c_str());
    delete canvas;
    delete cutValueWithPurity;
    delete cutValueWithCompleteness;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Probability(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/saved_results/mc_ccqel.root");
    //TFile *f2 = TFile::Open("/usera/jjd49/new_LAr/CondorUtilities/saved_results/probability_comparison/probability_1.root");
    //TFile *f3 = TFile::Open("/usera/jjd49/new_LAr/CondorUtilities/saved_results/probability_comparison/probability_2.root");
    
    TTree *t1 = (TTree*)f1->Get("PFO");
    //TTree *t2 = (TTree*)f2->Get("Tree");
    //TTree *t3 = (TTree*)f3->Get("Tree");
    //
    
    int nBins(80);

    //--------------------------------------------------------------------------------------------

    TH1F *pMCForwardsHistogram = new TH1F("pMCDownardsHistogram","", 2, 0, 2);
    TH2F *MinChiSquaredWithN = new TH2F("MinChiSquaredWithN", "", 50, 0.0, 5.0, 50, 0.0, 400);
    TH2F *MinWithDelta = new TH2F("MinWithDelta", "", 50, 0.0, 5.0, 50, -15.0, 15.0);

    int MCForwards;
    float MinChiSquaredPerHit;
    float UpDownDeltaChiSquaredPerHit;
    int NumberHits;

    t1->SetBranchAddress("MCForwards", &MCForwards);
    t1->SetBranchAddress("MinChiSquaredPerHit", &MinChiSquaredPerHit);
    t1->SetBranchAddress("UpDownDeltaChiSquaredPerHit", &UpDownDeltaChiSquaredPerHit);
    t1->SetBranchAddress("NumberHits", &NumberHits);

    int nentries = (Int_t)t1->GetEntries();

    for (int i = 0; i < nentries; i++)
    {
        t1->GetEntry(i);
        pMCForwardsHistogram->Fill(MCForwards); 
        MinChiSquaredWithN->Fill(MinChiSquaredPerHit, NumberHits);
        MinWithDelta->Fill(MinChiSquaredPerHit, UpDownDeltaChiSquaredPerHit);
    }

    pMCForwardsHistogram->Scale(1/(pMCForwardsHistogram->GetEntries()));

    //Draw(pMCForwardsHistogram, kBlue, "mcdownwards", "MC Direction (0 = upwards, 1 = downwards)", "Fraction of Events", "Fraction of Backwards and Forwards muon in BNB + Cosmic sample");
    //Draw(MinChiSquaredWithN, "MinChiSquaredWithN", "#chi^{2}_{min}/N", "Number of Cluster Hits N", "Cluster Hits N with #chi^{2}_{min}/N");
    //Draw(MinWithDelta, "MinWithDelta", "#chi^{2}_{min}/N", "#Delta#chi^{2}/N", "#Delta#chi^{2}/N with #chi^{2}_{min}/N");

    std::string vanillaPlotTitle("Forwards Probability Distribution For All BNB + Cosmic Events;#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N;Forwards Probability");
    SingleProbabilityCurve(t1, "Single_probability_curve", 0, "MCForwards", "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", 15.0, nBins, 2.0, vanillaPlotTitle);

    //CreatePurityCompletenessCurve(t1, "MCForwards", "DeltaChiSquaredPerHit", "MinChiSquaredPerHit", 10.0, 0.25, "lesser", "#chi^{2}_{min}/N");

    //--------------------------------------------------------------------------------------------

    //std::string yFacePlotTitle("Forwards Probability Distribution Comparison Between All Entries and Top Face Intersecting Tracks;#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N;Forwards Probability");
    //ProbabilityYFaceComparison(t1, "YFace_comparison", 2, "MCForwards", "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", 15.0, nBins, 2.0, yFacePlotTitle);

    //--------------------------------------------------------------------------------------------

    /*
    std::vector<std::pair<float, float>> phiCutValues = {std::pair<float, float>(0.52, 2.62), std::pair<float, float>(1.05, 2.09)};
    std::string phiPlotTitle("Forwards Probability Distribution Variation with #phi (angle of track with X axis);#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N;Forwards Probability");
    ProbabilityComparison(phiCutValues, "#phi", t1, "phi_comparison", 2, 0, "MCForwards", "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", 15.0, nBins, 2.0, "MCPhi", phiPlotTitle);

    std::vector<std::pair<float, float>> thetaCutValues = {std::pair<float, float>(1.4, 3.14)};
    std::string thetaPlotTitle("Forwards Probability Distribution Variation with #theta (angle of track with Y axis);#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N;Forwards Probability");
    ProbabilityComparison(thetaCutValues, "#theta", t1, "theta_comparison", 2, 0, "MCForwards", "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", 15.0, nBins, 2.0, "MCTheta", thetaPlotTitle);

    std::vector<std::pair<float, float>> NCutValues = {std::pair<float, float>(0, 100), std::pair<float, float>(100, 200), std::pair<float, float>(200, 300), std::pair<float, float>(300, 1e6)};
    std::string NPlotTitle("Forwards Probability Distribution Variation with N (number of W cluster hits);#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N;Forwards Probability");
    ProbabilityComparison(NCutValues, "N", t1, "N_comparison", 2, 0, "MCForwards", "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", 15.0, nBins, 2.0, "NumberHits", NPlotTitle);

    std::vector<std::pair<float, float>> minChiCutValues = {std::pair<float, float>(0.0, 1.0), std::pair<float, float>(1.0, 2.0), std::pair<float, float>(2.0, 3.0), std::pair<float, float>(3.0, 1e6)};
    std::string minChiPlotTitle("Forwards Probability Distribution Variation with #chi^{2}_{min}/N (Best fit #chi^{2}/N in event);#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N;Forwards Probability");
    ProbabilityComparison(minChiCutValues, "#chi^{2}_{min}/N", t1, "minChi_comparison", 2, 0, "MCForwards", "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", 15.0, nBins, 2.0, "MinChiSquaredPerHit", minChiPlotTitle);

    std::vector<std::pair<float, float>> recoLengthCutValues = {std::pair<float, float>(0.0, 50.0), std::pair<float, float>(50.0, 100.0), std::pair<float, float>(100.0, 150.0), std::pair<float, float>(150.0, 200.0), std::pair<float, float>(200.0, 1e6)};
    std::string recoLengthPlotTitle("Forwards Probability Distribution Variation with L (Reconstructed 3D Track Length);#Delta#chi^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N;Forwards Probability");
    ProbabilityComparison(recoLengthCutValues, "L", t1, "recoLength_comparison", 2, 0, "MCForwards", "ForwardsChiSquaredPerHit", "BackwardsChiSquaredPerHit", 15.0, nBins, 2.0, "recoLength", recoLengthPlotTitle);
    */

    //--------------------------------------------------------------------------------------------
    
    /*
    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    TH1F *pCosmicHistogram = new TH1F("pCosmicHistogram","", 100, 0.0, 100.0);
    TH1F *pNeutrinoHistogram = new TH1F("pNeutrinoHistogram","", 100, 0.0, 100.0);
    
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "MCPhi", 3.14, "phi_distribution", "Angle with X axis #phi", "Number of Entries", "Normalised Distributions of Track Angle With X Axis #phi for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "MCTheta", 3.14, "theta_distribution", "Angle with Y axis #theta", "Number of Entries", "Normalised Distributions of Track Angle With Y Axis #theta for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "NumberHits", 300, "N_distribution", "Number of Cluster Hits N", "Number of Entries", "Normalised Distributions of Number of Cluster Hits N for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "MinChiSquaredPerHit", 4.0, "minchisquared_distribution", "#chi^{2}_{min}/N", "Number of Entries", "Normalised Distributions of #chi^{2}_{min}/N for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "recoLength", 300, "recolength_distribution", "Reconstructed endpoint-to-endpoint track length L", "Number of Entries", "Reconstructed endpoint-to-endpoint track length L");
    */

    //--------------------------------------------------------------------------------------------
}
