#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2DErrors.h"
#include "TF2.h"
#include "TF3.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TF1.h"
#include "TRandom.h"
#include "TText.h"
#include "TGaxis.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iomanip> // setprecision
#include <sstream> // stringstream

//------------------------------------------------------------------------------------------------------------------------------------------

//std::string globalMode("mc_cosmic");
std::string globalMode("data");

float p0_par1(0.f), p0_par2(0.f), p0_par3(0.f);
float alpha_par1(0.f), alpha_par2(0.f), alpha_par3(0.f);
float beta_par1(0.f), beta_par2(0.f), beta_par3(0.f);

float p0_par1_simple(0.f), p0_par2_simple(0.f), p0_par3_simple(0.f);
float alpha_par1_simple(0.f), alpha_par2_simple(0.f), alpha_par3_simple(0.f);
float beta_par1_simple(0.f), beta_par2_simple(0.f), beta_par3_simple(0.f);

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
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/distributions/underlying/" + histogramName + ".png");
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

    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/distributions/"+ histogramName + ".png");
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
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/distributions/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ExtendTails(TH1F* pHistogram, float tailDefinition)
{
    float tailStart(pHistogram->FindBin(tailDefinition));
    float leftAverageBinContent(0.f), rightAverageBinContent(0.f);
    int NLeftBins(0), NRightBins(0);

    for (int i = 1; i <= pHistogram->GetNbinsX(); i++) 
    {    
        if (i >= tailStart)
        {    
            rightAverageBinContent += pHistogram->GetBinContent(i);
            NRightBins++;
        }    
        if (i <= -tailStart)
        {    
            leftAverageBinContent += pHistogram->GetBinContent(i);
            NLeftBins++;
        }    
    }    

    leftAverageBinContent /= NLeftBins;
    rightAverageBinContent /= NRightBins;

    for (int i = 1; i <= pHistogram->GetNbinsX(); i++) 
    {    
        if (i >= tailStart)
            pHistogram->SetBinContent(i, rightAverageBinContent);
        if (i <= -tailStart)
            pHistogram->SetBinContent(i, leftAverageBinContent);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateProbabilityPlot(TTree* tree, std::string name, int neutrinoMode, int yFaceOnly, std::string direction, std::string variableOne, TGraphAsymmErrors* ProbabilityPlot, float histogramRange, int &nBins, float tailDefinition, bool enableFiltering, std::string filterVariable, float lowerFilterBound, float upperFilterBound, bool chiSquaredFiltering, float chsLowerBound, float chsUpperBound, std::string xTitle1, std::string title1, std::string xTitle2, std::string title2, bool drawUnderlying)
{
    //TFile f("probability.root","RECREATE");

    TH1F *pHistogramOne = new TH1F("pHistogramOne","", nBins, -histogramRange, histogramRange);

    int   Direction;
    float VariableOne;
    float FilterVariable;
    int   NeutrinoInduced;
    int   IntersectsYFace;
    int   NumberHits;
    float recoLength;
    float MinChiSquaredPerHit;
    int FileIdentifier;
    int EventNumber;

    tree->SetBranchAddress(direction.c_str(), &Direction);
    tree->SetBranchAddress(variableOne.c_str(), &VariableOne);
    if (filterVariable != "NumberHits")
        tree->SetBranchAddress(filterVariable.c_str(), &FilterVariable);
    tree->SetBranchAddress("NeutrinoInduced", &NeutrinoInduced);
    tree->SetBranchAddress("IntersectsYFace", &IntersectsYFace);
    tree->SetBranchAddress("NumberHits", &NumberHits);
    tree->SetBranchAddress("recoLength", &recoLength);
    tree->SetBranchAddress("MinChiSquaredPerHit", &MinChiSquaredPerHit);
    tree->SetBranchAddress("FileIdentifier", &FileIdentifier);
    tree->SetBranchAddress("EventNumber", &EventNumber);

    int nentries = (Int_t)tree->GetEntries();
    //std::cout << "There are " << nentries << " entries in the tree" << std::endl;

    float notAlongDirectionPositive(0.f), notAlongDirectionNegative(0.f), totalNotAlongDirection(0.f), alongDirectionPositive(0.f), alongDirectionNegative(0.f), totalAlongDirection(0.f);
    
    for (int i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        if (chiSquaredFiltering && (MinChiSquaredPerHit < chsLowerBound || MinChiSquaredPerHit > chsUpperBound))
            continue;

        if (enableFiltering && filterVariable == "recoLength")
            FilterVariable = recoLength;

        if (enableFiltering && filterVariable == "NumberHits")
            FilterVariable = NumberHits;

        if (enableFiltering && filterVariable == "MinChiSquaredPerHit")
            FilterVariable = MinChiSquaredPerHit;

        if (enableFiltering && !(FilterVariable >= lowerFilterBound && FilterVariable <= upperFilterBound))
            continue;

        if ((neutrinoMode == 0 && NeutrinoInduced != 1) || (neutrinoMode == 1 && NeutrinoInduced != 0))
            continue;

        if (yFaceOnly == 1 && IntersectsYFace != 1)
            continue;

        if (Direction == 1)
        {
            pHistogramOne->Fill(VariableOne);

            totalAlongDirection += 1.f;

            if (VariableOne <= 0.f)
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
    }

    ExtendTails(pHistogramOne, tailDefinition);

    if (drawUnderlying)
    {
        Draw(pHistogramOne, kBlue, name+variableOne+"_extended", (xTitle1 + " with Extended Tails").c_str(), "Number of Entries", (title1 + " with Extended Tails").c_str());
    }

    //float previousProbability(0.5), previousError(0.f);

    for (int i = 1; i <= nBins + 1; i++)
    {
        //Normalised distributions
        float negativeBinEntry(pHistogramOne->GetBinContent(i));
        float binCenter(pHistogramOne->GetBinCenter(i));

        if (tailDefinition > 0 && binCenter > 0)
            continue;

        if (tailDefinition < 0 && binCenter < 0)
            continue;

        float mirrorBinNumber(pHistogramOne->GetXaxis()->FindBin(-1.0 * binCenter));
        float mirrorBinEntry(pHistogramOne->GetBinContent(mirrorBinNumber));

        int N(negativeBinEntry + mirrorBinEntry);
        float probability(negativeBinEntry/(negativeBinEntry + mirrorBinEntry));

        if (probability == 1.0)
            probability = 0.98;

        float pointError(std::sqrt(probability * (1 - probability)/N));

        float upperError(probability + pointError > 1.0 ? 1.0 - probability : pointError);
        ProbabilityPlot->SetPointError(i-1, 0.f, 0.f, pointError, upperError);
        ProbabilityPlot->SetPoint(i-1, std::abs(binCenter), probability);

        /*
        if (negativeBinEntry == 0 || mirrorBinEntry == 0)
        {
            ProbabilityPlot->SetPoint(i-1, std::abs(binCenter), previousProbability);
            ProbabilityPlot->SetPointError(i-1, 0.f, previousError);
            continue;
        }
        else
        {
            ProbabilityPlot->SetPoint(i-1, std::abs(binCenter), probability);
            ProbabilityPlot->SetPointError(i-1, 0.f, pointError);
        }

        previousProbability = probability;
        previousError = pointError;
        */
    }

    delete pHistogramOne;

    //std::cout << "----------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProbabilityComparison(std::vector<std::pair<float, float>> cutValues1, bool chiSquaredFiltering, float chsLowerBound, float chsUpperBound, std::string tokenValue, TTree* tree, std::string name, int neutrinoMode, int yFaceOnly, std::string direction, std::string variableOne, float histogramRange, int &nBins, float tailDefinition, std::string filterVariable, std::string plotTitle)
{
    std::vector<EColor> colourVector = {kRed, kBlue, kBlack, kOrange, kGreen, kYellow, kCyan, kGray};

    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    auto legend = new TLegend(0.55,0.175,0.875,0.375);
    legend->SetHeader("Colour Legend"); 

    TGraphAsymmErrors *Probability = new TGraphAsymmErrors();
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(colourVector.at(0));
    Probability->SetLineColor(colourVector.at(0));
    Probability->SetMinimum(0.0);
    Probability->SetMaximum(1.0);

    std::stringstream stream3;
    stream3 << fixed << setprecision(0) << chsLowerBound;
    string chs_lb_str = stream3.str();

    std::stringstream stream4;
    stream4 << fixed << setprecision(0) << chsUpperBound;
    string chs_ub_str = stream4.str();
   
    if (chsUpperBound == 1e6)
        chs_ub_str = "#infty";

    if (chiSquaredFiltering)
    {
        std::string insertTitle = " for " + chs_lb_str + " #leq #chi^{2}_{min}/N #leq " + chs_ub_str;
        plotTitle.insert(plotTitle.find(";"), insertTitle);
    }

    Probability->SetTitle(plotTitle.c_str());

    CreateProbabilityPlot(tree, "unmodified_comparison_", neutrinoMode, yFaceOnly, direction, variableOne, Probability, histogramRange, nBins, tailDefinition, false, filterVariable, 0.0, 1.0, chiSquaredFiltering, chsLowerBound, chsUpperBound, "xTitle1", "title1", "xTitle2", "title2", false);
    legend->AddEntry(Probability,"Cosmics: All Entries","p");
    Probability->Draw("APL");
    Probability->GetXaxis()->SetLimits(-1.0, histogramRange);
    Probability->Draw("APL");

    int counter(1);

    for (auto &pair : cutValues1)
    {
        TGraphAsymmErrors *Probability_Filtered = new TGraphAsymmErrors();
        Probability_Filtered->SetMarkerStyle(6);
        Probability_Filtered->SetMarkerStyle(6);
        Probability_Filtered->SetLineColor(colourVector.at(counter));
        Probability_Filtered->SetMarkerColor(colourVector.at(counter));

        CreateProbabilityPlot(tree, name+to_string(counter)+"_comparison_", neutrinoMode, yFaceOnly, direction, variableOne, Probability_Filtered, histogramRange, nBins, tailDefinition, true, filterVariable, pair.first, pair.second, chiSquaredFiltering, chsLowerBound, chsUpperBound, "xTitle1", "title1", "xTitle2", "title2", false);

        Probability_Filtered->Draw("PLSame");

        std::stringstream stream;
        if (filterVariable == "NumberHits" || filterVariable == "recoLength" || filterVariable == "MinChiSquaredPerHit")
            stream << fixed << setprecision(0) << pair.first;
        else
            stream << fixed << setprecision(2) << pair.first;
        std::string lower_bound = stream.str();

        std::stringstream stream2;
        if (filterVariable == "NumberHits" || filterVariable == "recoLength" || filterVariable == "MinChiSquaredPerHit")
            stream2 << fixed << setprecision(0) << pair.second;
        else
            stream2 << fixed << setprecision(2) << pair.second;
        std::string upper_bound = stream2.str();
        if (pair.second == 1e6)
            upper_bound = "#infty";

        std::string legendEntry("Cosmics: "+ lower_bound +" #leq "+tokenValue+" #leq " + upper_bound);
        legend->AddEntry(Probability_Filtered, legendEntry.c_str(),"p");

        counter++;
    }

    legend->Draw("same");
    std::string saveName("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/comparisons/Probability_comparison_" + name);


    if (chiSquaredFiltering)
        saveName = saveName + "_chisquared_filtered_" + chs_lb_str + "_" + chs_ub_str;

    saveName += + ".png";

    c0->SaveAs(saveName.c_str());

    delete Probability;
    delete c0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredFilteredFilteredProbabilityComparison(std::vector<std::pair<float, float>> cutValues1, std::vector<std::pair<float, float>> cutValues2, std::string tokenValue, TTree* tree, std::string name, int neutrinoMode, int yFaceOnly, std::string direction, std::string variableOne, float histogramRange, int &nBins, float tailDefinition, std::string filterVariable, std::string plotTitle)
{
    for (auto &pair : cutValues2)
    {
        ProbabilityComparison(cutValues1, true, pair.first, pair.second, tokenValue, tree, name, neutrinoMode, yFaceOnly, direction, variableOne, histogramRange, nBins, tailDefinition, filterVariable, plotTitle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProbabilityYFaceComparison(TTree* tree, std::string name, int neutrinoMode, std::string direction, std::string variableOne, float histogramRange, int &nBins, float tailDefinition, std::string plotTitle)
{
    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    auto legend = new TLegend(0.55,0.175,0.875,0.375);
    legend->SetHeader("Colour Legend"); 

    TGraphAsymmErrors *Probability = new TGraphAsymmErrors();
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(kRed);
    Probability->SetLineColor(kRed);
    Probability->SetMinimum(0.0);
    Probability->SetMaximum(1.0);
    Probability->SetTitle(plotTitle.c_str());

    CreateProbabilityPlot(tree, "_all_", neutrinoMode, 0, direction, variableOne, Probability, histogramRange, nBins, tailDefinition, false, "MCPhi", 0.0, 1.0, false, 0.0, 1.0, "xTitle1", "title1", "xTitle2", "title2", false);
    legend->AddEntry(Probability,"Cosmics: All Entries","p");
    Probability->Draw("APL");
    Probability->GetXaxis()->SetLimits(-1.0, histogramRange);
    Probability->Draw("APL");

    TGraphAsymmErrors *Probability_YFace = new TGraphAsymmErrors();
    Probability_YFace->SetMarkerStyle(6);
    Probability_YFace->SetMarkerColor(kBlue);
    Probability_YFace->SetLineColor(kBlue);

    CreateProbabilityPlot(tree, "_yface_", neutrinoMode, 1, direction, variableOne, Probability_YFace, histogramRange, nBins, tailDefinition, false, "MCPhi", 0.0, 1.0, false, 0.0, 1.0, "xTitle1", "title1", "xTitle2", "title2", false);

    legend->AddEntry(Probability_YFace,"Cosmics: Tracks Intersect Top Face","p");
    Probability_YFace->Draw("PLSame");

    legend->Draw("same");
    std::string saveName("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/comparisons/Probability_comparison_" + name + ".png");
    c0->SaveAs(saveName.c_str());

    delete Probability;
    delete c0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SingleProbabilityCurve(TTree* tree, std::string name, int neutrinoMode, std::string direction, std::string variableOne, float histogramRange, int &nBins, float tailDefinition, std::string plotTitle, bool drawUnderlying)
{
    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    auto legend = new TLegend(0.65,0.175,0.85,0.375);
    legend->SetHeader("Colour Legend"); 

    TGraphAsymmErrors *Probability = new TGraphAsymmErrors();
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(kBlack);
    Probability->SetLineColor(kBlack);
    Probability->SetMinimum(0.0);
    Probability->SetMaximum(1.0);
    Probability->SetTitle(plotTitle.c_str());

    CreateProbabilityPlot(tree, name, neutrinoMode, 0, direction, variableOne, Probability, histogramRange, nBins, tailDefinition, false, "MCPhi", 0.0, 1.0, false, 0.0, 1.0, "#Delta#chi_{DU}^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N", "Distribution of #Delta#chi_{DU}^{2}/N", "xTitle2", "title2", drawUnderlying);
    legend->AddEntry(Probability,"Cosmics: All Entries","p");
    Probability->Draw("APL");
    Probability->GetXaxis()->SetLimits(-1.0, histogramRange);
    Probability->Draw("APL");
    //legend->Draw("same");
    
    std::string saveName("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/comparisons/Single_curve_" + name + ".png");
    c0->SaveAs(saveName.c_str());

    delete Probability;
    delete c0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SingleProbabilityCurve(TTree* tree, std::string name, int neutrinoMode, std::string direction, std::string variableOne, float histogramRange, int &nBins, float tailDefinition, bool enableFitting, std::vector<float> &parameterVector, std::vector<float> &errorVector, bool enableNFiltering, float NLowerBound, float NUpperBound, bool enableChiSquaredMinFiltering, float chiSquaredMinLowerBound, float chiSquaredMinUpperBound, bool simpleFit, bool constrainedFit, std::string plotTitle, bool drawUnderlying, bool saveUnderlying)
{
    TCanvas *c0 = new TCanvas("Probability", "Probability", 900, 600);
    auto legend = new TLegend(0.55,0.175,0.85,0.375);
    legend->SetHeader("Colour Legend"); 

    TGraphAsymmErrors *Probability = new TGraphAsymmErrors();
    Probability->SetMarkerStyle(6);
    Probability->SetMarkerColor(kBlack);
    Probability->SetLineColor(kBlack);
    Probability->SetMinimum(0.0);
    Probability->SetMaximum(1.0);
    Probability->SetTitle(plotTitle.c_str());

    CreateProbabilityPlot(tree, name, neutrinoMode, 0, direction, variableOne, Probability, histogramRange, nBins, tailDefinition, enableNFiltering, "NumberHits", NLowerBound, NUpperBound, enableChiSquaredMinFiltering, chiSquaredMinLowerBound, chiSquaredMinUpperBound, "#Delta#chi_{DU}^{2}/N #equiv #chi^{2}_{D}/N - #chi^{2}_{U}/N", "Distribution of #Delta#chi_{DU}^{2}/N", "xTitle2", "title2", drawUnderlying);
    Probability->Draw("APL");
    Probability->GetXaxis()->SetLimits(-1.0, histogramRange);
    Probability->Draw("APL");

    if (enableFitting)
    {
        if (simpleFit)
        {
            TF1  *f0 = new TF1("f0", "0.5 + ([0] - 0.5) * (1 - exp(-[1] * x)) * exp(-[2] * x)", 0, histogramRange);
            f0->SetParName(0, "height");
            f0->SetParName(1, "alpha");
            f0->SetParName(2, "beta");
            f0->SetParameter(0, 0.01);
            f0->SetParameter(1, 0.01);
            f0->SetParameter(2, 0.01);
            Probability->Fit("f0", "qr");

            TF1 *fit = Probability->GetFunction("f0");
            float parameter0(fit->GetParameter(0)), parameter1(fit->GetParameter(1)), parameter2(fit->GetParameter(2));
            float error0(fit->GetParError(0)), error1(fit->GetParError(1)), error2(fit->GetParError(2));

            parameterVector.push_back(parameter0);
            parameterVector.push_back(parameter1);
            parameterVector.push_back(parameter2);

            //std::cout << "Simple fit P0: " << parameter0 << std::endl;
            //std::cout << "Simple fit alpha: " << parameter1 << std::endl;
            //std::cout << "Simple fit beta: " << parameter2 << std::endl;

            errorVector.push_back(error0);
            errorVector.push_back(error1);
            errorVector.push_back(error2);

            fit->SetLineColor(kRed);
            fit->Draw("same");
            legend->AddEntry(fit, "Simple Direct Fit","l");
        }
        if (constrainedFit) 
        {
            //Fixed through 1.0
            //TF1  *f1 = new TF1("f1", "0.5 + (((([1] + ([1] + [2])*(([1] + [2])/[2])^([2]/[1])))/(2*[1])) - 0.5) * (1 - exp(-[1] * x)) * exp(-[2] * x)", 0, histogramRange);
            //TF1  *f1 = new TF1("f1", "0.5 + (((2*[1]*[0]*(([1] + [2])/[2])^([2]/[1]) + 2*[2]*[0]*(([1] + [2])/[2])^([2]/[1]) - [1]*(([1] + [2])/[2])^([2]/[1]) - [2]*(([1] + [2])/[2])^([2]/[1]) + [1])/(2*[1])) - 0.5) * (1 - exp(-[1] * x)) * exp(-[2] * x)", 0, histogramRange);
            //TF1  *f1 = new TF1("f1", "0.5 + ((([0] - 0.5) * (( ([1] + [2]) * ((([1] + [2])/[1])^([2]/[1]))  )/[2])) * (1 - exp(-[1] * x)) * (exp(-[2] * x)))", 0, histogramRange);
            TF1  *f1 = new TF1("f1", "0.5 + ( ( ([0] - 0.5) * ( (  ([1] + [2]) * (( ([1] + [2])/([2]) )^([2]/[1]))  )/([1])  ) ) * ( 1 - exp(-[1] * x) ) * ( exp(-[2] * x) ) )", 0, histogramRange);
            f1->SetParName(0, "height");
            f1->SetParName(1, "alpha");
            f1->SetParName(2, "beta");
            f1->SetParameter(0, 0.51);
            f1->SetParLimits(0, 0.5, 1.0);
            f1->SetParameter(1, 0.01);
            f1->SetParameter(2, 0.01);
            Probability->Fit("f1", "qr");

            TF1 *fit = Probability->GetFunction("f1");
            float parameter0(fit->GetParameter(0)), parameter1(fit->GetParameter(1)), parameter2(fit->GetParameter(2));
            float error0(fit->GetParError(0)), error1(fit->GetParError(1)), error2(fit->GetParError(2));

            parameterVector.push_back(parameter0);
            parameterVector.push_back(parameter1);
            parameterVector.push_back(parameter2);

            //std::cout << "Constrained fit P_max: " << parameter0 << std::endl;
            //std::cout << "Constrained fit alpha: " << parameter1 << std::endl;
            //std::cout << "Constrained fit beta: " << parameter2 << std::endl;

            errorVector.push_back(error0);
            errorVector.push_back(error1);
            errorVector.push_back(error2);

            fit->SetLineColor(kBlue);
            fit->Draw("same");
            legend->AddEntry(fit, "Constrained Direct Fit","l");
        }

        if (p0_par1 != 0 && alpha_par1 != 0 && beta_par1 != 0)
        {
            //TF1  *f2 = new TF1("f2", "0.5 + (((2*[1]*[0]*(([1] + [2])/[2])^([2]/[1]) + 2*[2]*[0]*(([1] + [2])/[2])^([2]/[1]) - [1]*(([1] + [2])/[2])^([2]/[1]) - [2]*(([1] + [2])/[2])^([2]/[1]) + [1])/(2*[1])) - 0.5) * (1 - exp(-[1] * x)) * exp(-[2] * x)", 0, histogramRange);
            TF1  *f2 = new TF1("f2", "0.5 + ((([0] - 0.5) * (( ([1] + [2]) * ((([1] + [2])/[2])^([2]/[1]))  )/[1])) * (1 - exp(-[1] * x)) * (exp(-[2] * x)))", 0, histogramRange);
            float N_mean(NLowerBound + ((NUpperBound - NLowerBound)/2)), chi_mean(chiSquaredMinLowerBound + ((chiSquaredMinUpperBound - chiSquaredMinLowerBound)/2));
            float p0(std::max(p0_par1 + p0_par3 * chi_mean + p0_par2 * N_mean, 0.5f)), alpha(std::max(alpha_par1  + alpha_par3 * chi_mean + alpha_par2 * N_mean, 0.1f)), beta(std::max(beta_par1 + beta_par3 * chi_mean + beta_par2 * N_mean, 0.01f));

            f2->FixParameter(0, p0);
            f2->FixParameter(1, alpha);
            f2->FixParameter(2, beta);
     
            f2->SetLineColor(kMagenta);

            if (enableNFiltering && enableChiSquaredMinFiltering)    
                f2->Draw("same");

            legend->AddEntry(f2, "Constrained Fit From Parameter Surface","l");
        }
        if (p0_par1_simple != 0 && alpha_par1_simple != 0 && beta_par1_simple != 0)
        {
            TF1  *f3 = new TF1("f3", "0.5 + ([0] - 0.5) * (1 - exp(-[1] * x)) * exp(-[2] * x)", 0, histogramRange);
            float N_mean(NLowerBound + ((NUpperBound - NLowerBound)/2)), chi_mean(chiSquaredMinLowerBound + ((chiSquaredMinUpperBound - chiSquaredMinLowerBound)/2));
            float p0(std::max(p0_par1_simple + p0_par3_simple * chi_mean + p0_par2_simple * N_mean, 0.5f)), alpha(std::max(alpha_par1_simple  + alpha_par3_simple * chi_mean + alpha_par2_simple * N_mean, 0.1f)), beta(std::max(beta_par1_simple + beta_par3_simple * chi_mean + beta_par2_simple * N_mean, 0.01f));

            f3->FixParameter(0, p0);
            f3->FixParameter(1, alpha);
            f3->FixParameter(2, beta);
     
            f3->SetLineColor(kMagenta);

            //legend->AddEntry(f3, "Simple Fit From Parameter Surface","l");

            //if (enableNFiltering && enableChiSquaredMinFiltering)    
            //    f3->Draw("same");
        }
    }

    if ((p0_par1 != 0 && alpha_par1 != 0 && beta_par1 != 0) && (p0_par1 != 0 && alpha_par1 != 0 && beta_par1 != 0))
        legend->Draw("same");

    if (saveUnderlying)
    {
        std::stringstream stream;
        stream << fixed << setprecision(0) << NLowerBound;
        string N_low = stream.str();

        std::stringstream stream2;
        stream2 << fixed << setprecision(0) << NUpperBound;
        string N_high = stream2.str();

        std::stringstream stream3;
        stream3 << fixed << setprecision(1) << chiSquaredMinLowerBound;
        string chsm_low = stream3.str();

        std::stringstream stream4;
        stream4 << fixed << setprecision(1) << chiSquaredMinUpperBound;
        string chsm_high = stream4.str();

        std::string saveName("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/fits/underlying/" + name + "_N_" + stream.str() + "_" + stream2.str() + "_" + "chisquaredmin" + "_" + stream3.str() + "_" + stream4.str() + ".png");

        if (enableNFiltering && !enableChiSquaredMinFiltering)    
            saveName = ("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/fits/underlying/" + name + "_N_" + stream.str() + "_" + stream2.str() + ".png");

        if (!enableNFiltering && enableChiSquaredMinFiltering)    
            saveName = ("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/fits/underlying/" + name + "_chisquaredmin_" + stream3.str() + "_" + stream4.str() + ".png");
    
        std::string simpleString(simpleFit ? "simplefit" : ""), constrainedString(constrainedFit ? "constrainedfit" : "");

        if (!enableNFiltering && !enableChiSquaredMinFiltering)    
            saveName = ("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/fits/" + name + "_all_" + simpleString + "_" + constrainedString + ".png");

        c0->SaveAs(saveName.c_str());
    }

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

    TGraphErrors *cutValueWithPurity = new TGraphErrors(maxCut/cutStep + 1);
    TGraphErrors *cutValueWithCompleteness = new TGraphErrors(maxCut/cutStep + 1);

    cutValueWithPurity->SetMarkerStyle(6);
    cutValueWithPurity->SetMarkerColor(kBlue);
    std::string title("Purity and Completeness as a function of a cut on " + tokenValue +";Cut T on " + tokenValue + ": " + tokenValue + "#geq T" + ";Cut Purity/Efficiency");
    if (comparison == "lesser")
        title = ("Purity and Completeness as a function of a cut on " + tokenValue +";Cut T on " + tokenValue + ": " + tokenValue + "#leq T" + ";Cut Purity/Efficiency");
    cutValueWithPurity->SetTitle(title.c_str());


    cutValueWithCompleteness->SetMarkerStyle(6);
    cutValueWithCompleteness->SetMarkerColor(kRed);

    int n(1);

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
    cutValueWithPurity->Draw("APL");
    cutValueWithCompleteness->Draw("Psame");
    legend->Draw("same");

    std::string saveName("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/Purity_completeness_" + variableName + ".png");
    canvas->SaveAs(saveName.c_str());
    delete canvas;
    delete cutValueWithPurity;
    delete cutValueWithCompleteness;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DrawParameterCurve(TGraphErrors* outputGraph, int fittingMode, std::vector<std::pair<float, float>> cutValues, std::string plotName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    auto firstPair(*(cutValues.begin()));
    int nBinSize(firstPair.second - firstPair.first);
    int nSteps(cutValues.size());

    TCanvas *canvas = new TCanvas(plotName.c_str(), plotName.c_str(), 900, 600);
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/fits/Parameter_curve_" + plotName + ".png");
    outputGraph->SetMarkerStyle(7);
    outputGraph->SetMarkerSize(1.5);
    outputGraph->SetTitle(plotTitle.c_str());
    outputGraph->GetXaxis()->SetTitle(xTitle.c_str());
    outputGraph->GetYaxis()->SetTitle(yTitle.c_str());
    //outputGraph->GetYaxis()->SetRangeUser(0, std::ceil(map.at(nBinSize)));

    TAxis *ax = outputGraph->GetHistogram()->GetXaxis();
    Double_t x1 = ax->GetBinLowEdge(1);
    Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
    outputGraph->GetHistogram()->GetXaxis()->Set(nSteps+1,x1,x2);

    int k(1), precision(fittingMode == 0 ? 0 : 1);
    for (auto &pair : cutValues)
    {
        std::stringstream stream1;
        stream1 << fixed << setprecision(precision) << pair.first;
        string pair_first = stream1.str();

        std::stringstream stream2;
        stream2 << fixed << setprecision(precision) << pair.second;
        string pair_second = stream2.str();

        if (pair.second == 1e6)
            pair_second = "#infty";

        std::string binName = pair_first + "-" + pair_second;
        outputGraph->GetHistogram()->GetXaxis()->SetBinLabel(k, binName.c_str());
        ++k;
    }

    outputGraph->Draw("ALP");
    TF1 *f1 = new TF1("f1","[0] + [1]*sqrt(x) + [2]*x",0,1000);
    f1->SetParameters(1, 1, 1);
    TF1 *f2 = new TF1("f2","[0] + [1]*x",0,1000);
    f2->SetParameters(0, 0);
    f2->SetLineColor(kBlue);
    outputGraph->Fit("f1");
    outputGraph->Fit("f2");
    f1->Draw("same");
    f2->Draw("same");
    canvas->SaveAs(savePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateParameterCurves(TTree* tree, std::string name, int neutrinoMode, std::string direction, std::string variableOne, float histogramRange, int &nBins, float tailDefinition, int fittingMode, std::vector<std::pair<float, float>> cutValues, std::vector<float> parameter0values, std::vector<float> parameter1values, std::vector<float> parameter2values, std::string plotTitle, bool drawUnderlying, bool saveUnderlying)
{
    //fittingMode: 0 is N and 1 is chi swuared min
    TGraphErrors* par0Graph = new TGraphErrors();
    TGraphErrors* par1Graph = new TGraphErrors();
    TGraphErrors* par2Graph = new TGraphErrors();
    
    bool enableNFiltering(fittingMode == 0 ? true : false), enableChiSquaredMinFiltering(fittingMode == 1 ? true : false);
    std::string tokenType(fittingMode == 0 ? "N" : "#chi^{2}_{min}/N");
    std::string tokenType2(fittingMode == 0 ? "N" : "minchisquared");
    int precision(fittingMode == 0 ? 0 : 1);

    int n(0);
    for (auto &pair : cutValues)
    {
        std::stringstream stream1;
        stream1 << fixed << setprecision(precision) << pair.first;
        string pair_first = stream1.str();

        std::stringstream stream2;
        stream2 << fixed << setprecision(precision) << pair.second;
        string pair_second = stream2.str();

        std::vector<float> fitParameters;
        std::vector<float> fitErrors;

        std::string fitPlotTitle("Fitted Probability Distribution For " + pair_first + " #leq " + tokenType + " #leq " + pair_second + ";|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
        SingleProbabilityCurve(tree, name, neutrinoMode, direction, variableOne, histogramRange, nBins, tailDefinition, true, fitParameters, fitErrors, enableNFiltering, pair.first, pair.second, enableChiSquaredMinFiltering, pair.first, pair.second, false, true, fitPlotTitle, drawUnderlying, saveUnderlying);

        par0Graph->SetPoint(n, n, fitParameters.at(0));
        par1Graph->SetPoint(n, n, fitParameters.at(1));
        par2Graph->SetPoint(n, n, fitParameters.at(2));

        par0Graph->SetPointError(n, 0, fitErrors.at(0));
        par1Graph->SetPointError(n, 0, fitErrors.at(1));
        par2Graph->SetPointError(n, 0, fitErrors.at(2));

        parameter0values.push_back(fitParameters.at(0));
        parameter1values.push_back(fitParameters.at(1));
        parameter2values.push_back(fitParameters.at(2));

        ++n;
    }

    DrawParameterCurve(par0Graph, fittingMode, cutValues, "Parameter_0_Values_for_Bins_in_"+tokenType2, "Bin in "+tokenType, "P_{max} (function maximum)", "P_{max} as a function of "+tokenType);
    DrawParameterCurve(par1Graph, fittingMode, cutValues, "Parameter_1_Values_for_Bins_in_"+tokenType2, "Bin in "+tokenType, "#alpha (growth factor)", "#alpha as a function of "+tokenType);
    DrawParameterCurve(par2Graph, fittingMode, cutValues, "Parameter_2_Values_for_Bins_in_"+tokenType2, "Bin in "+tokenType, "#beta (decay factor)", "#beta as a function of "+tokenType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SetAxisLabels(TAxis* ax, std::vector<std::pair<float, float>> cutValues, int precision)
{
    auto firstPair(*(cutValues.begin()));
    int nBinSize(firstPair.second - firstPair.first);
    int nSteps(cutValues.size());

    //Double_t x1 = ax->GetBinLowEdge(1);
    //Double_t x2 = ax->GetBinUpEdge(cutValues.size());
    ax->Set(nSteps,0,nSteps);

    int k(1);
    for (auto &pair : cutValues)
    {
        std::stringstream stream1;
        stream1 << fixed << setprecision(precision) << pair.first;
        string pair_first = stream1.str();

        std::stringstream stream2;
        stream2 << fixed << setprecision(precision) << pair.second;
        string pair_second = stream2.str();

        if (pair.second == 1e6)
            pair_second = "#infty";

        std::string binName = pair_first + "-" + pair_second;
        ax->SetBinLabel(k, binName.c_str());
        ++k;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

Double_t func(Double_t *x, Double_t *par)
{
    return par[0] + par[1]*x[0] + par[2]*x[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConfigureAxisLabels(TAxis* pAxis, std::vector<std::pair<float, float>> &cutValues, int &labelNumber, int precision)
{
    for (const auto &pair : cutValues)
    {
        std::stringstream stream;
        stream << fixed << setprecision(precision) << pair.first;
        string lower = stream.str();

        std::stringstream stream2;
        stream2 << fixed << setprecision(precision) << pair.second;
        string upper = stream2.str();

        string label = lower + "-" + upper;
        pAxis->ChangeLabel(labelNumber,-1,-1,-1,-1,-1,label.c_str());
        ++labelNumber;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateParameterSurface(TTree* t1, std::vector<std::pair<float, float>> parameterNCutValues, std::vector<std::pair<float, float>> parameterChiCutValues, int parameterNumber, float histogramRange, int &nBins, float tailDefinition, bool simpleFit, bool drawSurface, std::string plotName, bool drawUnderlying)
{
    double x[parameterNCutValues.size() * parameterChiCutValues.size()], y[parameterNCutValues.size() * parameterChiCutValues.size()], z[parameterNCutValues.size() * parameterChiCutValues.size()];
    double ex[parameterNCutValues.size() * parameterChiCutValues.size()], ey[parameterNCutValues.size() * parameterChiCutValues.size()], ez[parameterNCutValues.size() * parameterChiCutValues.size()];

    float xStep((*(parameterNCutValues.begin())).second - (*(parameterNCutValues.begin())).first);
    float yStep((*(parameterChiCutValues.begin())).second - (*(parameterChiCutValues.begin())).first);

    int i(0);
    float j(xStep/2 + (*(parameterNCutValues.begin())).first);
    float maxX(0.f);

    for (auto &N_pair : parameterNCutValues)
    {
        float k(yStep/2 + (*(parameterChiCutValues.begin())).first);

        for (auto &chi_pair : parameterChiCutValues)
        {
            std::stringstream stream;
            stream << fixed << setprecision(0) << N_pair.first;
            string N_low = stream.str();

            std::stringstream stream2;
            stream2 << fixed << setprecision(0) << N_pair.second;
            string N_high = stream2.str();

            std::stringstream stream3;
            stream3 << fixed << setprecision(0) << chi_pair.first;
            string chsm_low = stream3.str();

            std::stringstream stream4;
            stream4 << fixed << setprecision(0) << chi_pair.second;
            string chsm_high = stream4.str();

            std::vector<float> fitParameters, fitErrors;

            std::string fitPlotTitle1("Fitted Probability Distribution For " + chsm_low + " #leq #chi^{2}_{min}/N #leq " + chsm_high + " and " + N_low + " #leq N #leq " + N_high + ";|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
            SingleProbabilityCurve(t1, "fitted_N_chi_range", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, nBins, tailDefinition, true, fitParameters, fitErrors, true, N_pair.first, N_pair.second, true, chi_pair.first, chi_pair.second, simpleFit, !simpleFit, fitPlotTitle1, false, drawUnderlying);

            Double_t parameterValue(fitParameters.at(parameterNumber)), parameterError(fitErrors.at(parameterNumber));
            x[i] = j;
            y[i] = k;
            z[i] = parameterValue;

            ex[i] = 0;
            ey[i] = 0;
            ez[i] = parameterError;
            //outputGraph->SetPoint(i, j, k, parameterValue);
            //outputGraph->SetPointError(i, j, k, parameterError);

            if (x[i] > maxX)
                maxX = x[i];

            ++i;
            k += yStep;
        }
        
        j += xStep; 
    }

    TCanvas *canvas = new TCanvas("par0", "par0", 900, 600);
    std::string surfaceDrawnString(drawSurface ? "with" : "without");
    std::string savePath("/usera/jjd49/pandora_direction/CondorUtilities/figures/probability/" + globalMode + "/fits/Parameter_surface_" + plotName + "_" + surfaceDrawnString + "_surface" + ".png");

    TGraph2DErrors *g2 = new TGraph2DErrors(parameterNCutValues.size() * parameterChiCutValues.size(),x,y,z,ex,ey,ez);

    std::string graphTitle("P_{max} Parameter Surface");
    if (parameterNumber == 1 && !simpleFit)
        graphTitle = ("#alpha Parameter Surface");
    if (parameterNumber == 2 && !simpleFit)
        graphTitle = ("#beta Parameter Surface");
    if (parameterNumber == 0 && simpleFit)
        graphTitle = ("P_{0} Parameter Surface (simple fit)");
    if (parameterNumber == 1 && simpleFit)
        graphTitle = ("#alpha Parameter Surface (simple fit)");
    if (parameterNumber == 2 && simpleFit)
        graphTitle = ("#beta Parameter Surface (simple fit)");

    TF2 *f = new TF2("f",func,-1000,1000,-1000,1000,3);
    f->SetParameters(0.1,0.1,0.1);
    //f->SetMinimum(0.0);

    g2->Fit(f);
    TF2 *fresult = (TF2*)g2->FindObject("f");

    fresult->SetTitle(graphTitle.c_str());
    fresult->GetXaxis()->SetNdivisions(-(parameterNCutValues.size() - 1));   
    //fresult->GetXaxis()->SetLimits(0.0, 1.2*maxX);
    fresult->GetYaxis()->SetNdivisions(-(parameterChiCutValues.size() - 1));   
    //fresult->GetYaxis()->SetLimits(0.0, 1.2*maxY);

    if (drawSurface)
    {
        //X axis is N
        int labelNumber(1);
        for (const auto &pair : parameterNCutValues)
        {
            std::stringstream stream;
            stream << fixed << setprecision(0) << pair.first;
            string lower = stream.str();

            std::stringstream stream2;
            stream2 << fixed << setprecision(0) << pair.second;
            string upper = stream2.str();

            string label = lower + "-" + upper;
            fresult->GetXaxis()->ChangeLabel(labelNumber,-1,-1,-1,-1,-1,label.c_str());
            ++labelNumber;
        }

        //Y axis is min chi squared 
        labelNumber = 1;
        for (const auto &pair : parameterChiCutValues)
        {
            std::stringstream stream;
            stream << fixed << setprecision(1) << pair.first;
            string lower = stream.str();

            std::stringstream stream2;
            stream2 << fixed << setprecision(1) << pair.second;
            string upper = stream2.str();

            string label = lower + "-" + upper;
            fresult->GetYaxis()->ChangeLabel(labelNumber,-1,-1,-1,-1,-1,label.c_str());
            ++labelNumber;
        }

        fresult->GetXaxis()->SetTitleOffset(2.5);
        fresult->GetYaxis()->SetTitleOffset(2.5);
        fresult->GetZaxis()->SetTitleOffset(1.5);

        fresult->GetXaxis()->CenterTitle();
        fresult->GetYaxis()->CenterTitle();
        fresult->GetZaxis()->CenterTitle();

        fresult->GetXaxis()->SetTitle("Bins in N");
        fresult->GetYaxis()->SetTitle("Bins in #chi^{2}_{min}/N");
        fresult->GetZaxis()->SetTitle("Parameter Value");

        fresult->Draw("surf0");
        g2->Draw("same err p0");
    }
    else
    {
        g2->SetTitle(graphTitle.c_str());
        g2->GetXaxis()->SetNdivisions(-(parameterNCutValues.size() - 1));   
        g2->GetYaxis()->SetNdivisions(-(parameterChiCutValues.size() - 1));   

        //X axis is N
        int labelNumber(1);
        for (const auto &pair : parameterNCutValues)
        {
            std::stringstream stream;
            stream << fixed << setprecision(0) << pair.first;
            string lower = stream.str();

            std::stringstream stream2;
            stream2 << fixed << setprecision(0) << pair.second;
            string upper = stream2.str();

            string label = lower + "-" + upper;
            g2->GetXaxis()->ChangeLabel(labelNumber,-1,-1,-1,-1,-1,label.c_str());
            ++labelNumber;
        }

        //Y axis is min chi squared 
        labelNumber = 1;
        for (const auto &pair : parameterChiCutValues)
        {
            std::stringstream stream;
            stream << fixed << setprecision(1) << pair.first;
            string lower = stream.str();

            std::stringstream stream2;
            stream2 << fixed << setprecision(1) << pair.second;
            string upper = stream2.str();

            string label = lower + "-" + upper;
            g2->GetYaxis()->ChangeLabel(labelNumber,-1,-1,-1,-1,-1,label.c_str());
            ++labelNumber;
        }

        g2->GetXaxis()->SetTitleOffset(2.5);
        g2->GetYaxis()->SetTitleOffset(2.5);
        g2->GetZaxis()->SetTitleOffset(1.5);

        g2->GetXaxis()->CenterTitle();
        g2->GetYaxis()->CenterTitle();
        g2->GetZaxis()->CenterTitle();

        g2->GetXaxis()->SetTitle("Bins in N");
        g2->GetYaxis()->SetTitle("Bins in #chi^{2}_{min}/N");
        g2->GetZaxis()->SetTitle("Parameter Value");

        g2->Draw("err p0");
    }

    Double_t p0 = fresult->GetParameter(0);
    Double_t p1 = fresult->GetParameter(1);
    Double_t p2 = fresult->GetParameter(2);

    if (simpleFit)
    {
        if (parameterNumber == 0) 
        {
            p0_par1_simple = p0;
            p0_par2_simple = p1;
            p0_par3_simple = p2;    
        }

        if (parameterNumber == 1) 
        {
            alpha_par1_simple = p0;
            alpha_par2_simple = p1;
            alpha_par3_simple = p2;    
        }

        if (parameterNumber == 2) 
        {
            beta_par1_simple = p0;
            beta_par2_simple = p1;
            beta_par3_simple = p2;    
        }
    }
    else
    {
        if (parameterNumber == 0) 
        {
            p0_par1 = p0;
            p0_par2 = p1;
            p0_par3 = p2;    
        }

        if (parameterNumber == 1) 
        {
            alpha_par1 = p0;
            alpha_par2 = p1;
            alpha_par3 = p2;    
        }

        if (parameterNumber == 2) 
        {
            beta_par1 = p0;
            beta_par2 = p1;
            beta_par3 = p2;    
        }
    }

    //g2->Draw("err p0");
    //fresult->Draw("same surf0");

    canvas->SaveAs(savePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------


void Probability_Cosmic(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/mc_cosmic.root");

    if (globalMode == "data")
        f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/data.root");

    TTree *t1 = (TTree*)f1->Get("PFO");

    float histogramRange(15.0), tailDefinition(15.0);

    if (globalMode == "data")
    {
        histogramRange = 25.0;
        tailDefinition = 15.0;
    }
   
    int nBins(80);

    //--------------------------------------------------------------------------------------------
   
    //Single curves
    /*
    std::string singlePlotTitle("Probability Distribution For All Cosmics Events;|#Delta#chi_{DU}^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    SingleProbabilityCurve(t1, "all", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, nBins, tailDefinition, singlePlotTitle, true);
    */

    /*
    std::vector<float> fitParameters, fitErrors;
    bool simpleFit(true), constrainedFit(true);
    bool drawUnderlying(true), saveUnderlying(true);
    std::string fitPlotTitle("Fitted Probability Distribution For All Cosmics Events;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    SingleProbabilityCurve(t1, "all_fitted", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, nBins, tailDefinition, true, fitParameters, fitErrors, false, 100, 200, false, 0.0, 1.0, simpleFit, constrainedFit, fitPlotTitle, drawUnderlying, saveUnderlying);
    */

    /*
    std::string fitPlotTitle1("Fitted Probability Distribution For 0.0 #leq #chi^{2}_{min}/N #leq 1.0;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    SingleProbabilityCurve(t1, "fitted_N_100_200_minchisquared_leq_1.0", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, nBins, tailDefinition, true, fitParameters, fitErrors, false, 100, 200, true, 0.0, 1.0, fitPlotTitle1, true, true);
    */

    /*
    TH1F *pNumberHitsHistogram = new TH1F("pNumberHitsHistogram","", 75, 0, 750);
    TH1F *pMinChiSquaredPerHitHistogram = new TH1F("pMinChiSquaredPerHitHistogram","", 100, 0.0, 5.0);

    float MinChiSquaredPerHit;
    int NumberHits;
    int NeutrinoInduced;

    t1->SetBranchAddress("MinChiSquaredPerHit", &MinChiSquaredPerHit);
    t1->SetBranchAddress("NumberHits", &NumberHits);
    t1->SetBranchAddress("NeutrinoInduced", &NeutrinoInduced);

    for (int i = 0; i < t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if (NeutrinoInduced != 0)
            continue;

        pMinChiSquaredPerHitHistogram->Fill(MinChiSquaredPerHit);
        pNumberHitsHistogram->Fill(NumberHits);
    }

    Draw(pMinChiSquaredPerHitHistogram, kBlue, "Cosmic_minchisquaredperhit", "#chi^{2}_{min}/N", "Number of Entries", "Distribution of #chi^{2}_{min}/N");
    Draw(pNumberHitsHistogram, kBlue, "Cosmic_N", "Number of hits N", "Number of Entries", "Distribution of Number of Hits N");
    */

    int newNBins(30);

    //N, L, chisquaredmin comparison
    std::vector<std::pair<float, float>> NCutValues = {std::pair<float, float>(0, 100), std::pair<float, float>(100, 200), std::pair<float, float>(200, 1e6)};
    std::string NPlotTitle("Downwards Probability P_{D} Distribution Variation with N;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    ProbabilityComparison(NCutValues, false, 0.0, 1.0, "N", t1, "N_comparison", 1, 0, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, "NumberHits", NPlotTitle);

    std::vector<std::pair<float, float>> minChiCutValues = {std::pair<float, float>(0.0, 1.0), std::pair<float, float>(1.0, 2.0), std::pair<float, float>(2.0, 3.0), std::pair<float, float>(3.0, 1e6)};
    std::string minChiPlotTitle("Downwards Probability P_{D} Distribution Variation with #chi^{2}_{min}/N (Best fit #chi^{2}/N in event);|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    ProbabilityComparison(minChiCutValues, false, 0.0, 1.0, "#chi^{2}_{min}/N", t1, "minChi_comparison", 1, 0, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, "MinChiSquaredPerHit", minChiPlotTitle);

    std::vector<std::pair<float, float>> LCutValues = {std::pair<float, float>(0.0, 100.0), std::pair<float, float>(100.0, 200.0), std::pair<float, float>(200.0, 1e6)};
    std::string LPlotTitle("Downwards Probability P_{D} Distribution Variation with L (Reconstructed 3D Track Length);|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    ProbabilityComparison(LCutValues, false, 0.0, 1.0, "L", t1, "L_comparison", 1, 0, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, "recoLength", LPlotTitle);

    //Angles comparison
    std::vector<std::pair<float, float>> phiCutValues = {std::pair<float, float>(0.52, 2.62), std::pair<float, float>(1.05, 2.09)};
    std::string phiPlotTitle("Downwards Probability P_{D} Distribution Variation with #phi (angle of track with X axis);|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    ProbabilityComparison(phiCutValues, false, 0.0, 1.0, "#phi", t1, "phi_comparison", 1, 0, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, "MCPhi", phiPlotTitle);

    std::vector<std::pair<float, float>> thetaCutValues = {std::pair<float, float>(1.4, 3.14), std::pair<float, float>(2.3, 2.7), std::pair<float, float>(2.1, 2.9)};
    std::string thetaPlotTitle("Downwards Probability P_{D} Distribution Variation with #theta (angle of track with Y axis);|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    ProbabilityComparison(thetaCutValues, false, 0.0, 1.0, "#theta", t1, "theta_comparison", 1, 0, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, "MCTheta", thetaPlotTitle);

    //Comparison for fixed chi squared min

    for (auto &pair : minChiCutValues)
    {
        ProbabilityComparison(NCutValues, true, pair.first, pair.second, "N", t1, "N_comparison", 1, 0, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, "NumberHits", NPlotTitle); 
        //ProbabilityComparison(LCutValues, true, pair.first, pair.second, "L", t1, "L_comparison", 1, 0, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, "recoLength", LPlotTitle);
    }

    //Parameter curve/surface cut values
    std::vector<std::pair<float, float>> parameterNCutValues = {std::pair<float, float>(0, 100), std::pair<float, float>(100, 200), std::pair<float, float>(200, 400)};
    std::vector<std::pair<float, float>> parameterChiCutValues = {std::pair<float, float>(0.0, 1.0), std::pair<float, float>(1.0, 2.0), std::pair<float, float>(2.0, 3.0)};

    /*
    //Parameter curves
    std::vector<float> parameter0ValuesN, parameter1ValuesN, parameter2ValuesN, parameter0ValuesChi, parameter1ValuesChi, parameter2ValuesChi;

    std::string parameterPlotTitle("Downwards Probability P_{D} Distribution For All Cosmics Events;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    CreateParameterCurves(t1, "fitted_curve", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, nBins, tailDefinition, 0, parameterNCutValues, parameter0ValuesN, parameter1ValuesN, parameter2ValuesN, parameterPlotTitle, false, false);

    std::string parameterPlotTitle2("Downwards Probability P_{D} Distribution For All Cosmics Events;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    CreateParameterCurves(t1, "fitted_curve", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", 15.0, nBins, tailDefinition, 1, parameterChiCutValues, parameter0ValuesChi, parameter1ValuesChi, parameter2ValuesChi, parameterPlotTitle2, false, false);
    */

    //Parameter surfaces
    
    //without surfaces
    CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 0, histogramRange, newNBins, tailDefinition, false, false, "par0", false);
    CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 1, histogramRange, newNBins, tailDefinition, false, false, "par1", false);
    CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 2, histogramRange, newNBins, tailDefinition, false, false, "par2", false);

    //with surfaces
    CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 0, histogramRange, newNBins, tailDefinition, false, true, "par0", false);
    CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 1, histogramRange, newNBins, tailDefinition, false, true, "par1", false);
    CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 2, histogramRange, newNBins, tailDefinition, false, true, "par2", false);

    //CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 0, histogramRange, newNBins, tailDefinition, true, false, "par0_simple", false);
    //CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 1, histogramRange, newNBins, tailDefinition, true, false, "par1_simple", false);
    //CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 2, histogramRange, newNBins, tailDefinition, true, false, "par2_simple", false);

    //with underlying probability curves
    CreateParameterSurface(t1, parameterNCutValues, parameterChiCutValues, 2, histogramRange, newNBins, tailDefinition, false, false, "par2", true);

    /*
    //Example blue fit curve outside of N and chi range
    //IMPORTANT: will only draw blue curves (which are created from the Parameter Surfaces) if you actually create them with CreateParameterSurface above
    std::string singleFitPlotTitle1("Fitted Probability Distribution For 0.0 #leq #chi^{2}_{min}/N #leq 0.5 and 200 #leq N #leq 275;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    std::vector<float> fitParameters1, fitErrors1;
    SingleProbabilityCurve(t1, "fitted_N_chi_range", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, true, fitParameters1, fitErrors1, true, 200, 275, true, 0.0, 0.5, false, true, singleFitPlotTitle1, false, true);

    std::string singleFitPlotTitle2("Fitted Probability Distribution For 0.0 #leq #chi^{2}_{min}/N #leq 0.5 and 200 #leq N #leq 275;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    std::vector<float> fitParameters2, fitErrors2;
    SingleProbabilityCurve(t1, "fitted_N_chi_range", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, true, fitParameters2, fitErrors2, true, 200, 275, true, 0.5, 1.0, false, true, singleFitPlotTitle2, false, true);

    std::string singleFitPlotTitle3("Fitted Probability Distribution For 1.0 #leq #chi^{2}_{min}/N #leq 1.5 and 200 #leq N #leq 275;|#Delta#chi^{2}/N| #equiv |#chi^{2}_{D}/N - #chi^{2}_{U}/N|;Downwards Probability P_{D}");
    std::vector<float> fitParameters3, fitErrors3;
    SingleProbabilityCurve(t1, "fitted_N_chi_range", 2, "MCDownwards", "UpDownDeltaChiSquaredPerHit", histogramRange, newNBins, tailDefinition, true, fitParameters3, fitErrors3, true, 200, 275, true, 1.0, 1.5, false, true, singleFitPlotTitle3, false, true);
    */

    //Distributions of parameters
    /*
    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    TH1F *pCosmicHistogram = new TH1F("pCosmicHistogram","", 100, 0.0, 100.0);
    TH1F *pNeutrinoHistogram = new TH1F("pNeutrinoHistogram","", 100, 0.0, 100.0);
    
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "MCPhi", 3.14, "phi_distribution", "Angle with X axis #phi", "Number of Entries", "Normalised Distributions of Track Angle With X Axis #phi for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "MCTheta", 3.14, "theta_distribution", "Angle with Y axis #theta", "Number of Entries", "Normalised Distributions of Track Angle With Y Axis #theta for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "NumberHits", 500, "N_distribution", "Number of Cluster Hits N", "Number of Entries", "Normalised Distributions of Number of Cluster Hits N for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "MinChiSquaredPerHit", 4.0, "minchisquared_distribution", "#chi^{2}_{min}/N", "Number of Entries", "Normalised Distributions of #chi^{2}_{min}/N for Neutrino and Cosmic Events");
    DrawNormalisedHistograms(t1, canvas, pCosmicHistogram, pNeutrinoHistogram, "recoLength", 300, "recolength_distribution", "Reconstructed endpoint-to-endpoint track length L", "Number of Entries", "Reconstructed endpoint-to-endpoint track length L");
    */


    //--------------------------------------------------------------------------------------------
}
