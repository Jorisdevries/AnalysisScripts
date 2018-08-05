#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLegend.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>

//------------------------------------------------------------------------------------------------------------------------------------------


void Sigma(void)
{
    TFile *f = TFile::Open("/usera/jjd49/lardirection_pandora/CondorUtilities/catroots/0.root"); 
    TTree *t1 = (TTree*)f->Get("Hit");
    int nentries = (Int_t)t1->GetEntries();
   
    int forwardsNumberBins(80);
    float forwardsLowerBound(0.0), forwardsUpperBound(4.0);
 
    TH1F *forwardsHitChargeValues = new TH1F("forwardsHitChargeValues","", forwardsNumberBins, 0.0, 4.0); 

    int MCDirection;
    float ChargeOverWidth;
    float HitCharge;
    float HitWidth;
    float ForwardsFitCharge;
    float ForwardsSigma;
    float ForwardsDelta;
    float ForwardsChiSquared;
    float BackwardsFitCharge;
    float BackwardsSigma;
    float BackwardsDelta;
    float BackwardsChiSquared;

    t1->SetBranchAddress("MCDirection", &MCDirection);
    t1->SetBranchAddress("ChargeOverWidth", &ChargeOverWidth);
    t1->SetBranchAddress("HitCharge", &HitCharge); 
    t1->SetBranchAddress("HitWidth", &HitWidth); 
    t1->SetBranchAddress("ForwardsFitCharge", &ForwardsFitCharge); 
    t1->SetBranchAddress("ForwardsSigma", &ForwardsSigma); 
    t1->SetBranchAddress("ForwardsDelta", &ForwardsDelta); 
    t1->SetBranchAddress("ForwardsChiSquared", &ForwardsChiSquared); 
    t1->SetBranchAddress("BackwardsFitCharge", &BackwardsFitCharge); 
    t1->SetBranchAddress("BackwardsSigma", &BackwardsSigma); 
    t1->SetBranchAddress("BackwardsDelta", &BackwardsDelta); 
    t1->SetBranchAddress("BackwardsChiSquared", &BackwardsChiSquared); 

    //--------------------------------------------------------------------------------------------------------------------------------------

    std::cout << "There are " << nentries << " entries in the tree" << std::endl;

    bool truncate(false);

    for (int i = 1; i < nentries; i++)
    {
        t1->GetEntry(i);

        if (truncate && i > 1e6)
            break;

        if (MCDirection == 1)
        {
            forwardsHitChargeValues->Fill(ForwardsFitCharge); 
        }
        //if (MCDirection == 0)
        //{
        //}
    }

    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 900, 600);
    canvas0->cd();
    forwardsHitChargeValues->SetMarkerStyle(6);
    forwardsHitChargeValues->SetMarkerColor(kGreen+1);
    forwardsHitChargeValues->SetTitle("Distribution of Q_{fit} values; Q_{fit}; Number of entries"); 
    //forwardsHitChargeValues->SetTitleOffset(1.4);
    forwardsHitChargeValues->Draw();
    canvas0->SaveAs("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/sigmascaling/Qfit.png");

    TGraphErrors *ChargeWithRMS90 = new TGraphErrors(forwardsNumberBins);    
    TGraphErrors *ChargeWithRMS60 = new TGraphErrors(forwardsNumberBins);    
    TGraphErrors *ChargeWithRMS80 = new TGraphErrors(forwardsNumberBins);    
    TGraphErrors *ChargeWithGaussianSigma = new TGraphErrors(forwardsNumberBins);    

    for (int j = 0; j < forwardsNumberBins; j++)
    {
        float binLowerEdge(j * (forwardsUpperBound - forwardsLowerBound)/forwardsNumberBins), binUpperEdge((j + 1) * (forwardsUpperBound - forwardsLowerBound)/forwardsNumberBins);
        float binCenter((binLowerEdge + binUpperEdge)/2);

        if (binCenter <= 0.7 || binCenter > 2.)
            continue;

        std::cout << "Computing bin center: " << binCenter << std::endl;

        TH1F *binDeltas = new TH1F(Form("x%d",j),"", 100 , -0.80, 0.80);

        std::vector<float> deltasVector;

        for (int i = 1; i < nentries; i++)
        {
            t1->GetEntry(i);
            
            if (truncate && i > 1e6)
                break;

            if (ForwardsFitCharge >= binLowerEdge && ForwardsFitCharge < binUpperEdge)
            {
                deltasVector.push_back(ForwardsDelta);
                binDeltas->Fill(ForwardsDelta);
            }
        }

        if (deltasVector.size() == 0 || binDeltas->GetEntries() == 0)
            continue;

        float mean(binDeltas->GetBinCenter(binDeltas->GetMaximumBin()));
        //float mean(0.f);
        std::vector<float> distancesToMeanVector;

        for (unsigned int i = 0; i < deltasVector.size(); i++)
        {
            float ForwardsDelta(deltasVector.at(i));
            float distanceToMean(std::abs(mean - ForwardsDelta));
            distancesToMeanVector.push_back(distanceToMean);
        }

        std::sort(distancesToMeanVector.begin(), distancesToMeanVector.end());

        float delineator90(distancesToMeanVector.at(0.9 * distancesToMeanVector.size()));
        float delineator80(distancesToMeanVector.at(0.80 * distancesToMeanVector.size()));
        float delineator60(distancesToMeanVector.at(0.60 * distancesToMeanVector.size()));

        float meanSummedSquares90(0.f), meanSummedSquares60(0.f), meanSummedSquares80(0.f);
        int nEntries90(0), nEntries60(0), nEntries80(0);

        float lowest80(1.0e7), highest80(0.f);
        float lowest60(1.0e7), highest60(0.f);
        float lowest90(1.0e7), highest90(0.f);

        for (unsigned int i = 0; i < deltasVector.size(); i++)
        {
            float ForwardsDelta(deltasVector.at(i));
            
            if (std::abs(mean - ForwardsDelta) <= delineator90)
            {
                if (ForwardsDelta < lowest90)
                    lowest90 = ForwardsDelta;

                if (ForwardsDelta > highest90)
                    highest90 = ForwardsDelta;

                meanSummedSquares90 += ForwardsDelta * ForwardsDelta;
                nEntries90++;            
            }

            if (std::abs(mean - ForwardsDelta) <= delineator60)
            {
                if (ForwardsDelta < lowest60)
                    lowest60 = ForwardsDelta;

                if (ForwardsDelta > highest60)
                    highest60 = ForwardsDelta;

                meanSummedSquares60 += ForwardsDelta * ForwardsDelta;
                nEntries60++;
            }

            if (std::abs(mean - ForwardsDelta) <= delineator80)
            {
                if (ForwardsDelta < lowest80)
                    lowest80 = ForwardsDelta;

                if (ForwardsDelta > highest80)
                    highest80 = ForwardsDelta;

                meanSummedSquares80 += ForwardsDelta * ForwardsDelta;
                nEntries80++;
            }
        }

        meanSummedSquares90 /= nEntries90;
        meanSummedSquares60 /= nEntries60;
        meanSummedSquares80 /= nEntries80;

        float rms90(std::sqrt(meanSummedSquares90));
        float rms60(std::sqrt(meanSummedSquares60));
        float rms80(std::sqrt(meanSummedSquares80));

        TF1 *gaussianfit = new TF1("gaussianfit","gaus", 0.0, 4.0);
        //gaussianfit->FixParameter(0, mean);
        binDeltas->Fit("gaussianfit", "", "", lowest80, highest80);
        TF1 *fit = binDeltas->GetFunction("gaussianfit");
        fit->SetLineColor(kBlue);
        float gaussianConstant(fit->GetParameter(0));
        float gaussianMean(fit->GetParameter(1));
        float gaussianSigma(fit->GetParameter(2));
        float gaussianSigmaError(fit->GetParError(2));

        TF1 *gaussianfitdraw1 = new TF1("gaussianfitdraw1","gaus", -4.0, 4.0);
        gaussianfitdraw1->SetLineColor(kRed);
        gaussianfitdraw1->FixParameter(0, gaussianConstant);
        gaussianfitdraw1->FixParameter(1, gaussianMean);
        gaussianfitdraw1->FixParameter(2, gaussianSigma);

        TF1 *gaussianfitdraw2 = new TF1("gaussianfitdraw2","gaus", lowest80, highest80);
        gaussianfitdraw2->SetLineColor(kBlue);
        gaussianfitdraw2->FixParameter(0, gaussianConstant);
        gaussianfitdraw2->FixParameter(1, gaussianMean);
        gaussianfitdraw2->FixParameter(2, gaussianSigma);

        float errorOnError90(rms90/(  std::sqrt( (float)(2*(nEntries90 - 1))   )  ));
        float errorOnError60(rms60/(  std::sqrt( (float)(2*(nEntries60 - 1))   )  ));
        float errorOnError80(rms80/(  std::sqrt( (float)(2*(nEntries80 - 1))   )  ));

        ChargeWithRMS90->SetPoint(j, binCenter, rms90);
        ChargeWithRMS90->SetPointError(j, 0.0, errorOnError90); 
        ChargeWithRMS60->SetPoint(j, binCenter, rms60);
        ChargeWithRMS60->SetPointError(j, 0.0, errorOnError60); 
        ChargeWithRMS80->SetPoint(j, binCenter, rms80);
        ChargeWithRMS80->SetPointError(j, 0.0, errorOnError80); 
        ChargeWithGaussianSigma->SetPoint(j, binCenter, gaussianSigma);
        ChargeWithGaussianSigma->SetPointError(j, 0.0, gaussianSigmaError); 
        

        stringstream ss;//create a stringstream
        ss << ((4.0/forwardsNumberBins) * (j + 0.5));//add number to the stream
        //ss.str()

        //TCanvas *c0 = new TCanvas("binDeltas", "binDeltas", 900, 600);
        binDeltas->SetXTitle("#Delta = (Q/w)_{obs} - (Q/w)_{fit}");
        binDeltas->SetYTitle("Number of entries");
        std::string title = "True Forwards Muons: Delta Distribution with Bin Centre Q_{fit} = " + ss.str();
        binDeltas->SetTitle(title.c_str());
        binDeltas->GetYaxis()->SetTitleOffset(1.4);
        //binDeltas->GetYaxis()->SetRangeUser(0.0, 1000.0);
        binDeltas->SetLineColor(kBlue);
        //binDeltas->Draw();

        TH1F *c1 = (TH1F*)binDeltas->Clone(); 
        TH1F *c2 = (TH1F*)binDeltas->Clone(); 
        TH1F *c3 = (TH1F*)binDeltas->Clone(); 

        int xlower90(c3->GetXaxis()->FindBin(lowest90)), xupper90(c3->GetXaxis()->FindBin(highest90));
        c3->SetFillColor(36);
        c3->GetXaxis()->SetRange(xlower90,xupper90);
        //c3->Draw("same");

        int xlower80(c2->GetXaxis()->FindBin(lowest80)), xupper80(c2->GetXaxis()->FindBin(highest80));
        c2->SetFillColor(26);
        c2->GetXaxis()->SetRange(xlower80,xupper80);
        //c2->Draw("same");

        int xlower(c1->GetXaxis()->FindBin(lowest60)), xupper(c1->GetXaxis()->FindBin(highest60));
        c1->SetFillColor(46);
        c1->GetXaxis()->SetRange(xlower,xupper);
        //c1->Draw("same");

        //gaussianfitdraw1->Draw("same");
        //gaussianfitdraw2->Draw("same");

        auto legend = new TLegend(0.15,0.65,0.38,0.80);
        legend->SetHeader("Colour Legend"); 
        legend->AddEntry(c1,"Inner 60% of Entries","f");
        legend->AddEntry(c2,"Inner 80% of Entries","f");
        legend->AddEntry(c3,"Inner 90% of Entries","f");
        legend->AddEntry(gaussianfitdraw2,"Gaussian Fit on Inner 80%","l");
        legend->AddEntry(gaussianfitdraw1,"Extrapolation of Gaussian Fit","l");
        //legend->Draw("same");


        //std::string name = ("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/sigmascaling/peakDeltasFit_" + ss.str() + ".png");
        //c0->SaveAs(name.c_str());
        //delete c0;      
    }

    TF1 *sqrtfit80 = new TF1("sqrtfit80","[0] *  sqrt(x)", 0.0, 3.0);
    sqrtfit80->SetParameter(0, 0.0938242);
    sqrtfit80->SetLineColor(kRed);
    ChargeWithRMS80->Fit("sqrtfit80");
    TF1 *sqrtfit80result = ChargeWithRMS80->GetFunction("sqrtfit80");
    float sqrt_par80 = sqrtfit80result->GetParameter(0);

    TF1 *linearfit80 = new TF1("linearfit80","[0] *  x", 0.0, 3.0);
    linearfit80->SetParameter(0, 0.0946206);
    linearfit80->SetLineColor(kBlue);
    ChargeWithRMS80->Fit("linearfit80");
    TF1 *linearfit80result = ChargeWithRMS80->GetFunction("linearfit80");
    float linear_par80 = linearfit80result->GetParameter(0);

    TF1 *combined_fit80 = new TF1("combined_fit80","sqrt(([0] * x * x) + ([1] * x))", 0.0, 3.0);
    combined_fit80->SetParLimits(0, 0.0, 1.0);
    combined_fit80->SetParLimits(1, 0.0, 1.0);
    combined_fit80->SetParameter(0, 0.001);
    combined_fit80->SetParameter(1, 0.001);
    combined_fit80->SetLineColor(kMagenta);
    ChargeWithRMS80->Fit("combined_fit80");
    TF1 *combined_fit80result = ChargeWithRMS80->GetFunction("combined_fit80");
    float combined_par80_1 = combined_fit80result->GetParameter(0);
    float combined_par80_2 = combined_fit80result->GetParameter(1);

    std::cout << "Sqrt 80% fit parameter: " << sqrt_par80 << std::endl;
    std::cout << "Linear 80% fit parameter: " << linear_par80 << std::endl;
    std::cout << "Combined 80% fit parameter 1 (quadratic): " << combined_par80_1 << std::endl;
    std::cout << "Combined 80% fit parameter 2 (linear): " << combined_par80_2 << std::endl;

    //---------------------------------------------------------------------------------------

    TF1 *sqrtfit60 = new TF1("sqrtfit60","[0] *  sqrt(x)", 0.0, 4.0);
    sqrtfit60->SetParameter(0, 1.0/3.5);
    sqrtfit60->SetLineColor(kRed);
    ChargeWithRMS60->Fit("sqrtfit60");
    TF1 *sqrtfit60result = ChargeWithRMS60->GetFunction("sqrtfit60");
    float sqrt_par60 = sqrtfit60result->GetParameter(0);

    TF1 *linearfit60 = new TF1("linearfit60","[0] *  x", 0.0, 4.0);
    linearfit60->SetParameter(0, 1.0/3.5);
    linearfit60->SetLineColor(kBlue);
    ChargeWithRMS60->Fit("linearfit60");
    TF1 *linearfit60result = ChargeWithRMS60->GetFunction("linearfit60");
    float linear_par60 = linearfit60result->GetParameter(0);

    TF1 *combined_fit60 = new TF1("combined_fit60","sqrt(([0] * x * x) + ([1] * x))", 0.0, 4.0);
    combined_fit60->SetParLimits(0, 0.0, 1.0);
    combined_fit60->SetParLimits(1, 0.0, 1.0);
    combined_fit60->SetParameter(0, 1.0/3.5);
    combined_fit60->SetParameter(1, 1.0/3.5);
    combined_fit60->SetLineColor(kMagenta);
    ChargeWithRMS60->Fit("combined_fit60");
    TF1 *combined_fit60result = ChargeWithRMS60->GetFunction("combined_fit60");
    float combined_par60_1 = combined_fit60result->GetParameter(0);
    float combined_par60_2 = combined_fit60result->GetParameter(1);

    std::cout << "Sqrt 60% fit60 par60ameter: " << sqrt_par60 << std::endl;
    std::cout << "Linear 60% fit60 par60ameter: " << linear_par60 << std::endl;
    std::cout << "Combined 60% fit60 par60ameter 1 (quadratic): " << combined_par60_1 << std::endl;
    std::cout << "Combined 60% fit60 par60ameter 2 (linear): " << combined_par60_2 << std::endl;

    //---------------------------------------------------------------------------------------

    TF1 *sqrtfitsigma = new TF1("sqrtfitsigma","[0] *  sqrt(x)", 0.0, 3.0);
    sqrtfitsigma->SetParameter(0, 0.144306);
    sqrtfitsigma->SetLineColor(kRed);
    ChargeWithGaussianSigma->Fit("sqrtfitsigma");
    TF1 *sqrtfitsigmaresult = ChargeWithGaussianSigma->GetFunction("sqrtfitsigma");
    float sqrt_parsigma = sqrtfitsigmaresult->GetParameter(0);

    TF1 *linearfitsigma = new TF1("linearfitsigma","[0] *  x", 0.0, 3.0);
    linearfitsigma->SetParameter(0, 0.117233);
    linearfitsigma->SetLineColor(kBlue);
    ChargeWithGaussianSigma->Fit("linearfitsigma");
    TF1 *linearfitsigmaresult = ChargeWithGaussianSigma->GetFunction("linearfitsigma");
    float linear_parsigma = linearfitsigmaresult->GetParameter(0);

    TF1 *combined_fitsigma = new TF1("combined_fitsigma","sqrt(([0] * x * x) + ([1] * x))", 0.0, 3.0);
    combined_fitsigma->SetParLimits(0, 0.0, 1.0);
    combined_fitsigma->SetParLimits(1, 0.0, 1.0);
    combined_fitsigma->SetParameter(0, 0.001);
    combined_fitsigma->SetParameter(1, 0.001);
    combined_fitsigma->SetLineColor(kMagenta);
    ChargeWithGaussianSigma->Fit("combined_fitsigma");
    TF1 *combined_fitsigmaresult = ChargeWithGaussianSigma->GetFunction("combined_fitsigma");
    float combined_parsigma1 = combined_fitsigmaresult->GetParameter(0);
    float combined_parsigma2 = combined_fitsigmaresult->GetParameter(1);


    std::cout << "******************************************" << std::endl;
    std::cout << "Sqrt Gaussian core fit parameter: " << sqrt_parsigma << std::endl;
    std::cout << "Linear Gaussian core fit parameter: " << linear_parsigma << std::endl;
    std::cout << "Combined Gaussian core fit parameter 1 (quadratic): " << combined_parsigma1 << std::endl;
    std::cout << "Combined Gaussian core fit parameter 2 (linear): " << combined_parsigma2 << std::endl;
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    TCanvas *canvas90 = new TCanvas("canvas", "canvas", 900, 600);
    canvas90->cd();
    ChargeWithRMS90->SetMarkerStyle(6);
    ChargeWithRMS90->SetMarkerColor(kRed);
    ChargeWithRMS90->SetTitle("RMS of inner 90% of Deltas; Q_{fit} Bin Center; RMS90 of Deltas"); //Bethe-Bloch Theory Fit
    //ChargeWithRMS90->SetTitleOffset(1.4);
    ChargeWithRMS90->Draw("AP");
    canvas90->SaveAs("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/sigmascaling/ChargeWithRMS90.png");

    TCanvas *canvas80 = new TCanvas("canvas80", "canvas80", 900, 600);
    canvas80->cd();
    ChargeWithRMS80->SetMarkerStyle(6);
    ChargeWithRMS80->SetMarkerColor(kBlue);
    ChargeWithRMS80->SetTitle("RMS of inner 80% of Deltas; Q_{fit} Bin Center; RMS80 of Deltas"); //Bethe-Bloch Theory Fit
    //ChargeWithRMS80->SetTitleOffset(1.4);
    ChargeWithRMS80->Draw("AP");
    combined_fit80->Draw("same");
    linearfit80->Draw("same");
    sqrtfit80->Draw("same");
    canvas80->SaveAs("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/sigmascaling/ChargeWithRMS80.png");

    TCanvas *canvassigma = new TCanvas("canvassigma", "canvassigma", 900, 600);
    canvassigma->cd();
    ChargeWithGaussianSigma->SetMarkerStyle(6);
    ChargeWithGaussianSigma->SetMarkerColor(kBlue);
    ChargeWithGaussianSigma->SetTitle("#sigma of Gaussian Fit to Delta Core; Q_{fit} Bin Center; Gaussian #sigma"); //Bethe-Bloch Theory Fit
    //ChargeWithGaussianSigma->SetTitleOffset(1.4);
    ChargeWithGaussianSigma->Draw("AP");
    combined_fitsigma->Draw("same");
    linearfitsigma->Draw("same");
    sqrtfitsigma->Draw("same");
    canvassigma->SaveAs("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/sigmascaling/ChargeWithGaussianSigma.png");

    TCanvas *canvas60 = new TCanvas("canvas60", "canvas60", 900, 600);
    canvas60->cd();
    ChargeWithRMS60->SetMarkerStyle(6);
    ChargeWithRMS60->SetMarkerColor(kGreen+1);
    ChargeWithRMS60->SetTitle("RMS of inner 60% of Deltas; Q_{fit} Bin Center; RMS60 of Deltas"); //Bethe-Bloch Theory Fit
    //ChargeWithRMS60->SetTitleOffset(1.4);
    ChargeWithRMS60->Draw("AP");
    combined_fit60->Draw("same");
    linearfit60->Draw("same");
    sqrtfit60->Draw("same");
    //doublefit->Draw("same");
    canvas60->SaveAs("/usera/jjd49/lardirection_pandora/CondorUtilities/figures/sigmascaling/ChargeWithRMS60.png");
}
