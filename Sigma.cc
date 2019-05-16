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
#include <iomanip>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>

//------------------------------------------------------------------------------------------------------------------------------------------

void Sigma(void)
{
    TFile *f = TFile::Open("/usera/jjd49/pandora_direction/Scripts/SigmaScaling.root"); 
    TTree *t1 = (TTree*)f->Get("DirectionHitTree");
    int nentries = (Int_t)t1->GetEntries();
   
    int forwardsNumberBins(30);
    float forwardsLowerBound(0.0), forwardsUpperBound(3.0);
 
    TH1F *chargeOverWidth = new TH1F("#tilde{Q} Distribution","", 100, 0.5, 2.5); 
    TH1F *forwardsHitChargeValues = new TH1F("Q_{fit Distribution}","", forwardsNumberBins, forwardsLowerBound, forwardsUpperBound); 
    TH1F *forwardsHitChargeValuesUnbinned = new TH1F("Q_{fit Distribution} (100 bins)","", 100, 0.25, 1.25); 

    float HitChargeOverWidth;
    float HitWidth;
    float FitCharge;
    float Sigma;
    float Delta;

    t1->SetBranchAddress("HitChargeOverWidth", &HitChargeOverWidth); 
    t1->SetBranchAddress("HitWidth", &HitWidth); 
    t1->SetBranchAddress("FitCharge", &FitCharge); //this is already multiplied by hit width! 
    t1->SetBranchAddress("Sigma", &Sigma); 
    t1->SetBranchAddress("Delta", &Delta); 

    //--------------------------------------------------------------------------------------------------------------------------------------

    std::cout << "There are " << nentries << " entries in the tree" << std::endl;

    bool truncate(false);

    for (int i = 1; i < nentries; i++)
    {
        t1->GetEntry(i);

        if (truncate && i > 1e6)
            break;

        chargeOverWidth->Fill(HitChargeOverWidth);
        forwardsHitChargeValues->Fill(FitCharge); 
        forwardsHitChargeValuesUnbinned->Fill(FitCharge);
    }

    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 900, 600);
    canvas0->cd();
    forwardsHitChargeValues->SetTitle("Distribution of Q_{fit} values; Q_{fit}; Number of entries"); 
    //forwardsHitChargeValues->SetTitleOffset(1.4);
    forwardsHitChargeValues->Draw();
    canvas0->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/Q_fit_Distribution.pdf");

    TCanvas *canvas01010 = new TCanvas("canvas01010", "canvas01010", 900, 600);
    canvas01010->cd();
    forwardsHitChargeValuesUnbinned->SetTitle("Distribution of Q_{fit} values; Q_{fit}; Number of entries"); 
    //forwardsHitChargeValuesUnbinned->SetTitleOffset(1.4);
    forwardsHitChargeValuesUnbinned->Draw();
    canvas01010->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/Q_fit_Distribution_Unbinned.pdf");

    TCanvas *canvas1111 = new TCanvas("canvas1111", "canvas1111", 900, 600);
    canvas1111->cd();
    chargeOverWidth->SetTitle("Distribution of #tilde{Q} values; #tilde{Q}; Number of entries"); 
    //chargeOverWidth->SetTitleOffset(1.4);
    chargeOverWidth->Draw();
    canvas1111->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/Q_over_w_Distribution_Unbinned.pdf");

    TGraphErrors *ChargeWithRMS90 = new TGraphErrors(forwardsNumberBins);    
    TGraphErrors *ChargeWithRMS50 = new TGraphErrors(forwardsNumberBins);    
    TGraphErrors *ChargeWithRMS70 = new TGraphErrors(forwardsNumberBins);    
    TGraphErrors *ChargeWithGaussianSigma = new TGraphErrors(forwardsNumberBins);    

    for (int j = 0; j < forwardsNumberBins; j++)
    {
        float binLowerEdge(j * (forwardsUpperBound - forwardsLowerBound)/forwardsNumberBins), binUpperEdge((j + 1) * (forwardsUpperBound - forwardsLowerBound)/forwardsNumberBins);
        float binCentre((binLowerEdge + binUpperEdge)/2);

        if (binCentre <= 0.55)
            continue;

        std::cout << "Computing bin center: " << binCentre << std::endl;

        stringstream ss;
        ss << fixed << setprecision(2) << binCentre;

        std::string binDeltasTitle("Bin Centre " + ss.str());
        TH1F *binDeltas = new TH1F(binDeltasTitle.c_str(),"", 100 , -0.70, 0.70);

        std::vector<float> deltasVector;

        for (int i = 1; i < nentries; i++)
        {
            t1->GetEntry(i);
            
            if (truncate && i > 1e6)
                break;

            if (FitCharge >= binLowerEdge && FitCharge < binUpperEdge)
            {
                deltasVector.push_back(Delta);
                binDeltas->Fill(Delta);
            }
        }

        if (deltasVector.size() == 0 || binDeltas->GetEntries() == 0)
            continue;

        float mean(binDeltas->GetBinCenter(binDeltas->GetMaximumBin()));
        //float mean(0.f);
        std::vector<float> distancesToMeanVector;

        for (unsigned int i = 0; i < deltasVector.size(); i++)
        {
            float Delta(deltasVector.at(i));
            float distanceToMean(std::abs(mean - Delta));
            distancesToMeanVector.push_back(distanceToMean);
        }

        std::sort(distancesToMeanVector.begin(), distancesToMeanVector.end());

        float delineator90(distancesToMeanVector.at(0.9 * distancesToMeanVector.size()));
        float delineator70(distancesToMeanVector.at(0.70 * distancesToMeanVector.size()));
        float delineator50(distancesToMeanVector.at(0.50 * distancesToMeanVector.size()));

        float meanSummedSquares90(0.f), meanSummedSquares50(0.f), meanSummedSquares70(0.f);
        int nEntries90(0), nEntries50(0), nEntries70(0);

        float lowest70(1.0e7), highest70(0.f);
        float lowest50(1.0e7), highest50(0.f);
        float lowest90(1.0e7), highest90(0.f);

        for (unsigned int i = 0; i < deltasVector.size(); i++)
        {
            float Delta(deltasVector.at(i));
            
            if (std::abs(mean - Delta) <= delineator90)
            {
                if (Delta < lowest90)
                    lowest90 = Delta;

                if (Delta > highest90)
                    highest90 = Delta;

                meanSummedSquares90 += Delta * Delta;
                nEntries90++;            
            }

            if (std::abs(mean - Delta) <= delineator50)
            {
                if (Delta < lowest50)
                    lowest50 = Delta;

                if (Delta > highest50)
                    highest50 = Delta;

                meanSummedSquares50 += Delta * Delta;
                nEntries50++;
            }

            if (std::abs(mean - Delta) <= delineator70)
            {
                if (Delta < lowest70)
                    lowest70 = Delta;

                if (Delta > highest70)
                    highest70 = Delta;

                meanSummedSquares70 += Delta * Delta;
                nEntries70++;
            }
        }

        meanSummedSquares90 /= nEntries90;
        meanSummedSquares50 /= nEntries50;
        meanSummedSquares70 /= nEntries70;

        float rms90(std::sqrt(meanSummedSquares90));
        float rms50(std::sqrt(meanSummedSquares50));
        float rms70(std::sqrt(meanSummedSquares70));

        TF1 *gaussianfit = new TF1("gaussianfit","gaus", 0.0, 4.0);
        //gaussianfit->FixParameter(0, mean);
        binDeltas->Fit("gaussianfit", "", "", lowest70, highest70);
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

        TF1 *gaussianfitdraw2 = new TF1("gaussianfitdraw2","gaus", lowest70, highest70);
        gaussianfitdraw2->SetLineColor(kBlue);
        gaussianfitdraw2->FixParameter(0, gaussianConstant);
        gaussianfitdraw2->FixParameter(1, gaussianMean);
        gaussianfitdraw2->FixParameter(2, gaussianSigma);

        float errorOnError90(rms90/(  std::sqrt( (float)(2*(nEntries90 - 1))   )  ));
        float errorOnError50(rms50/(  std::sqrt( (float)(2*(nEntries50 - 1))   )  ));
        float errorOnError70(rms70/(  std::sqrt( (float)(2*(nEntries70 - 1))   )  ));

        ChargeWithRMS90->SetPoint(j, binCentre, rms90);
        ChargeWithRMS90->SetPointError(j, 0.0, errorOnError90); 
        ChargeWithRMS50->SetPoint(j, binCentre, rms50);
        ChargeWithRMS50->SetPointError(j, 0.0, errorOnError50); 
        ChargeWithRMS70->SetPoint(j, binCentre, rms70);
        ChargeWithRMS70->SetPointError(j, 0.0, errorOnError70); 
        ChargeWithGaussianSigma->SetPoint(j, binCentre, gaussianSigma);
        ChargeWithGaussianSigma->SetPointError(j, 0.0, gaussianSigmaError); 

        TCanvas *c0 = new TCanvas("binDeltas", "binDeltas", 900, 600);
        binDeltas->SetXTitle("#Delta = #tilde{Q} - #tilde{Q}_{fit}");
        binDeltas->SetYTitle("Number of entries");
        std::string title = "True Forwards Muons: Distribution of #Delta at Bin Centre Q_{fit} = " + ss.str();
        binDeltas->SetTitle(title.c_str());
        binDeltas->GetYaxis()->SetTitleOffset(1.4);
        //binDeltas->GetYaxis()->SetRangeUser(0.0, 1000.0);
        binDeltas->SetLineColor(kBlue);
        binDeltas->Draw();

        TH1F *c1 = (TH1F*)binDeltas->Clone(); 
        TH1F *c2 = (TH1F*)binDeltas->Clone(); 
        TH1F *c3 = (TH1F*)binDeltas->Clone(); 

        int xlower90(c3->GetXaxis()->FindBin(lowest90)), xupper90(c3->GetXaxis()->FindBin(highest90));
        c3->SetFillColor(36);
        c3->GetXaxis()->SetRange(xlower90,xupper90);
        c3->Draw("same");

        int xlower70(c2->GetXaxis()->FindBin(lowest70)), xupper70(c2->GetXaxis()->FindBin(highest70));
        c2->SetFillColor(26);
        c2->GetXaxis()->SetRange(xlower70,xupper70);
        c2->Draw("same");

        int xlower(c1->GetXaxis()->FindBin(lowest50)), xupper(c1->GetXaxis()->FindBin(highest50));
        c1->SetFillColor(46);
        c1->GetXaxis()->SetRange(xlower,xupper);
        c1->Draw("same");

        gaussianfitdraw1->Draw("same");
        gaussianfitdraw2->Draw("same");

        auto legend = new TLegend(0.15,0.6,0.38,0.85);
        legend->SetHeader("Colour Legend"); 
        legend->AddEntry(c1,"Inner 50% of Entries","f");
        legend->AddEntry(c2,"Inner 70% of Entries","f");
        legend->AddEntry(c3,"Inner 90% of Entries","f");
        legend->AddEntry(gaussianfitdraw2,"Gaussian Fit on Inner 70%","l");
        legend->AddEntry(gaussianfitdraw1,"Extrapolation of Gaussian Fit","l");
        legend->Draw("same");

        std::string name = ("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/peakDeltasFit_" + ss.str() + ".pdf");
        c0->SaveAs(name.c_str());
        delete c0;      
    }

    TF1 *sqrtfit70 = new TF1("sqrtfit70","[0] *  sqrt(x)", 0.0, forwardsUpperBound);
    sqrtfit70->SetParameter(0, 0.0938242);
    sqrtfit70->SetLineColor(kRed);
    ChargeWithRMS70->Fit("sqrtfit70", "R");
    TF1 *sqrtfit70result = ChargeWithRMS70->GetFunction("sqrtfit70");
    float sqrt_par70 = sqrtfit70result->GetParameter(0);

    TF1 *linearfit70 = new TF1("linearfit70","[0] *  x", 0.0, forwardsUpperBound);
    linearfit70->SetParameter(0, 0.0946206);
    linearfit70->SetLineColor(kBlue);
    ChargeWithRMS70->Fit("linearfit70", "R");
    TF1 *linearfit70result = ChargeWithRMS70->GetFunction("linearfit70");
    float linear_par70 = linearfit70result->GetParameter(0);

    TF1 *combined_fit70 = new TF1("combined_fit70","sqrt(([0] * x * x) + ([1] * x))", 0.0, forwardsUpperBound);
    combined_fit70->SetParLimits(0, 0.0, 1.0);
    combined_fit70->SetParLimits(1, 0.0, 1.0);
    combined_fit70->SetParameter(0, 0.001);
    combined_fit70->SetParameter(1, 0.001);
    combined_fit70->SetLineColor(kMagenta);
    ChargeWithRMS70->Fit("combined_fit70", "R");
    TF1 *combined_fit70result = ChargeWithRMS70->GetFunction("combined_fit70");
    float combined_par70_1 = combined_fit70result->GetParameter(0);
    float combined_par70_2 = combined_fit70result->GetParameter(1);

    std::cout << "Sqrt 70% fit parameter: " << sqrt_par70 << std::endl;
    std::cout << "Linear 70% fit parameter: " << linear_par70 << std::endl;
    std::cout << "Combined 70% fit parameter 1 (quadratic): " << combined_par70_1 << std::endl;
    std::cout << "Combined 70% fit parameter 2 (linear): " << combined_par70_2 << std::endl;

    //---------------------------------------------------------------------------------------

    TF1 *sqrtfit50 = new TF1("sqrtfit60","[0] *  sqrt(x)", 0.0, forwardsUpperBound);
    sqrtfit50->SetParameter(0, 1.0/3.5);
    sqrtfit50->SetLineColor(kRed);
    ChargeWithRMS50->Fit("sqrtfit60", "R");
    TF1 *sqrtfit60result = ChargeWithRMS50->GetFunction("sqrtfit60");
    float sqrt_par50 = sqrtfit60result->GetParameter(0);

    TF1 *linearfit50 = new TF1("linearfit60","[0] *  x", 0.0, forwardsUpperBound);
    linearfit50->SetParameter(0, 1.0/3.5);
    linearfit50->SetLineColor(kBlue);
    ChargeWithRMS50->Fit("linearfit60", "R");
    TF1 *linearfit60result = ChargeWithRMS50->GetFunction("linearfit60");
    float linear_par50 = linearfit60result->GetParameter(0);

    TF1 *combined_fit50 = new TF1("combined_fit60","sqrt(([0] * x * x) + ([1] * x))", 0.0, forwardsUpperBound);
    combined_fit50->SetParLimits(0, 0.0, 1.0);
    combined_fit50->SetParLimits(1, 0.0, 1.0);
    combined_fit50->SetParameter(0, 1.0/3.5);
    combined_fit50->SetParameter(1, 1.0/3.5);
    combined_fit50->SetLineColor(kMagenta);
    ChargeWithRMS50->Fit("combined_fit60", "R");
    TF1 *combined_fit60result = ChargeWithRMS50->GetFunction("combined_fit60");
    float combined_par50_1 = combined_fit60result->GetParameter(0);
    float combined_par50_2 = combined_fit60result->GetParameter(1);

    std::cout << "Sqrt 50% fit50 par60ameter: " << sqrt_par50 << std::endl;
    std::cout << "Linear 50% fit50 par60ameter: " << linear_par50 << std::endl;
    std::cout << "Combined 50% fit50 par60ameter 1 (quadratic): " << combined_par50_1 << std::endl;
    std::cout << "Combined 50% fit50 par60ameter 2 (linear): " << combined_par50_2 << std::endl;

    //---------------------------------------------------------------------------------------

    TF1 *sqrtfitsigma = new TF1("sqrtfitsigma","[0] *  sqrt(x)", 0.0, forwardsUpperBound);
    sqrtfitsigma->SetParameter(0, 0.144306);
    sqrtfitsigma->SetLineColor(kRed);
    ChargeWithGaussianSigma->Fit("sqrtfitsigma", "R");
    TF1 *sqrtfitsigmaresult = ChargeWithGaussianSigma->GetFunction("sqrtfitsigma");
    float sqrt_parsigma = sqrtfitsigmaresult->GetParameter(0);

    TF1 *linearfitsigma = new TF1("linearfitsigma","[0] *  x", 0.0, forwardsUpperBound);
    linearfitsigma->SetParameter(0, 0.117233);
    linearfitsigma->SetLineColor(kBlue);
    ChargeWithGaussianSigma->Fit("linearfitsigma", "R");
    TF1 *linearfitsigmaresult = ChargeWithGaussianSigma->GetFunction("linearfitsigma");
    float linear_parsigma = linearfitsigmaresult->GetParameter(0);

    TF1 *combined_fitsigma = new TF1("combined_fitsigma","sqrt(([0] * x * x) + ([1] * x))", 0.0, forwardsUpperBound);
    combined_fitsigma->SetParLimits(0, 0.001, 1.0);
    combined_fitsigma->SetParLimits(1, 0.0, 1.0);
    combined_fitsigma->SetParameter(0, 0.001);
    combined_fitsigma->SetParameter(1, 0.001);
    combined_fitsigma->SetLineColor(kMagenta);
    ChargeWithGaussianSigma->Fit("combined_fitsigma", "R");
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
    ChargeWithRMS90->SetTitle("RMSE of inner 90% of #Delta Entries; Q_{fit} Bin Centre; RMSE of inner 90% of #Delta Entries"); //Bethe-Bloch Theory Fit
    //ChargeWithRMS90->SetTitleOffset(1.4);
    ChargeWithRMS90->Draw("AP");
    canvas90->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/ChargeWithRMS90.pdf");

    TCanvas *canvas70 = new TCanvas("canvas70", "canvas70", 900, 600);
    canvas70->cd();
    ChargeWithRMS70->SetMarkerStyle(6);
    ChargeWithRMS70->SetMarkerColor(kBlue);
    ChargeWithRMS70->SetTitle("RMSE of inner 70% of #Delta Entries; Q_{fit} Bin Centre; RMSE of inner 70% of #Delta Entries"); //Bethe-Bloch Theory Fit
    //ChargeWithRMS70->SetTitleOffset(1.4);
    ChargeWithRMS70->Draw("AP");
    combined_fit70->Draw("same");
    linearfit70->Draw("same");
    sqrtfit70->Draw("same");
    canvas70->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/ChargeWithRMS70.pdf");

    TCanvas *canvassigma = new TCanvas("canvassigma", "canvassigma", 900, 600);
    canvassigma->cd();
    ChargeWithGaussianSigma->SetMarkerStyle(6);
    ChargeWithGaussianSigma->SetMarkerColor(kBlue);
    ChargeWithGaussianSigma->SetTitle("#sigma of Gaussian Fit to Inner 70% of #Delta Entries; Q_{fit} Bin Centre; Gaussian Fit #sigma"); //Bethe-Bloch Theory Fit
    //ChargeWithGaussianSigma->SetTitleOffset(1.4);
    ChargeWithGaussianSigma->Draw("AP");
    combined_fitsigma->Draw("same");
    linearfitsigma->Draw("same");
    sqrtfitsigma->Draw("same");
    canvassigma->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/ChargeWithGaussianSigma.pdf");

    TCanvas *canvas50 = new TCanvas("canvas60", "canvas60", 900, 600);
    canvas50->cd();
    ChargeWithRMS50->SetMarkerStyle(6);
    ChargeWithRMS50->SetMarkerColor(kGreen+1);
    ChargeWithRMS50->SetTitle("RMSE of inner 50% of #Delta Entries; Q_{fit} Bin Centre; RMSE of inner 50% of #Delta Entries"); //Bethe-Bloch Theory Fit
    //ChargeWithRMS50->SetTitleOffset(1.4);
    ChargeWithRMS50->Draw("AP");
    combined_fit50->Draw("same");
    linearfit50->Draw("same");
    sqrtfit50->Draw("same");
    //doublefit->Draw("same");
    canvas50->SaveAs("/usera/jjd49/pandora_direction/Scripts/Figures/Sigma/ChargeWithRMS50.pdf");
}
