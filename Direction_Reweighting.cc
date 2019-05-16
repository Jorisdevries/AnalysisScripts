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
#include "TRandom.h"
#include "TGaxis.h"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "Drawing_Helper_Functions.h"

//------------------------------------------------------------------------------------------------------------------------------------------

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/catroots/0.root");
TTree *t1 = (TTree*)f1->Get("Validation");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/catroots/1.root");
TTree *t2 = (TTree*)f2->Get("Validation");

//------------------------------------------------------------------------------------------------------------------------------------------

void Direction_Reweighting(void)
{
    int fileIdentifier, eventNumber;

    t1->SetBranchAddress("fileIdentifier", &fileIdentifier); 
    t1->SetBranchAddress("eventNumber", &eventNumber); 
    
    std::vector<int> unmodifiedFiles;

    for (int i = 0; i < t1->GetEntries(); i++)
    {    
        t1->GetEntry(i);

        if (std::find(unmodifiedFiles.begin(), unmodifiedFiles.end(), fileIdentifier) == unmodifiedFiles.end())
            unmodifiedFiles.push_back(fileIdentifier);
    }

    std::cout << "Loop 1 done." << std::endl;

    //////

    int fileIdentifier_direction, eventNumber_direction;

    t2->SetBranchAddress("fileIdentifier", &fileIdentifier_direction); 
    t2->SetBranchAddress("eventNumber", &eventNumber_direction); 
    
    std::vector<int> directionFiles;

    for (int i = 0; i < t2->GetEntries(); i++)
    {    
        t2->GetEntry(i);

        if (std::find(directionFiles.begin(), directionFiles.end(), fileIdentifier_direction) == directionFiles.end())
            directionFiles.push_back(fileIdentifier_direction);
    }

    std::cout << "Loop 2 done." << std::endl;

    std::vector<int> diff; //in unmodified but not in direction

    std::set_difference(unmodifiedFiles.begin(), unmodifiedFiles.end(), directionFiles.begin(), directionFiles.end(), std::inserter(diff, diff.begin()));

    for (const auto i : diff)
        std::cout << "Missing file: " << i << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
