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

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

//------------------------------------------------------------------------------------------------------------------------------------------

//SETTINGS
bool applyTrackFilter(false);

TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/0.root");
TTree *t1 = (TTree*)f1->Get("EventSelection");

TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/1.root");
TTree *t2 = (TTree*)f2->Get("EventSelection");

TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/2.root");
TTree *t3 = (TTree*)f3->Get("EventSelection");

TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/3.root");
TTree *t4 = (TTree*)f4->Get("EventSelection");

//------------------------------------------------------------------------------------------------------------------------------------------

enum InteractionType : int
{
    CCQEL_MU,
    CCQEL_MU_P,
    CCQEL_MU_P_P,
    CCQEL_MU_P_P_P,
    CCQEL_MU_P_P_P_P,
    CCQEL_MU_P_P_P_P_P,
    CCQEL_E,
    CCQEL_E_P,
    CCQEL_E_P_P,
    CCQEL_E_P_P_P,
    CCQEL_E_P_P_P_P,
    CCQEL_E_P_P_P_P_P,
    NCQEL_P,
    NCQEL_P_P,
    NCQEL_P_P_P,
    NCQEL_P_P_P_P,
    NCQEL_P_P_P_P_P,
    CCRES_MU,
    CCRES_MU_P,
    CCRES_MU_P_P,
    CCRES_MU_P_P_P,
    CCRES_MU_P_P_P_P,
    CCRES_MU_P_P_P_P_P,
    CCRES_MU_PIPLUS,
    CCRES_MU_P_PIPLUS,
    CCRES_MU_P_P_PIPLUS,
    CCRES_MU_P_P_P_PIPLUS,
    CCRES_MU_P_P_P_P_PIPLUS,
    CCRES_MU_P_P_P_P_P_PIPLUS,
    CCRES_MU_PHOTON,
    CCRES_MU_P_PHOTON,
    CCRES_MU_P_P_PHOTON,
    CCRES_MU_P_P_P_PHOTON,
    CCRES_MU_P_P_P_P_PHOTON,
    CCRES_MU_P_P_P_P_P_PHOTON,
    CCRES_MU_PIZERO,
    CCRES_MU_P_PIZERO,
    CCRES_MU_P_P_PIZERO,
    CCRES_MU_P_P_P_PIZERO,
    CCRES_MU_P_P_P_P_PIZERO,
    CCRES_MU_P_P_P_P_P_PIZERO,
    CCRES_E,
    CCRES_E_P,
    CCRES_E_P_P,
    CCRES_E_P_P_P,
    CCRES_E_P_P_P_P,
    CCRES_E_P_P_P_P_P,
    CCRES_E_PIPLUS,
    CCRES_E_P_PIPLUS,
    CCRES_E_P_P_PIPLUS,
    CCRES_E_P_P_P_PIPLUS,
    CCRES_E_P_P_P_P_PIPLUS,
    CCRES_E_P_P_P_P_P_PIPLUS,
    CCRES_E_PHOTON,
    CCRES_E_P_PHOTON,
    CCRES_E_P_P_PHOTON,
    CCRES_E_P_P_P_PHOTON,
    CCRES_E_P_P_P_P_PHOTON,
    CCRES_E_P_P_P_P_P_PHOTON,
    CCRES_E_PIZERO,
    CCRES_E_P_PIZERO,
    CCRES_E_P_P_PIZERO,
    CCRES_E_P_P_P_PIZERO,
    CCRES_E_P_P_P_P_PIZERO,
    CCRES_E_P_P_P_P_P_PIZERO,
    NCRES_P,
    NCRES_P_P,
    NCRES_P_P_P,
    NCRES_P_P_P_P,
    NCRES_P_P_P_P_P,
    NCRES_PIPLUS,
    NCRES_P_PIPLUS,
    NCRES_P_P_PIPLUS,
    NCRES_P_P_P_PIPLUS,
    NCRES_P_P_P_P_PIPLUS,
    NCRES_P_P_P_P_P_PIPLUS,
    NCRES_PIMINUS,
    NCRES_P_PIMINUS,
    NCRES_P_P_PIMINUS,
    NCRES_P_P_P_PIMINUS,
    NCRES_P_P_P_P_PIMINUS,
    NCRES_P_P_P_P_P_PIMINUS,
    NCRES_PHOTON,
    NCRES_P_PHOTON,
    NCRES_P_P_PHOTON,
    NCRES_P_P_P_PHOTON,
    NCRES_P_P_P_P_PHOTON,
    NCRES_P_P_P_P_P_PHOTON,
    NCRES_PIZERO,
    NCRES_P_PIZERO,
    NCRES_P_P_PIZERO,
    NCRES_P_P_P_PIZERO,
    NCRES_P_P_P_P_PIZERO,
    NCRES_P_P_P_P_P_PIZERO,
    CCDIS_MU,
    CCDIS_MU_P,
    CCDIS_MU_P_P,
    CCDIS_MU_P_P_P,
    CCDIS_MU_P_P_P_P,
    CCDIS_MU_P_P_P_P_P,
    CCDIS_MU_PIPLUS,
    CCDIS_MU_P_PIPLUS,
    CCDIS_MU_P_P_PIPLUS,
    CCDIS_MU_P_P_P_PIPLUS,
    CCDIS_MU_P_P_P_P_PIPLUS,
    CCDIS_MU_P_P_P_P_P_PIPLUS,
    CCDIS_MU_PHOTON,
    CCDIS_MU_P_PHOTON,
    CCDIS_MU_P_P_PHOTON,
    CCDIS_MU_P_P_P_PHOTON,
    CCDIS_MU_P_P_P_P_PHOTON,
    CCDIS_MU_P_P_P_P_P_PHOTON,
    CCDIS_MU_PIZERO,
    CCDIS_MU_P_PIZERO,
    CCDIS_MU_P_P_PIZERO,
    CCDIS_MU_P_P_P_PIZERO,
    CCDIS_MU_P_P_P_P_PIZERO,
    CCDIS_MU_P_P_P_P_P_PIZERO,
    NCDIS_P,
    NCDIS_P_P,
    NCDIS_P_P_P,
    NCDIS_P_P_P_P,
    NCDIS_P_P_P_P_P,
    NCDIS_PIPLUS,
    NCDIS_P_PIPLUS,
    NCDIS_P_P_PIPLUS,
    NCDIS_P_P_P_PIPLUS,
    NCDIS_P_P_P_P_PIPLUS,
    NCDIS_P_P_P_P_P_PIPLUS,
    NCDIS_PIMINUS,
    NCDIS_P_PIMINUS,
    NCDIS_P_P_PIMINUS,
    NCDIS_P_P_P_PIMINUS,
    NCDIS_P_P_P_P_PIMINUS,
    NCDIS_P_P_P_P_P_PIMINUS,
    NCDIS_PHOTON,
    NCDIS_P_PHOTON,
    NCDIS_P_P_PHOTON,
    NCDIS_P_P_P_PHOTON,
    NCDIS_P_P_P_P_PHOTON,
    NCDIS_P_P_P_P_P_PHOTON,
    NCDIS_PIZERO,
    NCDIS_P_PIZERO,
    NCDIS_P_P_PIZERO,
    NCDIS_P_P_P_PIZERO,
    NCDIS_P_P_P_P_PIZERO,
    NCDIS_P_P_P_P_P_PIZERO,
    CCCOH,
    NCCOH,
    COSMIC_RAY_MU,
    COSMIC_RAY_P,
    COSMIC_RAY_E,
    COSMIC_RAY_PHOTON,
    COSMIC_RAY_OTHER,
    BEAM_PARTICLE_MU,
    BEAM_PARTICLE_P,
    BEAM_PARTICLE_E,
    BEAM_PARTICLE_PHOTON,
    BEAM_PARTICLE_PI_PLUS,
    BEAM_PARTICLE_PI_MINUS,
    BEAM_PARTICLE_KAON_PLUS,
    BEAM_PARTICLE_KAON_MINUS,
    BEAM_PARTICLE_OTHER,
    OTHER_INTERACTION,
    ALL_INTERACTIONS,
    COSMIC_RAY,
    NO_RECONSTRUCTABLE
};

//------------------------------------------------------------------------------------------------------------------------------------------

std::string ToString(const InteractionType interactionType)
{
    if (static_cast<int>(interactionType) == 166)
        return "NO_RECONSTRUCTABLE";

    if (static_cast<int>(interactionType) == 165)
        return "COSMIC_RAY";

    switch (interactionType)
    {
    case CCQEL_MU: return "CCQEL_MU";
    case CCQEL_MU_P: return "CCQEL_MU_P";
    case CCQEL_MU_P_P: return "CCQEL_MU_P_P";
    case CCQEL_MU_P_P_P: return "CCQEL_MU_P_P_P";
    case CCQEL_MU_P_P_P_P: return "CCQEL_MU_P_P_P_P";
    case CCQEL_MU_P_P_P_P_P: return "CCQEL_MU_P_P_P_P_P";
    case CCQEL_E: return "CCQEL_E";
    case CCQEL_E_P: return "CCQEL_E_P";
    case CCQEL_E_P_P: return "CCQEL_E_P_P";
    case CCQEL_E_P_P_P: return "CCQEL_E_P_P_P";
    case CCQEL_E_P_P_P_P: return "CCQEL_E_P_P_P_P";
    case CCQEL_E_P_P_P_P_P: return "CCQEL_E_P_P_P_P_P";
    case NCQEL_P: return "NCQEL_P";
    case NCQEL_P_P: return "NCQEL_P_P";
    case NCQEL_P_P_P: return "NCQEL_P_P_P";
    case NCQEL_P_P_P_P: return "NCQEL_P_P_P_P";
    case NCQEL_P_P_P_P_P: return "NCQEL_P_P_P_P_P";
    case CCRES_MU: return "CCRES_MU";
    case CCRES_MU_P: return "CCRES_MU_P";
    case CCRES_MU_P_P: return "CCRES_MU_P_P";
    case CCRES_MU_P_P_P: return "CCRES_MU_P_P_P";
    case CCRES_MU_P_P_P_P: return "CCRES_MU_P_P_P_P";
    case CCRES_MU_P_P_P_P_P: return "CCRES_MU_P_P_P_P_P";
    case CCRES_MU_PIPLUS: return "CCRES_MU_PIPLUS";
    case CCRES_MU_P_PIPLUS: return "CCRES_MU_P_PIPLUS";
    case CCRES_MU_P_P_PIPLUS: return "CCRES_MU_P_P_PIPLUS";
    case CCRES_MU_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_P_PIPLUS";
    case CCRES_MU_PHOTON: return "CCRES_MU_PHOTON";
    case CCRES_MU_P_PHOTON: return "CCRES_MU_P_PHOTON";
    case CCRES_MU_P_P_PHOTON: return "CCRES_MU_P_P_PHOTON";
    case CCRES_MU_P_P_P_PHOTON: return "CCRES_MU_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_P_PHOTON";
    case CCRES_MU_PIZERO: return "CCRES_MU_PIZERO";
    case CCRES_MU_P_PIZERO: return "CCRES_MU_P_PIZERO";
    case CCRES_MU_P_P_PIZERO: return "CCRES_MU_P_P_PIZERO";
    case CCRES_MU_P_P_P_PIZERO: return "CCRES_MU_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_P_PIZERO";
    case CCRES_E: return "CCRES_E";
    case CCRES_E_P: return "CCRES_E_P";
    case CCRES_E_P_P: return "CCRES_E_P_P";
    case CCRES_E_P_P_P: return "CCRES_E_P_P_P";
    case CCRES_E_P_P_P_P: return "CCRES_E_P_P_P_P";
    case CCRES_E_P_P_P_P_P: return "CCRES_E_P_P_P_P_P";
    case CCRES_E_PIPLUS: return "CCRES_E_PIPLUS";
    case CCRES_E_P_PIPLUS: return "CCRES_E_P_PIPLUS";
    case CCRES_E_P_P_PIPLUS: return "CCRES_E_P_P_PIPLUS";
    case CCRES_E_P_P_P_PIPLUS: return "CCRES_E_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_P_PIPLUS";
    case CCRES_E_PHOTON: return "CCRES_E_PHOTON";
    case CCRES_E_P_PHOTON: return "CCRES_E_P_PHOTON";
    case CCRES_E_P_P_PHOTON: return "CCRES_E_P_P_PHOTON";
    case CCRES_E_P_P_P_PHOTON: return "CCRES_E_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_P_PHOTON";
    case CCRES_E_PIZERO: return "CCRES_E_PIZERO";
    case CCRES_E_P_PIZERO: return "CCRES_E_P_PIZERO";
    case CCRES_E_P_P_PIZERO: return "CCRES_E_P_P_PIZERO";
    case CCRES_E_P_P_P_PIZERO: return "CCRES_E_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_P_PIZERO";
    case NCRES_P: return "NCRES_P";
    case NCRES_P_P: return "NCRES_P_P";
    case NCRES_P_P_P: return "NCRES_P_P_P";
    case NCRES_P_P_P_P: return "NCRES_P_P_P_P";
    case NCRES_P_P_P_P_P: return "NCRES_P_P_P_P_P";
    case NCRES_PIPLUS: return "NCRES_PIPLUS";
    case NCRES_P_PIPLUS: return "NCRES_P_PIPLUS";
    case NCRES_P_P_PIPLUS: return "NCRES_P_P_PIPLUS";
    case NCRES_P_P_P_PIPLUS: return "NCRES_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_P_PIPLUS";
    case NCRES_PIMINUS: return "NCRES_PIMINUS";
    case NCRES_P_PIMINUS: return "NCRES_P_PIMINUS";
    case NCRES_P_P_PIMINUS: return "NCRES_P_P_PIMINUS";
    case NCRES_P_P_P_PIMINUS: return "NCRES_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_P_PIMINUS";
    case NCRES_PHOTON: return "NCRES_PHOTON";
    case NCRES_P_PHOTON: return "NCRES_P_PHOTON";
    case NCRES_P_P_PHOTON: return "NCRES_P_P_PHOTON";
    case NCRES_P_P_P_PHOTON: return "NCRES_P_P_P_PHOTON";
    case NCRES_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_PHOTON";
    case NCRES_P_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_P_PHOTON";
    case NCRES_PIZERO: return "NCRES_PIZERO";
    case NCRES_P_PIZERO: return "NCRES_P_PIZERO";
    case NCRES_P_P_PIZERO: return "NCRES_P_P_PIZERO";
    case NCRES_P_P_P_PIZERO: return "NCRES_P_P_P_PIZERO";
    case NCRES_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_PIZERO";
    case NCRES_P_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_P_PIZERO";
    case CCDIS_MU: return "CCDIS_MU";
    case CCDIS_MU_P: return "CCDIS_MU_P";
    case CCDIS_MU_P_P: return "CCDIS_MU_P_P";
    case CCDIS_MU_P_P_P: return "CCDIS_MU_P_P_P";
    case CCDIS_MU_P_P_P_P: return "CCDIS_MU_P_P_P_P";
    case CCDIS_MU_P_P_P_P_P: return "CCDIS_MU_P_P_P_P_P";
    case CCDIS_MU_PIPLUS: return "CCDIS_MU_PIPLUS";
    case CCDIS_MU_P_PIPLUS: return "CCDIS_MU_P_PIPLUS";
    case CCDIS_MU_P_P_PIPLUS: return "CCDIS_MU_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_P_P_PIPLUS";
    case CCDIS_MU_PHOTON: return "CCDIS_MU_PHOTON";
    case CCDIS_MU_P_PHOTON: return "CCDIS_MU_P_PHOTON";
    case CCDIS_MU_P_P_PHOTON: return "CCDIS_MU_P_P_PHOTON";
    case CCDIS_MU_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_PHOTON";
    case CCDIS_MU_P_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_P_PHOTON";
    case CCDIS_MU_P_P_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_P_P_PHOTON";
    case CCDIS_MU_PIZERO: return "CCDIS_MU_PIZERO";
    case CCDIS_MU_P_PIZERO: return "CCDIS_MU_P_PIZERO";
    case CCDIS_MU_P_P_PIZERO: return "CCDIS_MU_P_P_PIZERO";
    case CCDIS_MU_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_PIZERO";
    case CCDIS_MU_P_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_P_PIZERO";
    case CCDIS_MU_P_P_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_P_P_PIZERO";
    case NCDIS_P: return "NCDIS_P";
    case NCDIS_P_P: return "NCDIS_P_P";
    case NCDIS_P_P_P: return "NCDIS_P_P_P";
    case NCDIS_P_P_P_P: return "NCDIS_P_P_P_P";
    case NCDIS_P_P_P_P_P: return "NCDIS_P_P_P_P_P";
    case NCDIS_PIPLUS: return "NCDIS_PIPLUS";
    case NCDIS_P_PIPLUS: return "NCDIS_P_PIPLUS";
    case NCDIS_P_P_PIPLUS: return "NCDIS_P_P_PIPLUS";
    case NCDIS_P_P_P_PIPLUS: return "NCDIS_P_P_P_PIPLUS";
    case NCDIS_P_P_P_P_PIPLUS: return "NCDIS_P_P_P_P_PIPLUS";
    case NCDIS_P_P_P_P_P_PIPLUS: return "NCDIS_P_P_P_P_P_PIPLUS";
    case NCDIS_PIMINUS: return "NCDIS_PIMINUS";
    case NCDIS_P_PIMINUS: return "NCDIS_P_PIMINUS";
    case NCDIS_P_P_PIMINUS: return "NCDIS_P_P_PIMINUS";
    case NCDIS_P_P_P_PIMINUS: return "NCDIS_P_P_P_PIMINUS";
    case NCDIS_P_P_P_P_PIMINUS: return "NCDIS_P_P_P_P_PIMINUS";
    case NCDIS_P_P_P_P_P_PIMINUS: return "NCDIS_P_P_P_P_P_PIMINUS";
    case NCDIS_PHOTON: return "NCDIS_PHOTON";
    case NCDIS_P_PHOTON: return "NCDIS_P_PHOTON";
    case NCDIS_P_P_PHOTON: return "NCDIS_P_P_PHOTON";
    case NCDIS_P_P_P_PHOTON: return "NCDIS_P_P_P_PHOTON";
    case NCDIS_P_P_P_P_PHOTON: return "NCDIS_P_P_P_P_PHOTON";
    case NCDIS_P_P_P_P_P_PHOTON: return "NCDIS_P_P_P_P_P_PHOTON";
    case NCDIS_PIZERO: return "NCDIS_PIZERO";
    case NCDIS_P_PIZERO: return "NCDIS_P_PIZERO";
    case NCDIS_P_P_PIZERO: return "NCDIS_P_P_PIZERO";
    case NCDIS_P_P_P_PIZERO: return "NCDIS_P_P_P_PIZERO";
    case NCDIS_P_P_P_P_PIZERO: return "NCDIS_P_P_P_P_PIZERO";
    case NCDIS_P_P_P_P_P_PIZERO: return "NCDIS_P_P_P_P_P_PIZERO";
    case CCCOH: return "CCCOH";
    case NCCOH: return "NCCOH";
    case COSMIC_RAY_MU: return "COSMIC_RAY_MU";
    case COSMIC_RAY_P: return "COSMIC_RAY_P";
    case COSMIC_RAY_E: return "COSMIC_RAY_E";
    case COSMIC_RAY_PHOTON: return "COSMIC_RAY_PHOTON";
    case COSMIC_RAY_OTHER: return "COSMIC_RAY_OTHER";
    case BEAM_PARTICLE_MU: return "BEAM_PARTICLE_MU";
    case BEAM_PARTICLE_P: return "BEAM_PARTICLE_P";
    case BEAM_PARTICLE_E: return "BEAM_PARTICLE_E";
    case BEAM_PARTICLE_PHOTON: return "BEAM_PARTICLE_PHOTON";
    case BEAM_PARTICLE_PI_PLUS: return "BEAM_PARTICLE_PI_PLUS";
    case BEAM_PARTICLE_PI_MINUS: return "BEAM_PARTICLE_PI_MINUS";
    case BEAM_PARTICLE_KAON_PLUS: return "BEAM_PARTICLE_KAON_PLUS";
    case BEAM_PARTICLE_KAON_MINUS: return "BEAM_PARTICLE_KAON_MINUS";
    case BEAM_PARTICLE_OTHER: return "BEAM_PARTICLE_OTHER";
    case OTHER_INTERACTION: return "OTHER_INTERACTION";
    case ALL_INTERACTIONS: return "ALL_INTERACTIONS";
    default: return "NO_RECONSTRUCTABLE";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

InteractionType FromString(std::string &interactionType)
{
    if (interactionType == "NO_RECONSTRUCTABLE")
        return NO_RECONSTRUCTABLE;

    if (interactionType == "COSMIC_RAY")
        return COSMIC_RAY;

    if (interactionType =="CCQEL_MU") return CCQEL_MU;
    else if (interactionType =="CCQEL_MU_P") return CCQEL_MU_P;
    else if (interactionType =="CCQEL_MU_P_P") return CCQEL_MU_P_P;
    else if (interactionType =="CCQEL_MU_P_P_P") return CCQEL_MU_P_P_P;
    else if (interactionType =="CCQEL_MU_P_P_P_P") return CCQEL_MU_P_P_P_P;
    else if (interactionType =="CCQEL_MU_P_P_P_P_P") return CCQEL_MU_P_P_P_P_P;
    else if (interactionType =="CCQEL_E") return CCQEL_E;
    else if (interactionType =="CCQEL_E_P") return CCQEL_E_P;
    else if (interactionType =="CCQEL_E_P_P") return CCQEL_E_P_P;
    else if (interactionType =="CCQEL_E_P_P_P") return CCQEL_E_P_P_P;
    else if (interactionType =="CCQEL_E_P_P_P_P") return CCQEL_E_P_P_P_P;
    else if (interactionType =="CCQEL_E_P_P_P_P_P") return CCQEL_E_P_P_P_P_P;
    else if (interactionType =="NCQEL_P") return NCQEL_P;
    else if (interactionType =="NCQEL_P_P") return NCQEL_P_P;
    else if (interactionType =="NCQEL_P_P_P") return NCQEL_P_P_P;
    else if (interactionType =="NCQEL_P_P_P_P") return NCQEL_P_P_P_P;
    else if (interactionType =="NCQEL_P_P_P_P_P") return NCQEL_P_P_P_P_P;
    else if (interactionType =="CCRES_MU") return CCRES_MU;
    else if (interactionType =="CCRES_MU_P") return CCRES_MU_P;
    else if (interactionType =="CCRES_MU_P_P") return CCRES_MU_P_P;
    else if (interactionType =="CCRES_MU_P_P_P") return CCRES_MU_P_P_P;
    else if (interactionType =="CCRES_MU_P_P_P_P") return CCRES_MU_P_P_P_P;
    else if (interactionType =="CCRES_MU_P_P_P_P_P") return CCRES_MU_P_P_P_P_P;
    else if (interactionType =="CCRES_MU_PIPLUS") return CCRES_MU_PIPLUS;
    else if (interactionType =="CCRES_MU_P_PIPLUS") return CCRES_MU_P_PIPLUS;
    else if (interactionType =="CCRES_MU_P_P_PIPLUS") return CCRES_MU_P_P_PIPLUS;
    else if (interactionType =="CCRES_MU_P_P_P_PIPLUS") return CCRES_MU_P_P_P_PIPLUS;
    else if (interactionType =="CCRES_MU_P_P_P_P_PIPLUS") return CCRES_MU_P_P_P_P_PIPLUS;
    else if (interactionType =="CCRES_MU_P_P_P_P_P_PIPLUS") return CCRES_MU_P_P_P_P_P_PIPLUS;
    else if (interactionType =="CCRES_MU_PHOTON") return CCRES_MU_PHOTON;
    else if (interactionType =="CCRES_MU_P_PHOTON") return CCRES_MU_P_PHOTON;
    else if (interactionType =="CCRES_MU_P_P_PHOTON") return CCRES_MU_P_P_PHOTON;
    else if (interactionType =="CCRES_MU_P_P_P_PHOTON") return CCRES_MU_P_P_P_PHOTON;
    else if (interactionType =="CCRES_MU_P_P_P_P_PHOTON") return CCRES_MU_P_P_P_P_PHOTON;
    else if (interactionType =="CCRES_MU_P_P_P_P_P_PHOTON") return CCRES_MU_P_P_P_P_P_PHOTON;
    else if (interactionType =="CCRES_MU_PIZERO") return CCRES_MU_PIZERO;
    else if (interactionType =="CCRES_MU_P_PIZERO") return CCRES_MU_P_PIZERO;
    else if (interactionType =="CCRES_MU_P_P_PIZERO") return CCRES_MU_P_P_PIZERO;
    else if (interactionType =="CCRES_MU_P_P_P_PIZERO") return CCRES_MU_P_P_P_PIZERO;
    else if (interactionType =="CCRES_MU_P_P_P_P_PIZERO") return CCRES_MU_P_P_P_P_PIZERO;
    else if (interactionType =="CCRES_MU_P_P_P_P_P_PIZERO") return CCRES_MU_P_P_P_P_P_PIZERO;
    else if (interactionType =="CCRES_E") return CCRES_E;
    else if (interactionType =="CCRES_E_P") return CCRES_E_P;
    else if (interactionType =="CCRES_E_P_P") return CCRES_E_P_P;
    else if (interactionType =="CCRES_E_P_P_P") return CCRES_E_P_P_P;
    else if (interactionType =="CCRES_E_P_P_P_P") return CCRES_E_P_P_P_P;
    else if (interactionType =="CCRES_E_P_P_P_P_P") return CCRES_E_P_P_P_P_P;
    else if (interactionType =="CCRES_E_PIPLUS") return CCRES_E_PIPLUS;
    else if (interactionType =="CCRES_E_P_PIPLUS") return CCRES_E_P_PIPLUS;
    else if (interactionType =="CCRES_E_P_P_PIPLUS") return CCRES_E_P_P_PIPLUS;
    else if (interactionType =="CCRES_E_P_P_P_PIPLUS") return CCRES_E_P_P_P_PIPLUS;
    else if (interactionType =="CCRES_E_P_P_P_P_PIPLUS") return CCRES_E_P_P_P_P_PIPLUS;
    else if (interactionType =="CCRES_E_P_P_P_P_P_PIPLUS") return CCRES_E_P_P_P_P_P_PIPLUS;
    else if (interactionType =="CCRES_E_PHOTON") return CCRES_E_PHOTON;
    else if (interactionType =="CCRES_E_P_PHOTON") return CCRES_E_P_PHOTON;
    else if (interactionType =="CCRES_E_P_P_PHOTON") return CCRES_E_P_P_PHOTON;
    else if (interactionType =="CCRES_E_P_P_P_PHOTON") return CCRES_E_P_P_P_PHOTON;
    else if (interactionType =="CCRES_E_P_P_P_P_PHOTON") return CCRES_E_P_P_P_P_PHOTON;
    else if (interactionType =="CCRES_E_P_P_P_P_P_PHOTON") return CCRES_E_P_P_P_P_P_PHOTON;
    else if (interactionType =="CCRES_E_PIZERO") return CCRES_E_PIZERO;
    else if (interactionType =="CCRES_E_P_PIZERO") return CCRES_E_P_PIZERO;
    else if (interactionType =="CCRES_E_P_P_PIZERO") return CCRES_E_P_P_PIZERO;
    else if (interactionType =="CCRES_E_P_P_P_PIZERO") return CCRES_E_P_P_P_PIZERO;
    else if (interactionType =="CCRES_E_P_P_P_P_PIZERO") return CCRES_E_P_P_P_P_PIZERO;
    else if (interactionType =="CCRES_E_P_P_P_P_P_PIZERO") return CCRES_E_P_P_P_P_P_PIZERO;
    else if (interactionType =="NCRES_P") return NCRES_P;
    else if (interactionType =="NCRES_P_P") return NCRES_P_P;
    else if (interactionType =="NCRES_P_P_P") return NCRES_P_P_P;
    else if (interactionType =="NCRES_P_P_P_P") return NCRES_P_P_P_P;
    else if (interactionType =="NCRES_P_P_P_P_P") return NCRES_P_P_P_P_P;
    else if (interactionType =="NCRES_PIPLUS") return NCRES_PIPLUS;
    else if (interactionType =="NCRES_P_PIPLUS") return NCRES_P_PIPLUS;
    else if (interactionType =="NCRES_P_P_PIPLUS") return NCRES_P_P_PIPLUS;
    else if (interactionType =="NCRES_P_P_P_PIPLUS") return NCRES_P_P_P_PIPLUS;
    else if (interactionType =="NCRES_P_P_P_P_PIPLUS") return NCRES_P_P_P_P_PIPLUS;
    else if (interactionType =="NCRES_P_P_P_P_P_PIPLUS") return NCRES_P_P_P_P_P_PIPLUS;
    else if (interactionType =="NCRES_PIMINUS") return NCRES_PIMINUS;
    else if (interactionType =="NCRES_P_PIMINUS") return NCRES_P_PIMINUS;
    else if (interactionType =="NCRES_P_P_PIMINUS") return NCRES_P_P_PIMINUS;
    else if (interactionType =="NCRES_P_P_P_PIMINUS") return NCRES_P_P_P_PIMINUS;
    else if (interactionType =="NCRES_P_P_P_P_PIMINUS") return NCRES_P_P_P_P_PIMINUS;
    else if (interactionType =="NCRES_P_P_P_P_P_PIMINUS") return NCRES_P_P_P_P_P_PIMINUS;
    else if (interactionType =="NCRES_PHOTON") return NCRES_PHOTON;
    else if (interactionType =="NCRES_P_PHOTON") return NCRES_P_PHOTON;
    else if (interactionType =="NCRES_P_P_PHOTON") return NCRES_P_P_PHOTON;
    else if (interactionType =="NCRES_P_P_P_PHOTON") return NCRES_P_P_P_PHOTON;
    else if (interactionType =="NCRES_P_P_P_P_PHOTON") return NCRES_P_P_P_P_PHOTON;
    else if (interactionType =="NCRES_P_P_P_P_P_PHOTON") return NCRES_P_P_P_P_P_PHOTON;
    else if (interactionType =="NCRES_PIZERO") return NCRES_PIZERO;
    else if (interactionType =="NCRES_P_PIZERO") return NCRES_P_PIZERO;
    else if (interactionType =="NCRES_P_P_PIZERO") return NCRES_P_P_PIZERO;
    else if (interactionType =="NCRES_P_P_P_PIZERO") return NCRES_P_P_P_PIZERO;
    else if (interactionType =="NCRES_P_P_P_P_PIZERO") return NCRES_P_P_P_P_PIZERO;
    else if (interactionType =="NCRES_P_P_P_P_P_PIZERO") return NCRES_P_P_P_P_P_PIZERO;
    else if (interactionType =="CCDIS_MU") return CCDIS_MU;
    else if (interactionType =="CCDIS_MU_P") return CCDIS_MU_P;
    else if (interactionType =="CCDIS_MU_P_P") return CCDIS_MU_P_P;
    else if (interactionType =="CCDIS_MU_P_P_P") return CCDIS_MU_P_P_P;
    else if (interactionType =="CCDIS_MU_P_P_P_P") return CCDIS_MU_P_P_P_P;
    else if (interactionType =="CCDIS_MU_P_P_P_P_P") return CCDIS_MU_P_P_P_P_P;
    else if (interactionType =="CCDIS_MU_PIPLUS") return CCDIS_MU_PIPLUS;
    else if (interactionType =="CCDIS_MU_P_PIPLUS") return CCDIS_MU_P_PIPLUS;
    else if (interactionType =="CCDIS_MU_P_P_PIPLUS") return CCDIS_MU_P_P_PIPLUS;
    else if (interactionType =="CCDIS_MU_P_P_P_PIPLUS") return CCDIS_MU_P_P_P_PIPLUS;
    else if (interactionType =="CCDIS_MU_P_P_P_P_PIPLUS") return CCDIS_MU_P_P_P_P_PIPLUS;
    else if (interactionType =="CCDIS_MU_P_P_P_P_P_PIPLUS") return CCDIS_MU_P_P_P_P_P_PIPLUS;
    else if (interactionType =="CCDIS_MU_PHOTON") return CCDIS_MU_PHOTON;
    else if (interactionType =="CCDIS_MU_P_PHOTON") return CCDIS_MU_P_PHOTON;
    else if (interactionType =="CCDIS_MU_P_P_PHOTON") return CCDIS_MU_P_P_PHOTON;
    else if (interactionType =="CCDIS_MU_P_P_P_PHOTON") return CCDIS_MU_P_P_P_PHOTON;
    else if (interactionType =="CCDIS_MU_P_P_P_P_PHOTON") return CCDIS_MU_P_P_P_P_PHOTON;
    else if (interactionType =="CCDIS_MU_P_P_P_P_P_PHOTON") return CCDIS_MU_P_P_P_P_P_PHOTON;
    else if (interactionType =="CCDIS_MU_PIZERO") return CCDIS_MU_PIZERO;
    else if (interactionType =="CCDIS_MU_P_PIZERO") return CCDIS_MU_P_PIZERO;
    else if (interactionType =="CCDIS_MU_P_P_PIZERO") return CCDIS_MU_P_P_PIZERO;
    else if (interactionType =="CCDIS_MU_P_P_P_PIZERO") return CCDIS_MU_P_P_P_PIZERO;
    else if (interactionType =="CCDIS_MU_P_P_P_P_PIZERO") return CCDIS_MU_P_P_P_P_PIZERO;
    else if (interactionType =="CCDIS_MU_P_P_P_P_P_PIZERO") return CCDIS_MU_P_P_P_P_P_PIZERO;
    else if (interactionType =="NCDIS_P") return NCDIS_P;
    else if (interactionType =="NCDIS_P_P") return NCDIS_P_P;
    else if (interactionType =="NCDIS_P_P_P") return NCDIS_P_P_P;
    else if (interactionType =="NCDIS_P_P_P_P") return NCDIS_P_P_P_P;
    else if (interactionType =="NCDIS_P_P_P_P_P") return NCDIS_P_P_P_P_P;
    else if (interactionType =="NCDIS_PIPLUS") return NCDIS_PIPLUS;
    else if (interactionType =="NCDIS_P_PIPLUS") return NCDIS_P_PIPLUS;
    else if (interactionType =="NCDIS_P_P_PIPLUS") return NCDIS_P_P_PIPLUS;
    else if (interactionType =="NCDIS_P_P_P_PIPLUS") return NCDIS_P_P_P_PIPLUS;
    else if (interactionType =="NCDIS_P_P_P_P_PIPLUS") return NCDIS_P_P_P_P_PIPLUS;
    else if (interactionType =="NCDIS_P_P_P_P_P_PIPLUS") return NCDIS_P_P_P_P_P_PIPLUS;
    else if (interactionType =="NCDIS_PIMINUS") return NCDIS_PIMINUS;
    else if (interactionType =="NCDIS_P_PIMINUS") return NCDIS_P_PIMINUS;
    else if (interactionType =="NCDIS_P_P_PIMINUS") return NCDIS_P_P_PIMINUS;
    else if (interactionType =="NCDIS_P_P_P_PIMINUS") return NCDIS_P_P_P_PIMINUS;
    else if (interactionType =="NCDIS_P_P_P_P_PIMINUS") return NCDIS_P_P_P_P_PIMINUS;
    else if (interactionType =="NCDIS_P_P_P_P_P_PIMINUS") return NCDIS_P_P_P_P_P_PIMINUS;
    else if (interactionType =="NCDIS_PHOTON") return NCDIS_PHOTON;
    else if (interactionType =="NCDIS_P_PHOTON") return NCDIS_P_PHOTON;
    else if (interactionType =="NCDIS_P_P_PHOTON") return NCDIS_P_P_PHOTON;
    else if (interactionType =="NCDIS_P_P_P_PHOTON") return NCDIS_P_P_P_PHOTON;
    else if (interactionType =="NCDIS_P_P_P_P_PHOTON") return NCDIS_P_P_P_P_PHOTON;
    else if (interactionType =="NCDIS_P_P_P_P_P_PHOTON") return NCDIS_P_P_P_P_P_PHOTON;
    else if (interactionType =="NCDIS_PIZERO") return NCDIS_PIZERO;
    else if (interactionType =="NCDIS_P_PIZERO") return NCDIS_P_PIZERO;
    else if (interactionType =="NCDIS_P_P_PIZERO") return NCDIS_P_P_PIZERO;
    else if (interactionType =="NCDIS_P_P_P_PIZERO") return NCDIS_P_P_P_PIZERO;
    else if (interactionType =="NCDIS_P_P_P_P_PIZERO") return NCDIS_P_P_P_P_PIZERO;
    else if (interactionType =="NCDIS_P_P_P_P_P_PIZERO") return NCDIS_P_P_P_P_P_PIZERO;
    else if (interactionType =="CCCOH") return CCCOH;
    else if (interactionType =="NCCOH") return NCCOH;
    else if (interactionType =="COSMIC_RAY_MU") return COSMIC_RAY_MU;
    else if (interactionType =="COSMIC_RAY_P") return COSMIC_RAY_P;
    else if (interactionType =="COSMIC_RAY_E") return COSMIC_RAY_E;
    else if (interactionType =="COSMIC_RAY_PHOTON") return COSMIC_RAY_PHOTON;
    else if (interactionType =="COSMIC_RAY_OTHER") return COSMIC_RAY_OTHER;
    else if (interactionType =="BEAM_PARTICLE_MU") return BEAM_PARTICLE_MU;
    else if (interactionType =="BEAM_PARTICLE_P") return BEAM_PARTICLE_P;
    else if (interactionType =="BEAM_PARTICLE_E") return BEAM_PARTICLE_E;
    else if (interactionType =="BEAM_PARTICLE_PHOTON") return BEAM_PARTICLE_PHOTON;
    else if (interactionType =="BEAM_PARTICLE_PI_PLUS") return BEAM_PARTICLE_PI_PLUS;
    else if (interactionType =="BEAM_PARTICLE_PI_MINUS") return BEAM_PARTICLE_PI_MINUS;
    else if (interactionType =="BEAM_PARTICLE_KAON_PLUS") return BEAM_PARTICLE_KAON_PLUS;
    else if (interactionType =="BEAM_PARTICLE_KAON_MINUS") return BEAM_PARTICLE_KAON_MINUS;
    else if (interactionType =="BEAM_PARTICLE_OTHER") return BEAM_PARTICLE_OTHER;
    else if (interactionType =="OTHER_INTERACTION") return OTHER_INTERACTION;
    else if (interactionType =="ALL_INTERACTIONS") return ALL_INTERACTIONS;
    else return NO_RECONSTRUCTABLE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH1F* pHistogram, EColor colour, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    float smallestBinEntry(1e6), largestBinEntry(0.f);
 
    for (int i = 1; i <= pHistogram->GetNbinsX(); i++) 
    {    
        if (pHistogram->GetBinContent(i) < smallestBinEntry) 
            smallestBinEntry = pHistogram->GetBinContent(i);

        if (pHistogram->GetBinContent(i) > largestBinEntry) 
            largestBinEntry = pHistogram->GetBinContent(i);
    }

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram->SetXTitle(xTitle.c_str());
    pHistogram->SetYTitle(yTitle.c_str());
    pHistogram->SetTitle(plotTitle.c_str());
    pHistogram->GetYaxis()->SetTitleOffset(1.3);
    pHistogram->GetYaxis()->SetRangeUser(0.95 * smallestBinEntry, 1.05 * largestBinEntry);
    pHistogram->SetLineColor(colour);
    pHistogram->Draw();
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Comparisons/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH1F* pHistogram1, TH1F* pHistogram2, EColor colour1, EColor colour2, bool normalise, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle, std::string legend1, std::string legend2)
{
    if (normalise)
    {
        pHistogram1->Scale(1/(pHistogram1->GetEntries()));
        pHistogram2->Scale(1/(pHistogram2->GetEntries()));
    }

    float largestBinEntry(0.f);
 
    for (int i = 1; i <= pHistogram1->GetNbinsX(); i++) 
    {    
        if (pHistogram1->GetBinContent(i) > largestBinEntry) 
            largestBinEntry = pHistogram1->GetBinContent(i);

        if (pHistogram2->GetBinContent(i) > largestBinEntry) 
            largestBinEntry = pHistogram2->GetBinContent(i);
    }   

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram1->SetXTitle(xTitle.c_str());
    pHistogram1->SetYTitle(yTitle.c_str());
    pHistogram1->SetTitle(plotTitle.c_str());
    pHistogram1->GetYaxis()->SetTitleOffset(1.3);
    pHistogram1->GetYaxis()->SetRangeUser(0.0, 1.05 * largestBinEntry);

    pHistogram1->SetLineColor(colour1);
    pHistogram1->Draw("HIST");
    pHistogram2->SetLineColor(colour2);

    auto legend = new TLegend(0.55,0.75,0.75,0.85);
    legend->SetHeader("Event Type"); 
    legend->AddEntry(pHistogram1, legend1.c_str(),"l");
    legend->AddEntry(pHistogram2, legend2.c_str(),"l");
    legend->Draw("same");

    pHistogram2->Draw("HISTsame");
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Comparisons/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(std::vector<TH1F*> histogramVector, std::vector<EColor> colourVector, std::vector<std::string> legendVector, bool normalise, bool scale, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    if (normalise)
    {
        for (auto pHistogram : histogramVector)
            pHistogram->Scale(1/(pHistogram->GetEntries()));
    }

    float largestBinEntry(0.f);
 
    for (const auto pHistogram : histogramVector)
    {
        for (int i = 1; i <= pHistogram->GetNbinsX(); i++) 
        {    
            if (pHistogram->GetBinContent(i) > largestBinEntry) 
                largestBinEntry = pHistogram->GetBinContent(i);
        }   
    }

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);

    histogramVector.front()->SetXTitle(xTitle.c_str());
    histogramVector.front()->SetYTitle(yTitle.c_str());
    histogramVector.front()->SetTitle(plotTitle.c_str());
    histogramVector.front()->GetYaxis()->SetTitleOffset(1.3);

    if (scale)
        histogramVector.front()->GetYaxis()->SetRangeUser(0.0, 1.05 * largestBinEntry);

    int colourIndex(0);
    for (const auto pHistogram : histogramVector)
    {
        pHistogram->SetLineColor(colourVector.at(colourIndex));
        ++colourIndex;
    }

    histogramVector.front()->Draw("HIST");

    auto legend = new TLegend(0.45,0.65,0.75,0.85);
    legend->SetHeader("Legend"); 

    int legendIndex(0);
    for (const auto pHistogram : histogramVector)
    {
        legend->AddEntry(pHistogram, legendVector.at(legendIndex).c_str(),"l");
        ++legendIndex;
    }

    legend->Draw("same");

    for (int histogramIndex = 1; histogramIndex < static_cast<int>(histogramVector.size()); ++histogramIndex)
        histogramVector.at(histogramIndex)->Draw("HISTsame");

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Comparisons/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());

    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string interactionTypeToUse, int targetInteractionType)
{
    int interactionType;
    int variable;

    pTree->SetBranchAddress(interactionTypeToUse.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType != targetInteractionType)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName)
{
    float variable;

    pTree->SetBranchAddress(variableName.c_str(), &variable);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);
        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string interactionTypeToUse, int targetInteractionType)
{
    int interactionType;
    float variable;

    pTree->SetBranchAddress(interactionTypeToUse.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType != targetInteractionType)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateAntiHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, int targetInteractionType)
{
    int interactionType;
    int variable;

    pTree->SetBranchAddress("TrueInteractionType", &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType == targetInteractionType)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateFilteredHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string interactionTypeName, int targetInteractionType)
{
    int interactionType;
    int variable;
    int filterVariable;
    int nParticles;
    int nTracks;

    pTree->SetBranchAddress(interactionTypeName.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedTracks", &nTracks);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType != targetInteractionType || filterVariable != targetValue)
            continue;

        if (applyTrackFilter && variable != nTracks)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string interactionTypeName, int targetInteractionType)
{
    float variable;
    int interactionType;
    int filterVariable;
    int nParticles;
    int nTracks;

    pTree->SetBranchAddress(interactionTypeName.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);

    if (filterVariableName != "RecoNeutrinoNumberAssociatedParticles")
        pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &nParticles);

    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedTracks", &nTracks);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType != targetInteractionType || filterVariable != targetValue)
            continue;

        if (applyTrackFilter && filterVariableName != "RecoNeutrinoNumberAssociatedParticles" && nParticles != nTracks)
            continue;

        if (applyTrackFilter && filterVariableName == "RecoNeutrinoNumberAssociatedParticles" && filterVariable != nTracks)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateParticlePDGMassHistogram(TTree* pTree, TH1F* pHistogram, int targetPdg, bool onlyContained)
{
    float longestPfoFitMass;
    int longestPfoPdg;
    int longestPfoIsContained;

    float shortestPfoFitMass;
    int shortestPfoPdg;
    int shortestPfoIsContained;

    pTree->SetBranchAddress("LongestPfoFitMass3D", &longestPfoFitMass);
    pTree->SetBranchAddress("LongestPfoMCPDG3D", &longestPfoPdg);
    pTree->SetBranchAddress("LongestPfoContained3D", &longestPfoIsContained);

    pTree->SetBranchAddress("ShortestPfoFitMass3D", &shortestPfoFitMass);
    pTree->SetBranchAddress("ShortestPfoMCPDG3D", &shortestPfoPdg);
    pTree->SetBranchAddress("ShortestPfoContained3D", &shortestPfoIsContained);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (onlyContained && longestPfoIsContained != 1)
            continue;

        if (longestPfoPdg != targetPdg)
            continue;

        pHistogram->Fill(longestPfoFitMass);
    }

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (shortestPfoFitMass == longestPfoFitMass) //no double counting
            continue;

        if (onlyContained && shortestPfoIsContained != 1)
            continue;

        if (shortestPfoPdg != targetPdg)
            continue;

        pHistogram->Fill(shortestPfoFitMass);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateParticlePDGMinChiHistogram(TTree* pTree, TH1F* pHistogram, int targetPdg, bool onlyContained)
{
    float longestPfoMinChi;
    int longestPfoPdg;
    int longestPfoIsContained;

    float shortestPfoMinChi;
    int shortestPfoPdg;
    int shortestPfoIsContained;

    pTree->SetBranchAddress("LongestPfoMinChiSquaredPerHit3D", &longestPfoMinChi);
    pTree->SetBranchAddress("LongestPfoMCPDG3D", &longestPfoPdg);
    pTree->SetBranchAddress("LongestPfoContained3D", &longestPfoIsContained);

    pTree->SetBranchAddress("ShortestPfoMinChiSquaredPerHit3D", &shortestPfoMinChi);
    pTree->SetBranchAddress("ShortestPfoMCPDG3D", &shortestPfoPdg);
    pTree->SetBranchAddress("ShortestPfoContained3D", &shortestPfoIsContained);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (onlyContained && longestPfoIsContained != 1)
            continue;

        if (longestPfoPdg != targetPdg)
            continue;

        pHistogram->Fill(longestPfoMinChi);
    }

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (shortestPfoMinChi == longestPfoMinChi) //no double counting
            continue;

        if (onlyContained && shortestPfoIsContained != 1)
            continue;

        if (shortestPfoPdg != targetPdg)
            continue;

        pHistogram->Fill(shortestPfoMinChi);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateDoubleFilteredHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string filterVariableName2, int targetValue2, std::string interactionTypeName, int targetInteractionType)
{
    int interactionType;
    int variable;
    int filterVariable;
    int filterVariable2;

    pTree->SetBranchAddress(interactionTypeName.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType != targetInteractionType || filterVariable != targetValue || filterVariable2 != targetValue2)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateDoubleFilteredHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string filterVariableName2, int targetValue2)
{
    int variable;
    int filterVariable;
    int filterVariable2;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);

    //std::map<int, int> variableCount;

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable != targetValue || filterVariable2 != targetValue2)
            continue;

        pHistogram->Fill(variable);
        //variableCount[variable]++;
    }

   /* 
    int totalCount(0);
    for (const auto &entry : variableCount)
    {
        std::cout << "PDG: " << entry.first << " with count: " << entry.second << std::endl;
        totalCount += entry.second;
    }

    std::cout << "Total: " << totalCount << std::endl;
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateDoubleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string filterVariableName2, int targetValue2)
{
    float variable;
    int filterVariable;
    int filterVariable2;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable != targetValue || filterVariable2 != targetValue2)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateTripleFilteredHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3)
{
    int variable;
    int filterVariable;
    int filterVariable2;
    int filterVariable3;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable != targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateTripleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3)
{
    float variable;
    int filterVariable;
    int filterVariable2;
    int filterVariable3;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable != targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateQuadrupleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3, std::string filterVariableName4, int targetValue4)
{
    float variable;
    int filterVariable;
    int filterVariable2;
    int filterVariable3;
    int filterVariable4;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);
    pTree->SetBranchAddress(filterVariableName4.c_str(), &filterVariable4);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable != targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3 || filterVariable4 != targetValue4)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateQuintupleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3, std::string filterVariableName4, int targetValue4, std::string filterVariableName5, int targetValue5)
{
    float variable;
    int filterVariable;
    int filterVariable2;
    int filterVariable3;
    int filterVariable4;
    int filterVariable5;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);
    pTree->SetBranchAddress(filterVariableName4.c_str(), &filterVariable4);
    pTree->SetBranchAddress(filterVariableName5.c_str(), &filterVariable5);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable != targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3 || filterVariable4 != targetValue4 || filterVariable5 != targetValue5)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateLowerBoundQuadrupleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, float targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3, std::string filterVariableName4, int targetValue4)
{
    float variable;
    float filterVariable;
    int filterVariable2;
    int filterVariable3;
    int filterVariable4;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);
    pTree->SetBranchAddress(filterVariableName4.c_str(), &filterVariable4);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable < targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3 || filterVariable4 != targetValue4)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateUpperBoundQuadrupleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, float targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3, std::string filterVariableName4, int targetValue4)
{
    float variable;
    float filterVariable;
    int filterVariable2;
    int filterVariable3;
    int filterVariable4;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);
    pTree->SetBranchAddress(filterVariableName4.c_str(), &filterVariable4);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable > targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3 || filterVariable4 != targetValue4)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateLowerBoundQuintupleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, float targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3, std::string filterVariableName4, int targetValue4, std::string filterVariableName5, int targetValue5)
{
    float variable;
    float filterVariable;
    int filterVariable2;
    int filterVariable3;
    int filterVariable4;
    int filterVariable5;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);
    pTree->SetBranchAddress(filterVariableName4.c_str(), &filterVariable4);
    pTree->SetBranchAddress(filterVariableName5.c_str(), &filterVariable5);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable < targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3 || filterVariable4 != targetValue4 || filterVariable5 != targetValue5)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateLowerBoundSextupleFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, float targetValue, std::string filterVariableName2, int targetValue2, std::string filterVariableName3, int targetValue3, std::string filterVariableName4, int targetValue4, std::string filterVariableName5, int targetValue5, std::string filterVariableName6, int targetValue6)
{
    float variable;
    float filterVariable;
    int filterVariable2;
    int filterVariable3;
    int filterVariable4;
    int filterVariable5;
    int filterVariable6;

    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);
    pTree->SetBranchAddress(filterVariableName2.c_str(), &filterVariable2);
    pTree->SetBranchAddress(filterVariableName3.c_str(), &filterVariable3);
    pTree->SetBranchAddress(filterVariableName4.c_str(), &filterVariable4);
    pTree->SetBranchAddress(filterVariableName5.c_str(), &filterVariable5);
    pTree->SetBranchAddress(filterVariableName6.c_str(), &filterVariable6);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable < targetValue || filterVariable2 != targetValue2 || filterVariable3 != targetValue3 || filterVariable4 != targetValue4 || filterVariable5 != targetValue5 || filterVariable6 != targetValue6)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStack(TTree* pTree, std::string variableName, std::vector<int> &targetInteractionTypes, int targetNumberParticles, int nBins, float minRange, float maxRange, std::string xTitle, std::string yTitle, std::string plotTitle, std::string typeName)
{
    THStack *pStack = new THStack("pStack","");
    std::vector<TH1F*> histogramVector;

    TLegend* pLegend = new TLegend(0.68,0.72,0.98,0.92);

    if (variableName.find("ChiSquared") != std::string::npos)
        minRange = -maxRange;

    for (const auto &interactionType2 : targetInteractionTypes)
    {
        std::string histogramName("pHistogram" + std::to_string(interactionType2));
        TH1F *pHistogram = new TH1F(histogramName.c_str(),"", nBins, minRange, maxRange);
        pHistogram->SetFillStyle(3004);

        if (typeName == "int")
            CreateFilteredHistogram(pTree, pHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetNumberParticles, "ModifiedInteractionType", interactionType2);
        else if (typeName == "float")
            CreateFilteredFloatHistogram(pTree, pHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetNumberParticles, "ModifiedInteractionType", interactionType2);
        else
            std::cout << "Unrecognised typename" << std::endl;

        pLegend->AddEntry(pHistogram, ToString(static_cast<InteractionType>(interactionType2)).c_str());
        histogramVector.push_back(pHistogram);
    }

    for (const auto pHistogram : histogramVector)
        pStack->Add(pHistogram);

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    std::string fullTitle(plotTitle + ";" + xTitle + ";" + yTitle);
    pStack->SetTitle(fullTitle.c_str());
    pStack->Draw("pfc plc");
    pLegend->Draw("same");

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Comparisons/N" + std::to_string(targetNumberParticles) + "_" + variableName + "_Stack.png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;

    for (const auto pHistogram : histogramVector)
       delete pHistogram; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateInteractionTypeTable(TTree* pTree, std::string variableName, int targetValue, std::string interactionTypeToReport, std::string cosmicRayDefinition)
{
    std::map<std::string, int> eventTypeCount;

    int interactionType;
    int variable;
    int nCosmicRays;
    int signal;
    int nParticles;
    //int nTracks;
    int fileIdentifier;
    int eventNumber;
    float bdtResponse;

    pTree->SetBranchAddress(interactionTypeToReport.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(cosmicRayDefinition.c_str(), &nCosmicRays);
    pTree->SetBranchAddress("Signal", &signal);
    //pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedTracks", &nTracks);
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier);
    pTree->SetBranchAddress("EventNumber", &eventNumber);

    static TString invalidBranch("BDT_response");
    TBranch* reponseBranch = static_cast<TBranch*>(pTree->GetListOfBranches()->FindObject(invalidBranch));

    if (reponseBranch) 
    pTree->SetBranchAddress("BDT_response", &bdtResponse);

    //ofstream myfile;
    //myfile.open ("my_events.txt");

    int signalCount(0);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);
    
        if (variable != targetValue)
            continue;

        //if (applyTrackFilter && variable != nTracks)
        //    continue;

        if (reponseBranch && bdtResponse < -0.02)
            continue;

        std::string nCosmicRaysString("");

        if (nCosmicRays == 1)
            nCosmicRaysString = "ONE_COSMIC_RAY";
        if (nCosmicRays == 2)
            nCosmicRaysString = "TWO_COSMIC_RAYS";
        if (nCosmicRays > 2)
            nCosmicRaysString = "MORE_THAN_TWO_COSMIC_RAYS";

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));
    
        //if (interactionType != 165 && isReconstructable == 0)
        //    interactionTypeString = "NO_RECONSTRUCTABLE";
    
        if (interactionType != 165 && nCosmicRaysString != "")
            interactionTypeString += "_" + nCosmicRaysString;

        if (signal == 1)
            ++signalCount;

        eventTypeCount[interactionTypeString]++;
    }

    //myfile.close();

    int totalCount(0), otherCount(0);

    for (const auto &pair : eventTypeCount)
        totalCount += pair.second;

    std::cout << "---------------------------------------------" << std::endl;

    // Declaring the type of Predicate that accepts 2 pairs and return a bool
    typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;
 
    // Defining a lambda function to compare two pairs. It will compare two pairs using second field
    Comparator compFunctor = [](std::pair<std::string, int> elem1 ,std::pair<std::string, int> elem2) { return elem1.second > elem2.second; };
 
    // Declaring a set that will store the pairs using above comparision logic
    std::set<std::pair<std::string, int>, Comparator> interactionTypeSet(eventTypeCount.begin(), eventTypeCount.end(), compFunctor);
 
    // Iterate over a set using range base for loop
    // It will display the items in sorted order of values
    for (std::pair<std::string, int> pair : interactionTypeSet)
    {
        if (100.0 * (float)pair.second/totalCount > 1)
            std::cout << pair.first << " (" << FromString(pair.first) << ") : " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;
        else
            otherCount += pair.second; 
    }

    std::cout << "............................................." << std::endl;
    std::cout << "SIGNAL: " << signalCount << " (" << 100.0 * (float)signalCount/totalCount << "%)" << std::endl;
    std::cout << "OTHER: " << otherCount << " (" << 100.0 * (float)otherCount/totalCount << "%)" << std::endl;
    std::cout << "TOTAL: " << totalCount << " (100%)" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << std::setprecision(6);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateCosmicInteractionTypeTable(TTree* pTree, std::string variableName, int targetValue, std::string interactionTypeToReport)
{
    std::map<std::string, int> eventTypeCount;

    int interactionType;
    int variable;
    int nCosmicRays;

    pTree->SetBranchAddress(interactionTypeToReport.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress("ChosenSliceNumberCosmicRays", &nCosmicRays);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (variable != targetValue)
            continue;

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));
    
        if (nCosmicRays == targetValue)
            eventTypeCount[interactionTypeString]++;
    }

    std::cout << "---------------------------------------------" << std::endl;
    int totalCount(0), otherCount(0);

    for (const auto &pair : eventTypeCount)
        totalCount += pair.second;

    for (const auto &pair : eventTypeCount)
    {
        if (100.0 * (float)pair.second/totalCount > 1)
            std::cout << pair.first << ": " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;
        else
            otherCount += pair.second;
    }

    std::cout << "OTHER: " << otherCount << " (" << 100.0 * (float)otherCount/totalCount << "%)" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateSlicingStack(TTree* pTree, int targetInteractionType)
{
    TLegend* pLegend = new TLegend(0.68,0.72,0.98,0.92);
    
    TH1F *CorrectSlice = new TH1F("CorrectSlice","", 5, 0, 5);
    CorrectSlice->SetFillStyle(3004);
    CreateFilteredHistogram(pTree, CorrectSlice, "RecoNeutrinoNumberAssociatedParticles", "ChosenSliceContainsTrueNeutrino", 1, "TrueInteractionType", targetInteractionType);
    pLegend->AddEntry(CorrectSlice, "Slice Contains #nu_{true}"); 

    TH1F *IncorrectSlice = new TH1F("IncorrectSlice","", 5, 0, 5);
    IncorrectSlice->SetFillStyle(3004);
    CreateFilteredHistogram(pTree, IncorrectSlice, "RecoNeutrinoNumberAssociatedParticles", "ChosenSliceContainsTrueNeutrino", 0, "TrueInteractionType", targetInteractionType);
    pLegend->AddEntry(IncorrectSlice, "Slice Does Not Contain #nu_{true}"); 

    THStack *pStack = new THStack("pStack","");
    pStack->Add(CorrectSlice);
    pStack->Add(IncorrectSlice);

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    std::string fullTitle("Slice does (not) contain #nu_{true} N_{p} and N_{t} for " + ToString(static_cast<InteractionType>(targetInteractionType)) + "; Number of Particles; Number of Events");
    pStack->SetTitle(fullTitle.c_str());
    pStack->Draw("pfc plc");
    pLegend->Draw("same");

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Comparisons/CorrectSlicingStack_" + std::to_string(targetInteractionType) + "_Stack.png");
    canvas->SaveAs(savePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicAnalysis(TTree* pTree, int numberNuRecoParticles)
{
    //Create interaction type table of underlying true neutrino interaction types (these interactions are not selected as nu_reco because there are only CRs in nu_reco)
    /*
    std::cout << "Reconstructed neutrino with " << numberNuRecoParticles << " particle(s) (FULL RECO) COSMIC interaction type table:" << std::endl;
    CreateCosmicInteractionTypeTable(pTree, "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "TrueInteractionType"); 

    //--------------------------------------------------------------

    //Cosmic Direction
    //CosmicDirectionAnalysis(pTree, 1);

    //--------------------------------------------------------------

    //Number of true nu particles in event
    TH1F* pTrueNuNumberParticlesHistogram = new TH1F("pTrueNuNumberParticlesHistogram", "", 5, 0, 5); 
    CreateDoubleFilteredHistogram(pTree, pTrueNuNumberParticlesHistogram, "TrueNeutrinoNumberAssociatedParticles", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", numberNuRecoParticles);
    Draw(pTrueNuNumberParticlesHistogram, kBlue, "pTrueNuNumberParticlesHistogram_CR_N_" + std::to_string(numberNuRecoParticles), "Number of #nu_{true} Particles", "Number of Entries", "Number of #nu_{true} Particles in event (N=" + std::to_string(numberNuRecoParticles) + " COSMIC_RAY(S))");

    //Number of tagging failures
    TH1F* pTaggingFailuresHistogram = new TH1F("pTaggingFailuresHistogram", "", 2, 0, 2); 
    CreateDoubleFilteredHistogram(pTree, pTaggingFailuresHistogram, "TaggingFailure", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", numberNuRecoParticles);
    Draw(pTaggingFailuresHistogram, kBlue, "pTaggingFailuresHistogram_CR_N_" + std::to_string(numberNuRecoParticles), "Tagging Failure (Boolean) ", "Number of Entries", "Whether CR tagging is correct (N=" + std::to_string(numberNuRecoParticles) + " COSMIC_RAY(S))");
    */

    //Opening angle between PFOs (identify split CRs for N=2)
    TH1F* pPhiHistogram = new TH1F("pPhiHistogram", "", 100, 0, 4.0); 
    CreateDoubleFilteredFloatHistogram(pTree, pPhiHistogram, "Phi", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "RecoNeutrinoNumberCosmicRays", numberNuRecoParticles);
    Draw(pPhiHistogram, kBlue, "pPhiHistogram_CR_N_" + std::to_string(numberNuRecoParticles), "Opening Angle #phi Between PFOs (rad) ", "Number of Entries", "Opening Angle #phi Between PFOs (rad) (N=" + std::to_string(numberNuRecoParticles) + " COSMIC_RAY(S))");

    //--------------------------------------------------------------

    /*
    //Number of true nu hits in event, event PFOs, nu reco
    //true nu
    TH1F* pTrueNuNumberHitsHistogram = new TH1F("pTrueNuNumberHitsHistogram", "", 100, 0, 100); 
    CreateDoubleFilteredHistogram(pTree, pTrueNuNumberHitsHistogram, "TrueNeutrinoNumberInducedHits", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", numberNuRecoParticles);

    //PFOs
    TH1F* pTrueNuNumberHitsInPFOsHistogram = new TH1F("pTrueNuNumberHitsInPFOsHistogram", "", 100, 0, 100); 
    CreateDoubleFilteredHistogram(pTree, pTrueNuNumberHitsInPFOsHistogram, "PfoAllTrueNeutrinoHits", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", numberNuRecoParticles);

    //nu reco
    TH1F* pTrueNuNumberHitsInRecoNuHistogram = new TH1F("pTrueNuNumberHitsInRecoNuHistogram", "", 100, 0, 100); 
    CreateDoubleFilteredHistogram(pTree, pTrueNuNumberHitsInRecoNuHistogram, "RecoNuAllTrueNeutrinoHits", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", numberNuRecoParticles);

    //Draw histograms
    std::vector<TH1F*> histogramVector = {pTrueNuNumberHitsHistogram, pTrueNuNumberHitsInPFOsHistogram, pTrueNuNumberHitsInRecoNuHistogram}; 
    std::vector<EColor> colourVector = {kRed, kBlue, kGreen, kOrange};
    std::vector<std::string> legendVector = {"Number of #nu_{true} hits in Event", "Number of #nu_{true} hits in PFOs     ", "Number of #nu_{true} hits in #nu_{reco}"};
    Draw(histogramVector, colourVector, legendVector, false, false, "true_nu_hits_distributions_CR_N_" + std::to_string(numberNuRecoParticles), "Number of #nu_{true} hits", "Number of Entries", "Distributions of number of #nu_{true} hits (N=" + std::to_string(numberNuRecoParticles) + " COSMIC_RAY(S))");
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VariableComparison(TTree* pTree, std::string variableName, int targetInteractionType, int nBins, float minRange, float maxRange, std::string xTitle, std::string yTitle, std::string plotTitle, std::string legend1, std::string legend2)
{
    TH1F* pVariableHistogram = new TH1F("pVariableHistogram", "", nBins, minRange, maxRange); 
    CreateHistogram(pTree, pVariableHistogram, variableName, "TrueInteractionType", targetInteractionType);

    TH1F* pVariableAntiHistogram = new TH1F("pVariableAntiHistogram", "", nBins, minRange, maxRange); 
    CreateAntiHistogram(pTree, pVariableAntiHistogram, variableName, targetInteractionType);

    Draw(pVariableHistogram, pVariableAntiHistogram, kBlue, kRed, true, "VariableComparison_" + std::to_string(targetInteractionType), xTitle.c_str(), yTitle.c_str(), plotTitle.c_str(), legend1.c_str(), legend2.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FindEvents(TTree* pTree, int targetMultiplicity, int targetInteractionType, int targetNumberCosmicRays)
{
    int interactionType;
    int multiplicity;
    int nCosmicRays;
    int fileIdentifier;
    int eventNumber;

    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &multiplicity);
    pTree->SetBranchAddress("TrueInteractionType", &interactionType);
    pTree->SetBranchAddress("ChosenSliceNumberCosmicRays", &nCosmicRays);
    pTree->SetBranchAddress("FileIdentifier", &fileIdentifier);
    pTree->SetBranchAddress("EventNumber", &eventNumber);

    int totalEvents(0);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (multiplicity != targetMultiplicity || interactionType != targetInteractionType || nCosmicRays != targetNumberCosmicRays)
            continue;

        ++totalEvents;
        
        if (eventNumber <= 1)
            std::cout << "Target event (multiplicity = " << multiplicity << ", interactionType = " << interactionType << ", nCosmicRays = "<< nCosmicRays << ") at " << fileIdentifier << ":" << eventNumber << std::endl;
    }

    std::cout << "Total number of events that qualify: " << totalEvents << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateInteractionTypeTables(int targetMultiplicity)
{
    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t1, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", "RecoNeutrinoNumberCosmicRays"); 

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t2, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", "RecoNeutrinoNumberCosmicRays"); 

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t3, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", "RecoNeutrinoNumberCosmicRays"); 

    std::cout << "Reconstructed neutrino with " << targetMultiplicity << " particle(s) (CHEATING: Slicing, Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t4, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", "RecoNeutrinoNumberCosmicRays"); 

    //std::cout << "Post-BDT reconstructed neutrino with " << targetMultiplicity << " particle(s) (FULL RECO) interaction type table:" << std::endl;
    //CreateInteractionTypeTable(t5, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", "RecoNeutrinoNumberCosmicRays"); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NonReconstructableAnalysis(TTree* pTree, int targetMultiplicity)
{
    TH1F* pNoReconstructablePDGHistogram = new TH1F("pNoReconstructablePDGHistogram", "", 3000, 0, 3000); 
    CreateDoubleFilteredHistogram(pTree, pNoReconstructablePDGHistogram, "FirstNeutrinoPfoMCPDG", "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", 166);
    Draw(pNoReconstructablePDGHistogram, kBlue, "NoReconstructablePDGs", "PDG Code", "Number of Events", "PDG codes for N=" + std::to_string(targetMultiplicity) + " fully cheated NO_RECONSTRUCTABLE");

    TH1F* pNoReconstructableParentPDGHistogram = new TH1F("pNoReconstructableParentPDGHistogram", "", 3000, 0, 3000); 
    CreateDoubleFilteredHistogram(pTree, pNoReconstructableParentPDGHistogram, "FirstNeutrinoPfoParentMCPDG", "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", 166);
    Draw(pNoReconstructableParentPDGHistogram, kBlue, "NoReconstructableParentPDGs", "Parent PDG Code", "Number of Events", "Parent PDG codes for N=" + std::to_string(targetMultiplicity) + " fully cheated NO_RECONSTRUCTABLE");

    TH1F* pNoReconstructableNumberHitsHistogram = new TH1F("pNoReconstructableNumberHitsHistogram", "", 3000, 0, 3000); 
    CreateDoubleFilteredHistogram(pTree, pNoReconstructableNumberHitsHistogram, "FirstNeutrinoPfoNumberHits", "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", 166);
    Draw(pNoReconstructableNumberHitsHistogram, kBlue, "NoReconstructableNumberHitss", "NumberHits Code", "Number of Events", "NumberHits codes for N=" + std::to_string(targetMultiplicity) + " fully cheated NO_RECONSTRUCTABLE");

    /* 
    //No_reconstructable category PDG codes
    TH1F* pNoReconstructablePDGHistogram = new TH1F("pNoReconstructablePDGHistogram", "", 3000, 0, 3000); 
    CreateDoubleFilteredHistogram(t6, pNoReconstructablePDGHistogram, "FirstNeutrinoPfoMCParentPDG", "RecoNeutrinoNumberAssociatedParticles", 1, "TrueInteractionType", 165);
    Draw(pNoReconstructablePDGHistogram, kBlue, "NoReconstructablePDGs", "PDG Code", "Number of Events", "PDG codes for N=1 fully cheated NO_RECONSTRUCTABLE");
    FindEvents(1, 165, 0);
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicDirectionAnalysis(TTree* pTree, int targetNumberParticles)
{
    TH1F* pCosmicDeltaChiSquaredHistogram = new TH1F("pCosmicDeltaChiSquaredHistogram", "", 100, -5, 5); 
    TH1F* pSignalDeltaChiSquaredHistogram = new TH1F("pSignalDeltaChiSquaredHistogram", "", 100, -5, 5); 

    if (targetNumberParticles == 1)
    {
        CreateFilteredFloatHistogram(pTree, pCosmicDeltaChiSquaredHistogram, "UpDownDeltaChiSquaredPerHit", "RecoNeutrinoNumberAssociatedParticles", 1, "ModifiedInteractionType", 165);
        CreateFilteredFloatHistogram(pTree, pSignalDeltaChiSquaredHistogram, "UpDownDeltaChiSquaredPerHit", "RecoNeutrinoNumberAssociatedParticles", 1, "ModifiedInteractionType", 0);
    }
    if (targetNumberParticles == 2)
    {
        CreateFilteredFloatHistogram(pTree, pCosmicDeltaChiSquaredHistogram, "UpDownDeltaChiSquaredPerHit", "RecoNeutrinoNumberAssociatedParticles", 2, "ModifiedInteractionType", 165);
        CreateFilteredFloatHistogram(pTree, pSignalDeltaChiSquaredHistogram, "UpDownDeltaChiSquaredPerHit", "RecoNeutrinoNumberAssociatedParticles", 2, "ModifiedInteractionType", 1);
    }

    Draw(pSignalDeltaChiSquaredHistogram, pCosmicDeltaChiSquaredHistogram, kBlue, kRed, true, "UpDownDeltaChiSquaredPerHit_Histogram", "#Delta#chi^{2}_{UD}/N", "Number of Events", "#Delta#chi^{2}_{UD}/N for CR Background and Signal", "Signal", "Cosmics");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStackedPlots(TTree* pTree, std::vector<int> &interactionTypes, int particleMultiplicity)
{
    gStyle->SetPalette(kPastel);

    std::vector<std::string> variableNames = {"TotalEventCharge", "LongestPfoCharge", "ShortestPfoCharge", "LongestPfoLength", "ShortestPfoLength", "Theta", "Phi", "NeutrinoMomentumX", "NeutrinoMomentumY", "NeutrinoMomentumZ", "FitParameterZero", "FitParameterOne", "FitParameterTwo", "FitParameterThree"};

    std::vector<float> upperBounds = {2e5, 1e5, 3e4, 400, 100, 3.14, 3.14, 1.0, 1.0, 1.5, 2e2, 10, 2e2, 1.5};

    for (int i = 0; i < static_cast<int>(variableNames.size()); ++i)
        CreateStack(pTree, variableNames.at(i), interactionTypes, particleMultiplicity, 100, 0, upperBounds.at(i), variableNames.at(i), "Number of Events", "Distribution of variable " + variableNames.at(i) + " for N=" + std::to_string(particleMultiplicity), "float");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateNumberMatchesPlots(TTree* pTree, int interactionType)
{
    TH1F* pTrueNuParticlesHistogram = new TH1F("pTrueNuParticlesHistogram", "", 5, 0, 5); 
    TH1F* pTrueNuTracksHistogram = new TH1F("pTrueNuTracksHistogram", "", 5, 0, 5); 
    CreateHistogram(pTree, pTrueNuParticlesHistogram, "TrueNeutrinoNumberAssociatedParticles", "TrueInteractionType", interactionType);
    CreateHistogram(pTree, pTrueNuTracksHistogram, "TrueNeutrinoNumberAssociatedTracks", "TrueInteractionType", interactionType);

    Draw(pTrueNuParticlesHistogram, pTrueNuTracksHistogram, kBlue, kRed, false, "TrueNuPTHistogram_InteractionType_" + std::to_string(interactionType), "Number of Particles", "Number of Events", "Number of reconstructed particles associated to the true neutrino for event type " + ToString(static_cast<InteractionType>(interactionType)), "#nu_{true} Particles", "#nu_{true} Tracks");

    TH1F* pTrueMuonNumberMatches = new TH1F("pTrueMuonNumberMatches", "", 5, 0, 5); 
    TH1F* pTrueProtonNumberMatches = new TH1F("pTrueProtonNumberMatches", "", 5, 0, 5); 
    CreateHistogram(pTree, pTrueMuonNumberMatches, "TrueMuonMuonmberAssociatedParticles", "TrueInteractionType", interactionType);
    CreateHistogram(pTree, pTrueProtonNumberMatches, "TrueProtonProtonmberAssociatedParticles", "TrueInteractionType", interactionType);

    Draw(pTrueMuonNumberMatches, pTrueProtonNumberMatches, kViolet, kOrange, false, "MuonProtonMatches_InteractionType_" + std::to_string(interactionType), "Number of Particles", "Number of Events", "Number of reconstructed particles associated to the true Muon/Proton for event type " + ToString(static_cast<InteractionType>(interactionType)), "Muon Particle Matches", "Proton Particle Matches");

    TH1F* pRecoNuParticlesHistogram = new TH1F("pRecoNuParticlesHistogram", "", 5, 0, 5); 
    TH1F* pRecoNuTracksHistogram = new TH1F("pRecoNuTracksHistogram", "", 5, 0, 5); 
    CreateHistogram(pTree, pRecoNuParticlesHistogram, "RecoNeutrinoNumberAssociatedParticles", "TrueInteractionType", interactionType);
    CreateHistogram(pTree, pRecoNuTracksHistogram, "RecoNeutrinoNumberAssociatedTracks", "TrueInteractionType", interactionType);

    Draw(pRecoNuParticlesHistogram, pRecoNuTracksHistogram, kBlue, kRed, false, "RecoNuPTHistogram_InteractionType_" + std::to_string(interactionType), "Number of Particles", "Number of Events", "Number of reconstructed particles associated to the reconstructed neutrino for event type " + ToString(static_cast<InteractionType>(interactionType)), "#nu_{reco} Particles", "#nu_{reco} Tracks");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareDistributions(TTree* pTree, int targetMultiplicity, std::string variableName, int firstInteractionType, int secondInteractionType, int nBins, float lowerBound, float upperBound)
{
    TH1F* pFirstInteractionTypeHistogram = new TH1F("pFirstInteractionTypeHistogram", "", nBins, lowerBound, upperBound); 
    TH1F* pSecondInteractionTypeHistogram = new TH1F("pSecondInteractionTypeHistogram", "", nBins, lowerBound, upperBound); 

    CreateFilteredFloatHistogram(pTree, pFirstInteractionTypeHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", firstInteractionType);
    CreateFilteredFloatHistogram(pTree, pSecondInteractionTypeHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", secondInteractionType);

    Draw(pFirstInteractionTypeHistogram, pSecondInteractionTypeHistogram, kBlue, kRed, true, variableName + "_Distribution_Comparison", variableName, "Fraction of Events", variableName + " Normalised Variable Distribution Comparison for Event Types " + ToString(static_cast<InteractionType>(firstInteractionType)) + " and " + ToString(static_cast<InteractionType>(secondInteractionType)), "Event Type " + std::to_string(firstInteractionType), "Event Type " + std::to_string(secondInteractionType));

    delete pFirstInteractionTypeHistogram;
    delete pSecondInteractionTypeHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ComparePDGFitMassDistributions(TTree* pTree, bool onlyContained, int targetPdg1, int targetPdg2, int nBins, float lowerBound, float upperBound)
{
    TH1F* pFirstPDGHistogram = new TH1F("pFirstPDGHistogram", "", nBins, lowerBound, upperBound); 
    TH1F* pSecondPDGHistogram = new TH1F("pSecondPDGHistogram", "", nBins, lowerBound, upperBound); 

    CreateParticlePDGMassHistogram(pTree, pFirstPDGHistogram, targetPdg1, onlyContained);
    CreateParticlePDGMassHistogram(pTree, pSecondPDGHistogram, targetPdg2, onlyContained);

    Draw(pFirstPDGHistogram, pSecondPDGHistogram, kBlue, kRed, true, "FitMass_Distribution_Comparison_" + std::to_string(targetPdg1) + "_" + std::to_string(targetPdg2), "Fit Mass", "Fraction of Events", "Normalised Fit Mass Distributions PDG Codes "  + std::to_string(targetPdg1) + " and " + std::to_string(targetPdg2), "PDG " + std::to_string(targetPdg1), "PDG " + std::to_string(targetPdg2));

    delete pFirstPDGHistogram;
    delete pSecondPDGHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ComparePDGFitMinChiDistributions(TTree* pTree, bool onlyContained, int targetPdg1, int targetPdg2, int nBins, float lowerBound, float upperBound)
{
    TH1F* pFirstPDGHistogram = new TH1F("pFirstPDGHistogram", "", nBins, lowerBound, upperBound); 
    TH1F* pSecondPDGHistogram = new TH1F("pSecondPDGHistogram", "", nBins, lowerBound, upperBound); 

    CreateParticlePDGMinChiHistogram(pTree, pFirstPDGHistogram, targetPdg1, onlyContained);
    CreateParticlePDGMinChiHistogram(pTree, pSecondPDGHistogram, targetPdg2, onlyContained);

    Draw(pFirstPDGHistogram, pSecondPDGHistogram, kBlue, kRed, true, "MinChi_Distribution_Comparison_" + std::to_string(targetPdg1) + "_" + std::to_string(targetPdg2), "Fit #chi^{2}_{min}", "Fraction of Events", "Normalised #chi^{2}_{min}/N Distributions For PDG Codes "  + std::to_string(targetPdg1) + " and " + std::to_string(targetPdg2), "PDG " + std::to_string(targetPdg1), "PDG " + std::to_string(targetPdg2));

    delete pFirstPDGHistogram;
    delete pSecondPDGHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareYFaceIntersectingDistributions(TTree* pTree, int targetMultiplicity, std::string variableName, int firstInteractionType, int secondInteractionType, int nBins, float lowerBound, float upperBound)
{
    TH1F* pFirstInteractionTypeHistogram = new TH1F("pFirstInteractionTypeHistogram", "", nBins, lowerBound, upperBound); 
    TH1F* pSecondInteractionTypeHistogram = new TH1F("pSecondInteractionTypeHistogram", "", nBins, lowerBound, upperBound); 

    CreateTripleFilteredHistogram(pTree, pFirstInteractionTypeHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", firstInteractionType, "LongestPfoIntersectsYFace", 1);
    CreateTripleFilteredHistogram(pTree, pSecondInteractionTypeHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "ModifiedInteractionType", secondInteractionType, "LongestPfoIntersectsYFace", 1);

    Draw(pFirstInteractionTypeHistogram, pSecondInteractionTypeHistogram, kBlue, kRed, true, variableName + "_Distribution_Comparison", variableName, "Fraction of Events", variableName + " Normalised Variable Distribution Comparison for Event Types " + ToString(static_cast<InteractionType>(firstInteractionType)) + " and " + ToString(static_cast<InteractionType>(secondInteractionType)), "Event Type " + std::to_string(firstInteractionType), "Event Type " + std::to_string(secondInteractionType));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CompareDistributions(TTree* pTree, int targetMultiplicity, std::string variableName, int nBins, float lowerBound, float upperBound)
{
    TH1F* pSignalHistogram = new TH1F("pSignalHistogram", "", nBins, lowerBound, upperBound); 
    TH1F* pBackgroundHistogram = new TH1F("pBackgroundHistogram", "", nBins, lowerBound, upperBound); 

    CreateFilteredFloatHistogram(pTree, pSignalHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "Signal", 1);
    CreateFilteredFloatHistogram(pTree, pBackgroundHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "Signal", 0);

    Draw(pSignalHistogram, pBackgroundHistogram, kBlue, kRed, true, variableName + "_Distribution_SB_Comparison", variableName, "Fraction of Events", variableName + " Normalised Variable Distribution Comparison for Signal and Background", "Signal", "Background");

    delete pSignalHistogram;
    delete pBackgroundHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateVariableSBComparisons(TTree* pTree, int targetMultiplicity)
{
    /*
    std::vector<std::string> variableNamesVector = 
    {
     "TotalEventCharge", "OpeningAngle", "NeutrinoMomentumX", "NeutrinoMomentumY", "NeutrinoMomentumZ",
     "LongestPfoCharge", "LongestPfoLength", "LongestPfoTrackProbability", "LongestPfoTheta", "LongestPfoPhi",
     "LongestPfoDaughterHitFraction", "LongestPfoMinX", "LongestPfoMinY", "LongestPfoMinZ", "LongestPfoMaxX",
     "LongestPfoMaxY", "LongestPfoMaxZ", 
     "ShortestPfoCharge", "ShortestPfoLength", "ShortestPfoTrackProbability", "ShortestPfoTheta", "ShortestPfoPhi",
     "ShortestPfoDaughterHitFraction", "ShortestPfoMinX", "ShortestPfoMinY", "ShortestPfoMinZ", "ShortestPfoMaxX",
     "ShortestPfoMaxY", "ShortestPfoMaxZ", 
     "LongestPfoDirectionFitSuccesful", "LongestPfoFitParameterZero", "LongestPfoFitParameterOne", "LongestPfoFitParameterTwo", "LongestPfoUpDownDeltaChiSquaredPerHit",
     "LongestPfoDeltaChiSquaredPerHit", "LongestPfoMinChiSquaredPerHit", "LongestPfoSplitApplied", "LongestPfoFitChargeRangeOverLength", "LongestPfoAngleCorrectedFitChargeRangeOverLength",
     "LongestPfoCosmicProbability", "LongestPfoMinSumQW", "LongestPfoMaxSumQW", "LongestPfoSumQWRatio", "LongestPfoStraightnessRatio",
     "ShortestPfoDirectionFitSuccesful", "ShortestPfoFitParameterZero", "ShortestPfoFitParameterOne", "ShortestPfoFitParameterTwo", "ShortestPfoUpDownDeltaChiSquaredPerHit",
     "ShortestPfoDeltaChiSquaredPerHit", "ShortestPfoMinChiSquaredPerHit", "ShortestPfoSplitApplied", "ShortestPfoFitChargeRangeOverLength", "ShortestPfoAngleCorrectedFitChargeRangeOverLength",
     "ShortestPfoCosmicProbability", "ShortestPfoMinSumQW", "ShortestPfoMaxSumQW", "ShortestPfoSumQWRatio", "ShortestPfoStraightnessRatio"
    };

    std::vector<int> numberBinsVector = 
    {
    100, 100, 100, 100, 100,
    100, 100, 100, 100, 100,
    100, 100, 100, 100, 100,
    100, 100, 
    100, 100, 100, 100, 100,
    100, 100, 100, 100, 100,
    100, 100, 
    2, 100, 100, 100, 100,
    100, 100, 100, 100, 100,
    100, 100, 100, 100, 100,
    2, 100, 100, 100, 100,
    100, 100, 100, 100, 100,
    100, 100, 100, 100, 100
    };

    std::vector<float> lowerBoundsVector =
    {
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, -300, 0, 0,
    -300, 0, 
    0, 0, 0, 0, 0,
    0, 0, -300, 0, 0,
    -300, 0, 
    0, 0, 0, 0, -15,
    -15, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, -15,
    -15, 0, 0, 0, 0,
    0, 0, 0, 0, 0
    };

    std::vector<float> upperBoundsVector = 
    {
    250000, 3.14, 1.1, 1.1, 1.1,
    2000, 600, 1.1, 3.14, 3.14,
    0.2, 300.0, 300.0, 1100.0, 300.0,
    300.0, 1100.0, 
    2000, 100, 1.1, 3.14, 3.14,
    0.2, 300.0, 300.0, 1100.0, 300.0,
    300.0, 1100.0, 
    1, 1100.0, 2.0, 1100.0, 5.0,
    5.0, 5.0, 1.1, 0.2, 0.2,
    1.1, 10.0, 10.0, 10.0, 1.3,
    1, 1100.0, 2.0, 1100.0, 5.0,
    5.0, 5.0, 1.1, 0.2, 0.2,
    1.1, 10.0, 10.0, 10.0, 1.3
    };
    */

    /*
    std::vector<std::string> variableNamesVector = {
    "TotalEventCharge", 
    "OpeningAngle", 
    "NeutrinoMomentumX", 
    "NeutrinoMomentumY", 
    "NeutrinoMomentumZ", 
    "LongestPfoCharge",
    "LongestPfoLength",
    "LongestPfoTrackProbability",
    "LongestPfoPolarAngle",
    "LongestPfoAzimuthalAngle",
    "LongestPfoDaughterHitFraction",
    "LongestPfoMinX",
    "LongestPfoMaxX",
    "LongestPfoMinY",
    "LongestPfoMaxY",
    "LongestPfoMinZ",
    "LongestPfoMaxZ",
    "ShortestPfoCharge",
    "ShortestPfoLength",
    "ShortestPfoTrackProbability",
    "ShortestPfoPolarAngle",
    "ShortestPfoAzimuthalAngle",
    "ShortestPfoDaughterHitFraction",
    "ShortestPfoMinX",
    "ShortestPfoMaxX",
    "ShortestPfoMinY",
    "ShortestPfoMaxY",
    "ShortestPfoMinZ",
    "ShortestPfoMaxZ",
    "LongestPfoUpDownDeltaChiSquaredPerHit",
    "LongestPfoDeltaChiSquaredPerHit",
    "LongestPfoMinChiSquaredPerHit",
    "LongestPfoFitChargeRangeOverLength",
    "LongestPfoCosmicProbability",
    "LongestPfoMinSumQW",
    "LongestPfoMaxSumQW",
    "LongestPfoSumQWRatio",
    "LongestPfoMeanQW",
    "LongestPfoMeanHitCharge",
    "LongestPfoStraightnessRatio",
    "ShortestPfoUpDownDeltaChiSquaredPerHit",
    "ShortestPfoDeltaChiSquaredPerHit",
    "ShortestPfoMinChiSquaredPerHit",
    "ShortestPfoFitChargeRangeOverLength",
    "ShortestPfoCosmicProbability",
    "ShortestPfoMinSumQW",
    "ShortestPfoMaxSumQW",
    "ShortestPfoSumQWRatio",
    "ShortestPfoMeanQW",
    "ShortestPfoMeanHitCharge",
    "ShortestPfoStraightnessRatio",
    };
    std::vector<int> numberBinsVector = {
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    100, 
    };
    std::vector<float> lowerBoundsVector = {
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    -150.0, 
    -150.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    -150.0, 
    -150.0, 
    0.0, 
    0.0, 
    -15.0,
    -15.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.9,
    -15.0,
    -15.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.9,
    };
    std::vector<float> upperBoundsVector = {
    500000.0, 
    3.14, 
    1.0, 
    1.0, 
    1.0, 
    500000.0, 
    500.0, 
    1.4, 
    3.14, 
    3.14, 
    0.1, 
    300.0, 
    300.0, 
    150.0, 
    150.0, 
    1100.0, 
    1100.0, 
    500000.0, 
    500.0, 
    1.4, 
    3.14, 
    3.14, 
    0.1, 
    300.0, 
    300.0, 
    150.0, 
    150.0, 
    1100.0, 
    1100.0, 
    15.0,
    15.0,
    5.0,
    0.1,
    1.4,
    25.0,
    25.0,
    5.0,
    5.0,
    5.0,
    1.1,
    15.0,
    15.0,
    5.0,
    0.1,
    1.4,
    25.0,
    25.0,
    5.0,
    5.0,
    5.0,
    1.1,
    };
    */

    std::vector<std::string> variableNamesVector = {
    "LongestPfoCosmicProbability", "ShortestPfoCosmicProbability"}; 
    std::vector<int> numberBinsVector = {
    100, 
    100}; 
    std::vector<float> lowerBoundsVector = {
    0.5, 
    0.5}; 
    std::vector<float> upperBoundsVector = {
    1.0, 
    1.0}; 

    for (int i = 0; i < static_cast<int>(variableNamesVector.size()); ++i)
        CompareDistributions(pTree, targetMultiplicity, variableNamesVector.at(i), numberBinsVector.at(i), lowerBoundsVector.at(i), upperBoundsVector.at(i));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OtherInteractionAnalysis(TTree* pTree, int targetMultiplicity)
{
    std::cout << "---------------------------------------------" << std::endl;

    int nuanceCode;
    int particleMultiplicity;
    int interactionType;

    pTree->SetBranchAddress("NeutrinoNuanceCode", &nuanceCode);
    pTree->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &particleMultiplicity);
    pTree->SetBranchAddress("ModifiedInteractionType", &interactionType);

    std::map<int, int> nuanceCodeCount;

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (particleMultiplicity != targetMultiplicity || interactionType != 163)
            continue;

        nuanceCodeCount[nuanceCode]++;
    }

    int totalCount(0);

    for (const auto &pair : nuanceCodeCount)
        totalCount += pair.second;

    typedef std::function<bool(std::pair<int, int>, std::pair<int, int>)> Comparator;
    Comparator compFunctor = [](std::pair<int, int> elem1 ,std::pair<int, int> elem2) { return elem1.second > elem2.second; };
    std::set<std::pair<int, int>, Comparator> interactionTypeSet(nuanceCodeCount.begin(), nuanceCodeCount.end(), compFunctor);

    for (std::pair<int, int> pair : interactionTypeSet)
        std::cout << "Nuance " << pair.first << ": " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;

    std::cout << "Total: " << totalCount << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    TH1F* pNumberMuonsHistogram = new TH1F("pNumberMuonsHistogram", "", 5, 0, 5); 
    CreateDoubleFilteredHistogram(t1, pNumberMuonsHistogram, "NumberMuons", "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "NeutrinoNuanceCode", 1000, "ModifiedInteractionType", 163);
    Draw(pNumberMuonsHistogram, kBlue, "NumberMuonsHistogram", "Number of Muons", "Number of Events", "Number of True Muons for N=" + std::to_string(targetMultiplicity) + " OTHER_INTERACTION Nuance Code 1000 Events");

    TH1F* pNumberProtonsHistogram = new TH1F("pNumberProtonsHistogram", "", 5, 0, 5); 
    CreateDoubleFilteredHistogram(t1, pNumberProtonsHistogram, "NumberProtons", "RecoNeutrinoNumberAssociatedParticles", targetMultiplicity, "NeutrinoNuanceCode", 1000, "ModifiedInteractionType", 163);
    Draw(pNumberProtonsHistogram, kBlue, "NumberProtonsHistogram", "Number of Protons", "Number of Events", "Number of True Protons for N=" + std::to_string(targetMultiplicity) + " OTHER_INTERACTION Nuance Code 1000 Events");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicDirectionAnalysis()
{
    //Find events that are neutrino-induced and look like they're going downwards
    float udDeltaChiSquaredPerHit;
    float deltaChiSquaredPerHit;
    float mcDeltaY;
    float mcVertexY;
    float recoDeltaY;
    float recoVertexY;
    int fiducialLowY;
    int recoFiducialLowY;
    int mcCosmicRay;
    int mcPDG;
    int intersectsYFace;
    int recoIntersectsYFace;
    int contained;
    int forwards;
    int downwards;
    float maxY;
    float minY;
    float braggCharge;
    int modifiedInteractionType;
    int eventNumber;
    int fileIdentifier;

    t1->SetBranchAddress("LongestPfoUpDownDeltaChiSquaredPerHit", &udDeltaChiSquaredPerHit);
    t1->SetBranchAddress("LongestPfoDeltaChiSquaredPerHit", &deltaChiSquaredPerHit);
    t1->SetBranchAddress("LongestPfoMCDeltaY", &mcDeltaY);
    t1->SetBranchAddress("LongestPfoMCVertexY", &mcVertexY);
    t1->SetBranchAddress("LongestPfoRecoDeltaY", &recoDeltaY);
    t1->SetBranchAddress("LongestPfoRecoVertexY", &recoVertexY);
    t1->SetBranchAddress("LongestPfoMCFiducialLowY", &fiducialLowY);
    t1->SetBranchAddress("LongestPfoRecoFiducialLowY", &recoFiducialLowY);
    t1->SetBranchAddress("LongestPfoMCCosmicRay", &mcCosmicRay);
    t1->SetBranchAddress("LongestPfoMCPDG", &mcPDG);
    t1->SetBranchAddress("LongestPfoMCIntersectsYFace", &intersectsYFace);
    t1->SetBranchAddress("LongestPfoRecoIntersectsYFace", &recoIntersectsYFace);
    t1->SetBranchAddress("LongestPfoMCContained", &contained);
    t1->SetBranchAddress("LongestPfoMCForwards", &forwards);
    t1->SetBranchAddress("LongestPfoMCDownwards", &downwards);
    t1->SetBranchAddress("LongestPfoMaxY", &maxY);
    t1->SetBranchAddress("LongestPfoMinY", &minY);
    t1->SetBranchAddress("LongestPfoRecoLowestTenCmTotalCharge", &braggCharge);
    t1->SetBranchAddress("ModifiedInteractionType", &modifiedInteractionType);
    t1->SetBranchAddress("FileIdentifier", &fileIdentifier);
    t1->SetBranchAddress("EventNumber", &eventNumber);

    TH1F* pCosmicDeltaChiHistogram = new TH1F("pCosmicDeltaChiHistogram", "", 100, -5.0, 5.0); 
    TH1F* pNeutrinoDeltaChiHistogram = new TH1F("pNeutrinoDeltaChiHistogram", "", 100, -5.0, 5.0); 

    int largeBraggSmallDeltaChiCosmics(0), smallBraggSmallDeltaChiCosmics(0), smallBraggLargeDeltaChiCosmics(0), largeBraggLargeDeltaChiCosmics(0);
    int largeBraggSmallDeltaChiTotal(0), smallBraggSmallDeltaChiTotal(0), smallBraggLargeDeltaChiTotal(0), largeBraggLargeDeltaChiTotal(0);

    for (int i = 0; i < t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if (mcCosmicRay == 1 && mcPDG == 13 && recoFiducialLowY == 1 && braggCharge >= 3400 && maxY >= 100.0 && recoVertexY <= 100.0)
            pCosmicDeltaChiHistogram->Fill(udDeltaChiSquaredPerHit);

        if (mcCosmicRay == 0 && mcPDG == 13 && recoFiducialLowY == 1 && braggCharge <= 3400 && maxY >= 100.0 && recoVertexY <= 100.0)
            pNeutrinoDeltaChiHistogram->Fill(udDeltaChiSquaredPerHit);

        if (mcPDG == 13 && recoFiducialLowY == 1 && braggCharge >= 3400 && udDeltaChiSquaredPerHit <= -1.0  && maxY >= 100.0)
        {
            ++largeBraggSmallDeltaChiTotal;
            if (mcCosmicRay == 1) ++largeBraggSmallDeltaChiCosmics;
        }

        if (mcPDG == 13 && recoFiducialLowY == 1 && braggCharge < 3400 && udDeltaChiSquaredPerHit <= -1.0  && maxY >= 100.0)
        {
            ++smallBraggSmallDeltaChiTotal;
            if (mcCosmicRay == 1) ++smallBraggSmallDeltaChiCosmics;
        }

        if (mcPDG == 13 && recoFiducialLowY == 1 && braggCharge < 3400 && udDeltaChiSquaredPerHit > 1.0  && maxY >= 100.0)
        {
            ++smallBraggLargeDeltaChiTotal;
            if (mcCosmicRay == 1) ++smallBraggLargeDeltaChiCosmics;
        }

        if (mcPDG == 13 && recoFiducialLowY == 1 && braggCharge >= 3400 && udDeltaChiSquaredPerHit > 1.0  && maxY >= 100.0)
        {
            ++largeBraggLargeDeltaChiTotal;
            if (mcCosmicRay == 1) ++largeBraggLargeDeltaChiCosmics;
        }
        
    }

    std::cout << "Fraction cosmics with Bragg charge < 3400 and udDeltaChiSquaredPerHit <= -1.0: " << smallBraggSmallDeltaChiCosmics << " ("  << static_cast<float>(smallBraggSmallDeltaChiCosmics)/smallBraggSmallDeltaChiTotal << ")" << std::endl;
    std::cout << "Fraction cosmics with Bragg charge < 3400 and udDeltaChiSquaredPerHit > 1.0: " << smallBraggLargeDeltaChiCosmics << " ("  << static_cast<float>(smallBraggLargeDeltaChiCosmics)/smallBraggLargeDeltaChiTotal << ")"  << std::endl;
    std::cout << "Fraction cosmics with Bragg charge >= 3400 and udDeltaChiSquaredPerHit <= -1.0: " << largeBraggSmallDeltaChiCosmics << " (" << static_cast<float>(largeBraggSmallDeltaChiCosmics)/largeBraggSmallDeltaChiTotal << ")" << std::endl;
    std::cout << "Fraction cosmics with Bragg charge >= 3400 and udDeltaChiSquaredPerHit > 1.0: " << largeBraggLargeDeltaChiCosmics << "(" << static_cast<float>(largeBraggLargeDeltaChiCosmics)/largeBraggLargeDeltaChiTotal << ")" << std::endl;

    Draw(pCosmicDeltaChiHistogram, pNeutrinoDeltaChiHistogram, kBlue, kRed, true, "du_delta_chi_comparison_histogram_nocheat", "#Delta#chi^{2}_{DU}/N", "Number of Entries", "#Delta#chi^{2}_{DU}/N (No Cheating)", "Cosmic", "Neutrino"); 

    /*
    //Lowest 10cm of track total charge in ADC (Bragg charge)
    TH1F* pMCNeutrinoBraggChargeHistogram = new TH1F("pMCNeutrinoBraggChargeHistogram", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoBraggChargeHistogram, "LongestPfoMCLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    TH1F* pMCCosmicBraggChargeHistogram = new TH1F("pMCCosmicBraggChargeHistogram", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pMCCosmicBraggChargeHistogram, "LongestPfoMCLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoBraggChargeHistogram, pMCCosmicBraggChargeHistogram, kBlue, kRed, true, "MC_BraggCharge_Distribution_Comparison_Neutrino_Cosmic", "Bragg Charge Lowest 10 CM (ADC)", "Fraction of Events", "Normalised Bragg Charge Distributions: Cosmic And Neutrino-Induced Muons (Cheating Using HitToMCMap)", "Neutrino-Induced", "Cosmic");

    //Reco Bragg charge
    TH1F* pRecoNeutrinoBraggChargeHistogram = new TH1F("pRecoNeutrinoBraggChargeHistogram", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pRecoNeutrinoBraggChargeHistogram, "LongestPfoRecoLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    TH1F* pRecoCosmicBraggChargeHistogram = new TH1F("pRecoCosmicBraggChargeHistogram", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pRecoCosmicBraggChargeHistogram, "LongestPfoRecoLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pRecoNeutrinoBraggChargeHistogram, pRecoCosmicBraggChargeHistogram, kBlue, kRed, true, "Reco_BraggCharge_Distribution_Comparison_Neutrino_Cosmic", "Bragg Charge Lowest 10 CM (ADC)", "Fraction of Events", "Normalised Bragg Charge Distributions: Cosmic And Neutrino-Induced Muons (Reconstructed)", "Neutrino-Induced", "Cosmic");
    */

    //Bragg hit distance to NN distributions
    /*
    TH1F* pBraggHitDistanceHistogram = new TH1F("pBraggHitDistanceHistogram", "", 500.0, 0.0, 200.0); 
    TH1F* pNonBraggHitDistanceHistogram = new TH1F("pNonBraggHitDistanceHistogram", "", 500.0, 0.0, 200.0); 
    CreateDoubleFilteredFloatHistogram(t0, pBraggHitDistanceHistogram, "HitNNDistance", "MCPure", 1, "MCFiducialLowY", 1);
    CreateDoubleFilteredFloatHistogram(t0, pNonBraggHitDistanceHistogram, "HitNNDistance", "MCPure", 0, "MCFiducialLowY", 1);
    Draw(pBraggHitDistanceHistogram, pNonBraggHitDistanceHistogram, kBlue, kRed, true, "BraggHit_NN_Distance", "Bragg Normalised Distance to 10 NN", "Fraction of Events", "Normalised Bragg Distance to 10 NN: True Bragg and Non-Bragg Hits", "Bragg", "Non-Bragg");
    */
    
    /*
    //FilteredBragg charge
    TH1F* pFilteredNeutrinoBraggChargeHistogram = new TH1F("pFilteredNeutrinoBraggChargeHistogram", "", 100, 0.0, 10000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pFilteredNeutrinoBraggChargeHistogram, "LongestPfoFilteredLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 0, "LongestPfoRecoFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    TH1F* pFilteredCosmicBraggChargeHistogram = new TH1F("pFilteredCosmicBraggChargeHistogram", "", 100, 0.0, 10000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pFilteredCosmicBraggChargeHistogram, "LongestPfoFilteredLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 1, "LongestPfoRecoFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pFilteredNeutrinoBraggChargeHistogram, pFilteredCosmicBraggChargeHistogram, kBlue, kRed, true, "Filtered_BraggCharge_Distribution_Comparison_Neutrino_Cosmic", "Bragg Charge Lowest 10 CM (ADC)", "Fraction of Events", "Normalised Bragg Charge Distributions: Cosmic And Neutrino-Induced Muons (Filterednstructed)", "Neutrino-Induced", "Cosmic");
    */

    //Characterise true upwards neutrino-induced muons in terms of DeltaY and MC vertex poaition
    TH1F* pMCNeutrinoMCDeltaYHistogram = new TH1F("pMCNeutrinoMCDeltaYHistogram", "", 100, 0.0, 240.0); 
    TH1F* pMCNeutrinoMCVertexYHistogram = new TH1F("pMCNeutrinoMCVertexYHistogram", "", 100, -120.0, 120.0); 

    //CreateTripleFilteredFloatHistogram(t1, pMCNeutrinoMCDeltaYHistogram, "LongestPfoMCDeltaY", "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 
    //CreateTripleFilteredFloatHistogram(t1, pMCNeutrinoMCVertexYHistogram, "LongestPfoMCVertexY", "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoMCDeltaYHistogram, "LongestPfoMCDeltaY", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoMCVertexYHistogram, "LongestPfoMCVertexY", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    Draw(pMCNeutrinoMCDeltaYHistogram, kBlue, "pMCNeutrinoMCDeltaYHistogram", "MC Delta Y (cm)", "Number of Events", "MC Delta Y Distribution");
    Draw(pMCNeutrinoMCVertexYHistogram, kBlue, "pMCNeutrinoMCVertexYHistogram", "MC Vertex Y (cm)", "Number of Events", "MC Vertex Y Distribution");

    /*
    //Delta chi squared distributions for cosmic and neutrino-induced muons
    TH1F* pMCNeutrinoUpDownDeltaChiHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCNeutrinoUpDownDeltaChiHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiHistogram", "", 100, -5.0, 5.0); 

    CreateQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoUpDownDeltaChiHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateQuadrupleFilteredFloatHistogram(t1, pMCCosmicUpDownDeltaChiHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoUpDownDeltaChiHistogram, pMCCosmicUpDownDeltaChiHistogram, kBlue, kRed, true, "MC_UpDownDeltaChi_Distribution_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised #Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");

    //The same but max Y > 100cm
    TH1F* pMCNeutrinoUpDownDeltaChiHighYHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiHighYHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCCosmicUpDownDeltaChiHighYHistogram = new TH1F("pMCCosmicUpDownDeltaChiHighYHistogram", "", 100, -5.0, 5.0); 

    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCNeutrinoUpDownDeltaChiHighYHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCCosmicUpDownDeltaChiHighYHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoUpDownDeltaChiHighYHistogram, pMCCosmicUpDownDeltaChiHighYHistogram, kBlue, kRed, true, "MC_UpDownDeltaChi_Distribution_HighY_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And Neutrino-Induced Muons (Max Y > 100cm)", "Neutrino-Induced", "Cosmic");

    //The same but max Y > 100cm and no direction cheating
    TH1F* pMCNeutrinoUpDownDeltaChiHighYNoCheatHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiHighYNoCheatHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCCosmicUpDownDeltaChiHighYNoCheatHistogram = new TH1F("pMCCosmicUpDownDeltaChiHighYNoCheatHistogram", "", 100, -5.0, 5.0); 

    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoUpDownDeltaChiHighYNoCheatHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13);
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCCosmicUpDownDeltaChiHighYNoCheatHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13); 

    Draw(pMCNeutrinoUpDownDeltaChiHighYNoCheatHistogram, pMCCosmicUpDownDeltaChiHighYNoCheatHistogram, kBlue, kRed, true, "MC_UpDownDeltaChi_Distribution_HighYNoCheat_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And Neutrino-Induced Muons (Max Y > 100cm, No Cheating)", "Neutrino-Induced", "Cosmic");
    */

    /*
    //The same but by modified interaction type max Y > 100cm and no direction cheating
    TH1F* pMCNeutrinoUpDownDeltaChiModHighYNoCheatHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiModHighYNoCheatHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCCosmicUpDownDeltaChiModHighYNoCheatHistogram = new TH1F("pMCCosmicUpDownDeltaChiModHighYNoCheatHistogram", "", 100, -5.0, 5.0); 

    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoUpDownDeltaChiModHighYNoCheatHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "ModifiedInteractionType", 0);
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCCosmicUpDownDeltaChiModHighYNoCheatHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "ModifiedInteractionType", 165); 

    Draw(pMCNeutrinoUpDownDeltaChiModHighYNoCheatHistogram, pMCCosmicUpDownDeltaChiModHighYNoCheatHistogram, kBlue, kRed, true, "MC_UpDownDeltaChiMod_Distribution_HighYNoCheat_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And CCQE1M0P Neutrino-Induced Muons (Max Y > 100cm, No Cheating)", "Neutrino-Induced", "Cosmic");
    */

    /*
    //Cosmic probability
    TH1F* pMCNeutrinoCosmicProbabilityHistogram = new TH1F("pMCNeutrinoCosmicProbabilityHistogram", "", 100, 0.5, 1.0); 
    TH1F* pMCCosmicCosmicProbabilityHistogram = new TH1F("pMCCosmicCosmicProbabilityHistogram", "", 100, 0.5, 1.0); 

    CreateQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoCosmicProbabilityHistogram, "LongestPfoCosmicProbability", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateQuadrupleFilteredFloatHistogram(t1, pMCCosmicCosmicProbabilityHistogram, "LongestPfoCosmicProbability", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoCosmicProbabilityHistogram, pMCCosmicCosmicProbabilityHistogram, kBlue, kRed, true, "MC_CosmicProbability_Distribution_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised #Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");

    //Cosmic probability max Y > 100cm
    TH1F* pMCNeutrinoCosmicProbabilityHighYHistogram = new TH1F("pMCNeutrinoCosmicProbabilityHighYHistogram", "", 100, 0.5, 1.0); 
    TH1F* pMCCosmicCosmicProbabilityHighYHistogram = new TH1F("pMCCosmicCosmicProbabilityHighYHistogram", "", 100, 0.5, 1.0); 

    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCNeutrinoCosmicProbabilityHighYHistogram, "LongestPfoCosmicProbability", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCCosmicCosmicProbabilityHighYHistogram, "LongestPfoCosmicProbability", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoCosmicProbabilityHighYHistogram, pMCCosmicCosmicProbabilityHighYHistogram, kBlue, kRed, true, "MC_CosmicProbability_Distribution_HighY_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And Neutrino-Induced Muons (Max Y > 100cm)", "Neutrino-Induced", "Cosmic");
    */

    //Compare neutrino-induced forwards/backwards muons
    /*
    TH1F* pMCForwardsFBDeltaChiHistogram = new TH1F("pMCForwardsFBDeltaChiHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCBackwardsFBDeltaChiHistogram = new TH1F("pMCBackwardsFBDeltaChiHistogram", "", 100, -5.0, 5.0); 

    CreateQuadrupleFilteredFloatHistogram(t1, pMCForwardsFBDeltaChiHistogram, "LongestPfoDeltaChiSquaredPerHit", "LongestPfoMCCosmicRay", 0, "LongestPfoMCContained", 1, "LongestPfoMCPDG", 13, "LongestPfoMCForwards", 1); 
    CreateQuadrupleFilteredFloatHistogram(t1, pMCBackwardsFBDeltaChiHistogram, "LongestPfoDeltaChiSquaredPerHit", "LongestPfoMCCosmicRay", 0, "LongestPfoMCContained", 1, "LongestPfoMCPDG", 13, "LongestPfoMCForwards", 0); 

    Draw(pMCForwardsFBDeltaChiHistogram, pMCBackwardsFBDeltaChiHistogram, kBlue, kRed, true, "MC_FBDeltaChi_Distribution_Comparison_Neutrino_NewSplit", "#Delta#chi^{2}_{FB}/N", "Fraction of Events", "Normalised #Delta#chi^{2}_{FB}/N Distributions: Neutrino-Induced Muons", "MC Forwards", "MC Backwards");
    */

    /*
    //Hit width and input energy distributions
    TH1F* pHitWidthHistogram = new TH1F("pHitWidthHistogram", "", 100, 0.0, 3.0); 
    CreateFloatHistogram(t0, pHitWidthHistogram, "HitWidth");
    Draw(pHitWidthHistogram, kBlue, "hitwidth_distribution", "Hit Width (cm)", "Number of Events", "Hit Width Distribution");

    TH1F* pHitInputEnergyHistogram = new TH1F("pHitInputEnergyHistogram", "", 100, 0.0, 1000.0); 
    CreateFloatHistogram(t0, pHitInputEnergyHistogram, "HitInputEnergy");
    Draw(pHitInputEnergyHistogram, kBlue, "hitinputenergy_distribution", "Hit InputEnergy (cm)", "Number of Events", "Hit InputEnergy Distribution");
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DetectorCalibration()
{
    TH1F* pDriftCoordinateHistogram = new TH1F("pDriftCoordinateHistogram", "", 50, 0.0, 255.0); 
    TH1F* pYZCoordinateHistogram = new TH1F("pYZCoordinateHistogram", "", 100, 0.0, 1500.0); 

    pDriftCoordinateHistogram->SetBinContent(1, 1.0175); 
    pDriftCoordinateHistogram->SetBinContent(2, 1.025); 
    pDriftCoordinateHistogram->SetBinContent(3, 1.028); 
    pDriftCoordinateHistogram->SetBinContent(4, 1.033); 
    pDriftCoordinateHistogram->SetBinContent(5, 1.036); 
    pDriftCoordinateHistogram->SetBinContent(6, 1.034); 
    pDriftCoordinateHistogram->SetBinContent(7, 1.0375); 
    pDriftCoordinateHistogram->SetBinContent(8, 1.0375); 
    pDriftCoordinateHistogram->SetBinContent(9, 1.038); 
    pDriftCoordinateHistogram->SetBinContent(10, 1.0375); 
    pDriftCoordinateHistogram->SetBinContent(11, 1.03725); 
    pDriftCoordinateHistogram->SetBinContent(12, 1.0375); 
    pDriftCoordinateHistogram->SetBinContent(13, 1.035); 
    pDriftCoordinateHistogram->SetBinContent(14, 1.035); 
    pDriftCoordinateHistogram->SetBinContent(15, 1.0325); 
    pDriftCoordinateHistogram->SetBinContent(16, 1.032); 
    pDriftCoordinateHistogram->SetBinContent(17, 1.03); 
    pDriftCoordinateHistogram->SetBinContent(18, 1.028); 
    pDriftCoordinateHistogram->SetBinContent(19, 1.025); 
    pDriftCoordinateHistogram->SetBinContent(20, 1.024); 
    pDriftCoordinateHistogram->SetBinContent(21, 1.023); 
    pDriftCoordinateHistogram->SetBinContent(22, 1.02); 
    pDriftCoordinateHistogram->SetBinContent(23, 1.0175); 
    pDriftCoordinateHistogram->SetBinContent(24, 1.0125); 
    pDriftCoordinateHistogram->SetBinContent(25, 1.01); 
    pDriftCoordinateHistogram->SetBinContent(26, 1.008); 
    pDriftCoordinateHistogram->SetBinContent(27, 1.006); 
    pDriftCoordinateHistogram->SetBinContent(28, 1.005); 
    pDriftCoordinateHistogram->SetBinContent(29, 1.0); 
    pDriftCoordinateHistogram->SetBinContent(30, 0.998); 
    pDriftCoordinateHistogram->SetBinContent(31, 0.994); 
    pDriftCoordinateHistogram->SetBinContent(32, 0.993); 
    pDriftCoordinateHistogram->SetBinContent(33, 0.988); 
    pDriftCoordinateHistogram->SetBinContent(34, 0.986); 
    pDriftCoordinateHistogram->SetBinContent(35, 0.983); 
    pDriftCoordinateHistogram->SetBinContent(36, 0.98); 
    pDriftCoordinateHistogram->SetBinContent(37, 0.976); 
    pDriftCoordinateHistogram->SetBinContent(38, 0.972); 
    pDriftCoordinateHistogram->SetBinContent(39, 0.968); 
    pDriftCoordinateHistogram->SetBinContent(40, 0.966); 
    pDriftCoordinateHistogram->SetBinContent(41, 0.962); 
    pDriftCoordinateHistogram->SetBinContent(42, 0.9605); 
    pDriftCoordinateHistogram->SetBinContent(43, 0.957); 
    pDriftCoordinateHistogram->SetBinContent(44, 0.954); 
    pDriftCoordinateHistogram->SetBinContent(45, 0.953); 
    pDriftCoordinateHistogram->SetBinContent(46, 0.95); 
    pDriftCoordinateHistogram->SetBinContent(47, 0.949); 
    pDriftCoordinateHistogram->SetBinContent(48, 0.947); 
    pDriftCoordinateHistogram->SetBinContent(49, 0.946); 
    pDriftCoordinateHistogram->SetBinContent(50, 0.948); 

    //TF1 *fit = new TF1("fit", langaufun, 0, 250, 4);
    //fit->SetParameters(2.80976e+02, 2.34936e+01, 1.62744e+03, 6.66806e+01);
    //fit->SetParNames("Width","MP","Area","GSigma");
    TF1 *fit = new TF1("fit", "[0] + [1] * sin([2] * x - [3])", 0, 100);
    fit->SetParameters(0.01, 0.01, 0.01, 0.01);
    pDriftCoordinateHistogram->Fit("fit", "R");

    TF1 *fit2 = new TF1("fit2", "[0] + [1] * x", 100, 240);
    fit2->SetParameters(1.0, 1.0);
    pDriftCoordinateHistogram->Fit("fit2", "R+");

    TF1 *fit3 = new TF1("fit3", "[0] + [1] * (x - [2]) + [3] * sin([4] * x - [5])", 0, 255);
    fit3->SetParameters(0.01, 0.01, 0.01, 0.01, 0.01, 0.01);
    fit3->SetLineColor(kBlue);
    pDriftCoordinateHistogram->Fit("fit3", "R");

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pDriftCoordinateHistogram->SetXTitle("X Coordinate (cm)");
    pDriftCoordinateHistogram->SetYTitle("dQ/dx Correction Factor");
    pDriftCoordinateHistogram->SetTitle("W Plane dQ/dx Correction Factor: Drift Coordinate");
    pDriftCoordinateHistogram->GetYaxis()->SetTitleOffset(1.3);
    pDriftCoordinateHistogram->GetYaxis()->SetRangeUser(0.935, 1.04);
    pDriftCoordinateHistogram->SetLineColor(kBlue);
    pDriftCoordinateHistogram->Draw();
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/Comparisons/pDriftCoordinateCorrectionHistogram.png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;

    //Draw(pDriftCoordinateHistogram, kBlue, "pDriftCoordinateCorrectionHistogram", "X Coordinate (cm)", "dQ/dx Correction Factor", "W Plane dQ/dx Correction Factor: Drift Coordinate");

    delete pDriftCoordinateHistogram;
    delete pYZCoordinateHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BraggChargeTreatment()
{
    //Lowest 10cm of track total charge in ADC (Bragg charge)
    TH1F* pMCNeutrinoBraggChargeHistogram = new TH1F("pMCNeutrinoBraggChargeHistogram", "", 100, 0.0, 21000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoBraggChargeHistogram, "LongestPfoMCLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    TH1F* pMCCosmicBraggChargeHistogram = new TH1F("pMCCosmicBraggChargeHistogram", "", 100, 0.0, 21000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pMCCosmicBraggChargeHistogram, "LongestPfoMCLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoBraggChargeHistogram, pMCCosmicBraggChargeHistogram, kBlue, kRed, true, "MC_BraggCharge_Distribution_Comparison_Neutrino_Cosmic", "Bragg Charge Lowest 10 CM (ADC)", "Fraction of Events", "Normalised Bragg Charge Distributions: Cosmic And Neutrino-Induced Muons (Cheating Using HitToMCMap)", "Neutrino-Induced", "Cosmic");

    //Reco Bragg charge W parent hits only
    TH1F* pRecoNeutrinoBraggChargeHistogramW = new TH1F("pRecoNeutrinoBraggChargeHistogramW", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pRecoNeutrinoBraggChargeHistogramW, "LongestPfoRecoLowestTenCmTotalChargeWView", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    TH1F* pRecoCosmicBraggChargeHistogramW = new TH1F("pRecoCosmicBraggChargeHistogramW", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pRecoCosmicBraggChargeHistogramW, "LongestPfoRecoLowestTenCmTotalChargeWView", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pRecoNeutrinoBraggChargeHistogramW, pRecoCosmicBraggChargeHistogramW, kBlue, kRed, true, "Reco_BraggCharge_Distribution_Comparison_Neutrino_Cosmic_W", "Bragg Charge Lowest 10 CM (ADC)", "Fraction of Events", "Normalised Bragg Charge Distributions: Cosmic And Neutrino-Induced Muons (Reconstructed, W parent hits only)", "Neutrino-Induced", "Cosmic");

    //Reco Bragg charge all 3D hits
    TH1F* pRecoNeutrinoBraggChargeHistogram = new TH1F("pRecoNeutrinoBraggChargeHistogram", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pRecoNeutrinoBraggChargeHistogram, "LongestPfoRecoLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    TH1F* pRecoCosmicBraggChargeHistogram = new TH1F("pRecoCosmicBraggChargeHistogram", "", 100, 0.0, 7000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pRecoCosmicBraggChargeHistogram, "LongestPfoRecoLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pRecoNeutrinoBraggChargeHistogram, pRecoCosmicBraggChargeHistogram, kBlue, kRed, true, "Reco_BraggCharge_Distribution_Comparison_Neutrino_Cosmic", "Bragg Charge Lowest 10 CM (ADC)", "Fraction of Events", "Normalised Bragg Charge Distributions: Cosmic And Neutrino-Induced Muons (Reconstructed)", "Neutrino-Induced", "Cosmic");


    /*
    //FilteredBragg charge
    TH1F* pFilteredNeutrinoBraggChargeHistogram = new TH1F("pFilteredNeutrinoBraggChargeHistogram", "", 100, 0.0, 10000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pFilteredNeutrinoBraggChargeHistogram, "LongestPfoFilteredLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 0, "LongestPfoRecoFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    TH1F* pFilteredCosmicBraggChargeHistogram = new TH1F("pFilteredCosmicBraggChargeHistogram", "", 100, 0.0, 10000.0); 
    CreateQuadrupleFilteredFloatHistogram(t1, pFilteredCosmicBraggChargeHistogram, "LongestPfoFilteredLowestTenCmTotalCharge", "LongestPfoMCCosmicRay", 1, "LongestPfoRecoFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pFilteredNeutrinoBraggChargeHistogram, pFilteredCosmicBraggChargeHistogram, kBlue, kRed, true, "Filtered_BraggCharge_Distribution_Comparison_Neutrino_Cosmic", "Bragg Charge Lowest 10 CM (ADC)", "Fraction of Events", "Normalised Bragg Charge Distributions: Cosmic And Neutrino-Induced Muons (Filterednstructed)", "Neutrino-Induced", "Cosmic");
    */

    /*
    TH1F* pPureNNDistanceHistogram = new TH1F("pPureNNDistanceHistogram", "", 100, 0.0, 100.0); 
    TH1F* pImpureNNDistanceHistogram = new TH1F("pImpureNNDistanceHistogram", "", 100, 0.0, 100.0); 

    CreateTripleFilteredFloatHistogram(t0, pPureNNDistanceHistogram, "HitNNDistance", "MCBraggPeak", 1, "MCFiducialLowY", 1, "MCPure", 1);
    CreateTripleFilteredFloatHistogram(t0, pImpureNNDistanceHistogram, "HitNNDistance", "MCBraggPeak", 1, "MCFiducialLowY", 1, "MCPure", 0);

    std::cout << "Pure overflow: " << pPureNNDistanceHistogram->GetBinContent(101) << std::endl;
    std::cout << "Impure overflow: " << pImpureNNDistanceHistogram->GetBinContent(101) << std::endl;

    Draw(pPureNNDistanceHistogram, pImpureNNDistanceHistogram, kBlue, kRed, true, "Pure_Impure_NN_Distance_Comparison", "Summed NN Distance (5)", "Fraction of Events", "Normalised NN Distance Distributions", "Pure Hits", "Impure Hits");

    delete pPureNNDistanceHistogram;
    delete pImpureNNDistanceHistogram;
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Test()
{
    //std::vector<int> interactionTypes = {165, 1, 163, 12, 0};

    //std::string variableName("FitChargeRangeOverLength");

    //CreateStack(pTestTree, variableName, interactionTypes, 1, 100, 0.0, 0.04, variableName, "Number of Events", "Distribution of variable " + variableName + " for N=1", "float");

    //CompareDistributions(t1, 1, "LongestPfoFitParameterTwo", 0, 12, 100, 0, 1100);

    //TH1F* pNumberMuonsHistogram = new TH1F("pNumberMuonsHistogram", "", 5, 0, 5); 
    //CreateFilteredHistogram(t1, pNumberMuonsHistogram, "NumberMuons", "RecoNeutrinoNumberAssociatedParticles", 1, "ModifiedInteractionType", 0);
    //Draw(pNumberMuonsHistogram, kBlue, "NumberMuonsHistogram", "Number of Muons", "Number of Events", "Number of True Muons for N=1 OTHER_INTERACTION Nuance Code 1000 Events");

    /*
    std::vector<std::string> variableNamesVector = {"LongestPfoMinSumQW", "LongestPfoMaxSumQW", "LongestPfoSumQWRatio", "LongestPfoMeanQW", "LongestPfoFitChargeRangeOverLength"};
    std::vector<int> numberBinsVector = {100, 100, 100, 100, 100};
    std::vector<float> lowerBoundsVector = {0.0, 0.0, 0.5, 1.0, 0.0};
    std::vector<float> upperBoundsVector = {50.0, 50.0, 3.0, 5.0, 10.0};

    for (int i = 0; i < static_cast<int>(variableNamesVector.size()); ++i)
        CompareDistributions(t1, 1, variableNamesVector.at(i), 0, 12, numberBinsVector.at(i), lowerBoundsVector.at(i), upperBoundsVector.at(i));
    */

    /*
    std::vector<std::string> variableNamesVector = {"LongestPfoFitMass3D"};
    std::vector<int> numberBinsVector = {100};
    std::vector<float> lowerBoundsVector = {0.0};
    std::vector<float> upperBoundsVector = {1000.0};

    for (int i = 0; i < static_cast<int>(variableNamesVector.size()); ++i)
        ComparePDGFitMassDistributions(t1, 1, variableNamesVector.at(i), numberBinsVector.at(i), lowerBoundsVector.at(i), upperBoundsVector.at(i));
    */

    //ComparePDGFitMassDistributions(t1, false, 13, 2212, 100, 0.0, 1500.0);
    //ComparePDGFitMassDistributions(t1, false, 211, 2112, 100, 0.0, 1500.0);

    /*
    TH1F* pFirstPDGHistogram = new TH1F("pFirstPDGHistogram", "", 100, 0.0, 1500.0); 
    TH1F* pSecondPDGHistogram = new TH1F("pSecondPDGHistogram", "", 100, 0.0, 1500.0); 

    CreateParticlePDGMassHistogram(t1, pFirstPDGHistogram, 211, false);
    CreateParticlePDGMassHistogram(t1, pSecondPDGHistogram, 2212, false);

    Draw(pFirstPDGHistogram, pSecondPDGHistogram, kBlue, kRed, true, "FitMass_Distribution_Comparison_" + std::to_string(211) + "_" + std::to_string(2212), "Fit Mass", "Fraction of Events", "Normalised Fit Mass Distributions PDG Codes " + std::to_string(211) + " and " + std::to_string(2212), "PDG " + std::to_string(211), "PDG " + std::to_string(2212));

    delete pFirstPDGHistogram;
    delete pSecondPDGHistogram;
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Event_Selection(void)
{
    //CreateInteractionTypeTables(1);

    //Test();

    //CosmicDirectionAnalysis();

    //ComparePDGFitMassDistributions(t1, false, 13, 2212, 100, 0.0, 1500.0);

    //ComparePDGFitMinChiDistributions(t1, false, 13, 2212, 100, 0.0, 80.0);

    //BraggChargeTreatment();

    //DetectorCalibration();

    //CompareYFaceIntersectingDistributions(t1, 1, "LongestPfoCosmicProbability", 0, 165, 100, 0.5, 1.0);

    //CreateVariableSBComparisons(t1, 2);

    //OtherInteractionAnalysis(t1, 2);

    //CosmicAnalysis(t1, 2);

    //std::vector<int> interactionTypesOne = {165, 1, 163, 12, 0};
    //CreateStackedPlots(t1, interactionTypesOne, 1);

    NonReconstructableAnalysis(t1, 1);

    //CreateSlicingStack(t1, 0);
    //CreateSlicingStack(t1, 1);

    //FindEvents(t1, 1, 165, 1);

    //CompareDistributions(t1, 1, "NuRecoContainmentFraction", 0, 165, 100, 0, 1.1);

    //Find weird cosmics that are entirely contained
    /*
    int interactionType;
    int multiplicity;
    float containmentFraction;
    float recoBraggCharge;
    int fileIdentifier;
    int eventNumber;

    t1->SetBranchAddress("RecoNeutrinoNumberAssociatedParticles", &multiplicity);
    t1->SetBranchAddress("ModifiedInteractionType", &interactionType);
    t1->SetBranchAddress("NuRecoContainmentFraction", &containmentFraction);
    t1->SetBranchAddress("LongestPfoRecoLowestTenCmTotalCharge", &recoBraggCharge);
    t1->SetBranchAddress("FileIdentifier", &fileIdentifier);
    t1->SetBranchAddress("EventNumber", &eventNumber);

    for (int i = 0; i < t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if (multiplicity != 1 || interactionType != 0 || recoBraggCharge != 0.0)
            continue;

        if (eventNumber <= 5)
            std::cout << "Fully contained cosmic ray at: " << fileIdentifier << ":" << eventNumber << std::endl;
    }
    */

    //-----------------------------------------------------------

    /*
    TH1F* pMCDownwardsNeutrinoPolarAngleHistogram = new TH1F("pMCDownwardsNeutrinoPolarAngleHistogram", "", 100, 0.0, 3.14); 
    TH1F* pMCUpwardsNeutrinoPolarAngleHistogram = new TH1F("pMCUpwardsNeutrinoPolarAngleHistogram", "", 100, 0.0, 3.14); 
    TH1F* pMCCosmicPolarAngleHistogram = new TH1F("pMCCosmicPolarAngleHistogram", "", 100, 0.0, 3.14); 

    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCDownwardsNeutrinoPolarAngleHistogram, "LongestPfoPolarAngle", "LongestPfoMaxY", -100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCUpwardsNeutrinoPolarAngleHistogram, "LongestPfoPolarAngle", "LongestPfoMaxY", -100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCCosmicPolarAngleHistogram, "LongestPfoPolarAngle", "LongestPfoMaxY", -100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    std::vector<TH1F*> histogramVector = {pMCDownwardsNeutrinoPolarAngleHistogram, pMCUpwardsNeutrinoPolarAngleHistogram, pMCCosmicPolarAngleHistogram};
    std::vector<EColor> colourVector = {kRed, kBlue, kGreen};
    std::vector<std::string> legendVector = {"Downwards Neutrino-Induced", "Upwards Neutrino-Induced", "Cosmic"};

    Draw(histogramVector, colourVector, legendVector, true, true, "pMCCosmicPolarAngleHistogram", "Polar Angle #theta (rad)", "Number of Events", "Cosmic/Neutrino Muon Polar Angle Distribution");
    */

    //-----------------------------------------------------------

    TH1F* pMCNeutrinoUpDownDeltaChiNoCheatHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiNoCheatHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCCosmicUpDownDeltaChiNoCheatHistogram = new TH1F("pMCCosmicUpDownDeltaChiNoCheatHistogram", "", 100, -5.0, 5.0); 

    float downUpChiSquaredPerHit;
    int isCosmicRay;
    int mcPDG;
    int mcDownwards;
    int fiducialLowY;
    float polarAngle;
    float maxY;
    int nuRecoContained;

    t1->SetBranchAddress("LongestPfoUpDownDeltaChiSquaredPerHit", &downUpChiSquaredPerHit);
    t1->SetBranchAddress("LongestPfoMCCosmicRay", &isCosmicRay);
    t1->SetBranchAddress("LongestPfoMCPDG", &mcPDG);
    t1->SetBranchAddress("LongestPfoMCDownwards", &mcDownwards);
    t1->SetBranchAddress("LongestPfoMCFiducialLowY", &fiducialLowY);
    t1->SetBranchAddress("LongestPfoPolarAngle", &polarAngle);
    t1->SetBranchAddress("LongestPfoMaxY", &maxY);
    t1->SetBranchAddress("NuRecoContained", &nuRecoContained);

    for (int i = 0; i < t1->GetEntries(); i++)
    {
        t1->GetEntry(i);

        if (isCosmicRay == 0 && mcPDG == 13 && fiducialLowY == 1 && polarAngle < 3.14/2.0 && nuRecoContained == 1)
            pMCNeutrinoUpDownDeltaChiNoCheatHistogram->Fill(downUpChiSquaredPerHit);

        if (isCosmicRay == 1 && mcPDG == 13 && fiducialLowY == 1 && polarAngle < 3.14/2.0 && nuRecoContained == 1)
            pMCCosmicUpDownDeltaChiNoCheatHistogram->Fill(downUpChiSquaredPerHit);
    }

    Draw(pMCNeutrinoUpDownDeltaChiNoCheatHistogram, pMCCosmicUpDownDeltaChiNoCheatHistogram, kBlue, kRed, true, "MC_UpDownDeltaChi_Distribution_NoCheat_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And Neutrino-Induced Muons (Max Y > 100cm) and #nu_{reco} Contained", "Neutrino-Induced", "Cosmic");

    //-----------------------------------------------------------

    /*
    //Cosmic probability
    TH1F* pMCNeutrinoCosmicProbabilityHistogram = new TH1F("pMCNeutrinoCosmicProbabilityHistogram", "", 100, 0.5, 1.0); 
    TH1F* pMCCosmicCosmicProbabilityHistogram = new TH1F("pMCCosmicCosmicProbabilityHistogram", "", 100, 0.5, 1.0); 

    CreateQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoCosmicProbabilityHistogram, "LongestPfoCosmicProbability", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateQuadrupleFilteredFloatHistogram(t1, pMCCosmicCosmicProbabilityHistogram, "LongestPfoCosmicProbability", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoCosmicProbabilityHistogram, pMCCosmicCosmicProbabilityHistogram, kBlue, kRed, true, "MC_CosmicProbability_Distribution_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised #Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");

    //Cosmic probability max Y > 100cm
    TH1F* pMCNeutrinoCosmicProbabilityHighYHistogram = new TH1F("pMCNeutrinoCosmicProbabilityHighYHistogram", "", 100, 0.5, 1.0); 
    TH1F* pMCCosmicCosmicProbabilityHighYHistogram = new TH1F("pMCCosmicCosmicProbabilityHighYHistogram", "", 100, 0.5, 1.0); 

    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCNeutrinoCosmicProbabilityHighYHistogram, "LongestPfoCosmicProbability", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCCosmicCosmicProbabilityHighYHistogram, "LongestPfoCosmicProbability", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoCosmicProbabilityHighYHistogram, pMCCosmicCosmicProbabilityHighYHistogram, kBlue, kRed, true, "MC_CosmicProbability_Distribution_HighY_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And Neutrino-Induced Muons (Max Y > 100cm)", "Neutrino-Induced", "Cosmic");

    //Cosmic probability max Y > 100cm and contained
    TH1F* pMCNeutrinoCosmicProbabilityHighYContainedHistogram = new TH1F("pMCNeutrinoCosmicProbabilityHighYContainedHistogram", "", 100, 0.5, 1.0); 
    TH1F* pMCCosmicCosmicProbabilityHighYContainedHistogram = new TH1F("pMCCosmicCosmicProbabilityHighYContainedHistogram", "", 100, 0.5, 1.0); 

    CreateLowerBoundSextupleFilteredFloatHistogram(t1, pMCNeutrinoCosmicProbabilityHighYContainedHistogram, "LongestPfoCosmicProbability", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0, "NuRecoContained", 1);
    CreateLowerBoundSextupleFilteredFloatHistogram(t1, pMCCosmicCosmicProbabilityHighYContainedHistogram, "LongestPfoCosmicProbability", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1, "NuRecoContained", 1); 

    Draw(pMCNeutrinoCosmicProbabilityHighYContainedHistogram, pMCCosmicCosmicProbabilityHighYContainedHistogram, kBlue, kRed, true, "MC_CosmicProbability_Distribution_HighYContained_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And Neutrino-Induced Muons (Max Y > 100cm)", "Neutrino-Induced", "Cosmic");
    */

    /*
    //If polar angle > pi/2, look at LongestPfoMCDeltaY and LongestPfoMaxY for cosmics and neutrino-induced muons
    //TH1F* pMCDownwardsNeutrinoMCDeltaYHistogram = new TH1F("pMCDownwardsNeutrinoMCDeltaYHistogram", "", 100, 0.0, 235.0); 
    //TH1F* pMCDownwardsCosmicMCDeltaYHistogram = new TH1F("pMCDownwardsCosmicMCDeltaYHistogram", "", 100, 0.0, 235.0); 

    TH1F* pMCDownwardsNeutrinoMaxYHistogram = new TH1F("pMCDownwardsNeutrinoMaxYHistogram", "", 100, -116.5, 116.5); 
    TH1F* pMCDownwardsCosmicMaxYHistogram = new TH1F("pMCDownwardsCosmicMaxYHistogram", "", 100, -116.5, 116.5); 

    //CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCDownwardsNeutrinoMCDeltaYHistogram, "LongestPfoMCDeltaY", "LongestPfoPolarAngle", 3.1415/2.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    //CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCDownwardsCosmicMCDeltaYHistogram, "LongestPfoMCDeltaY", "LongestPfoPolarAngle", 3.1415/2.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCDownwardsNeutrinoMaxYHistogram, "LongestPfoMaxY", "LongestPfoPolarAngle", 3.1415/2.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCDownwardsCosmicMaxYHistogram, "LongestPfoMaxY", "LongestPfoPolarAngle", 3.1415/2.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    
    //Draw(pMCDownwardsNeutrinoMCDeltaYHistogram, pMCDownwardsCosmicMCDeltaYHistogram, kBlue, kRed, true, "pMCNeutrinoMCDeltaYHistogram", "MC Delta Y (cm)", "Number of Events", "Neutrino-Induced Muon MC Delta Y Distribution", "Cosmic", "Neutrino");
    //Draw(pMCDownwardsNeutrinoMaxYHistogram, pMCDownwardsCosmicMaxYHistogram, kBlue, kRed, true, "pMCNeutrinoMCVertexYHistogram", "Max Y (cm)", "Number of Events", "Neutrino-Induced Muon Max Y Distribution", "Cosmic", "Neutrino");

    //Characterise true upwards neutrino-induced muons in terms of DeltaY and MC vertex position
    //i.e. aim to get rid of downwards neutrino-induced muons
    //TH1F* pMCNeutrinoMCDeltaYHistogram = new TH1F("pMCNeutrinoMCDeltaYHistogram", "", 100, 0.0, 240.0); 
    //TH1F* pMCNeutrinoMCVertexYHistogram = new TH1F("pMCNeutrinoMCVertexYHistogram", "", 100, -120.0, 120.0); 
    TH1F* pMCDownwardsNeutrinoPolarAngleHistogram = new TH1F("pMCDownwardsNeutrinoPolarAngleHistogram", "", 100, 0.0, 3.14); 
    TH1F* pMCUpwardsNeutrinoPolarAngleHistogram = new TH1F("pMCUpwardsNeutrinoPolarAngleHistogram", "", 100, 0.0, 3.14); 

    //CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoMCDeltaYHistogram, "LongestPfoMCDeltaY", "LongestPfoMaxY", 100.0, ""LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    //CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoMCVertexYHistogram, "LongestPfoMCVertexY", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCDownwardsNeutrinoPolarAngleHistogram, "LongestPfoPolarAngle", "LongestPfoMaxY", -100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCUpwardsNeutrinoPolarAngleHistogram, "LongestPfoPolarAngle", "LongestPfoMaxY", -100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 

    //Draw(pMCNeutrinoMCDeltaYHistogram, kBlue, "pMCNeutrinoMCDeltaYHistogram", "MC Delta Y (cm)", "Number of Events", "Neutrino-Induced Muon MC Delta Y Distribution");
    //Draw(pMCNeutrinoMCVertexYHistogram, kBlue, "pMCNeutrinoMCVertexYHistogram", "MC Vertex Y (cm)", "Number of Events", "Neutrino-Induced Muon MC Vertex Y Distribution");

    //TH1F* pMCCosmicMCDeltaYHistogram = new TH1F("pMCCosmicMCDeltaYHistogram", "", 100, 0.0, 240.0); 
    //TH1F* pMCCosmicMCVertexYHistogram = new TH1F("pMCCosmicMCVertexYHistogram", "", 100, -120.0, 120.0); 
    TH1F* pMCCosmicPolarAngleHistogram = new TH1F("pMCCosmicPolarAngleHistogram", "", 100, 0.0, 3.14); 

    //CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCCosmicMCDeltaYHistogram, "LongestPfoMCDeltaY", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    //CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCCosmicMCVertexYHistogram, "LongestPfoMCVertexY", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 
    CreateLowerBoundQuadrupleFilteredFloatHistogram(t1, pMCCosmicPolarAngleHistogram, "LongestPfoPolarAngle", "LongestPfoMaxY", -100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    //Draw(pMCNeutrinoMCDeltaYHistogram, pMCCosmicMCDeltaYHistogram, kBlue, kRed, true, "pMCCosmicMCDeltaYHistogram", "MC Delta Y (cm)", "Number of Events", "Cosmic MC Delta Y Distribution", "Neutrino", "Cosmic");
    //Draw(pMCNeutrinoMCDeltaYHistogram, pMCCosmicMCVertexYHistogram, kBlue, kRed, true, "pMCCosmicMCVertexYHistogram", "MC Vertex Y (cm)", "Number of Events", "Cosmic MC Vertex Y Distribution", "Neutrino", "Cosmic");
    std::vector<TH1F*> histogramVector = {pMCDownwardsNeutrinoPolarAngleHistogram, pMCUpwardsNeutrinoPolarAngleHistogram, pMCCosmicPolarAngleHistogram};
    std::vector<EColor> colourVector = {kRed, kBlue, kGreen};
    std::vector<std::string> legendVector = {"Downwards Neutrino-Induced", "Upwards Neutrino-Induced", "Cosmic"};

    //Draw(histogramVector, colourVector, legendVector, true, true, "pMCCosmicPolarAngleHistogram", "Polar Angle #theta (rad)", "Number of Events", "Cosmic/Neutrino Muon Polar Angle Distribution");
    */

    /*
    //Delta chi squared distributions for cosmic and neutrino-induced muons
    TH1F* pMCNeutrinoUpDownDeltaChiHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCCosmicUpDownDeltaChiHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiHistogram", "", 100, -5.0, 5.0); 

    CreateQuadrupleFilteredFloatHistogram(t1, pMCNeutrinoUpDownDeltaChiHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateQuadrupleFilteredFloatHistogram(t1, pMCCosmicUpDownDeltaChiHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoUpDownDeltaChiHistogram, pMCCosmicUpDownDeltaChiHistogram, kBlue, kRed, true, "MC_UpDownDeltaChi_Distribution_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised #Delta#chi^{2}_{DU}/N Distributions: Cosmic And Neutrino-Induced Muons", "Neutrino-Induced", "Cosmic");

    //The same but max Y > 100cm
    TH1F* pMCNeutrinoUpDownDeltaChiHighYHistogram = new TH1F("pMCNeutrinoUpDownDeltaChiHighYHistogram", "", 100, -5.0, 5.0); 
    TH1F* pMCCosmicUpDownDeltaChiHighYHistogram = new TH1F("pMCCosmicUpDownDeltaChiHighYHistogram", "", 100, -5.0, 5.0); 

    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCNeutrinoUpDownDeltaChiHighYHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 0, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0);
    CreateLowerBoundQuintupleFilteredFloatHistogram(t1, pMCCosmicUpDownDeltaChiHighYHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMaxY", 100.0, "LongestPfoMCCosmicRay", 1, "LongestPfoMCFiducialLowY", 1, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 1); 

    Draw(pMCNeutrinoUpDownDeltaChiHighYHistogram, pMCCosmicUpDownDeltaChiHighYHistogram, kBlue, kRed, true, "MC_UpDownDeltaChi_Distribution_HighY_Comparison_Neutrino_Cosmic", "#Delta#chi^{2}_{DU}/N", "Fraction of Events", "Normalised Cosmic Probability Distributions: Cosmic And Neutrino-Induced Muons (Max Y > 100cm)", "Neutrino-Induced", "Cosmic");
    */
    //TH1F* pMCUpwardsNeutrinoPolarAngleHistogram = new TH1F("pMCUpwardsNeutrinoPolarAngleHistogram", "", 100, 0.0, 3.14); 
    //CreateTripleFilteredFloatHistogram(t1, pMCUpwardsNeutrinoPolarAngleHistogram, "LongestPfoPolarAngle", "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 
    //Draw(pMCUpwardsNeutrinoPolarAngleHistogram, kBlue, "Test", "Polar Angle #phi", "Number of Entries", "Polar Angle #phi Downwards Neutrino-Induced Muons");

    //TH1F* pMCUpwardsNeutrinoDUCHSHistogram = new TH1F("pMCUpwardsNeutrinoDUCHSHistogram", "", 100, -5.0, 5.0); 
    //CreateTripleFilteredFloatHistogram(t1, pMCUpwardsNeutrinoDUCHSHistogram, "LongestPfoUpDownDeltaChiSquaredPerHit", "LongestPfoMCCosmicRay", 0, "LongestPfoMCPDG", 13, "LongestPfoMCDownwards", 0); 
    //Draw(pMCUpwardsNeutrinoDUCHSHistogram, kBlue, "Test2", "#Delta#chi^{2}_{DU}", "Number of Entries", "#Delta#chi^{2}_{DU}/N Dowwnards Neutrino-Induced Muons");


}

//------------------------------------------------------------------------------------------------------------------------------------------
