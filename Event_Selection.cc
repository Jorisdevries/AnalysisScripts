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
    ALL_INTERACTIONS
};

//------------------------------------------------------------------------------------------------------------------------------------------

std::string ToString(const InteractionType interactionType)
{
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
    default: return "EMPTY";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<int, std::string> interactionTypeToStringMap = 
{
{0, "#nu_{#mu} + n #rightarrow #mu^{-}"},
{1, "#nu_{#mu} + n #rightarrow #mu^{-}p"},
{2, "#nu_{#mu} + n #rightarrow #mu^{-}pp"},
{3, "#nu_{#mu} + n #rightarrow #mu^{-}ppp"},
};

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(TH1F* pHistogram, EColor colour, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pHistogram->SetXTitle(xTitle.c_str());
    pHistogram->SetYTitle(yTitle.c_str());
    pHistogram->SetTitle(plotTitle.c_str());
    pHistogram->GetYaxis()->SetTitleOffset(1.3);
    pHistogram->SetLineColor(colour);
    pHistogram->Draw();
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/" + histogramName + ".png");
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
 
    for (int i = 1; i < pHistogram1->GetNbinsX(); i++) 
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
    pHistogram1->Draw();
    pHistogram2->SetLineColor(colour2);

    auto legend = new TLegend(0.55,0.75,0.75,0.85);
    legend->SetHeader("Event Type"); 
    legend->AddEntry(pHistogram1, legend1.c_str(),"l");
    legend->AddEntry(pHistogram2, legend2.c_str(),"l");
    legend->Draw("same");

    pHistogram2->Draw("same");
    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());
    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Draw(std::vector<TH1F*> histogramVector, std::vector<EColor> colourVector, std::vector<std::string> legendVector, bool normalise, std::string histogramName, std::string xTitle, std::string yTitle, std::string plotTitle)
{
    if (normalise)
    {
        for (const auto pHistogram : histogramVector)
            pHistogram->Scale(1/(pHistogram->GetEntries()));
    }

    float largestBinEntry(0.f);
 
    for (const auto pHistogram : histogramVector)
    {
        for (int i = 1; i < pHistogram->GetNbinsX(); i++) 
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
    histogramVector.front()->GetYaxis()->SetRangeUser(0.0, 1.05 * largestBinEntry);

    int colourIndex(0);
    for (const auto pHistogram : histogramVector)
    {
        pHistogram->SetLineColor(colourVector.at(colourIndex));
        ++colourIndex;
    }

    histogramVector.front()->Draw();

    auto legend = new TLegend(0.55,0.75,0.75,0.85);
    legend->SetHeader("Legend"); 

    int legendIndex(0);
    for (const auto pHistogram : histogramVector)
    {
        legend->AddEntry(pHistogram, legendVector.at(legendIndex).c_str(),"l");
        ++legendIndex;
    }

    legend->Draw("same");

    for (int histogramIndex = 1; histogramIndex < static_cast<int>(histogramVector.size()); ++histogramIndex)
        histogramVector.at(histogramIndex)->Draw("same");

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/" + histogramName + ".png");
    canvas->SaveAs(savePath.c_str());

    delete canvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, int targetInteractionType)
{
    int interactionType;
    int variable;

    pTree->SetBranchAddress("TrueInteractionType", &interactionType);
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

    pTree->SetBranchAddress(interactionTypeName.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType != targetInteractionType || filterVariable != targetValue)
            continue;

        pHistogram->Fill(variable);
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateFilteredFloatHistogram(TTree* pTree, TH1F* pHistogram, std::string variableName, std::string filterVariableName, int targetValue, std::string interactionTypeName, int targetInteractionType)
{
    int interactionType;
    float variable;
    int filterVariable;

    pTree->SetBranchAddress(interactionTypeName.c_str(), &interactionType);
    pTree->SetBranchAddress(variableName.c_str(), &variable);
    pTree->SetBranchAddress(filterVariableName.c_str(), &filterVariable);

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (interactionType != targetInteractionType || filterVariable != targetValue)
            continue;

        pHistogram->Fill(variable);
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

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (filterVariable != targetValue || filterVariable2 != targetValue2)
            continue;

        pHistogram->Fill(variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStack(TTree* pTree, std::string variableName, std::vector<int> &targetInteractionTypes, int nBins, float minRange, float maxRange)
{
    THStack *pStack = new THStack("pStack","");
    std::vector<TH1F*> histogramVector;

    //TLegend* pLegend = new TLegend(0.68,0.72,0.98,0.92);

    for (const auto &interactionType : targetInteractionTypes)
    {
        std::string histogramName("pHistogram" + std::to_string(interactionType));
        TH1F *pHistogram = new TH1F(histogramName.c_str(),"", nBins, minRange, maxRange);
        //pHistogram->SetFillStyle(3004);
        CreateHistogram(pTree, pHistogram, variableName, interactionType);

        //pLegend->AddEntry(pHistogram, interactionTypeToStringMap.at(interactionType).c_str());
        //pLegend->AddEntry(pHistogram, ToString(static_cast<InteractionType>(interactionType)).c_str());
        //histogramVector.push_back((TH1F*)pHistogram->Clone());
        pStack->Add((TH1F*)pHistogram->Clone());
    }

    //for (const auto pHistogram : histogramVector)
    //    pStack->Add(pHistogram);

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    pStack->Draw("pfc plc");
    //pLegend->Draw("same");

    //std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/" + variableName + "_Stack.png");
    //canvas->SaveAs(savePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateStack(std::string variableName, std::vector<int> &targetInteractionTypes, int targetNumberParticles, int nBins, float minRange, float maxRange, std::string xTitle, std::string yTitle, std::string plotTitle, std::string typeName)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/full_reco_isolatedhits.root");
    TTree *t1 = (TTree*)f1->Get("EventSelection");

    THStack *pStack = new THStack("pStack","");

    std::vector<TH1F*> histogramVector;

    TLegend* pLegend = new TLegend(0.68,0.72,0.98,0.92);

    for (const auto &interactionType2 : targetInteractionTypes)
    {
        std::string histogramName("pHistogram" + std::to_string(interactionType2));
        TH1F *pHistogram = new TH1F(histogramName.c_str(),"", nBins, minRange, maxRange);
        pHistogram->SetFillStyle(3004);

        if (typeName == "int")
            CreateFilteredHistogram(t1, pHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetNumberParticles, "TrueInteractionType", interactionType2);
        else if (typeName == "float")
            CreateFilteredFloatHistogram(t1, pHistogram, variableName, "RecoNeutrinoNumberAssociatedParticles", targetNumberParticles, "TrueInteractionType", interactionType2);
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

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/" + variableName + "_" + std::to_string(targetNumberParticles) + "_Stack.png");
    canvas->SaveAs(savePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateInteractionTypeTable(TTree* pTree, std::string variableName, int targetValue, std::string interactionTypeToReport)
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

        std::string nCosmicRaysString("");

        if (interactionTypeToReport == "TrueInteractionType" && nCosmicRays == 1)
            nCosmicRaysString = "ONE_COSMIC_RAY";
        if (interactionTypeToReport == "TrueInteractionType" && nCosmicRays == 2)
            nCosmicRaysString = "TWO_COSMIC_RAYS";
        if (interactionTypeToReport == "TrueInteractionType" && nCosmicRays > 2)
            nCosmicRaysString = "MORE_THAN_TWO_COSMIC_RAYS";

        std::string interactionTypeString(ToString(static_cast<InteractionType>(interactionType)));
    
        if (nCosmicRays != targetValue && nCosmicRaysString != "")
            interactionTypeString += "_" + nCosmicRaysString;
        else if (nCosmicRays == targetValue && nCosmicRaysString != "")
            interactionTypeString = nCosmicRaysString;

        eventTypeCount[interactionTypeString]++;
    }

    int totalCount(0), otherCount(0);

    for (const auto &pair : eventTypeCount)
        totalCount += pair.second;

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Total event count: " << totalCount << std::endl;

    std::vector<float> percentagesVector;

    for (const auto &pair : eventTypeCount)
    {
        if (100.0 * (float)pair.second/totalCount > 1)
        {
            std::cout << pair.first << ": " << pair.second << " (" << 100.0 * (float)pair.second/totalCount << "%)" << std::endl;
            percentagesVector.push_back(100.0 * (float)pair.second/totalCount);
        }
        else
            otherCount += pair.second; 
    }

    percentagesVector.push_back((float)otherCount/totalCount );

    std::cout << "OTHER: " << otherCount << " (" << 100.0 * (float)otherCount/totalCount << "%)" << std::endl;

    std::sort(percentagesVector.begin(), percentagesVector.end(), std::greater<float>());
    std::cout << std::setprecision(4);
    std::cout << "ORDER: ";

    for (const auto &entry : percentagesVector)
        std::cout << entry << " ";

    std::cout << std::endl;
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

void CreateSlicingStack(int targetInteractionType)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/particle_counting.root");
    TTree *t1 = (TTree*)f1->Get("EventDescription");

    TLegend* pLegend = new TLegend(0.68,0.72,0.98,0.92);
    
    TH1F *CorrectSlice = new TH1F("CorrectSlice","", 5, 0, 5);
    CorrectSlice->SetFillStyle(3004);
    CreateFilteredHistogram(t1, CorrectSlice, "RecoNeutrinoNumberAssociatedParticles", "ChosenSliceContainsTrueNeutrino", 1, "TrueInteractionType", targetInteractionType);
    pLegend->AddEntry(CorrectSlice, "Slice Contains #nu_{true}"); 

    TH1F *IncorrectSlice = new TH1F("IncorrectSlice","", 5, 0, 5);
    IncorrectSlice->SetFillStyle(3004);
    CreateFilteredHistogram(t1, IncorrectSlice, "RecoNeutrinoNumberAssociatedParticles", "ChosenSliceContainsTrueNeutrino", 0, "TrueInteractionType", targetInteractionType);
    pLegend->AddEntry(IncorrectSlice, "Slice Does Not Contain #nu_{true}"); 

    THStack *pStack = new THStack("pStack","");
    pStack->Add(CorrectSlice);
    pStack->Add(IncorrectSlice);

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 900, 600);
    std::string fullTitle("Slice does (not) contain #nu_{true} N_{p} and N_{t} for " + ToString(static_cast<InteractionType>(targetInteractionType)) + "; Number of Particles; Number of Events");
    pStack->SetTitle(fullTitle.c_str());
    pStack->Draw("pfc plc");
    pLegend->Draw("same");

    std::string savePath("/usera/jjd49/pandora_direction/Scripts/Figures/Event_Selection/CorrectSlicingStack_" + std::to_string(targetInteractionType) + "_Stack.png");
    canvas->SaveAs(savePath.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicAnalysis(int numberNuRecoParticles)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/cheating_II/full_reco.root");
    TTree *t1 = (TTree*)f1->Get("Variables");

    TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/cheating_II/cheating_id.root");
    TTree *t2 = (TTree*)f2->Get("Variables");

    TTree* chosenTree(t1);

    std::cout << "Reconstructed neutrino with " << numberNuRecoParticles << " particle(s) (FULL RECO) COSMIC interaction type table:" << std::endl;
    CreateCosmicInteractionTypeTable(chosenTree, "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "TrueInteractionType"); 

    TH1F* pTrueNuNumberParticlesHistogram = new TH1F("pTrueNuNumberParticlesHistogram", "", 5, 0, 5); 
    CreateDoubleFilteredHistogram(chosenTree, pTrueNuNumberParticlesHistogram, "TrueNeutrinoNumberAssociatedParticles", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "TrueInteractionType", -1);
    Draw(pTrueNuNumberParticlesHistogram, kBlue, "pTrueNuNumberParticlesHistogram_CR_N_" + std::to_string(numberNuRecoParticles), "Number of #nu_{true} Particles", "Number of Entries", "Number of #nu_{true} Particles in event (N=" + std::to_string(numberNuRecoParticles) + " COSMIC_RAY(S))");

    TH1F* pTrueNuNumberHitsHistogram = new TH1F("pTrueNuNumberHitsHistogram", "", 105, -10, 200); 
    CreateDoubleFilteredHistogram(chosenTree, pTrueNuNumberHitsHistogram, "TrueNeutrinoNumberInducedHits", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", 1);

    TH1F* pTrueNuNumberHitsInPFOsHistogram = new TH1F("pTrueNuNumberHitsInPFOsHistogram", "", 105, -10, 200); 
    CreateDoubleFilteredHistogram(chosenTree, pTrueNuNumberHitsInPFOsHistogram, "TrueNeutrinoNumberInducedHitsInPfos", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", 1);

    TH1F* pTrueNuNumberHitsInRecoNuHistogram = new TH1F("pTrueNuNumberHitsInRecoNuHistogram", "", 105, -10, 200); 
    CreateDoubleFilteredHistogram(chosenTree, pTrueNuNumberHitsInRecoNuHistogram, "TrueNeutrinoNumberInducedHitsInRecoNeutrino", "RecoNeutrinoNumberAssociatedParticles", numberNuRecoParticles, "ChosenSliceNumberCosmicRays", 1);

    std::vector<TH1F*> histogramVector = {pTrueNuNumberHitsHistogram, pTrueNuNumberHitsInPFOsHistogram}; 
    std::vector<EColor> colourVector = {kRed, kBlue};
    std::vector<std::string> legendVector = {"Number of #nu_{true} hits in Event", "Number of #nu_{true} hits in PFOs     "};
    Draw(histogramVector, colourVector, legendVector, false, "true_nu_hits_distributions_CR_N_" + std::to_string(numberNuRecoParticles), "Number of #nu_{true} hits", "Number of Entries", "Distributions of number of #nu_{true} hits (N=" + std::to_string(numberNuRecoParticles) + " COSMIC_RAY(S))");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VariableComparison(std::string variableName, int targetInteractionType, int nBins, float minRange, float maxRange, std::string xTitle, std::string yTitle, std::string plotTitle, std::string legend1, std::string legend2)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/full_reco_isolatedhits.root");
    TTree *t1 = (TTree*)f1->Get("Variables");

    TH1F* pVariableHistogram = new TH1F("pVariableHistogram", "", nBins, minRange, maxRange); 
    CreateHistogram(t1, pVariableHistogram, variableName, targetInteractionType);

    TH1F* pVariableAntiHistogram = new TH1F("pVariableAntiHistogram", "", nBins, minRange, maxRange); 
    CreateAntiHistogram(t1, pVariableAntiHistogram, variableName, targetInteractionType);

    Draw(pVariableHistogram, pVariableAntiHistogram, kBlue, kRed, true, "VariableComparison_" + std::to_string(targetInteractionType), xTitle.c_str(), yTitle.c_str(), plotTitle.c_str(), legend1.c_str(), legend2.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FindEvents(int targetMultiplicity, int targetInteractionType, int targetNumberCosmicRays)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/full_reco_isolatedhits.root");
    TTree *t1 = (TTree*)f1->Get("EventSelection");

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

    for (int i = 0; i < pTree->GetEntries(); i++)
    {
        pTree->GetEntry(i);

        if (multiplicity != targetMultiplicity || interactionType != targetInteractionType || nCosmicRays != targetNumberCosmicRays)
            continue;
        
        if (eventNumber <= 5)
            std::cout << "Target event " << fileIdentifier << ":" << eventNumber << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Event_Selection(void)
{
    TFile *f1 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/particle_counting.root");
    TTree *t1 = (TTree*)f1->Get("EventDescription");

    TFile *f2 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/validation_with_counting.root");
    TTree *t2 = (TTree*)f2->Get("EventSelection");

    TFile *f3 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/cheating/Cheating_Slicing_Tagging_ID.root");
    TTree *t3 = (TTree*)f3->Get("EventDescription");

    TFile *f4 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/cheating/Cheating_Tagging_ID.root");
    TTree *t4 = (TTree*)f4->Get("EventDescription");

    TFile *f5 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/cheating/Cheating_ID.root");
    TTree *t5 = (TTree*)f5->Get("EventDescription");

    TFile *f6 = TFile::Open("/usera/jjd49/pandora_direction/CondorUtilities/saved_results/event_selection/with_variables/full_reco.root");
    TTree *t6 = (TTree*)f6->Get("Variables");

    //std::cout << "Reconstructed neutrino with 1 particle interaction type table:" << std::endl;
    //CreateInteractionTypeTable(t1, "RecoNeutrinoNumberAssociatedParticles", 1, "TrueInteractionType"); 

    //std::cout << "Reconstructed neutrino with 2 particles interaction type table:" << std::endl;
    //CreateInteractionTypeTable(t1, "RecoNeutrinoNumberAssociatedParticles", (FULL RE 2, "TrueInteractionType"); 

    /*

    std::cout << "Reconstructed neutrino with 2 particle (FULL RECO) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t1, "RecoNeutrinoNumberAssociatedParticles", 2, "TrueInteractionType"); 

    std::cout << "Reconstructed neutrino with 2 particle (CHEATING: Slicing, Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t3, "RecoNeutrinoNumberAssociatedParticles", 2, "TrueInteractionType"); 

    std::cout << "Reconstructed neutrino with 2 particle (CHEATING: Tagging, Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t4, "RecoNeutrinoNumberAssociatedParticles", 2, "TrueInteractionType"); 

    std::cout << "Reconstructed neutrino with 2 particle (CHEATING: Neutrino ID) interaction type table:" << std::endl;
    CreateInteractionTypeTable(t5, "RecoNeutrinoNumberAssociatedParticles", 2, "TrueInteractionType"); 
    */

    /*
    int interactionType(1);

    TH1F* pTrueNuParticlesHistogram = new TH1F("pTrueNuParticlesHistogram", "", 5, 0, 5); 
    TH1F* pTrueNuTracksHistogram = new TH1F("pTrueNuTracksHistogram", "", 5, 0, 5); 
    CreateHistogram(t2, pTrueNuParticlesHistogram, "TrueNeutrinoNumberAssociatedParticles", interactionType);
    CreateHistogram(t2, pTrueNuTracksHistogram, "TrueNeutrinoNumberAssociatedTracks", interactionType);

    //Draw(pTrueNuParticlesHistogram, pTrueNuTracksHistogram, kBlue, kRed, false, "TrueNuPTHistogram_InteractionType_" + std::to_string(interactionType), "Number of Particles", "Number of Events", "Number of reconstructed particles associated to the true neutrino for event type " + ToString(static_cast<InteractionType>(interactionType)));

    TH1F* pTrueMuonNumberMatches = new TH1F("pTrueMuonNumberMatches", "", 5, 0, 5); 
    TH1F* pTrueProtonNumberMatches = new TH1F("pTrueProtonNumberMatches", "", 5, 0, 5); 
    CreateHistogram(t2, pTrueMuonNumberMatches, "TrueMuonMuonmberAssociatedParticles", interactionType);
    CreateHistogram(t2, pTrueProtonNumberMatches, "TrueProtonProtonmberAssociatedParticles", interactionType);

    //Draw(pTrueMuonNumberMatches, pTrueProtonNumberMatches, kViolet, kOrange, false, "MuonProtonMatches_InteractionType_" + std::to_string(interactionType), "Number of Particles", "Number of Events", "Number of reconstructed particles associated to the true Muon/Proton for event type " + ToString(static_cast<InteractionType>(interactionType)), "Muon Particle Matches", "Proton Particle Matches");

    TH1F* pRecoNuParticlesHistogram = new TH1F("pRecoNuParticlesHistogram", "", 5, 0, 5); 
    TH1F* pRecoNuTracksHistogram = new TH1F("pRecoNuTracksHistogram", "", 5, 0, 5); 
    CreateHistogram(t1, pRecoNuParticlesHistogram, "RecoNeutrinoNumberAssociatedParticles", interactionType);
    CreateHistogram(t1, pRecoNuTracksHistogram, "RecoNeutrinoNumberAssociatedTracks", interactionType);

    //Draw(pRecoNuParticlesHistogram, pRecoNuTracksHistogram, kBlue, kRed, false, "RecoNuPTHistogram_InteractionType_" + std::to_string(interactionType), "Number of Particles", "Number of Events", "Number of reconstructed particles associated to the reconstructed neutrino for event type " + ToString(static_cast<InteractionType>(interactionType)));

    gStyle->SetPalette(kPastel);

    std::vector<int> oneParticleInteractionTypes = {-1, 1, 163, 12, 0};
    //CreateStack("ChosenSliceNumberHits", oneParticleInteractionTypes, 1, 100, 0, 2000, "Number of #nu_{reco} Hits", "Number of Events", "Number of #nu_{reco} Hits for #nu_{reco} with 1 associated particle", "int");

    std::vector<int> twoParticleInteractionTypes = {-1, 0, 163, 100, 101, 88, 18, 1};
    //CreateStack("ChosenSliceNumberHits", twoParticleInteractionTypes, 2, 100, 0, 2000, "Number of #nu_{reco} Hits", "Number of Events", "Number of #nu_{reco} Hits for #nu_{reco} with 2 associated particles", "int");

    //CreateSlicingStack(0);
    //CreateSlicingStack(1);
    */


    FindEvents(1, -1, 1);

    //Cosmic breakdown
    //CosmicAnalysis(1);

    /*
    //Drawing variable stacked plots
    gStyle->SetPalette(kPastel);

    std::vector<int> oneParticleInteractionTypes = {-1, 1, 163, 12, 0};
    std::vector<std::string> variableNames = {"LongestPfoCharge"};
    int particleMultiplicity(1);

    for (const auto variableName : variableNames)
        CreateStack(variableName, oneParticleInteractionTypes, particleMultiplicity, 100, 0, 200000, variableName, "Number of Events", variableName + " for N=1", "float");
    */

    //--------------------------------------------------------------------------------------------------------------------------------------
}
