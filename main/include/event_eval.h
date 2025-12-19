#ifndef _test_def
#define _test_def

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TSystem.h"
#include "TString.h"
#include "TBranch.h"
#include "TBenchmark.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include <TChain.h>
#include "TMinuit.h"
#include <fstream>
#include <iostream>
#include <TCanvas.h>
#include <TColor.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TProfile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector3.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "math.h"
#include "TVector.h"
#include "TH2.h"
#include "TMatrixD.h"
#include <numeric>
#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;

using namespace std;

//functions
void chain_input(const char * filename);
void get_binning();
void reset();
void initiate();
Int_t get_spin();
Int_t get_init();
Int_t get_initlep();
Int_t get_initp();
Int_t get_mclepscat();
Int_t get_lepscat();
void get_kine(TLorentzVector* lscat, Double_t &Q2, Double_t &y, Double_t &xB, Double_t &W2);
void get_dis(void);
Int_t check_dis(void);
Int_t get_parton(void);
void fill_dis();
void get_hadrons(void);
void get_mchadrons(void);
Double_t get_fuu(Double_t x, Double_t q2, Double_t z);
Double_t get_deltaq(Double_t x, Double_t q2, Double_t z);
void fill_had(Double_t w, Double_t z);
void fill_mchad(Double_t w, Double_t z);

void write_file(const char * outdir);

#ifdef MULTI_GLOBAL
TChain* ch;
TTreeReader* tr;
Int_t PID;//input chosen PID to analyze
Int_t PARTON;
Int_t SPIN;
Int_t NPos;
Int_t NNeg;
Int_t Ntot;
//const PDF* pdf;
const Int_t NUMy=12;
const Int_t NUMq2=13;
const Int_t NUMxb=25;
const Int_t NUMw=10;
const Int_t NUMz=13;
const Int_t NUMpt=14;
const Double_t cross_angle=-0.0250;
const Double_t b1=0.0481;
const Double_t b2=0.6114;
const Double_t b3=-0.3509;
const Double_t b4=-0.4611;
const Double_t b5=0.7172;
const Double_t b6=-0.0317;
Double32_t y_bins[NUMy+1];
Double32_t q2_bins[NUMq2+1];
Double32_t xb_bins[NUMxb+1];
Double32_t w_bins[NUMw+1];
Double32_t z_bins[NUMz+1];
Double32_t pt_bins[NUMpt+1];

Double_t rQ2;
Double_t ry;
Double_t rxB;
Double_t rW2;

Double_t mcQ2;
Double_t mcy;
Double_t mcxB;
Double_t mcW2;

Double_t depol;

//void find_mcjpsi(void);
//Int_t find_recelec(void);
//void get_dis(void);
//void loop_hadrons(void);

//TLorentzVector* l;//incoming beam lepton
//TLorentzVector* p;//incoming beam proton
//TLorentzVector* lscatmc;//scattered beam lepton
TLorentzVector* lscatrec;//reconstructed beam lepton
//Double_t const mass_p=0.93827;

TLorentzVector* l;
TLorentzVector* p;
TLorentzVector* lscatmc;
TLorentzVector* lscatr;

Int_t mcid_escat;

//PDFs and FFs
//Number of DIS events
Double32_t nbr_disU[NUMq2][NUMxb];
//Number of hadrons
Double32_t nbr_hadU[NUMq2][NUMxb][NUMz];
Double32_t nbr_hadU_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_xbU[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2U[NUMq2][NUMxb][NUMz];
Double32_t nbr_yU[NUMq2][NUMxb][NUMz];
Double32_t nbr_depolU[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2_globU[NUMq2][NUMz];
Double32_t nbr_zU[NUMq2][NUMxb][NUMz];

//Number of DIS events
Double32_t nbr_disD[NUMq2][NUMxb];
//Number of hadrons
Double32_t nbr_hadD[NUMq2][NUMxb][NUMz];
Double32_t nbr_hadD_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_xbD[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2D[NUMq2][NUMxb][NUMz];
Double32_t nbr_yD[NUMq2][NUMxb][NUMz];
Double32_t nbr_depolD[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2_globD[NUMq2][NUMz];
Double32_t nbr_zD[NUMq2][NUMxb][NUMz];

//Number of DIS events
Double32_t nbr_mcdisU[NUMq2][NUMxb];
//Number of hadrons
Double32_t nbr_mchadU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mchadU_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcxbU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2U[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcyU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcdepolU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2_globU[NUMq2][NUMz];
Double32_t nbr_mczU[NUMq2][NUMxb][NUMz];

//Number of DIS events
Double32_t nbr_mcdisD[NUMq2][NUMxb];
//Number of hadrons
Double32_t nbr_mchadD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mchadD_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcxbD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2D[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcyD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcdepolD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2_globD[NUMq2][NUMz];
Double32_t nbr_mczD[NUMq2][NUMxb][NUMz];

Double_t Lumi;

//particle masses [GeV]
const double mass_e= 0.00051;
const double mass_mu=0.10566;
const double mass_pi= 0.13957;
const double mass_k= 0.49368;
const double mass_p=0.93827;
const double particle_mass[5] = {mass_e,mass_mu,mass_pi,mass_k,mass_p};

#else
extern TChain* ch;
extern TTreeReader* tr;
//extern event_tree* ev(NULL);
extern Int_t PID;//input chosen PID to analyze;
extern Int_t PARTON;
extern Int_t SPIN;
extern Int_t NPos;
extern Int_t NNeg;
extern Int_t Ntot;
//extern const PDF* pdf;
extern const Int_t NUMy=12;
extern const Int_t NUMq2=13;
extern const Int_t NUMxb=25;
extern const Int_t NUMw=10;
extern const Int_t NUMz=13;
extern const Int_t NUMpt=14;
extern const Double_t cross_angle=0.0250;
extern const Double_t b1=0.0481;
extern const Double_t b2=0.6114;
extern const Double_t b3=-0.3509;
extern const Double_t b4=-0.4611;
extern const Double_t b5=0.7172;
extern const Double_t b6=-0.0317;
extern Double32_t y_bins[NUMy+1];
extern Double32_t q2_bins[NUMq2+1];
extern Double32_t xb_bins[NUMxb+1];
extern Double32_t w_bins[NUMw+1];
extern Double32_t z_bins[NUMz+1];
extern Double32_t pt_bins[NUMpt+1];

extern Double_t Q2;
extern Double_t y;
extern Double_t xB;
extern Double_t W2;

extern Double_t mcQ2;
extern Double_t mcy;
extern Double_t mcxB;
extern Double_t mcW2;

extern Double_t depol;

extern Double_t hepQ2;
extern Double_t hepy;
extern Double_t hepxB;
extern Double_t hepW2;

extern TLorentzVector* l;
extern TLorentzVector* p;
extern TLorentzVector* lscatmc;
extern TLorentzVector* lscat;
extern TLorentzVector* lscatrec;//reconstructed beam lepton
extern Double_t rQ2;
extern Double_t ry;
extern Double_t rxB;
extern Double_t rW2;

extern Int_t mcid_escat;

//Number of DIS events
extern Double32_t nbr_disU[NUMq2][NUMxb];
//Number of hadrons
extern Double32_t nbr_hadU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_hadU_err[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_xbU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_q2U[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_yU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_depolU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_q2_globU[NUMq2][NUMz];
extern Double32_t nbr_zU[NUMq2][NUMxb][NUMz];

//Number of DIS events
extern Double32_t nbr_disD[NUMq2][NUMxb];
//Number of hadrons
extern Double32_t nbr_hadD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_hadD_err[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_xbD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_q2D[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_yD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_depolD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_q2_globD[NUMq2][NUMz];
extern Double32_t nbr_zD[NUMq2][NUMxb][NUMz];

//Number of DIS events
extern Double32_t nbr_mcdisU[NUMq2][NUMxb];
//Number of hadrons
extern Double32_t nbr_mchadU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mchadU_err[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcxbU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcq2U[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcyU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcdepolU[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcq2_globU[NUMq2][NUMz];
extern Double32_t nbr_mczU[NUMq2][NUMxb][NUMz];

//Number of DIS events
extern Double32_t nbr_mcdisD[NUMq2][NUMxb];
//Number of hadrons
extern Double32_t nbr_mchadD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mchadD_err[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcxbD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcq2D[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcyD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcdepolD[NUMq2][NUMxb][NUMz];
extern Double32_t nbr_mcq2_globD[NUMq2][NUMz];
extern Double32_t nbr_mczD[NUMq2][NUMxb][NUMz];

extern Double_t Lumi;

//particle masses [GeV]
extern const double mass_e= 0.00051;
extern const double mass_mu=0.10566;
extern const double mass_pi= 0.13957;
extern const double mass_k= 0.49368;
extern const double mass_p=0.93827;
extern const double particle_mass[5] = {mass_e,mass_mu,mass_pi,mass_k,mass_p};

#endif
#endif
