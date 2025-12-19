#define MULTI_GLOBAL
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "math.h"
#include "TChain.h"
#include "TTree.h"
#include "TVector.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCut.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TBranch.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TPostScript.h>
#include <numeric>

using namespace std;

const Int_t NUMy=12;
const Int_t NUMq2=13;
const Int_t NUMxb=25;
const Int_t NUMw=10;
const Int_t NUMz=13;
const Int_t NUMpt=14;

Int_t NPos;
Int_t NNeg;
Int_t Ntot;

Int_t PID;
char mc[20];
char en[50];

//Number of hadrons
Double32_t nbr_hadU[NUMq2][NUMxb][NUMz];
Double32_t nbr_hadU_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_xbU[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2U[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2_globU[NUMq2][NUMz];
Double32_t nbr_zU[NUMq2][NUMxb][NUMz];
Double32_t nbr_yU[NUMq2][NUMxb][NUMz];

//Number of DIS events
Double32_t nbr_disD[NUMq2][NUMxb];
//Number of hadrons
Double32_t nbr_hadD[NUMq2][NUMxb][NUMz];
Double32_t nbr_hadD_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_xbD[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2D[NUMq2][NUMxb][NUMz];
Double32_t nbr_q2_globD[NUMq2][NUMz];
Double32_t nbr_zD[NUMq2][NUMxb][NUMz];
Double32_t nbr_yD[NUMq2][NUMxb][NUMz];

char outDirName[400];

void reset();
void write_file();

Int_t main(int argc, char *argv[])
{
  char tmp1[40],tmp2[40],tmp3[40],tmp4[40],tmp5[40],tmp6[40],tmp7[40],tmp8[40],tmp9[40],tmp10[40],tmp11[40],tmp12[40];
  FILE* f;
  const char* inFileName = argv[1];
  strcpy(outDirName,argv[2]);
  PID=atoi(argv[3]);
  strcpy(mc,argv[4]);
  strcpy(en,argv[5]);

  printf("Processing \n");
  printf("input %s \n",inFileName);
  printf("ouput %s \n",outDirName);
  printf("PID %d \n",PID);
  printf("mc %s \n",mc);
  printf("energy %s \n",en);

  ifstream flist;
  char inputFile[300];
  reset();

  flist.open(inFileName);
  assert(flist);

  //read file containing list of root files
  while(flist >> inputFile)
    {
      printf("Processing %s \n",inputFile);
      //Write out to sum
      f=fopen(inputFile,"r");
      assert(f);
      fscanf(f, "%s %s %s \n",tmp1,tmp2,tmp3);
      Ntot+=atoi(tmp1);
      NPos+=atoi(tmp2);
      NNeg+=atoi(tmp3);

      for(Int_t binq=0;binq<NUMq2;binq++)
	for(Int_t binx=0;binx<NUMxb;binx++)
	  for(Int_t binz=0;binz<NUMz;binz++)
	    {
	      fscanf(f, "%s %s %s %s %s %s %s %s %s %s %s %s \n",tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12);
	      nbr_hadU[binq][binx][binz]+=atof(tmp1);
	      nbr_hadU_err[binq][binx][binz]+=atof(tmp2);
	      nbr_hadD[binq][binx][binz]+=atof(tmp3);
	      nbr_hadD_err[binq][binx][binz]+=atof(tmp4);
	      nbr_q2U[binq][binx][binz]+=atof(tmp5);
	      nbr_q2D[binq][binx][binz]+=atof(tmp6);
	      nbr_xbU[binq][binx][binz]+=atof(tmp7);
	      nbr_xbD[binq][binx][binz]+=atof(tmp8);
	      nbr_zU[binq][binx][binz]+=atof(tmp9);
	      nbr_zD[binq][binx][binz]+=atof(tmp10);
	      nbr_yU[binq][binx][binz]+=atof(tmp11);
	      nbr_yD[binq][binx][binz]+=atof(tmp12);
	    }
      fclose(f);
    }
  flist.close();

  write_file();
}

void reset()
{
  NPos=0.;
  NNeg=0.;
  Ntot=0.;
  
  //fprintf(f, "%d %d %d \n",Ntot,NPos,NNeg);
  memset(nbr_hadU,0,sizeof(nbr_hadU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_hadU_err,0,sizeof(nbr_hadU_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_xbU,0,sizeof(nbr_xbU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2U,0,sizeof(nbr_q2U[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2_globU,0,sizeof(nbr_q2_globU[0][0])*NUMq2*NUMz);
  memset(nbr_zU,0,sizeof(nbr_zU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_yU,0,sizeof(nbr_yU[0][0][0])*NUMq2*NUMxb*NUMz);

  memset(nbr_hadD,0,sizeof(nbr_hadD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_hadD_err,0,sizeof(nbr_hadD_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_xbD,0,sizeof(nbr_xbD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2D,0,sizeof(nbr_q2D[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2_globD,0,sizeof(nbr_q2_globD[0][0])*NUMq2*NUMz);
  memset(nbr_zD,0,sizeof(nbr_zD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_yD,0,sizeof(nbr_yD[0][0][0])*NUMq2*NUMxb*NUMz);
}

//write all away
void write_file()
{
  FILE *f;
  char outfile[700];
  char part_id[40];
  char charge[40];

  if(abs(PID)==211)
    sprintf(part_id,"pions");
  else if(abs(PID)==321)
    sprintf(part_id,"kaons");
  else if(abs(PID)==2212)
    sprintf(part_id,"protons");
  if(PID>0)
    sprintf(charge,"pos");
  else
    sprintf(charge,"neg");

  //Write out to sum
  if(strcmp("mc",mc)==0)
    sprintf(outfile,"%s/mc_had_%s_%s%s.out",outDirName,en,part_id,charge);
  else
    sprintf(outfile,"%s/had_%s_%s%s.out",outDirName,en,part_id,charge);
  f=fopen(outfile,"w");
  assert(f);
  fprintf(f, "%d %d %d \n",Ntot,NPos,NNeg);
  for(Int_t binq=0;binq<NUMq2;binq++)
    for(Int_t binx=0;binx<NUMxb;binx++)
      for(Int_t binz=0;binz<NUMz;binz++)
        fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f \n",nbr_hadU[binq][binx][binz],nbr_hadU_err[binq][binx][binz],nbr_hadD[binq][binx][binz],nbr_hadD_err[binq][binx][binz],nbr_q2U[binq][binx][binz],nbr_q2D[binq][binx][binz],nbr_xbU[binq][binx][binz],nbr_xbD[binq][binx][binz],nbr_zU[binq][binx][binz],nbr_zD[binq][binx][binz],nbr_yU[binq][binx][binz],nbr_yD[binq][binx][binz]);
  fclose(f);
}
 
