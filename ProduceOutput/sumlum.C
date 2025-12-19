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
#include "TChain.h"
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

using namespace std;

const Int_t NUMy=12;
const Int_t NUMq2=13;
const Int_t NUMxb=25;
const Int_t NUMw=10;
const Int_t NUMz=13;
const Int_t NUMpt=14;

Double32_t q2_bins[NUMq2+1];
Double32_t xb_bins[NUMxb+1];
Double32_t z_bins[NUMz+1];

//Rec
Float_t NPos;
Float_t NNeg;
Float_t Ntot;

Float_t scale_lum[4];

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

//mc
Float_t mcNPos;
Float_t mcNNeg;
Float_t mcNtot;

//Number of hadrons
Double32_t nbr_mchadU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mchadU_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcxbU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2U[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2_globU[NUMq2][NUMz];
Double32_t nbr_mczU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcyU[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcdepolU[NUMq2][NUMxb][NUMz];

//Number of DIS events
Double32_t nbr_mcdisD[NUMq2][NUMxb];
//Number of hadrons
Double32_t nbr_mchadD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mchadD_err[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcxbD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2D[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcq2_globD[NUMq2][NUMz];
Double32_t nbr_mczD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcyD[NUMq2][NUMxb][NUMz];
Double32_t nbr_mcdepolD[NUMq2][NUMxb][NUMz];

//average bin
Double_t av_q2[NUMq2][NUMz][NUMxb];
Double_t av_xb[NUMq2][NUMz][NUMxb];
Double_t av_z[NUMq2][NUMz][NUMxb];

Double_t mcav_q2[NUMq2][NUMz][NUMxb];
Double_t mcav_xb[NUMq2][NUMz][NUMxb];
Double_t mcav_z[NUMq2][NUMz][NUMxb];

Int_t pid;
char en[40];
char charge[20];
char part_id[20];
Float_t lumi1, lumi10, lumi100,lumi1000;

void get_binning();
void reset();
void init();
void sum_input();
void write_out();

void sumlum(Int_t const id,char const * enn)
{
  pid=id;
  snprintf(en,40,"%s",enn);

  printf("PID %d \n",pid);
  printf("Energy %s \n",en);
  get_binning();
  reset();
  init();
  sum_input();
  write_out();
}

void init()
{  
  //Make plot name
  if(abs(pid)==211)
    snprintf(part_id,20,"pions");
  else if(abs(pid)==321)
    snprintf(part_id,20,"kaons");
  else if(abs(pid)==2212)
    snprintf(part_id,20,"protons");
  if(pid>0)
    snprintf(charge,20,"pos");
  else
    snprintf(charge,20,"neg");
  
  //lumi value in pb^-1
  if(strcmp(en,"18x275")==0)
    {
      lumi1=6.73;
      lumi10=72.1;
      lumi100=1480;
      lumi1000=63000;
    }
  else if(strcmp(en,"5x41")==0)
    {
      lumi1=12.5;
      lumi10=246;
      lumi100=17400;
      lumi1000=0;
    }
  else
    {
      printf("Unknown lumi \n");
      exit(1);
    }

  scale_lum[0]=1.;
  scale_lum[1]=lumi1/(lumi10);
  scale_lum[2]=lumi1/(lumi100);
  //scale_lum[1]=lumi1/(lumi1+lumi10);
  //scale_lum[2]=lumi1/(lumi1+lumi10+lumi100);
  if(strcmp(en,"18x275")==0)
    scale_lum[3]=lumi1/(lumi1000);
  else
    scale_lum[3]=0.0;

  for(Int_t i=0;i<4;i++)
    printf("scale %d %f \n",i,scale_lum[i]);
}

void sum_input()
{
  char tmp1[50],tmp2[50],tmp3[50],tmp4[50],tmp5[50],tmp6[50],tmp7[50],tmp8[50],tmp9[50],tmp10[50],tmp11[50],tmp12[50],tmp13[50],tmp14[50];
  FILE* f;
  char inFileName[100];
  char inmcFileName[100];

  for(Int_t i=0;i<4;i++)
    {
      Int_t qr=pow(10,i);
      snprintf(inFileName,100,"input/input_withdepol/had_%sQ%d_%s%s.out",en,qr,part_id,charge);
      snprintf(inmcFileName,100,"input/input_withdepol/mc_had_%sQ%d_%s%s.out",en,qr,part_id,charge);

      if(strcmp(en,"5x41")==0 && i==3)
	continue;
      
      printf("data file %s \n", inFileName);
      printf("mc file %s \n", inmcFileName);
      //Open reconstructed data
      f=fopen(inFileName,"r");
      if(!f)
	printf("Problem opening file %s \n",inFileName);

      fscanf(f, "%s %s %s \n",tmp1,tmp2,tmp3);
      if(i==0)
	{
	  Ntot=atof(tmp1);
	  NPos=atof(tmp2);
	  NNeg=atof(tmp3);
	}
      printf("i %d Ntot etc %f %f %f \n",i,atof(tmp1),atof(tmp2),atof(tmp3));
      for(Int_t binq=0;binq<NUMq2;binq++)
        for(Int_t binx=0;binx<NUMxb;binx++)
          for(Int_t binz=0;binz<NUMz;binz++)
            {
	      Int_t qscal;
              fscanf(f, "%s %s %s %s %s %s %s %s %s %s %s %s \n",tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12);
	      //if(binq==6 && binx==18 && binz==0)
	      //if(binq==12 && binx==22 && binz==1)
	      //{
	      //  printf("bins q x z %d %d %d \n",binq,binx,binz);
	      //  printf("had %f %f %f %f %f %f %f %f %f %f %f %f \n",atof(tmp1),atof(tmp2),atof(tmp3),atof(tmp4),atof(tmp5),atof(tmp6),atof(tmp7),atof(tmp8),atof(tmp9),atof(tmp10),atof(tmp11),atof(tmp12));
	      //}
	      qscal=i;
	      //if(binq<4)
	      //qscal=0;
		//if(binq>=4 && binq<8)
		//qscal=1;
		//if(binq>=8 && binq<12)
		//qscal=2;
		//if(binq>=12)
		//qscal=3;
	      //printf("binq %d scale %f \n",binq,scale_lum[qscal]);
	      
              nbr_hadU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp1));
              nbr_hadU_err[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp2));
              nbr_hadD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp3));
              nbr_hadD_err[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp4));
              nbr_q2U[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp5));
	      nbr_q2D[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp6));
	      nbr_xbU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp7));
              nbr_xbD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp8));
              nbr_zU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp9));
              nbr_zD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp10));
	      nbr_yU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp11));
	      nbr_yD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp12));
            }
      
      //Open MC data
      f=fopen(inmcFileName,"r");
      if(!f)
	printf("Problem opening file %s \n",inmcFileName);
      fscanf(f, "%s %s %s \n",tmp1,tmp2,tmp3);
      if(i==0)
	{
	  mcNtot=atof(tmp1);
	  mcNPos=atof(tmp2);
	  mcNNeg=atof(tmp3);
	}
      for(Int_t binq=0;binq<NUMq2;binq++)
        for(Int_t binx=0;binx<NUMxb;binx++)
          for(Int_t binz=0;binz<NUMz;binz++)
            {
	      Int_t qscal=i;
	      //if(binq<4)
	      //qscal=0;
	      //if(binq>=4 && binq<8)
	      //qscal=1;
	      //if(binq>=8 && binq<12)
	      //qscal=2;
	      //if(binq>=12)
	      //qscal=3;
	      //printf("MC binq %d scale %f \n",binq,scale_lum[qscal]);

              fscanf(f, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s \n",tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14);
              nbr_mchadU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp1));
	      nbr_mchadU_err[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp2));
              nbr_mchadD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp3));
	      nbr_mchadD_err[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp4));
	      nbr_mcq2U[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp5));
              nbr_mcq2D[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp6));
	      nbr_mcxbU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp7));
	      nbr_mcxbD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp8));
              nbr_mczU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp9));
	      nbr_mczD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp10));
              nbr_mcyU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp11));
	      nbr_mcyD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp12));
              nbr_mcdepolU[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp13));
	      nbr_mcdepolD[binq][binx][binz]+=(scale_lum[qscal]*atof(tmp14));
	      //if(binq==6 && binx==18 && binz==0)
	      //if(binq==12 && binx==22 && binz==1)
	      //{
	      //  printf("bins q x z %d %d %d \n",binq,binx,binz);
	      //  printf("MC had %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",atof(tmp1),atof(tmp2),atof(tmp3),atof(tmp4),atof(tmp5),atof(tmp6),atof(tmp7),atof(tmp8),atof(tmp9),atof(tmp10),atof(tmp11),atof(tmp12),atof(tmp13),atof(tmp14));
	      //}
            }
      fclose(f);
    }
  
  //scale data and MC
  for(Int_t binq=0;binq<NUMq2;binq++)
    for(Int_t binx=0;binx<NUMxb;binx++)
      for(Int_t binz=0;binz<NUMz;binz++)
	{
	  Float_t fac;
	  
	  if(binq<4)
	    fac=1.;
	  if(binq>=4 && binq<8)
	    fac=2.;
	  if(binq>=8 && binq<12)
	    fac=3.;
	  if(binq>=12)
	    fac=4.;
	  //printf("binq %d scale %f \n",binq,scale_lum[qscal]);
	  
	  nbr_hadU[binq][binx][binz]/=fac;
	  nbr_hadU_err[binq][binx][binz]/=fac;
	  nbr_hadD[binq][binx][binz]/=fac;
	  nbr_hadD_err[binq][binx][binz]/=fac;
	  nbr_q2U[binq][binx][binz]/=fac;
	  nbr_q2D[binq][binx][binz]/=fac;
	  nbr_xbU[binq][binx][binz]/=fac;
	  nbr_xbD[binq][binx][binz]/=fac;
	  nbr_zU[binq][binx][binz]/=fac;
	  nbr_zD[binq][binx][binz]/=fac;
	  nbr_yU[binq][binx][binz]/=fac;
	  nbr_yD[binq][binx][binz]/=fac;

	  nbr_mchadU[binq][binx][binz]/=fac;
	  nbr_mchadU_err[binq][binx][binz]/=fac;
	  nbr_mchadD[binq][binx][binz]/=fac;
	  nbr_mchadD_err[binq][binx][binz]/=fac;
	  nbr_mcq2U[binq][binx][binz]/=fac;
	  nbr_mcq2D[binq][binx][binz]/=fac;
	  nbr_mcxbU[binq][binx][binz]/=fac;
	  nbr_mcxbD[binq][binx][binz]/=fac;
	  nbr_mczU[binq][binx][binz]/=fac;
	  nbr_mczD[binq][binx][binz]/=fac;
	  nbr_mcyU[binq][binx][binz]/=fac;
	  nbr_mcyD[binq][binx][binz]/=fac;
	  nbr_mcdepolU[binq][binx][binz]/=fac;
	  nbr_mcdepolD[binq][binx][binz]/=fac;
	  printf("nbr had U %f had D %f dpeol %f %f \n",nbr_hadU[binq][binx][binz],nbr_hadD[binq][binx][binz],nbr_mcdepolU[binq][binx][binz],nbr_mcdepolD[binq][binx][binz]);
	}
}

void write_out()
{
  FILE* f;
  char FileName[100];
  char mcFileName[100];
  Double_t LPos;
  Double_t LNeg;
  Double_t Ltot;

  printf("Ntot and more %f %f %f \n",Ntot,NPos,NNeg);
  printf("mcNtot and more %f %f %f \n",mcNtot,mcNPos,mcNNeg);
  
  Ltot=(Ntot/5000000)*lumi1;
  LPos=(NPos/5000000)*lumi1;
  LNeg=(NNeg/5000000)*lumi1;

  snprintf(FileName,100,"output/withdepol/scaledhad_%s_%s%s.out",en,part_id,charge);
  snprintf(mcFileName,100,"output/withdepol/scaledmc_had_%s_%s%s.out",en,part_id,charge);

  printf("data file %s \n", FileName);
  printf("mc file %s \n", mcFileName);

  //Open data
  f=fopen(FileName,"w");
  if(!f)
    printf("Problem opening file %s \n",FileName);

  //fprintf(f, "%f %f %f \n",Ntot,NPos,NNeg);
  fprintf(f, "%f %f %f \n",Ltot,LPos,LNeg);
  for(Int_t binq=0;binq<NUMq2;binq++)
    for(Int_t binx=0;binx<NUMxb;binx++)
      for(Int_t binz=0;binz<NUMz;binz++)
	{
	  //if(binq==6 && binx==18 && binz==0)
	  //  printf("written out data %f %f %f %f %f %f %f %f %f %f %f %f \n",nbr_hadU[binq][binx][binz],nbr_hadU_err[binq][binx][binz],nbr_hadD[binq][binx][binz],nbr_hadD_err[binq][binx][binz],nbr_q2U[binq][binx][binz],nbr_q2D[binq][binx][binz],nbr_xbU[binq][binx][binz],nbr_xbD[binq][binx][binz],nbr_zU[binq][binx][binz],nbr_zD[binq][binx][binz],nbr_yU[binq][binx][binz],nbr_yD[binq][binx][binz]);
	  
	  fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f \n",nbr_hadU[binq][binx][binz],nbr_hadU_err[binq][binx][binz],nbr_hadD[binq][binx][binz],nbr_hadD_err[binq][binx][binz],nbr_q2U[binq][binx][binz],nbr_q2D[binq][binx][binz],nbr_xbU[binq][binx][binz],nbr_xbD[binq][binx][binz],nbr_zU[binq][binx][binz],nbr_zD[binq][binx][binz],nbr_yU[binq][binx][binz],nbr_yD[binq][binx][binz]);
	}
  fclose(f);

  //Open MC
  f=fopen(mcFileName,"w");
  if(!f)
    printf("Problem opening file %s \n",mcFileName);

  //fprintf(f, "%f %f %f \n",mcNtot,mcNPos,mcNNeg);
  Ltot=(mcNtot/5000000)*lumi1;
  LPos=(mcNPos/5000000)*lumi1;
  LNeg=(mcNNeg/5000000)*lumi1;
  fprintf(f, "%f %f %f \n",Ltot,LPos,LNeg);
 
  for(Int_t binq=0;binq<NUMq2;binq++)
    for(Int_t binx=0;binx<NUMxb;binx++)
      for(Int_t binz=0;binz<NUMz;binz++)
	{
	  //if(binq==6 && binx==18 && binz==0)
	  //printf("written out MC %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",nbr_mchadU[binq][binx][binz],nbr_mchadU_err[binq][binx][binz],nbr_mchadD[binq][binx][binz],nbr_mchadD_err[binq][binx][binz],nbr_mcq2U[binq][binx][binz],nbr_mcq2D[binq][binx][binz],nbr_mcxbU[binq][binx][binz],nbr_mcxbD[binq][binx][binz],nbr_mczU[binq][binx][binz],nbr_mczD[binq][binx][binz],nbr_mcyU[binq][binx][binz],nbr_mcyD[binq][binx][binz],nbr_mcdepolU[binq][binx][binz],nbr_mcdepolD[binq][binx][binz]);
	  fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",nbr_mchadU[binq][binx][binz],nbr_mchadU_err[binq][binx][binz],nbr_mchadD[binq][binx][binz],nbr_mchadD_err[binq][binx][binz],nbr_mcq2U[binq][binx][binz],nbr_mcq2D[binq][binx][binz],nbr_mcxbU[binq][binx][binz],nbr_mcxbD[binq][binx][binz],nbr_mczU[binq][binx][binz],nbr_mczD[binq][binx][binz],nbr_mcyU[binq][binx][binz],nbr_mcyD[binq][binx][binz],nbr_mcdepolU[binq][binx][binz],nbr_mcdepolD[binq][binx][binz]);
	}
  fclose(f);
}

void reset()
{
  //Rec
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

  //MC
  mcNPos=0.;
  mcNNeg=0.;
  mcNtot=0.;

  //fprintf(f, "%d %d %d \n",Ntot,NPos,NNeg);
  memset(nbr_mchadU,0,sizeof(nbr_mchadU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mchadU_err,0,sizeof(nbr_mchadU_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcxbU,0,sizeof(nbr_mcxbU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2U,0,sizeof(nbr_mcq2U[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2_globU,0,sizeof(nbr_mcq2_globU[0][0])*NUMq2*NUMz);
  memset(nbr_mczU,0,sizeof(nbr_mczU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcyU,0,sizeof(nbr_mcyU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcdepolU,0,sizeof(nbr_mcdepolU[0][0][0])*NUMq2*NUMxb*NUMz);

  memset(nbr_mchadD,0,sizeof(nbr_mchadD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mchadD_err,0,sizeof(nbr_mchadD_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcxbD,0,sizeof(nbr_mcxbD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2D,0,sizeof(nbr_mcq2D[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2_globD,0,sizeof(nbr_mcq2_globD[0][0])*NUMq2*NUMz);
  memset(nbr_mczD,0,sizeof(nbr_mczD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcyD,0,sizeof(nbr_mcyD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcdepolD,0,sizeof(nbr_mcdepolD[0][0][0])*NUMq2*NUMxb*NUMz);
}



void get_binning()
{
  //get binning in x
  ifstream x_in("binning/x-binsunpol.txt");
  for(Int_t i=0;i<NUMxb+1;i++)
    x_in >> xb_bins[i];
  x_in.close();
  printf("x binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMxb+1;i++)
    printf("xbins %f \n",xb_bins[i]);
  printf("\n");

  //get binning in q2
  ifstream q2_in("binning/q2-binsunpol.txt");
  for(Int_t i=0;i<NUMq2+1;i++)
    q2_in >> q2_bins[i];
  q2_in.close();
  printf("q2 binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMq2+1;i++)
    printf("q2bins %f \n",q2_bins[i]);
  printf("\n");

  //get binning in z
  ifstream z_in("binning/z-binsunpol.txt");
  for(Int_t i=0;i<NUMz+1;i++)
    z_in >> z_bins[i];
  z_in.close();
  printf("z binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMz+1;i++)
    printf("zbins %f \n",z_bins[i]);
  printf("\n");
}
