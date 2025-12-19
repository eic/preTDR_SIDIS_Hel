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

Double32_t asym[NUMq2][NUMz][NUMxb];
Double32_t asymerr[NUMq2][NUMz][NUMxb];
Double32_t asymsys[NUMq2][NUMz][NUMxb];
Double32_t asymtoterr[NUMq2][NUMz][NUMxb];
//average bin
Double_t av_q2[NUMq2][NUMz][NUMxb];
Double_t av_xb[NUMq2][NUMz][NUMxb];
Double_t av_z[NUMq2][NUMz][NUMxb];

Double_t ratio_err[NUMq2][NUMz][NUMxb];

Int_t pid;
char en[40];

void get_input();
void init();
void get_binning();
void plot_graph();
void write_file();

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny, Float_t lMargin, Float_t rMargin,Float_t bMargin, Float_t tMargin);

char part_id[20];
char charge[20];

void plot_sysPub(Int_t const id,char const * enn)
{
  pid=id;
  snprintf(en,40,"%s",enn);
  printf("PID %d \n",pid);
  printf("Energy %s \n",en);
  get_binning();
  init();
  get_input();
  plot_graph();
  write_file();
}


void plot_graph()
{
  char plot_name[100];
  Color_t color;
  char sleg[60];
  TH2F* frame;
  TGaxis *H;
  Double_t err_var[NUMq2][NUMz][NUMxb];
  Double_t err0_var[NUMq2][NUMz][NUMxb];
  // Number of PADS
  const Int_t Nx =3;//first y bin is empty
  const Int_t Ny = 2;
  TPad *pad[Nx][Ny];

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
  
  snprintf(plot_name,100,"output_plots/sys_%s_%s%s.eps",en,charge,part_id);
  TPostScript *ps= new TPostScript(plot_name,113);
  printf("plotname %s \n",plot_name);

  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(kFALSE);
  gStyle->SetEndErrorSize(8);
  
  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1500,1000);

  //Margins
  Float_t lMargin = 0.12;
  Float_t rMargin = 0.12;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.07;

  ps->NewPage();
  //Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);//Because of the left and right margins, all the pads do not have the
  TLegend* leg2;

  Double_t ymax=-999.;
    
  for(Int_t binq=0;binq<NUMq2;binq++)
    for(Int_t binx=0;binx<NUMxb;binx++)
      for(Int_t binz=0;binz<NUMz;binz++)
	{
	  //printf("xbins %f %f \n",xb_bins[binx],xb_bins[binx+1]);
	  //printf("logs %f %f ",fabs(log(xb_bins[binx])),fabs(log(xb_bins[binx+1])));
	  err0_var[binq][binz][binx]=0.;
	  err_var[binq][binz][binx]=0.5*(fabs(log10(xb_bins[binx]))-fabs(log10(xb_bins[binx+1])));
	  if(ymax<asym[binq][binz][binx])
	    ymax=asym[binq][binz][binx];
	}
  ymax=1.2;
  //printf("YMAX %f \n",ymax);
  TGraphErrors* gr[NUMq2][NUMz+1];//create one dummy graph for legend pad
  TGraphErrors* grsys[NUMq2][NUMz+1];//create one dummy graph for legend pad
  for(Int_t binq=0;binq<NUMq2;binq++)
    {
      for(Int_t binz=0;binz<NUMz;binz++)
	{
	  gr[binq][binz]=new TGraphErrors(NUMxb,av_xb[binq][binz],asym[binq][binz],err0_var[binq][binz],asymerr[binq][binz]);
	  grsys[binq][binz]=new TGraphErrors(NUMxb,av_xb[binq][binz],asym[binq][binz],err_var[binq][binz],asymtoterr[binq][binz]);
	}
      Int_t binz=NUMz-1;
      gr[binq][binz+1]=new TGraphErrors(NUMxb,av_xb[binq][binz],asym[binq][binz],err_var[binq][binz],err_var[binq][binz]);
      grsys[binq][binz+1]=new TGraphErrors(NUMxb,av_xb[binq][binz],asym[binq][binz],err_var[binq][binz],asymtoterr[binq][binz]);

      int marker=20+binq;
      if(marker==28)
        marker+=4;
      if(marker==31)
        marker+=2;

      if(binq==0)
        color=kGreen+2;
      if(binq==1)
        color=kTeal-7;
      if(binq==2)
        color=kCyan+2;
      if(binq==3)
        color=kAzure-9;
      if(binq==4)
        color=kBlue;
      if(binq==5)
        color=kBlue+1;
      if(binq==6)
        color=kViolet+8;
      if(binq==7)
        color=kViolet-3;
      if(binq==8)
        color=kMagenta;
      if(binq==9)
        color=kPink;
      if(binq==10)
        color=kRed+1;
      if(binq==11)
        color=kOrange+7;
      if(binq==12)
        color=kBlack;

      for(Int_t binz=0;binz<NUMz+1;binz++)
        {
	  gr[binq][binz]->SetMarkerSize(1.4);
	  grsys[binq][binz]->SetMarkerSize(1.);

	  marker=20;
          gr[binq][binz]->SetMarkerStyle(marker);
          gr[binq][binz]->SetMarkerColor(color);
          gr[binq][binz]->SetLineColor(color);
	  gr[binq][binz]->SetLineWidth(1);
		  
	  //marker=25;
	  grsys[binq][binz]->SetMarkerStyle(marker);
          grsys[binq][binz]->SetMarkerColor(color);
          grsys[binq][binz]->SetLineColor(color);
	  grsys[binq][binz]->SetLineWidth(1.);
	  //grsys[binq][binz]->SetFillColor(color);
	  //grsys[binq][binz]->SetFillStyle(1001);
	  grsys[binq][binz]->SetFillStyle(3002);
	  grsys[binq][binz]->SetFillColorAlpha(color, 0.5);
        }
    }

  for (Int_t wx=0;wx<Nx;wx++)//can apparantly not just loop over z, cf. needs to follow order
    {
      for (Int_t wy=0;wy<Ny;wy++)
        {
	  //binz=0, 3, 6, 9, 12
          Int_t binz=9*(1-wy)+3*wx;
	  if(binz==0)
            binz=1;
          if(binz==12)
            binz=11;
          printf("BINZ %d \n",binz);
	  
	  //Int_t binz=4*(3-wy)+wx;
          C->cd(0);
          // Get the pads previosly created.
          char pname[20];
          snprintf(pname,20,"pad_%i_%i",wx,wy);
          pad[wx][wy] = (TPad*) gROOT->FindObject(pname);
          pad[wx][wy]->Draw();
          pad[wx][wy]->cd();
	  // Size factors
	  Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[wx][wy]->GetAbsWNDC();
	  Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[wx][wy]->GetAbsHNDC();
	  frame=new TH2F("","",200,log10(1e-05)+0.0001,log10(1)-0.0001,200,0.0001,ymax);
	  frame->Draw("");
	  frame->SetStats(kFALSE);
	  // Format for y axis
          frame->GetYaxis()->SetLabelFont(43);
          frame->GetYaxis()->SetLabelSize(28);
          frame->GetYaxis()->SetLabelOffset(0.02);
          frame->GetYaxis()->SetTitleFont(43);
          frame->GetYaxis()->SetTitleSize(28);
          frame->GetYaxis()->SetTitleOffset(1.7);
          frame->GetYaxis()->SetNdivisions(505);
          // TICKS Y Axis
          frame->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

          // Format for x axis
          frame->GetXaxis()->SetLabelFont(43);
          frame->GetXaxis()->SetLabelSize(28);
          frame->GetXaxis()->SetLabelOffset(0.02);
          frame->GetXaxis()->SetTitleFont(43);
          frame->GetXaxis()->SetTitleSize(28);
          frame->GetXaxis()->SetTitleOffset(1.7);
          frame->GetXaxis()->SetNdivisions(505);
          // TICKS X Axis
          frame->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);

          if((wx==Nx-1) && wy==0)
            frame->GetXaxis()->SetTitle("log(x_{B})");
          if(wx==0 && wy==Ny-1)
            frame->GetYaxis()->SetTitle("Uncertainty_{A_{1}}");
	  for(Int_t binq=0;binq<NUMq2;binq++)
	    {
	      if(binz<11)
		{
		  grsys[binq][binz]->Draw("2");//2
		  gr[binq][binz]->Draw("p");//||
		  //gr[binq][binz]->Draw("||");//||
		  //grsys[binq][binz]->Draw("a2");//2
		  //grsys[binq][binz]->SetLineWidth(-802);
		}
	    }
	  if(wx==Nx-1 && wy==0)
            {
              //TLegend* leg = new TLegend(-3.,-0.4,-1.,0.4);
	      //TLegend* leg = new TLegend(0.0001,0.45,0.7,0.98);
	      TLegend* leg = new TLegend(0.03,0.325,0.6,0.905);
              leg->SetTextFont(43);
              leg->SetTextSize(24);
	      
              for(Int_t binq=0;binq<NUMq2;binq++)
                {
		  if(binq==0)
                    snprintf(sleg,60,"1.00<Q^{2} [GeV^{2}]<1.78");
                  else if(binq==1)
                    snprintf(sleg,60,"1.78<Q^{2} [GeV^{2}]<3.16");
                  else if(binq==2)
                    snprintf(sleg,60,"3.16<Q^{2} [GeV^{2}]<5.62");
                  else if(binq==3)
                    snprintf(sleg,60,"5.62<Q^{2} [GeV^{2}]<10");
                  else if(binq==4)
                    snprintf(sleg,60,"10<Q^{2} [GeV^{2}]<17.8");
                  else if(binq==5)
                    snprintf(sleg,60,"17.8<Q^{2} [GeV^{2}]<31.6");
                  else if(binq==6)
                    snprintf(sleg,60,"31.6<Q^{2} [GeV^{2}]<56.2");
                  else if(binq==7)
                    snprintf(sleg,60,"56.2<Q^{2} [GeV^{2}]<100");
                  else if(binq==8)
                    snprintf(sleg,60,"100<Q^{2} [GeV^{2}]<177");
                  else if(binq==9)
                    snprintf(sleg,60,"177<Q^{2} [GeV^{2}]<316");
                  else if(binq==10)
                    snprintf(sleg,60,"316<Q^{2} [GeV^{2}]<562");
                  else if(binq==11)
                    snprintf(sleg,60,"562<Q^{2} [GeV^{2}]<1000");
                  else if(binq==12)
                    snprintf(sleg,60,"1000<Q^{2} [GeV^{2}]<10000");
                  if(binq==0 || binq==2 || binq==4 || binq==6 || binq==8 || binq==10 || binq==12)
                    leg->AddEntry(gr[binq][NUMz],sleg,"p");
		}
	      leg->SetFillStyle(0);
              leg->SetBorderSize(0);
              leg->Draw();
	    } 
	  if(wy==Ny-1)
	    {
	      TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 0.75,log10(1)-0.0001, 0.75,log10(1e-05)+0.0001, log10(1)-0.0001,100, "B-U");//"GS" for Log and respect thick size if defined
	      H->CenterTitle();
	      H->SetTitleFont(43);
	      H->SetTitleSize(28);
	      H->SetTitleOffset(1.5);
	      
	      if(wx==0)
		{
		  TLatex Text_ePIC;
		  Text_ePIC.SetTextSize(28);
		  Text_ePIC.SetTextFont(43);
		  Text_ePIC.DrawText(-4.6,1.05,"ePIC Performance");  // performance plot
		  TLatex l;
		  l.SetTextSize(28);
		  l.SetTextFont(43);
		  //l.DrawText(-4,0.64,"0.011 < z < 0.05");
		  //l.DrawText(-4,1.05,"0.011 < z < 0.05");
		  //l.DrawText(-4,1.05,"0.05 < z < 0.10");
 		  l.DrawText(-4,0.9,"0.05 < z < 0.10");
		}
	      if(wx==1)
		{
		  TLatex l;
		  l.SetTextSize(28);
		  l.SetTextFont(43);
		  //l.DrawText(-4,1.05,"0.05 < z < 0.10");
		  //l.DrawText(-4,1.05,"0.15 < z < 0.20");
		  l.DrawText(-4,0.9,"0.15 < z < 0.20");
		}
	      if(wx==2)
		{
		  TLatex l;
		  l.SetTextSize(28);
		  l.SetTextFont(43);
		  //l.DrawText(-4,1.05,"0.10 < z < 0.15");
		  //l.DrawText(-4,1.05,"0.30 < z < 0.40");
		  l.DrawText(-4,0.9,"0.30 < z < 0.40");
		}
	      if(wx==3)
		{
		  TLatex l;
		  l.SetTextSize(28);
		  l.SetTextFont(43);
		  //l.DrawText(-4,1.05,"0.15 < z < 0.20");
		}
	      H->Draw();
	    }
	  if(wy==2)
            {
              //TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 0.75, log10(1)-0.0001, 0.75,log10(1e-05)+0.0001, log10(1)-0.0001 , 100, "B+U");
	      TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 1.18, log10(1)-0.0001, 1.18,log10(1e-05)+0.0001, log10(1)-0.0001 , 100, "B+U");

              H->CenterTitle();
              H->SetTitleFont(43);
              H->SetTitleSize(28);
              H->SetTitleOffset(1.9);

              if(wx==0)
                H->SetTitle("0.20 < z < 0.25");
              if(wx==1)
                H->SetTitle("0.25 < z < 0.30");
              if(wx==2)
                H->SetTitle("0.30 < z < 0.40");
              if(wx==3)
                H->SetTitle("0.40 < z < 0.50");

              H->Draw();
            }
	  if(wy==1)
            {
	      //TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 0.75, log10(1)-0.0001, 0.75,log10(1e-05)+0.0001, log10(1)-0.0001 , 100, "B+U");
	      TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 1.18, log10(1)-0.0001, 1.18,log10(1e-05)+0.0001, log10(1)-0.0001 , 100, "B+U");
              H->CenterTitle();
              H->SetTitleFont(43);
              H->SetTitleSize(28);
              H->SetTitleOffset(1.9);
	      
              /*if(wx==0)
                H->SetTitle("0.50 < z < 0.60");
              if(wx==1)
                H->SetTitle("0.60 < z < 0.70");
              if(wx==2)
                H->SetTitle("0.70 < z < 0.80");
              if(wx==3)
	      H->SetTitle("0.80 < z < 0.90");*/
              H->Draw();
            }
	  if(wy==0)
            {
	      //TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 0.75, log10(1)-0.0001, 0.75,log10(1e-05)+0.0001, log10(1)-0.0001 , 100, "B+U");
	      //TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 1.2, log10(1)-0.0001, 1.2,log10(1e-05)+0.0001, log10(1)-0.0001 , 100, "B+U");
	      TGaxis *H=new TGaxis(log10(1e-05)+0.0001, 1.3, log10(1)-0.0001, 1.3,log10(1e-05)+0.0001, log10(1)-0.0001 , 100, "B+U");
              H->CenterTitle();
              H->SetTitleFont(43);
              H->SetTitleSize(28);
              H->SetTitleOffset(1.5);

	      if(wx==0)
                //H->SetTitle("0.60 < z < 0.70");
		{
		  TLatex l;
		  l.SetTextSize(28);
		  l.SetTextFont(43);
		  //l.DrawText(-4,1.05,"0.10 < z < 0.15");
		  //l.DrawText(-4,1.05,"0.30 < z < 0.40");
		  l.DrawText(-4,0.9,"0.60 < z < 0.70");
		}
	      
              //H->SetTitle("0.90 < z < 1.00");
              //if(wx==1)
	      //H->SetTitle("0.80 < z < 0.90");
              H->Draw();
            }

	  /*if(wx==Nx-1)
	    {
	      //TGaxis *H=new TGaxis(log10(1)-0.0001, 0.75, log10(1)-0.0002,-0.5, 0.75,-0.5, 100, "B-U");
	      TGaxis *H=new TGaxis(log10(1)-0.0001, 1.18, log10(1)-0.0002,-0.5, 1.18,-0.5, 100, "B-U");
	      
              H->CenterTitle();
              H->SetTitleFont(43);
              H->SetTitleSize(30);
              H->SetTitleOffset(3);

              H->Draw();
	      }*/
	
	  // if(wx==Nx-2 && wy==0)
	   if(wx==Nx-1 && wy==0)
            {
	      //TLegend* leg = new TLegend(0.01,0.93,0.3,0.98);
	      TLegend* leg = new TLegend(0.01,0.9,0.3,0.95);
              leg->SetTextFont(43);
              leg->SetTextSize(28);

              if(pid==211)
                snprintf(sleg,20,"#pi^{+}");
              else if(pid==-211)
                snprintf(sleg,20,"#pi^{-}");
              else if(pid==321)
                snprintf(sleg,20,"K^{+}");
              else if(pid==-321)
                snprintf(sleg,20,"K^{-}");
              else if(pid==2212)
                snprintf(sleg,20,"p");
              else if(pid==-2212)
                snprintf(sleg,20,"#bar{p}");
	      //To not point to anything. One can even use "l" option
              leg->AddEntry((TObject*)0, sleg, "");
              leg->SetFillStyle(0);
              leg->SetBorderSize(0);
              leg->Draw();
	    }
          if(wx==Nx-2 && wy==0)
            {
	      /*TMathText l;
              l.SetTextSize(28);
              l.SetTextFont(43);
              //l.DrawMathText(-4.6,0.74,"\\mathscr{L}=10 \\text{ fb}^{-1}");
 	      l.DrawMathText(-4.4,1.02,"\\mathscr{L}_{\\text{proj}}=10 \\text{ fb}^{-1}");
	      */
	      //TLegend* legen = new TLegend(0.01,0.84,0.3,0.88);
	      //TLegend* legen = new TLegend(0.01,0.88,0.3,0.90);
	      //TLegend* legen = new TLegend(0.01,0.88,0.3,0.90);
	      TLegend* legen = new TLegend(0.01,0.8,0.3,0.95);
              legen->SetTextFont(43);
              legen->SetTextSize(28);
	      //if(strcmp(en,"18x275")==0)
	      //snprintf(sleg,60,"ep, 18x275 GeV^{2}; 70#% pol.");
	      //else if(strcmp(en,"5x41")==0)
	      //snprintf(sleg,60,"ep, 5x41 GeV^{2}; 70#% pol.");
	      //To not point to anything. One can even use "l" option
              //legen->AddEntry((TObject*)0, sleg, "");
	      legen->AddEntry((TObject*)0, "ep, 18x275 GeV^{2}; 70% pol.", "");
	      snprintf(sleg,60,"\\mathscr{L}_{\\text{proj}} = 10 \\text{ fb}^{-1}");
	      legen->AddEntry((TObject*)0, sleg, "");
	      legen->SetFillStyle(0);
              legen->SetBorderSize(0);
              legen->Draw();

	      //TMathText l;
              //l.SetTextSize(28);
              //l.SetTextFont(43);
              //l.DrawMathText(-4.6,0.74,"\\mathscr{L}=10 \\text{ fb}^{-1}");
 	      //l.DrawMathText(-4.4,0.9,"\\mathscr{L}_{\\text{proj}}=10 \\text{ fb}^{-1}");

	      //TLegend* leg2 = new TLegend(0.0001,0.645,0.7,0.745);
	      //TLegend* leg2 = new TLegend(0.046,0.53,0.7,0.76);
	      //TLegend* leg2 = new TLegend(0.05,0.67,0.7,0.87);
	      TLegend* leg2 = new TLegend(0.05,0.55,0.7,0.79);
              leg2->SetTextFont(43);
              leg2->SetTextSize(28);
              leg2->AddEntry(gr[12][NUMz],"stat. unc.","ep");
              leg2->AddEntry(grsys[12][NUMz],"total unc.","f");
	      leg2->AddEntry((TObject*)0, "2% pol. unc.", "");
              leg2->SetFillStyle(0);
              leg2->SetBorderSize(0);
              leg2->Draw();
 	    }
	}
    }
  C->Update();
  C->cd();
  ps->Close();
  
  char plot_name_pdf[200];
  snprintf(plot_name_pdf,200,"output_plots/sysPub_%s_%s%s.pdf",en,charge,part_id);
  gSystem->Setenv("file_ps",plot_name);
  gSystem->Setenv("file_pdf",plot_name_pdf);
  gSystem->Exec("ps2pdf $file_ps tmp.pdf");
  gSystem->Exec("pdfcrop --noverbose tmp.pdf $file_pdf");
  gSystem->Exec("rm $file_ps");
  gSystem->Exec("rm tmp.pdf");
}

void get_input()
{
  char tmp1[50],tmp2[50],tmp3[50],tmp4[50],tmp5[50],tmp6[50],tmp7[50],tmp8[50],tmp9[50],tmp10[50];
  FILE* f;
  char FileName[100];

  if(pid==211)
    snprintf(FileName,100,"output/withdepol/asymunf_%s_pionspos.out",en);
  else if(pid==-211)
    snprintf(FileName,100,"output/withdepol/asymunf_%s_pionsneg.out",en);
  else if(pid==321)
    snprintf(FileName,100,"output/withdepol/asymunf_%s_kaonspos.out",en);
  else if(pid==-321)
    snprintf(FileName,100,"output/withdepol/asymunf_%s_kaonsneg.out",en);
  else
    {
      printf("Unknown particle type\n");
      exit(1);
    }
  f=fopen(FileName,"r");
  assert(f);
  if(!f)
    printf("Problem opening file %s \n",FileName);

  for(Int_t binq=0; binq<NUMq2;binq++)
    {
      Double_t as=(Double_t)2/(NUMq2+1)+(Double_t)binq/(NUMq2+1);
      for(Int_t binz=0; binz<NUMz;binz++)
	for(Int_t binx=0; binx<NUMxb;binx++)
	  {
	    fscanf(f, "%s %s %s %s %s %s \n",tmp1,tmp2,tmp3,tmp4,tmp5,tmp6);

	    if(atof(tmp1)!=-999.)
	      asym[binq][binz][binx]=as;
	    else
	      asym[binq][binz][binx]=-999.;
	    if(!isnan(atof(tmp2)))
	      asymerr[binq][binz][binx]=atof(tmp2);
	    else
	      asymerr[binq][binz][binx]=-999;
	    //asymsys[binq][binz][binx]=atof(tmp3);
	    av_q2[binq][binz][binx]=atof(tmp4);
	    av_z[binq][binz][binx]=atof(tmp5);
	    av_xb[binq][binz][binx]=atof(tmp6);
	  }
    }
  fclose(f);

  //Read in file for systematic uncertainty
  if(pid==211)
    snprintf(FileName,100,"output/withoutdepol/asymunf_%s_pionspos.out",en);
  else if(pid==-211)
    snprintf(FileName,100,"output/withoutdepol/asymunf_%s_pionsneg.out",en);
  else if(pid==321)
    snprintf(FileName,100,"output/withoutdepol/asymunf_%s_kaonspos.out",en);
  else if(pid==-321)
    snprintf(FileName,100,"output/withoutdepol/asymunf_%s_kaonsneg.out",en);
  else
    {
      printf("Unknown particle type\n");
      exit(1);
    }
  f=fopen(FileName,"r");
  assert(f);
  if(!f)
    printf("Problem opening file %s \n",FileName);
  
  for(Int_t binq=0; binq<NUMq2;binq++)
    {
      for(Int_t binz=0; binz<NUMz;binz++)
	for(Int_t binx=0; binx<NUMxb;binx++)
	  {
	    fscanf(f, "%s %s %s %s %s %s \n",tmp1,tmp2,tmp3,tmp4,tmp5,tmp6);
	    asymsys[binq][binz][binx]=atof(tmp3);
	    
	    if(asymerr[binq][binz][binx]==-999 || asymsys[binq][binz][binx]==-999)
	      {
		asym[binq][binz][binx]=-999;
		asymerr[binq][binz][binx]=-999;
		asymsys[binq][binz][binx]=-999;
		asymtoterr[binq][binz][binx]=0.;
	      }
	    else
	      asymtoterr[binq][binz][binx]=sqrt(asymerr[binq][binz][binx]*asymerr[binq][binz][binx]+asymsys[binq][binz][binx]*asymsys[binq][binz][binx]);

	  }
    }
  fclose(f);
}

void write_file()
{
  char tmp1[50],tmp2[50],tmp3[50],tmp4[50],tmp5[50],tmp6[50],tmp7[50],tmp8[50],tmp9[50],tmp10[50];
  FILE* f;
  char FileName[100];

  FILE* fout;
  char FileNameOut[100];

  if(pid==211)
    snprintf(FileNameOut,100,"output/Unc_asymunf_%s_pionspos.out",en);
  else if(pid==-211)
    snprintf(FileNameOut,100,"output/Unc_asymunf_%s_pionsneg.out",en);
  else if(pid==321)
    snprintf(FileNameOut,100,"output/Unc_asymunf_%s_kaonspos.out",en);
  else if(pid==-321)
      snprintf(FileNameOut,100,"output/Unc_asymunf_%s_kaonsneg.out",en);
  else
    {
      printf("Unknown particle type\n");
      exit(1);
    }
  fout=fopen(FileNameOut,"w");
  assert(fout);
  if(!fout)
    printf("Problem opening file %s \n",FileName);
  
  for(Int_t binq=0; binq<NUMq2;binq++)
    for(Int_t binz=0; binz<NUMz;binz++)
      {
        for(Int_t binx=0; binx<NUMxb;binx++)
          {
            fprintf(fout,"%f %f %f %f %f %f \n",asym[binq][binz][binx],asymerr[binq][binz][binx],asymsys[binq][binz][binx],av_q2[binq][binz][binx],av_z[binq][binz][binx],av_xb[binq][binz][binx]);
          }
        fprintf(fout,"\n");
      }
  fclose(fout);
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

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,Float_t lMargin, Float_t rMargin,Float_t bMargin, Float_t tMargin)
{
   if (!C) return;
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
   for (Int_t i=0;i<Nx;i++) {
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
      for (Int_t j=0;j<Ny;j++) {
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         }
	 else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         }
	 else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
	 C->cd(0);
         char name[16];
         snprintf(name,16,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
         pad->Draw();
      }
   }
}
