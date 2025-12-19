void get_hadrons(void)
{
  TLorentzVector q(l->Px()-lscatr->Px(),l->Py()-lscatr->Py(),l->Pz()-lscatr->Pz(),l->E()-lscatr->E());
  TLorentzVector qtr(l->Px()-lscatmc->Px(),l->Py()-lscatmc->Py(),l->Pz()-lscatmc->Pz(),l->E()-lscatmc->E());

  Double_t mass=-999;
  switch (abs(PID))
    {
    case 211: mass=mass_pi;
      break;
    case 321: mass=mass_k;
      break;
    case 2212: mass=mass_p;
      break;
    case 11: mass=mass_e;
      break;
    default: mass=mass_pi;
      break;
    }
  for(Int_t i=0;i<id->GetSize();i++)
    {
      Int_t imc=(*relsim)[i];
      if(PID==(*mcid)[imc])
	{
	  Double_t E=sqrt(((*px)[i])*((*px)[i])+((*py)[i])*((*py)[i])+((*pz)[i])*((*pz)[i])+mass*mass);
	  TLorentzVector h((*px)[i],(*py)[i],(*pz)[i],E);
	  Double_t z=((*p)*h)/((*p)*q);
	  Double_t weight=1.;
	  Double_t tmpweight=1.;
	  Double_t fuu;
	  Double_t dq;

	  //for weight
	  Double_t Etr=sqrt(((*mcpx)[imc])*((*mcpx)[imc])+((*mcpy)[imc])*((*mcpy)[imc])+((*mcpz)[imc])*((*mcpz)[imc])+mass*mass);
	  TLorentzVector htr((*mcpx)[imc],(*mcpy)[imc],(*mcpz)[imc],Etr);
	  Double_t ztr=((*p)*htr)/((*p)*qtr);
	  //printf("rec momentum %f %f %f %f \n",(*px)[i],(*py)[i],(*pz)[i],E);
	  //printf("Rec z %f \n",z);
	  if((abs(PARTON)==1 || abs(PARTON)==2 || abs(PARTON)==3  || PARTON==21) && z>0.011)
	    {
	      fuu=get_fuu(mcxB,mcQ2,ztr);
	      dq=get_deltaq(mcxB,mcQ2,ztr);

	      if(SPIN>0)
		tmpweight+=depol*(dq/fuu);
	      else if(SPIN<0)
		tmpweight-=depol*(dq/fuu);
	      if(tmpweight>1.999 || tmpweight<0.001)
		{
		  printf("weight %f \n",tmpweight);
		  printf("dq %f fuu %f \n",dq, fuu);
		}
	      else
		weight=tmpweight;
	      if(fuu<0)
		{
		  printf("¡¡WARNING!!: FUU<0: \n");
		  //printf("x z q2 %f %f %f;  process %d parton %d hadron %d; dq %f fuu %f \n",rxB,z,rQ2,ev->hepmcp_procid,PARTON,PID,dq, fuu);
		}
	    }
	  //printf("Rec weight %f \n",weight);
	  fill_had(weight,z);
	}
    }
}

void get_mchadrons()
{
  TLorentzVector q(l->Px()-lscatmc->Px(),l->Py()-lscatmc->Py(),l->Pz()-lscatmc->Pz(),l->E()-lscatmc->E());
  Double_t mass=-999;
  switch (abs(PID))
    {
    case 211: mass=mass_pi;
      break;
    case 321: mass=mass_k;
      break;
    case 2212: mass=mass_p;
      break;
    case 11: mass=mass_e;
      break;
    default: mass=mass_pi;
      break;
    }
  for(int imc=0;imc<mcid->GetSize();imc++)
    if((*mcid)[imc]==PID && (*mcgenstat)[imc]==1)//beam lepton
      {
        Double_t E=sqrt(((*mcpx)[imc])*((*mcpx)[imc])+((*mcpy)[imc])*((*mcpy)[imc])+((*mcpz)[imc])*((*mcpz)[imc])+mass*mass);
	TLorentzVector h((*mcpx)[imc],(*mcpy)[imc],(*mcpz)[imc],E);
	Double_t z=((*p)*h)/((*p)*q);
	Double_t weight=1.;
	Double_t tmpweight=1.;
	Double_t fuu;
	Double_t dq;

	//printf("MC momentum %f %f %f %f \n", (*mcpx)[imc],(*mcpy)[imc],(*mcpz)[imc],E);
	//printf("mc z %f \n",z);
	if((abs(PARTON)==1 || abs(PARTON)==2 || abs(PARTON)==3  || PARTON==21) && z>0.011)
	  {
	    fuu=get_fuu(mcxB,mcQ2,z);
	    dq=get_deltaq(mcxB,mcQ2,z);
	    if(SPIN>0)
	      tmpweight+=depol*(dq/fuu);
	    else if(SPIN<0)
	      tmpweight-=depol*(dq/fuu);
	    if(tmpweight>1.999 || tmpweight<0.001)
	      {
		printf("MC weight %f \n",tmpweight);
		printf("MC dq %f fuu %f \n",dq, fuu);
	      }
	    else
	      weight=tmpweight;
	    if(fuu<0)
	      {
		printf("¡¡MC WARNING!!: FUU<0: \n");
		//printf("x z q2 %f %f %f;  process %d parton %d hadron %d; dq %f fuu %f \n",mcxB,z,mcQ2,ev->hepmcp_procid,PARTON,PID,dq,fuu);
	      }
	  }
	//printf("MC weight %f \n",weight);
	fill_mchad(weight,z);
      }
}


extern "C" {
  void fuu_(int *ih,int *ic,double* x,double *q2, double* z,double *f1);
}

Double_t get_fuu(Double_t x, Double_t q2,Double_t z)
{
  //UNPOL FFs
  //int is=2;//use DSSV 14(pion) and 17(kaon)
  int ih=-999;
  int ic=-999;
  int nlo=1;
  switch (abs(PID))
    {
    case 211: ih=1;
      break;
    case 321: ih=2;
      break;
    case 2212: ih=3;
      break;
    }
  ic= PID<0? -1 : 1;
  double f1;
  fuu_(&ih,&ic,&x,&q2,&z,&f1);
  return f1;
}

extern "C" {
  void deltaq_(int *ih,int *ic,int *ipart,double* x,double *q2, double* z,double *delta);
}

Double_t get_deltaq(Double_t x, Double_t q2,Double_t z)
{
  //UNPOL FFs
  //int is=2;//use DSSV 14(pion) and 17(kaon)
  int ih=-999;
  int ic=-999;
  int nlo=1;
  switch (abs(PID))
    {
    case 211: ih=1;
      break;
    case 321: ih=2;
      break;
    case 2212: ih=3;
      break;
    }
  ic= PID<0? -1 : 1;
  double delta;
  deltaq_(&ih,&ic,&PARTON,&x,&q2,&z,&delta);
  return delta;
}

void fill_had(Double_t weight, Double_t z)
{
  for(Int_t binq=0;binq<NUMq2;binq++)
    if(rQ2>=q2_bins[binq] && rQ2<q2_bins[binq+1])
      for(Int_t binz=0;binz<NUMz;binz++)
        if(z>=z_bins[binz] && z<z_bins[binz+1])
          {
            if(SPIN>0)
              {
                nbr_q2_globU[binq][binz]+=(weight*rQ2);
              }
            else
              {
                nbr_q2_globD[binq][binz]+=(weight*rQ2);
              }
            for(Int_t binx=0;binx<NUMxb;binx++)
              if(rxB>=xb_bins[binx] && rxB<xb_bins[binx+1])
                {
                  if(SPIN>0)
                    {
                      nbr_hadU[binq][binx][binz]+=weight;
                      nbr_hadU_err[binq][binx][binz]+=(weight*weight);
                      nbr_xbU[binq][binx][binz]+=(weight*rxB);
                      nbr_q2U[binq][binx][binz]+=(weight*rQ2);
                      nbr_yU[binq][binx][binz]+=(weight*ry);
                      nbr_zU[binq][binx][binz]+=(weight*z);
                    }
                  else
                    {
                      nbr_hadD[binq][binx][binz]+=weight;
                      nbr_hadD_err[binq][binx][binz]+=(weight*weight);
                      nbr_xbD[binq][binx][binz]+=(weight*rxB);
                      nbr_q2D[binq][binx][binz]+=(weight*rQ2);
                      nbr_yD[binq][binx][binz]+=(weight*ry);
                      nbr_zD[binq][binx][binz]+=(weight*z);
                    }
                }
          }
}

void fill_mchad(Double_t weight, Double_t z)
{
  for(Int_t binq=0;binq<NUMq2;binq++)
    if(mcQ2>=q2_bins[binq] && mcQ2<q2_bins[binq+1])
      for(Int_t binz=0;binz<NUMz;binz++)
        if(z>=z_bins[binz] && z<z_bins[binz+1])
          {
            if(SPIN>0)
              {
                nbr_mcq2_globU[binq][binz]+=(weight*mcQ2);
              }
            else
              {
                nbr_mcq2_globD[binq][binz]+=(weight*mcQ2);
              }
            for(Int_t binx=0;binx<NUMxb;binx++)
              if(mcxB>=xb_bins[binx] && mcxB<xb_bins[binx+1])
                {
                  if(SPIN>0)
                    {
                      nbr_mchadU[binq][binx][binz]+=weight;
                      nbr_mchadU_err[binq][binx][binz]+=(weight*weight);
                      nbr_mcxbU[binq][binx][binz]+=(weight*mcxB);
                      nbr_mcq2U[binq][binx][binz]+=(weight*mcQ2);
                      nbr_mcyU[binq][binx][binz]+=(weight*mcy);
                      nbr_mczU[binq][binx][binz]+=(weight*z);
                      nbr_depolU[binq][binx][binz]+=(weight*depol);
                    }
                  else
                    {
                      nbr_mchadD[binq][binx][binz]+=weight;
                      nbr_mchadD_err[binq][binx][binz]+=(weight*weight);
                      nbr_mcxbD[binq][binx][binz]+=(weight*mcxB);
                      nbr_mcq2D[binq][binx][binz]+=(weight*mcQ2);
                      nbr_mcyD[binq][binx][binz]+=(weight*mcy);
                      nbr_mczD[binq][binx][binz]+=(weight*z);
                      nbr_depolD[binq][binx][binz]+=(weight*depol);
                    }
                }
          }
}
