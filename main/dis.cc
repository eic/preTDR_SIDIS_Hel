void fill_dis()
{
  for(Int_t binq=0;binq<NUMq2;binq++)
    if(rQ2>=q2_bins[binq] && rQ2<q2_bins[binq+1])
      for(Int_t binx=0;binx<NUMxb;binx++)
        if(rxB>=xb_bins[binx] && rxB<xb_bins[binx+1])
          {
            if(SPIN>0)
              nbr_disU[binq][binx]++;
            else
              nbr_disD[binq][binx]++;
          }
  for(Int_t binq=0;binq<NUMq2;binq++)
    if(mcQ2>=q2_bins[binq] && mcQ2<q2_bins[binq+1])
      for(Int_t binx=0;binx<NUMxb;binx++)
        if(mcxB>=xb_bins[binx] && mcxB<xb_bins[binx+1])
          {
            if(SPIN>0)
              nbr_mcdisU[binq][binx]++;
            else
              nbr_mcdisD[binq][binx]++;
          }
}

Int_t get_init()
{
  if(l)
    delete l;
  if(p)
    delete p;
  if(lscatmc)
    delete lscatmc;
  if(lscatr)
    delete lscatr;

  l=new TLorentzVector();
  p=new TLorentzVector();
  lscatmc=new TLorentzVector();
  lscatr=new TLorentzVector();

  if(get_initlep())
    return 1;
  if(get_initp())
    return 1;
  if(get_mclepscat())
    return 1;
  if(get_lepscat())
    return 1;
  return 0;
}

Int_t get_initlep()
{
  for(int imc=0;imc<mcid->GetSize();imc++)
    if((*mcid)[imc]==11 && (*mcgenstat)[imc]==4)//beam lepton
      {
	Double_t E=sqrt(((*mcpx)[imc])*((*mcpx)[imc])+((*mcpy)[imc])*((*mcpy)[imc])+((*mcpz)[imc])*((*mcpz)[imc])+mass_e*mass_e);
	l->SetPxPyPzE((*mcpx)[imc],(*mcpy)[imc],(*mcpz)[imc],E);
	return 0;
      }
  return 1;
}

Int_t get_initp()
{
  for(int imc=0;imc<mcid->GetSize();imc++)
    if((*mcid)[imc]==2212 && (*mcgenstat)[imc]==4)//beam lepton
      {
	Double_t E=sqrt(((*mcpx)[imc])*((*mcpx)[imc])+((*mcpy)[imc])*((*mcpy)[imc])+((*mcpz)[imc])*((*mcpz)[imc])+mass_e*mass_e);
	p->SetPxPyPzE((*mcpx)[imc],(*mcpy)[imc],(*mcpz)[imc],E);
	//p->Print();
	return 0;
      }
  return 1;
}

Int_t get_mclepscat()
{
  mcid_escat=-999;
  
  //printf("Searching for scattered lepton \n");
  for(Int_t imc=0;imc<mcid->GetSize();imc++)
    {
      if (mcid_escat>-999)//scattered lepton found
	break;
      if (((*mcid)[imc]==11) && ((*mcgenstat)[imc]==1))//found a final-state lepton
	{
	  //loop over parent particles
	  Int_t imcl=imc;
	  Int_t imcll=-999;
	  Int_t el=0;
	  while(1==1)//loop over (grand)parents until beam lepton parent is found.
	    {
	      if(mcid_escat>-999)//parent beam lepton found
		break;
	      //loop over parents of particle imcl
	      for(Int_t imcb=(*mcparb)[imcl];imcb<(*mcpare)[imcl];imcb++) {
		Int_t ind=(*mcpar)[imcb];//get index of parent
		if((*mcid)[ind]==11 && (*mcgenstat)[ind]==4)//check parent information
		  {
		    //printf("--> scattered lepton found \n");
		    mcid_escat=imc;
		    break;
		  }
		//check if parent at least is an electron. If so, we can check its parents
		else if((*mcid)[ind]==11) {
		  el=1;
		  imcll=ind;//keep parent electron ID to loop over its parents
		}
	      }
	      if(el==0)//no parent with PID 11 found
		break;
	      else//parent with PID 11 but generator status different from 4 found
		//—> assign imcl particle as parent particle and repeat loop
		{
		  imcl=imcll;
		  el=0;
		}
	    }
	}
    }//finished search for scattered beam lepton

  if(mcid_escat>-999)
    {
      Double_t E=sqrt(((*mcpx)[mcid_escat])*((*mcpx)[mcid_escat])+((*mcpy)[mcid_escat])*((*mcpy)[mcid_escat])+((*mcpz)[mcid_escat])*((*mcpz)[mcid_escat])+mass_e*mass_e);
      lscatmc->SetPxPyPzE((*mcpx)[mcid_escat],(*mcpy)[mcid_escat],(*mcpz)[mcid_escat],E);
      return 0;
    }
  return 1;
}

Int_t get_lepscat()
{
  //Search for reconstructed scattered electron:
  //for the moment, it is identified by navigating to the generated scattered lepton
  Int_t recid_e=-999;
  //Find the scattered mc ID and navigate to the reconstructed electron
  for(Int_t irel=0;irel<relrec->GetSize();irel++)
    {
      if((*relsim)[irel]==mcid_escat)
	{
	  recid_e=(*relrec)[irel];
	  break;
	}
    }
  if(recid_e>-999)
    {
      Double_t E=sqrt(((*px)[recid_e])*((*px)[recid_e])+((*py)[recid_e])*((*py)[recid_e])+((*pz)[recid_e])*((*pz)[recid_e])+mass_e*mass_e);
      lscatr->SetPxPyPzE((*px)[recid_e],(*py)[recid_e],(*pz)[recid_e],E);
      return 0;
    }
  return 1;
}

void get_dis(void)
{
  get_kine(lscatr,rQ2,ry,rxB,rW2);
  get_kine(lscatmc,mcQ2,mcy,mcxB,mcW2);

  Double_t theta=1+12*(mcQ2/(mcQ2+1))*(0.125*0.125/(0.125*0.125+mcxB*mcxB));
  Double_t R=(b1/log10(mcQ2/0.04))*theta+((b2/mcQ2)+b3/(mcQ2*mcQ2+0.3*0.3))*(1+b4*mcxB+b5*mcxB*mcxB)*pow(mcxB,b6);
  depol=(2*mcy-mcy*mcy)/(mcy*mcy+2+2*R-2*mcy-2*mcy*R);
  //depol=1.;
}

Int_t check_dis()
{
  if(mcQ2<1 || mcW2<10 || mcy<0.01 || mcy>0.95)
    return 1;

  if(rQ2<1 || rW2<10 || ry<0.01 || ry>0.95)
    return 1;

  if(depol<0.1)
    return 1;

  return 0;
}
