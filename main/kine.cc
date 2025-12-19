Int_t get_spin(void)
{
  //generate a random number in the range 0 to 99 in order to randomize orientaion of thrust axis
  //direction of thrust axis is well recontructed, but orientation is biased to forward direction
  Int_t random=rand()%100;

  if (random<50)
    return -1;
  else
    return 1;
  return 0;
}

void get_kine(TLorentzVector* lscat, Double_t &Q2, Double_t &y, Double_t &xB, Double_t &W2)
{
  TLorentzVector q(l->Px()-lscat->Px(),l->Py()-lscat->Py(),l->Pz()-lscat->Pz(),l->E()-lscat->E());
  Q2=-q*q;
  y=((*p)*q)/((*p)*(*l));
  xB=(Q2/(2*(*p)*q));
  W2=mass_p*mass_p+(1-xB)*Q2/xB;
}


Int_t get_parton(void)
{
  //Not available at the moment.
  //if(ev->hepmcp_procid!=99 && ev->hepmcp_procid!=131 && ev->hepmcp_procid!=132 && ev->hepmcp_procid!=133 && ev->hepmcp_procid!=134 && ev->hepmcp_procid!=135 && ev->hepmcp_procid!=136)
  //return -999;
 
  PARTON=-999;
  for(Int_t imc=0;imc<mcid->GetSize();imc++)
    {
      //      printf("IMC %d mcid %d mcgenstat %d mcsimstat %d gentyp %d \n",imc, (*mcid)[imc],(*mcgenstat)[imc],(*mcsimstat)[imc],(*gentyp)[imc]);
      if(((abs((*mcid)[imc])>=1 && abs((*mcid)[imc])<=6) || abs((*mcid)[imc])==21) && (*mcgenstat)[imc]==61)
	{
	  PARTON=(*mcid)[imc];
	  break;
	}
    }
  return PARTON;
}
