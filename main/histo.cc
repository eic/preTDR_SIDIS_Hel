void chain_input(const char * filename)
{
  ifstream flist;
  char inFileName[300];

  //Declare new Chain
  ch=new TChain("events");
  printf("filename %s \n",filename);
  //open file containing list of root files
  flist.open(filename);
  assert(flist);
  //read file containing list of root files
  while(flist >> inFileName)
    {
      printf("Processing %s \n",inFileName);
      ch->Add(inFileName);
    }
  flist.close();
 }

void get_binning()
{
  //get binning in x
  ifstream x_in("/eic/u/cvhulse/epic/analysis/binning/x-binsunpol.txt");
  for(Int_t i=0;i<NUMxb+1;i++)
    x_in >> xb_bins[i];
  x_in.close();
  printf("x binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMxb+1;i++)
    printf("xbins %f \n",xb_bins[i]);
  printf("\n");

  //get binning in y
  ifstream y_in("/eic/u/cvhulse/epic/analysis/binning/y-binsunpol.txt");
  for(Int_t i=0;i<NUMy+1;i++)
    y_in >> y_bins[i];
  y_in.close();
  printf("y binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMy+1;i++)
    printf("ybins %f \n",y_bins[i]);
  printf("\n");

  //get binning in W2
  ifstream w_in("/eic/u/cvhulse/epic/analysis/binning/w2-binsunpol.txt");
  for(Int_t i=0;i<NUMw+1;i++)
    w_in >> w_bins[i];
  w_in.close();
  printf("W2 binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMw+1;i++)
    printf("wbins %f \n",w_bins[i]);
  printf("\n");

  //get binning in q2
  ifstream q2_in("/eic/u/cvhulse/epic/analysis/binning/q2-binsunpol.txt");
  for(Int_t i=0;i<NUMq2+1;i++)
    q2_in >> q2_bins[i];
  q2_in.close();
  printf("q2 binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMq2+1;i++)
    printf("q2bins %f \n",q2_bins[i]);
  printf("\n");

  //get binning in z
  ifstream z_in("/eic/u/cvhulse/epic/analysis/binning/z-binsunpol.txt");
  for(Int_t i=0;i<NUMz+1;i++)
    z_in >> z_bins[i];
  z_in.close();
  printf("z binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMz+1;i++)
    printf("zbins %f \n",z_bins[i]);
  printf("\n");

  //get binning in pt
  ifstream pt_in("/eic/u/cvhulse/epic/analysis/binning/pt-binsunpol.txt");
  for(Int_t i=0;i<NUMpt+1;i++)
    pt_in >> pt_bins[i];
  pt_in.close();
  printf("pt binning \n");
  printf("========== \n");
  for(Int_t i=0;i<NUMpt+1;i++)
    printf("ptbins %f \n",pt_bins[i]);
  printf("\n");
}


void reset()
{
  memset(nbr_disU,0,sizeof(nbr_disU[0][0])*NUMq2*NUMxb);
  memset(nbr_hadU,0,sizeof(nbr_hadU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_hadU_err,0,sizeof(nbr_hadU_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_xbU,0,sizeof(nbr_xbU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2U,0,sizeof(nbr_q2U[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_yU,0,sizeof(nbr_yU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_depolU,0,sizeof(nbr_depolU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2_globU,0,sizeof(nbr_q2_globU[0][0])*NUMq2*NUMz);
  memset(nbr_zU,0,sizeof(nbr_zU[0][0][0])*NUMq2*NUMxb*NUMz);

  memset(nbr_disD,0,sizeof(nbr_disD[0][0])*NUMq2*NUMxb);
  memset(nbr_hadD,0,sizeof(nbr_hadD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_hadD_err,0,sizeof(nbr_hadD_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_xbD,0,sizeof(nbr_xbD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2D,0,sizeof(nbr_q2D[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_yD,0,sizeof(nbr_yD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_depolD,0,sizeof(nbr_depolD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_q2_globD,0,sizeof(nbr_q2_globD[0][0])*NUMq2*NUMz);
  memset(nbr_zD,0,sizeof(nbr_zD[0][0][0])*NUMq2*NUMxb*NUMz);

  memset(nbr_mcdisU,0,sizeof(nbr_mcdisU[0][0])*NUMq2*NUMxb);
  memset(nbr_mchadU,0,sizeof(nbr_mchadU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mchadU_err,0,sizeof(nbr_mchadU_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcxbU,0,sizeof(nbr_mcxbU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2U,0,sizeof(nbr_mcq2U[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcyU,0,sizeof(nbr_mcyU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcdepolU,0,sizeof(nbr_mcdepolU[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2_globU,0,sizeof(nbr_mcq2_globU[0][0])*NUMq2*NUMz);
  memset(nbr_mczU,0,sizeof(nbr_mczU[0][0][0])*NUMq2*NUMxb*NUMz);

  memset(nbr_mcdisD,0,sizeof(nbr_mcdisD[0][0])*NUMq2*NUMxb);
  memset(nbr_mchadD,0,sizeof(nbr_mchadD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mchadD_err,0,sizeof(nbr_mchadD_err[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcxbD,0,sizeof(nbr_mcxbD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2D,0,sizeof(nbr_mcq2D[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcyD,0,sizeof(nbr_mcyD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcdepolD,0,sizeof(nbr_mcdepolD[0][0][0])*NUMq2*NUMxb*NUMz);
  memset(nbr_mcq2_globD,0,sizeof(nbr_mcq2_globD[0][0])*NUMq2*NUMz);
  memset(nbr_mczD,0,sizeof(nbr_mczD[0][0][0])*NUMq2*NUMxb*NUMz);
}

void initiate(void)
{
  l=new TLorentzVector();
  p=new TLorentzVector();
  lscatmc=new TLorentzVector();
  lscatrec=new TLorentzVector();

  id=new TTreeReaderArray<int>((*tr),"ReconstructedChargedParticles.PDG");
  px=new TTreeReaderArray<float>((*tr),"ReconstructedChargedParticles.momentum.x");
  py=new TTreeReaderArray<float>((*tr),"ReconstructedChargedParticles.momentum.y");
  pz=new TTreeReaderArray<float>((*tr),"ReconstructedChargedParticles.momentum.z");
  m=new TTreeReaderArray<float>((*tr),"ReconstructedChargedParticles.mass");

  relsim=new TTreeReaderArray<UInt_t>((*tr),"ReconstructedChargedParticleAssociations.simID");
  relrec=new TTreeReaderArray<UInt_t>((*tr),"ReconstructedChargedParticleAssociations.recID");
  relall= new TTreeReaderArray<UInt_t>((*tr), "ReconstructedParticleAssociations.simID");

  genid=new TTreeReaderArray<int>((*tr),"GeneratedParticles.PDG");
  gentyp=new TTreeReaderArray<int>((*tr),"GeneratedParticles.type");
  genpx=new TTreeReaderArray<float>((*tr),"GeneratedParticles.momentum.x");
  genpy=new TTreeReaderArray<float>((*tr),"GeneratedParticles.momentum.y");
  genpz=new TTreeReaderArray<float>((*tr),"GeneratedParticles.momentum.z");
  genE=new TTreeReaderArray<float>((*tr),"GeneratedParticles.energy");
  genm=new TTreeReaderArray<float>((*tr),"GeneratedParticles.mass");

  mcid=new TTreeReaderArray<int>((*tr),"MCParticles.PDG");
  //status: 1: final-state part; 11-19: beam part; 21-29: part from hardest sub-process; 2: secondary part?
  //4: incoming particles
  mcgenstat=new TTreeReaderArray<int>((*tr),"MCParticles.generatorStatus");
  mcsimstat=new TTreeReaderArray<int>((*tr),"MCParticles.simulatorStatus");
  mccharge = new TTreeReaderArray<float>((*tr),"MCParticles.charge");
  mcm = new TTreeReaderArray<double>((*tr),"MCParticles.mass");
  mcbvx = new TTreeReaderArray<double>((*tr),"MCParticles.vertex.x");
  mcbvy = new TTreeReaderArray<double>((*tr),"MCParticles.vertex.y");
  mcbvz = new TTreeReaderArray<double>((*tr),"MCParticles.vertex.z");
  mcevx = new TTreeReaderArray<double>((*tr),"MCParticles.endpoint.x");
  mcevy = new TTreeReaderArray<double>((*tr),"MCParticles.endpoint.y");
  mcevz = new TTreeReaderArray<double>((*tr),"MCParticles.endpoint.z");
  mcpx = new TTreeReaderArray<double>((*tr),"MCParticles.momentum.x");
  mcpy = new TTreeReaderArray<double>((*tr),"MCParticles.momentum.y");
  mcpz = new TTreeReaderArray<double>((*tr),"MCParticles.momentum.z");
  mcparb= new TTreeReaderArray<UInt_t>((*tr),"MCParticles.parents_begin");
  mcpare= new TTreeReaderArray<UInt_t>((*tr),"MCParticles.parents_end");
  mcpar = new TTreeReaderArray<Int_t>((*tr),"_MCParticles_parents.index");
  mcdaub= new TTreeReaderArray<UInt_t>((*tr),"MCParticles.daughters_begin");
  mcdaue= new TTreeReaderArray<UInt_t>((*tr),"MCParticles.daughters_end");
  mcdau = new TTreeReaderArray<Int_t>((*tr),"_MCParticles_daughters.index");

  mcscatid= new TTreeReaderArray<Int_t>((*tr),"MCScatteredElectrons_objIdx.index");
  scatid= new TTreeReaderArray<Int_t>((*tr),"MCScatteredElectronAssociations_objIdx.index");

  tabmcQ2=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsTruth.Q2");
  tabmcxB=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsTruth.x");
  tabmcy=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsTruth.y");
  tabmcW=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsTruth.W");

  tabrecQ2=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsElectron.Q2");
  tabrecxB=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsElectron.x");
  tabrecy=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsElectron.y");
  tabrecW=new TTreeReaderArray<float>((*tr), "InclusiveKinematicsElectron.W");
}

void write_file(const char * outDirName)
{
  char outfile[300];
  FILE *f;

  //write all away
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

  Double_t error=0;

  //Write out to sum
  sprintf(outfile,"%s/had_%s%s.out",outDirName,part_id,charge);
  f=fopen(outfile,"w");
  assert(f);
  fprintf(f, "%d %d %d \n",Ntot,NPos,NNeg);
  for(Int_t binq=0;binq<NUMq2;binq++)
    for(Int_t binx=0;binx<NUMxb;binx++)
      for(Int_t binz=0;binz<NUMz;binz++)
        fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f \n",nbr_hadU[binq][binx][binz],nbr_hadU_err[binq][binx][binz],nbr_hadD[binq][binx][binz],nbr_hadD_err[binq][binx][binz],nbr_q2U[binq][binx][binz],nbr_q2D[binq][binx][binz],nbr_xbU[binq][binx][binz],nbr_xbD[binq][binx][binz],nbr_zU[binq][binx][binz],nbr_zD[binq][binx][binz],nbr_yU[binq][binx][binz],nbr_yD[binq][binx][binz]);
  fclose(f);

  sprintf(outfile,"%s/had_mc_%s%s.out",outDirName,part_id,charge);
  f=fopen(outfile,"w");
  assert(f);
  fprintf(f, "%d %d %d \n",Ntot,NPos,NNeg);
  for(Int_t binq=0;binq<NUMq2;binq++)
    for(Int_t binx=0;binx<NUMxb;binx++)
      for(Int_t binz=0;binz<NUMz;binz++)
        fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",nbr_mchadU[binq][binx][binz],nbr_mchadU_err[binq][binx][binz],nbr_mchadD[binq][binx][binz],nbr_mchadD_err[binq][binx][binz],nbr_mcq2U[binq][binx][binz],nbr_mcq2D[binq][binx][binz],nbr_mcxbU[binq][binx][binz],nbr_mcxbD[binq][binx][binz],nbr_mczU[binq][binx][binz],nbr_mczD[binq][binx][binz],nbr_mcyU[binq][binx][binz],nbr_mcyD[binq][binx][binz],nbr_depolU[binq][binx][binz],nbr_depolD[binq][binx][binz]);
  fclose(f);
}
