#ifndef _test2_def
#define _test2_def

#ifdef TTREE_GLOBAL
////Generated particles (Simply final-state MC particles. Redundant with MC particles)
TTreeReaderArray<float> *genpx;
TTreeReaderArray<float> *genpy;
TTreeReaderArray<float> *genpz;
TTreeReaderArray<float> *genE;
TTreeReaderArray<float> *genm;
TTreeReaderArray<int>   *genid;
TTreeReaderArray<int>   *gentyp;
//MC particles (I think After smearing etc.)
TTreeReaderArray<int> *mcid;
TTreeReaderArray<int> *mcgenstat;
TTreeReaderArray<int> *mcsimstat;
TTreeReaderArray<UInt_t> *mcparb;
TTreeReaderArray<UInt_t> *mcpare;
TTreeReaderArray<Int_t> *mcpar;
TTreeReaderArray<UInt_t> *mcdaub;
TTreeReaderArray<UInt_t> *mcdaue;
TTreeReaderArray<Int_t> *mcdau;
TTreeReaderArray<float> *mccharge;
TTreeReaderArray<double> *mcm;
TTreeReaderArray<double> *mcbvx;
TTreeReaderArray<double> *mcbvy;
TTreeReaderArray<double> *mcbvz;
TTreeReaderArray<double> *mcevx;
TTreeReaderArray<double> *mcevy;
TTreeReaderArray<double> *mcevz;
TTreeReaderArray<double> *mcpx;
TTreeReaderArray<double> *mcpy;
TTreeReaderArray<double> *mcpz;

//Reconstructed charged particles
TTreeReaderArray<float> *px;
TTreeReaderArray<float> *py;
TTreeReaderArray<float> *pz;
TTreeReaderArray<float> *m;
TTreeReaderArray<int> *id;
TTreeReaderArray<UInt_t> *relsim;
TTreeReaderArray<UInt_t> *relrec;
TTreeReaderArray<UInt_t> *relall;
TTreeReaderArray<int> *mcscatid;
TTreeReaderArray<int> *scatid;

TTreeReaderArray<float> *tabmcQ2;
TTreeReaderArray<float> *tabrecQ2;
TTreeReaderArray<float> *tabmcxB;
TTreeReaderArray<float> *tabrecxB;
TTreeReaderArray<float> *tabmcy;
TTreeReaderArray<float> *tabrecy;
TTreeReaderArray<float> *tabmcW;
TTreeReaderArray<float> *tabrecW;

#else
////Generated particles (Simply final-state MC particles. Redundant with MC particles)
extern TTreeReaderArray<float> *genpx;
extern TTreeReaderArray<float> *genpy;
extern TTreeReaderArray<float> *genpz;
extern TTreeReaderArray<float> *genE;
extern TTreeReaderArray<float> *genm;
extern TTreeReaderArray<int>   *genid;
extern TTreeReaderArray<int>   *gentyp;
//MC particles (I think After smearing etc.)
extern TTreeReaderArray<int> *mcid;
extern TTreeReaderArray<int> *mcgenstat;
extern TTreeReaderArray<int> *mcsimstat;
extern TTreeReaderArray<UInt_t> *mcparb;
extern TTreeReaderArray<UInt_t> *mcpare;
extern TTreeReaderArray<Int_t> *mcpar;
extern TTreeReaderArray<UInt_t> *mcdaub;
extern TTreeReaderArray<UInt_t> *mcdaue;
extern TTreeReaderArray<Int_t> *mcdau;
extern TTreeReaderArray<float> *mccharge;
extern TTreeReaderArray<double> *mcm;
extern TTreeReaderArray<double> *mcbvx;
extern TTreeReaderArray<double> *mcbvy;
extern TTreeReaderArray<double> *mcbvz;
extern TTreeReaderArray<double> *mcevx;
extern TTreeReaderArray<double> *mcevy;
extern TTreeReaderArray<double> *mcevz;
extern TTreeReaderArray<Float_t> *mcpx;
extern TTreeReaderArray<Float_t> *mcpy;
extern TTreeReaderArray<Float_t> *mcpz;

//Reconstructed charged particles
extern TTreeReaderArray<float> *px;
extern TTreeReaderArray<float> *py;
extern TTreeReaderArray<float> *pz;
extern TTreeReaderArray<float> *m;
extern TTreeReaderArray<int> *id;
extern TTreeReaderArray<UInt_t> *relsim;
extern TTreeReaderArray<UInt_t> *relrec;
extern TTreeReaderArray<UInt_t> *relall;
extern TTreeReaderArray<int> *mcscatid;
extern TTreeReaderArray<int> *scatid;

extern TTreeReaderArray<int> *tabmcQ2;
extern TTreeReaderArray<int> *tabrecQ2;
extern TTreeReaderArray<int> *tabmcxB;
extern TTreeReaderArray<int> *tabrecxB;
extern TTreeReaderArray<int> *tabmcy;
extern TTreeReaderArray<int> *tabrecy;
extern TTreeReaderArray<int> *tabmcW;
extern TTreeReaderArray<int> *tabrecW;
#endif
#endif
