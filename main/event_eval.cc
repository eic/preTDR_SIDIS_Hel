#define MULTI_GLOBAL
#define TTREE_GLOBAL
#include "event_eval.h"
#include "event_trees.h"
#include "histo.cc"
#include "kine.cc"
#include "dis.cc"
#include "hadrons.cc"

using namespace std;

extern "C" {
  void init_(void);
}

Int_t main(int argc, char* argv[])
{
  const char* inFileName = argv[1];
  const char* outDirName = argv[2];
  PID=atoi(argv[3]);

  printf("Processing PID %d \n",PID);
  printf("Outfile %s \n",outDirName);
  chain_input(inFileName);
  get_binning();
  
  init_();
  printf("Hello \n");

  NPos=0;
  NNeg=0;
  Ntot=0;

  //reset 3D array containing hadron counts
  reset();

  //get arrays
  int ev(0);
  tr=new TTreeReader(ch);
  initiate();
  
  while ((*tr).Next())
    {
      SPIN=get_spin();
      Ntot++;
      if(SPIN>0)
        NPos++;
      else
        NNeg++;
      PARTON=-999;

      if(get_init())
        continue;
      get_dis();

      if(check_dis())
        continue;
      //printf("Tab MC Q2 %f xB %f y %f W%f \n",(*tabmcQ2)[0],(*tabmcxB)[0],(*tabmcy)[0],(*tabmcW)[0]);
      //printf("My MC Q2 %f xB %f y %f W %f \n", mcQ2,mcxB,mcy,mcW2);
      //printf("My rec Q2 %f xB %f y %f W %f \n", rQ2,rxB,ry,rW2);
      //printf("Tab REC Q2 %f xB %f y %f W%f \n",(*tabrecQ2)[0],(*tabrecxB)[0],(*tabrecy)[0],(*tabrecW)[0]);
      //printf("\n" );

      PARTON=get_parton();
      fill_dis();
      get_hadrons();
      get_mchadrons();      
    }
  write_file(outDirName);
  printf("Processed %d events pos %d neg %d \n",Ntot,NPos, NNeg);
  return 0;
}
