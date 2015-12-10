#include "spectrum.h"
#include "batch.h"
#include "filters.h"

#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

template<class InputIterator> void fun1(InputIterator it, InputIterator endIt) {
  vector<float> dataV;
  //  cout<<(it==endIt)<<endl;
  
  it->serialize(dataV);
  cout<<dataV.size()<<endl;
}

int main(int argc, char **argv) {
  vector<Results_PA> pairsPA;
  vector<Results_ASP> pairsASP;
  int idx;

  cout<<"Load_resultsASPbin: "<<Load_resultsASPbin("testASP.bin",pairsASP)<<endl;
  cout<<"  -> "<<pairsASP.size()<<" pairs";
  idx = (int)pairsASP.size()-1;
  if(idx<0) cout<<endl; else cout<<", last pair = ("<<pairsASP[idx].spec1<<","<<pairsASP[idx].spec2<<","<<pairsASP[idx].shift1<<","<<pairsASP[idx].score1<<","<<pairsASP[idx].score2<<")\n";
  Save_resultsASPbin("testASPout1.bin",pairsASP);
  
  cout<<"Load_results_bin: "<<Load_results_bin("testASP.bin",pairsASP)<<endl;
  cout<<"  -> "<<pairsASP.size()<<" pairs";
  idx = (int)pairsASP.size()-1;
  if(idx<0) cout<<endl; else cout<<", last pair = ("<<pairsASP[idx].spec1<<","<<pairsASP[idx].spec2<<","<<pairsASP[idx].shift1<<","<<pairsASP[idx].score1<<","<<pairsASP[idx].score2<<")\n";
  Save_results_bin("testASPout2.bin",pairsASP.size(),pairsASP.begin());
}
