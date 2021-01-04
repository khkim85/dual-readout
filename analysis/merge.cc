#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "functions.h"

#include "GeoSvc.h"
#include "GridDRcalo.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TChain.h"
#include <iostream>
#include <string>
#include <exception>
#include <filesystem>
//./merge ../../box/ele_%.root 50 ../../eltest
//./merge ../../box/showers/pi_%d.root 1000 ../../pishower
namespace fs = std::filesystem;
int main(int , char* argv[]) {
  TString innameform = argv[1];
  int num_file = std::stoi(argv[2]);
  TString outname = argv[3];
  //TFile outfile(outname.Data(),"recreate");
  //TFile* file = TFile::Open(outname.Data(), "RECREATE");
  auto drsimchain = TChain("DRsim");
  auto recochain = TChain("Reco");
  TString inname;
  TFile* box;
  TList* keys;
  int count=0;
  for(int i=0; i<num_file;i++){
    inname.Form(innameform.Data(),i);
    if(!fs::exists(inname.Data()))continue;
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    if(keys->GetEntries()>1){
    if(strcmp(keys->At(0)->GetName(),"DRsim")==0 && strcmp(keys->At(1)->GetName(),"Reco")==0){
    printf("%d;\n",i);
    drsimchain.Add(inname);
    recochain.Add(inname);
    count++;
    box->Close();
    }
    }
    //mychain.Add(inname+"/DRsim");

  }
  drsimchain.Merge(outname+"drsim.root");
  printf("sim saved\n");
  recochain.Merge(outname+"reco.root");
  printf("reco saved\n");
  printf("count %d\n",count);
  //recochain.AddFriend("DRsim");
  //recochain.Merge(outname);
  //auto tree1 = recochain.GetTree();
  //auto tree2 = drsimchain.GetTree();

  //outfile.Write();
  //outfile.Close();
}
