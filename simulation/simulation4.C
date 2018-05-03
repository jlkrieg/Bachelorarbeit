#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TF2.h>

Double_t deg2rad = TMath::Pi()/180;
Double_t rad2deg = 180/TMath::Pi();

Double_t calc_az(Double_t el,Double_t az, Double_t a){
  return atan((tan(a)-cos(el)*tan(az))/(tan(a)*tan(az)+cos(el)));
  }

void simulation4(){
  TFile *file= TFile::Open("run281_cut.root");
  TNtuple *nt=(TNtuple*)file->Get("nt");
  Int_t N=nt->GetEntries();
  Double_t eld[N],azd[N],ela[N],aza[N];
  Double_t el[N],az[N],del[N],daz[N]; //radiant
  const Float_t* a = nt->GetArgs();
  for (int i=0; i<N; i++){
    nt->GetEntry(i);
    azd[i] = a[1];
    eld[i] = a[2];
    aza[i] = a[3];
    ela[i] = a[4];
    el[i] = ela[i]*deg2rad;
    az[i] = aza[i]*deg2rad;
    del[i]=(eld[i]-ela[i])*deg2rad;
    daz[i]=(azd[i]-aza[i])*deg2rad;
  }
  Double_t aa=0.265;
  Double_t chi=0;
  for(int i=0; i<N; i++){
    chi=chi+pow((calc_az(el[i],az[i],aa)),2);
  }
  std::cout<<chi<<std::endl;
}
