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

void simulation2(){

  Double_t deg2rad = TMath::Pi()/180;
  Double_t rad2deg = 180/TMath::Pi();
  TFile *file= TFile::Open("run281_cut.root");
  TNtuple *nt=(TNtuple*)file->Get("nt");
  Int_t N=nt->GetEntries();
  Double_t eld[N],azd[N],ela[N],aza[N];
  Double_t ela_rad[N],eld_rad[N],aza_rad[N],azd_rad[N];
  Double_t del[N],daz[N]; //radiant
  const Float_t* a = nt->GetArgs();
  for (int i=0; i<N; i++){
    nt->GetEntry(i);
    azd[i] = a[1];
    eld[i] = a[2];
    aza[i] = a[3];
    ela[i] = a[4];
    ela_rad[i] = ela[i]*deg2rad;
    eld_rad[i] = eld[i]*deg2rad;
    aza_rad[i] = aza[i]*deg2rad;
    azd_rad[i] = azd[i]*deg2rad;
    del[i]=eld_rad[i]-ela_rad[i];
    daz[i]=azd_rad[i]-aza_rad[i];
  }
  bool twoparam=true;

  Double_t alpha_eld,alpha_azd,alpha_ela,alpha_aza,phi_eld,phi_ela,phi_aza;
  TCanvas *can= new TCanvas("c","c",1200,600);
  gStyle->SetOptStat(1);
  can->Divide(2,2);
  can->cd(1);
  TGraph* g_eld=new TGraph(N,eld_rad,del);
  if (twoparam){
    TF1 *f_eld =new TF1("f_eld","x-asin([0]*sin(x)*cos([1])+cos(x)*sin([1]))+[1]");
    g_eld->Fit(f_eld);
    alpha_eld = acos(f_eld->GetParameter(0))*rad2deg;
    phi_eld = f_eld->GetParameter(1)*rad2deg;
  }
  else{
    TF1 *f_eld =new TF1("f_eld","x-asin([0]*sin(x))");
    g_eld->Fit(f_eld);
    alpha_eld = acos(f_eld->GetParameter(0))*rad2deg;
  }
  g_eld->SetTitle("#Delta Elevation (drive)");
  g_eld->GetXaxis()->SetTitle("elevation drive (rad)");
  g_eld->GetYaxis()->SetTitle("#Delta el (rad)");
  g_eld->Draw("A*");

  can->cd(2);
  can->Update();
  TGraph* g_azd=new TGraph(N,eld_rad,daz);
  TF1 *f_azd =new TF1("f_azd","-atan([0]/cos(x))");
  g_azd->Fit(f_azd);
  alpha_azd = atan(f_azd->GetParameter(0))*rad2deg;
  g_azd->SetTitle("#Delta Azemut (drive)");
  g_azd->GetXaxis()->SetTitle("elevation drive (rad)");
  g_azd->GetYaxis()->SetTitle("#Delta az (rad)");
  g_azd->Draw("A*");

  can->cd(3);
  can->Update();
  TGraph* g_ela=new TGraph(N,ela_rad,del);
  if (twoparam){
    TF1 *f_ela =new TF1("f_ela","asin(sin(x+[1])*[0])-x");
    g_ela->Fit(f_ela,"","",0,1.34);
    alpha_ela = acos(1/f_ela->GetParameter(0))*rad2deg;
    phi_ela = f_ela->GetParameter(1)*rad2deg;
  }
  else{
    TF1 *f_ela =new TF1("f_ela","asin(sin(x)*[0])-x");
    g_ela->Fit(f_ela,"","",0,1.34);
    alpha_ela = acos(1/f_ela->GetParameter(0))*rad2deg;
  }
  g_ela->SetTitle("#Delta Elevation (astro)");
  g_ela->GetXaxis()->SetTitle("elevation astro (rad)");
  g_ela->GetYaxis()->SetTitle("#Delta el (rad)");
  g_ela->Draw("A*");

  can->cd(4);
  can->Update();
  TGraph* g_aza=new TGraph(N,ela_rad,daz);
  // if (twoparam){
  //   TF1 *f_aza =new TF1("f_aza","-atan(sqrt((1-[0])/([0]+cos(2*(x+[1])))))");
  //   f_aza->SetParLimits(0,0.7,1);
  //   g_aza->Fit(f_aza);
  //   alpha_aza = acos(f_aza->GetParameter(0))/2*rad2deg;
  //   phi_aza = f_aza->GetParameter(1)*rad2deg;
  // }
  // else{
  TF1 *f_aza =new TF1("f_aza","-atan(sqrt((1-[0])/([0]+cos(2*x))))");
    f_aza->SetParLimits(0,0.7,1);
    g_aza->Fit(f_aza);
    alpha_aza = acos(f_aza->GetParameter(0))/2*rad2deg;
  // }
  g_aza->SetTitle("#Delta Azemut (astro)");
  g_aza->GetXaxis()->SetTitle("elevation astro (rad)");
  g_aza->GetYaxis()->SetTitle("#Delta az (rad)");
  g_aza->Draw("A*");

  can->Print("simulation2.png");

  std::cout<<"----------------------"<<std::endl;
  std::cout<<"alpha (eld) = "<<alpha_eld<<std::endl;
  std::cout<<"alpha (azd) = "<<alpha_azd<<std::endl;
  std::cout<<"alpha (ela) = "<<alpha_ela<<std::endl;
  std::cout<<"alpha (aza) = "<<alpha_aza<<std::endl;
  if (twoparam){
    std::cout<<"----------------------"<<std::endl;
    std::cout<<"phi (eld)   = "<<phi_eld<<std::endl;
    std::cout<<"phi (ela)   = "<<phi_ela<<std::endl;
    // std::cout<<"phi (aza)   = "<<phi_aza<<std::endl;
  }
}
