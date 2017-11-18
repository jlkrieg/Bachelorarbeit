#include <iostream>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>

void pedestal(){
  //linear function y=ax+b
  Double_t x[4]={0,7,14,21};
  Double_t a[4]={0.0355944,0.0762582,0.165461,0.369326};
  Double_t aerr[4]={0.00739128,0.00874256,0.0131825,0.0211513};
  Double_t b[4]={1.94881,1.87431,1.75895,1.46777};
  Double_t berr[4]={0.00990898,0.0108411,0.014623,0.0220131};
  TCanvas *can= new TCanvas("c","c",1200,600);
  gStyle->SetOptStat(1);
  can->Divide(2,1);
  can->cd(1);
  TGraphErrors* g1=new TGraphErrors(4,x,a,0,aerr);
  g1->SetTitle("dM/dT;gain;dM/dT");
  g1->SetLineColor(kRed);
  g1->Draw("AP");
  can->cd(2);
  can->Update();
  TGraphErrors* g2=new TGraphErrors(4,x,b,0,berr);
  g2->SetTitle("Pedestal;gain;pedestal");
  g2->SetLineColor(kRed);
  g2->Draw("AP");
  can->Print("plots/pedestal.pdf");
}
