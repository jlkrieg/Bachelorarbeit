#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

double wa(Double_t a[],Double_t b[],int c){
  double temp1=0,temp2=0,wa;
  for (int i=0;i<c;i++){
    temp1=temp1+a[i]/pow(b[i],2);
    temp2=temp2+1/pow(b[i],2);
  }
  wa=temp1/temp2;
  return wa;
}

double waerr(Double_t a[], Double_t b[],int c){
  double temp1=0,temp2=0,waerr;
  for (int i=0;i<c;i++){
    temp1=pow(a[i]/pow(b[i],2),2);
    temp2=temp2+1/pow(b[i],2);
  }
  waerr=sqrt(temp1)/temp2;
  return waerr;
}

void Korrelation3(){
  TFile *file= TFile::Open("imgdata.root");
  TTree *tree=(TTree*)file->Get("ntuple");
  Long64_t entries =tree->GetEntries();
  TLeaf*g =tree->GetLeaf("gain");
  TLeaf*e =tree->GetLeaf("exposure");
  TLeaf*m =tree->GetLeaf("median");
  TLeaf*w =tree->GetLeaf("width");
  double ve,vg;
  int e1g0=0,e1g7=0,e1g14=0,e1g21=0,e10g0=0,e10g7=0,e10g14=0,e10g21=0,e20g0=0,e20g7=0,e20g14=0,e20g21=0;
  Double_t e1g0m[entries],e1g0w[entries],e1g7m[entries],e1g7w[entries],e1g14m[entries],e1g14w[entries],e1g21m[entries],e1g21w[entries],e10g0m[entries],e10g0w[entries],e10g7m[entries],e10g7w[entries],e10g14m[entries],e10g14w[entries],e10g21m[entries],e10g21w[entries],e20g0m[entries],e20g0w[entries],e20g7m[entries],e20g7w[entries],e20g14m[entries],e20g14w[entries],e20g21m[entries],e20g21w[entries];
  for (int i=0; i<entries;i++){
    e->GetBranch()->GetEntry(i);
    ve=e->GetValue();
    if (ve==1){
      g->GetBranch()->GetEntry(i);
      vg=g->GetValue();
      if (vg==0){
	m->GetBranch()->GetEntry(i);
	e1g0m[e1g0]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e1g0w[e1g0]=w->GetValue();
	e1g0++;
      }
      if (vg==7){
	m->GetBranch()->GetEntry(i);
	e1g7m[e1g7]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e1g7w[e1g7]=w->GetValue();
	e1g7++;
      }
      if (vg==14){
	m->GetBranch()->GetEntry(i);
	e1g14m[e1g14]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e1g14w[e1g14]=w->GetValue();
	e1g14++;
      }
      if (vg==21){
	m->GetBranch()->GetEntry(i);
	e1g21m[e1g21]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e1g21w[e1g21]=w->GetValue();
	e1g21++;
      }
    }
    if (ve==10){
      g->GetBranch()->GetEntry(i);
      vg=g->GetValue();
      if (vg==0){
	m->GetBranch()->GetEntry(i);
	e10g0m[e10g0]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e10g0w[e10g0]=w->GetValue();
	e10g0++;
      }
      if (vg==7){
	m->GetBranch()->GetEntry(i);
	e10g7m[e10g7]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e10g7w[e10g7]=w->GetValue();
	e10g7++;
      }
      if (vg==14){
	m->GetBranch()->GetEntry(i);
	e10g14m[e10g14]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e10g14w[e10g14]=w->GetValue();
	e10g14++;
      }
      if (vg==21){
	m->GetBranch()->GetEntry(i);
	e10g21m[e10g21]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e10g21w[e10g21]=w->GetValue();
	e10g21++;
      }
    }
    if (ve==20){
      g->GetBranch()->GetEntry(i);
      vg=g->GetValue();
      if (vg==0){
	m->GetBranch()->GetEntry(i);
	e20g0m[e20g0]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e20g0w[e20g0]=w->GetValue();
	e20g0++;
      }
      if (vg==7){
	m->GetBranch()->GetEntry(i);
	e20g7m[e20g7]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e20g7w[e20g7]=w->GetValue();
	e20g7++;
      }
      if (vg==14){
	m->GetBranch()->GetEntry(i);
	e20g14m[e20g14]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e20g14w[e20g14]=w->GetValue();
	e20g14++;
      }
      if (vg==21){
	m->GetBranch()->GetEntry(i);
	e20g21m[e20g21]=m->GetValue();
	w->GetBranch()->GetEntry(i);
	e20g21w[e20g21]=w->GetValue();
	e20g21++;
      }
    }
  }
  Double_t x[3];
  x[0]=1;
  x[1]=10;
  x[2]=20;
  Double_t y_g0[3], yerr_g0[3], y_g7[3], yerr_g7[3], y_g14[3], yerr_g14[3], y_g21[3], yerr_g21[3];
  y_g0[0]=wa(e1g0m,e1g0w,e1g0);
  yerr_g0[0]=waerr(e1g0m,e1g0w,e1g0);
  y_g0[1]=wa(e10g0m,e10g0w,e10g0);
  yerr_g0[1]=waerr(e10g0m,e10g0w,e10g0);
  y_g0[2]=wa(e20g0m,e20g0w,e20g0);
  yerr_g0[2]=waerr(e20g0m,e20g0w,e20g0);
  gStyle->SetOptFit(0001);
  TGraphErrors* grg0=new TGraphErrors(3,x,y_g0,0,yerr_g0);
  //grg0->SetMarkerColor(3);
  //grg0->SetMarkerStyle(22);
  //grg0->Draw("ALP");

  y_g7[0]=wa(e1g7m,e1g7w,e1g7);
  yerr_g7[0]=waerr(e1g7m,e1g7w,e1g7);
  y_g7[1]=wa(e10g7m,e10g7w,e10g7);
  yerr_g7[1]=waerr(e10g7m,e10g7w,e10g7);
  y_g7[2]=wa(e20g7m,e20g7w,e20g7);
  yerr_g7[2]=waerr(e20g7m,e20g7w,e20g7);
  TGraphErrors* grg7=new TGraphErrors(3,x,y_g7,0,yerr_g7);
  //grg7->SetMarkerStyle(22);
  //grg7->SetMarkerColor(3);
  //grg7->Draw("ALP");

  y_g14[0]=wa(e1g14m,e1g14w,e1g14);
  yerr_g14[0]=waerr(e1g14m,e1g14w,e1g14);
  y_g14[1]=wa(e10g14m,e10g14w,e10g14);
  yerr_g14[1]=waerr(e10g14m,e10g14w,e10g14);
  y_g14[2]=wa(e20g14m,e20g14w,e20g14);
  yerr_g14[2]=waerr(e20g14m,e20g14w,e20g14);
  TGraphErrors* grg14=new TGraphErrors(3,x,y_g14,0,yerr_g14);
  //grg14->SetMarkerStyle(22);
  //grg14->SetMarkerColor(3);
  //grg14->Draw("ALP");

  y_g21[0]=wa(e1g21m,e1g21w,e1g21);
  yerr_g21[0]=waerr(e1g21m,e1g21w,e1g21);
  y_g21[1]=wa(e10g21m,e10g21w,e10g21);
  yerr_g21[1]=waerr(e10g21m,e10g21w,e10g21);
  y_g21[2]=wa(e20g21m,e20g21w,e20g21);
  yerr_g21[2]=waerr(e20g21m,e20g21w,e20g21);
  TGraphErrors* grg21=new TGraphErrors(3,x,y_g21,0,yerr_g21);
  //grg21->SetMarkerStyle(22);
  grg21->SetMarkerColor(kOrange);
  //grg14->Draw("ALP");

  TCanvas * c1=new TCanvas("c","c");
  TMultiGraph* mg=new TMultiGraph();
  mg->Add(grg0,"gain0");
  mg->Add(grg7,"gain7");
  mg->Add(grg14,"gain0");
  mg->Add(grg21,"gain7");
  mg->Draw("AL");
  grg0->SetLineColor(kBlue);
  grg0->Fit("pol1");
  grg7->SetLineColor(kGreen);
  grg7->Fit("pol1");
  grg14->SetLineColor(kRed);
  grg14->Fit("pol1");
  grg21->SetLineColor(kOrange);
  grg21->Fit("pol1");
  mg->SetTitle("exposure v median;exposure;median");
  TLegend *leg=new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry("grg21", "gain=21 (orange)","LE");
  leg->AddEntry("grg14","gain=14 (red)","LE");
  leg->AddEntry("grg7", "gain=7 (green)","LE");
  leg->AddEntry("grg0","gain=0 (blue)","LE");
  leg->Draw();
  std::cout<<"           | "<<"gain=0  | "<<"gain=7  | "<<"gain=14 | "<<"gain=21"<<std::endl;
  std::cout<<"-----------|---------|---------|---------|--------"<<std::endl;
  for (int i=0; i<3; i++){
    std::cout<<"exposure="<<x[i]<<" | "<<y_g0[i]<<" | "<<y_g7[i]<<" | "<<y_g14[i]<<" | "<<y_g21[i]<<std::endl;
  }
  //std::cout<<waerr(e1g0m,e1g0w,e1g0)<<std::endl;
  //std::cout<<wa(e10g0m,e10g0w,e10g0)<<std::endl;
  //std::cout<<waerr(e10g0m,e10g0w,e10g0)<<std::endl;
  //std::cout<<wa(e20g0m,e20g0w,e20g0)<<std::endl;
  //std::cout<<waerr(e20g0m,e20g0w,e20g0)<<std::endl;
}
