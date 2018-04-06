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

void Korrelation2(){
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
  Double_t x[4];
  x[0]=0;
  x[1]=7;
  x[2]=14;
  x[3]=21;
  Double_t y_e1[4], yerr_e1[4], y_e10[4], yerr_e10[4], y_e20[4], yerr_e20[4];
  y_e1[0]=wa(e1g0m,e1g0w,e1g0);
  yerr_e1[0]=waerr(e1g0m,e1g0w,e1g0);
  y_e1[1]=wa(e1g7m,e1g7w,e1g7);
  yerr_e1[1]=waerr(e1g7m,e1g7w,e1g7);
  y_e1[2]=wa(e1g14m,e1g14w,e1g14);
  yerr_e1[2]=waerr(e1g14m,e1g14w,e1g14);
  y_e1[3]=wa(e1g21m,e1g21w,e1g21);
  yerr_e1[3]=waerr(e1g21m,e1g21w,e1g21);
  gStyle->SetOptFit(0001);
  TGraphErrors* gre1=new TGraphErrors(4,x,y_e1,0,yerr_e1);
  //grg0->SetMarkerColor(3);
  //grg0->SetMarkerStyle(22);
  //grg0->Draw("ALP");

  y_e10[0]=wa(e10g0m,e10g0w,e10g0);
  yerr_e10[0]=waerr(e10g0m,e10g0w,e10g0);
  y_e10[1]=wa(e10g7m,e10g7w,e10g7);
  yerr_e10[1]=waerr(e10g7m,e10g7w,e10g7);
  y_e10[2]=wa(e10g14m,e10g14w,e10g14);
  yerr_e10[2]=waerr(e10g14m,e10g14w,e10g14);
  y_e10[3]=wa(e10g21m,e10g21w,e10g21);
  yerr_e10[3]=waerr(e10g21m,e10g21w,e10g21);
  gStyle->SetOptFit(0001);
  TGraphErrors* gre10=new TGraphErrors(4,x,y_e10,0,yerr_e10);
  //grg0->SetMarkerColor(3);
  //grg0->SetMarkerStyle(22);
  //grg0->Draw("ALP")

  y_e20[0]=wa(e20g0m,e20g0w,e20g0);
  yerr_e20[0]=waerr(e20g0m,e20g0w,e20g0);
  y_e20[1]=wa(e20g7m,e20g7w,e20g7);
  yerr_e20[1]=waerr(e20g7m,e20g7w,e20g7);
  y_e20[2]=wa(e20g14m,e20g14w,e20g14);
  yerr_e20[2]=waerr(e20g14m,e20g14w,e20g14);
  y_e20[3]=wa(e20g21m,e20g21w,e20g21);
  yerr_e20[3]=waerr(e20g21m,e20g21w,e20g21);
  gStyle->SetOptFit(0001);
  TGraphErrors* gre20=new TGraphErrors(4,x,y_e20,0,yerr_e20);
  //grg0->SetMarkerColor(3);
  //grg0->SetMarkerStyle(22);
  //grg0->Draw("ALP")

  TCanvas * c1=new TCanvas("c","c");
  TMultiGraph* mg=new TMultiGraph();
  mg->Add(gre1,"exosure1");
  mg->Add(gre10,"exposure10");
  mg->Add(gre20,"exposure20");
  mg->Draw("AL");
  gre1->SetLineColor(kBlue);
  gre1->SetMarkerSize(5);
  gre1->Fit("pol1");
  gre10->SetLineColor(kGreen);
  gre10->Fit("pol1");
  gre20->SetLineColor(kRed);
  gre20->Fit("pol1");
  mg->SetTitle("gain v median;gain;median");
  TLegend *leg=new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry("gre20", "exposure=20 (red)","LE");
  leg->AddEntry("gre10","exposure=10 (green)","LE");
  leg->AddEntry("gre1", "exposure=1 (blue)","LE");
  leg->Draw();


  std::cout<<"       |"<<"exposure=1| "<<"exposure=10|"<<"exposure=20"<<std::endl;
  std::cout<<"-------|----------|------------|-----------"<<std::endl;
  for (int i=0; i<4; i++){
    std::cout<<"gain="<<x[i]<<" |  "<<y_e1[i]<<"  |  "<<y_e10[i]<<"  |  "<<y_e20[i]<<std::endl;
  }
}
