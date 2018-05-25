
const Float_t f = TMath::Pi()/180.0;


Double_t fitfunc(Double_t* arg,Double_t* par){
  double el = arg[0];
  double az0 = par[0];
  double el0 = par[1];
  return asin(sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f))/f;
}

Double_t del(Double_t el,Double_t az0,Double_t el0){
  return asin(sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f))/f;
}

double AZ = 90;

Double_t funcphi2(Double_t* arg,Double_t* par){
  double el = arg[0];
  double phi = AZ;
  double phi0 = par[0];
  double el0 = par[1];
  double A = -sin(phi*f)*(cos(el*f)*cos(phi0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(phi*f)*sin(phi0*f)*cos(el0*f);
  double B = +cos(phi*f)*(cos(el*f)*cos(phi0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(phi*f)*sin(phi0*f)*cos(el0*f);
  return atan2(-A,B)/f;
}

Double_t daz(Double_t el,Double_t phi0,Double_t phi1){
  double phi = AZ;
  double el0 = phi1;
  double A = -sin(phi*f)*(cos(el*f)*cos(phi0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(phi*f)*sin(phi0*f)*cos(el0*f);
  double B = +cos(phi*f)*(cos(el*f)*cos(phi0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(phi*f)*sin(phi0*f)*cos(el0*f);
  return atan2(-A,B)/f;
}

void run281f(float _AZ=90){
  AZ = _AZ;
  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  TGraph* g1 = new TGraph();
  TGraph* g2= new TGraph();
  TGraph* g3 = new TGraph();
  TGraph* g4 = new TGraph();
  int k = 0;
  int kk = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Float_t azd = a[1];    
    Float_t eld = a[2];
    if( eld > 10 ){
    Float_t az = a[3];
    Float_t el = a[4];
    Float_t flag = a[5];
    Float_t u = a[0];
    g1->SetPoint(kk++,eld,el);
    if( fabs(a[1]-AZ)<2 ){
      g2->SetPoint(k,eld,az);
      //g3->SetPoint(k,az,el);
      k++;
    }
    }
  }
  g1->Sort();
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.80);
  TString tit("azimuth drive = ");
  tit += Int_t(AZ);
  g1->SetTitle(tit+" deg");
  g1->GetXaxis()->SetTitle("elevation drive (deg)");
  g1->GetYaxis()->SetTitle("elevation CCD (deg)");
  
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(0.80);
  g2->GetXaxis()->SetTitle("elevation drive (deg)");
  g2->GetYaxis()->SetTitle("azimuth CCD (deg)");

  
  g3->SetMarkerStyle(20);
  g3->SetMarkerSize(0.80);
  TCanvas* can = new TCanvas("can","was",0,0,800,600);
  can->Divide(2,2);
  can->cd(1);
  g1->Draw("AP");
  TF1*f = new TF1("f",fitfunc,0,90,2);
  g1->Fit("f","","",20,90);
  double az0 = g1->GetFunction("f")->GetParameter(0);
  double el0 = g1->GetFunction("f")->GetParameter(1);
  std::cout << "ZZ " << az0 << " " << el0 << std::endl;
      //(new TLine(10,10,80,80))->Draw();
  can->cd(2);
  g2->SetMarkerColor(kBlack);
  g2->Draw("AP");
  TF1*ff = new TF1("ff",funcphi2,0,90,2);
  g2->Fit("ff","","",20,90);
  Double_t phi0 = g2->GetFunction("ff")->GetParameter(0);
  Double_t phi1 = g2->GetFunction("ff")->GetParameter(1);
  std::cout << "ZYZ " << phi0 << " " << phi1 << std::endl;
  can->cd(0);
  TString nam("run281_az");
  nam += Int_t(AZ);
  nam += ".png";
  can->SaveAs(nam);
	std::cout<<daz(1,12.5,-1.5)<<std::endl;
  k = 0;
  kk = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Float_t azd = a[1];    
    Float_t eld = a[2];
    if( eld > 10 ){;
    Float_t az = a[3];
    Float_t el = a[4];
    Float_t flag = a[5];
    Float_t u = a[0];
    el=el-del(eld,az0,el0);
    g3->SetPoint(kk++,eld,el);
    if( fabs(a[1]-AZ)<2 ){
      az=az-daz(eld,az0,el0);
      g4->SetPoint(k,eld,az);
      k++;
    }
    } 
  }  
  can->cd(3);
  g3->Sort();
  g3->SetMarkerStyle(20);
  g3->SetMarkerSize(0.80);
  g3->Draw("AP");
  can->cd(4);
  g4->Sort();
  g4->SetMarkerStyle(20);
  g4->SetMarkerSize(0.80);
  g4->Draw("AP");
can->SaveAs(nam);
}


