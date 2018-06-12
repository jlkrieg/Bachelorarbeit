
const Float_t f = TMath::Pi()/180.0;
double AZ = 90;

Double_t fitEl(Double_t* arg,Double_t* par){
  double el = arg[0];
  double az0 = par[0];
  double el0 = par[1];
  double phi = AZ;
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  return asin(Z)/f;
}

Double_t funcEl(Double_t el,Double_t az0,Double_t el0){
  double phi = AZ;
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  return asin(Z)/f;
}

Double_t fitAz(Double_t* arg,Double_t* par){
  double el = arg[0];
  double az = AZ;
  double az0 = par[0];
  double el0 = par[1];
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  return atan2(-Y,X)/f;
}

Double_t funcAz(Double_t el,Double_t az0,Double_t el0){
  double az = AZ;
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  return atan2(-Y,X)/f;
}

void run281fit(float _AZ=90){
  double U =0;
  double DAZ =2;
  double cutEl=10;
  AZ = _AZ;
  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();

  TGraph* g1 = new TGraph();
  TGraph* g2= new TGraph();
  TGraph* g4 = new TGraph();
  TGraph* g0 = new TGraph();
  TGraph* g90 = new TGraph();
  TGraph* g180 = new TGraph();
  TGraph* g270 = new TGraph();
  TMultiGraph* g3 = new TMultiGraph();

  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.80);

  TString tit("azimuth drive = ");
  tit += Int_t(AZ);
  g2->SetTitle(tit+" deg");
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(0.80);

  g4->SetMarkerStyle(20);
  g4->SetMarkerSize(0.80);

  g0->SetMarkerStyle(20);
  g0->SetMarkerSize(0.80);
  g0->SetMarkerColor(kRed);

  g90->SetMarkerStyle(20);
  g90->SetMarkerSize(0.80);
  g90->SetMarkerColor(kBlue);

  g180->SetMarkerStyle(20);
  g180->SetMarkerSize(0.80);
  g180->SetMarkerColor(kGreen);

  g270->SetMarkerStyle(20);
  g270->SetMarkerSize(0.80);
  g270->SetMarkerColor(kOrange);

  TCanvas* can = new TCanvas("can","Plots",0,0,800,600);
  TString nam("run281_az");
  nam += Int_t(AZ);
  nam += ".png";
  can->Divide(2,2);

  int k = 0;
  int kk = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t u = a[0];
    //if( u == U ) continue; //evtl Hysterese
    Float_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Float_t azd = a[1];    
    Float_t eld = a[2];
    if( eld > cutEl ){ //grosse Ausreisser unter 10deg
      Float_t azd = a[1]; 
      Float_t az = a[3];
      Float_t el = a[4];
      Float_t flag = a[5];
      g1->SetPoint(kk++,eld,el);
      if( fabs(a[1]-AZ)<=DAZ ){
	g2->SetPoint(k,eld,az);
	k++;
      }
    }
  }

  g1->Sort();
  can->cd(1);
  g1->Draw("AP");
  g1->GetXaxis()->SetTitle("elevation drive (deg)");
  g1->GetYaxis()->SetTitle("elevation CCD (deg)");
  TF1*f = new TF1("fel",fitEl,0,90,2);
  g1->Fit("fel","","",20,90);
  Double_t az0el = g1->GetFunction("fel")->GetParameter(0);
  Double_t el0el = g1->GetFunction("fel")->GetParameter(1);
  std::cout << "el " << az0el << " " << el0el << std::endl;

  can->cd(2);
  g2->SetMarkerColor(kBlack);
  g2->Draw("AP");
  g2->GetXaxis()->SetTitle("elevation drive (deg)");
  g2->GetYaxis()->SetTitle("azimuth CCD (deg)");
  TF1*ff = new TF1("faz",fitAz,0,90,2);
  g2->Fit("faz","","",20,90);
  Double_t az0az = g2->GetFunction("faz")->GetParameter(0);
  Double_t el0az = g2->GetFunction("faz")->GetParameter(1);
  std::cout << "az " << az0az << " " << el0az << std::endl;

  Double_t del,daz;
  k = 0;
  kk = 0;
  int k0=0;
  int k90=0;
  int k180=0;
  int k270=0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t u = a[0];
    //if( u == U ) continue; //evtl Hysterese
    Float_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Float_t azd = a[1];    
    Float_t eld = a[2];
    if( eld > cutEl ){
      Float_t az = a[3];
      Float_t el = a[4];
      Float_t flag = a[5];
      del=el-funcEl(eld,az0el,el0el);
      if( fabs(a[1])<=DAZ ){
	g0->SetPoint(k0++,eld,del);
      }  
      if( fabs(a[1]-90)<=DAZ ){();
	g90->SetPoint(k90++,eld,del);
      }
      if( fabs(a[1]-180)<=DAZ ){
	g180->SetPoint(k180++,eld,del);
      }
      if( fabs(a[1]+91)<=DAZ ){
	g270->SetPoint(k270++,eld,del);
      }
      if( fabs(a[1]-AZ)<=DAZ ){
	g4->SetPoint(kk++,eld,el);
	daz=az-funcAz(eld,az0az,el0az);
	g4->SetPoint(k,eld,daz);
	k++;
      }
    } 
  }

  can->cd(3);
  g3->Add(g0);
  g3->Add(g90);
  g3->Add(g180);
  g3->Add(g270);
  g3->Draw("AP");
  g3->GetXaxis()->SetTitle("elevation drive (deg)");
  g3->GetYaxis()->SetTitle("elevation CCD (deg)");
  TLegend* leg = new TLegend(0.3,0.1,0.7,0.3);
  leg->AddEntry(g0,"az=0deg");
  leg->AddEntry(g90,"az=90deg");
  leg->AddEntry(g180,"az=180deg");
  leg->AddEntry(g270,"az=-91deg");
  leg->Draw();
  can->cd(4);
  g4->Sort();
  g4->Draw("AP");
  g4->GetXaxis()->SetTitle("elevation drive (deg)");
  g4->GetYaxis()->SetTitle("azimuth CCD (deg)");
can->SaveAs(nam);
}


