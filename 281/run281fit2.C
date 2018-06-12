
const Float_t f = TMath::Pi()/180.0;
double AZ = 90;

Double_t fitEl(Double_t* arg,Double_t* par){
  double el = arg[0];
  double az0 = par[0];
  double el0 = par[1];
  double alpha = par[2];
  double beta = par[3];
  double az = AZ;
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  return asin(-cos(alpha)*sin(beta)*X-sin(alpha)*Y+cos(alpha)*cos(beta)*Z)/f;
}

Double_t funcEl(Double_t el,Double_t az0,Double_t el0,Double_t alpha, Double_t beta){
  double az = AZ;
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  return asin(-cos(alpha)*sin(beta)*X-sin(alpha)*Y+cos(alpha)*cos(beta)*Z)/f;
}

Double_t fitAz(Double_t* arg,Double_t* par){
  double el = arg[0];
  double az = AZ;
  double az0 = par[0];
  double el0 = par[1];
  double alpha = par[2];
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  double temp = cos(alpha)*Y+sin(alpha)*Z;
  return atan2(-temp,X)/f;
}

Double_t funcAz(Double_t el,Double_t az0,Double_t el0,Double_t alpha){
  double az = AZ;
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  double temp = cos(alpha)*Y+sin(alpha)*Z;
  return atan2(-temp,X)/f;
}

void run281fit(float _AZ=90){
  double U =0;
  double DAZ =2;
  double cutEl=35;
  AZ = _AZ;
  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();

  TGraph* g1 = new TGraph();
  TGraph* g2= new TGraph();
  TGraph* g3 = new TGraph();
  TGraph* g4 = new TGraph();

  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.80);

  TString tit("azimuth drive = ");
  tit += Int_t(AZ);
  g2->SetTitle(tit+" deg");
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(0.80);

  g3->SetMarkerStyle(20);
  g3->SetMarkerSize(0.80);

  g4->SetMarkerStyle(20);
  g4->SetMarkerSize(0.80);



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
      if( fabs(a[1]-AZ)<=DAZ ){
	g1->SetPoint(kk++,eld,el);
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
  TF1*f = new TF1("fel",fitEl,0,90,4);
  g1->Fit("fel","","",20,90);
  Double_t az0el = g1->GetFunction("fel")->GetParameter(0);
  Double_t el0el = g1->GetFunction("fel")->GetParameter(1);
  Double_t elAlpha = g1->GetFunction("fel")->GetParameter(2);
  Double_t elBeta = g1->GetFunction("fel")->GetParameter(3);
  std::cout << "el " << az0el << " " << el0el << std::endl;

  can->cd(2);
  g2->SetMarkerColor(kBlack);
  g2->Draw("AP");
  g2->GetXaxis()->SetTitle("elevation drive (deg)");
  g2->GetYaxis()->SetTitle("azimuth CCD (deg)");
  TF1*ff = new TF1("faz",fitAz,0,90,3);
  g2->Fit("faz","","",20,90);
  Double_t az0az = g2->GetFunction("faz")->GetParameter(0);
  Double_t el0az = g2->GetFunction("faz")->GetParameter(1);
  Double_t azAlpha = g2->GetFunction("faz")->GetParameter(2);
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
      if( fabs(a[1]-AZ)<=DAZ ){
	daz=az-funcAz(eld,az0az,el0az,azAlpha);
	g4->SetPoint(k,eld,daz);
	del=el-funcEl(eld,az0el,el0el,elAlpha,elBeta);
	g3->SetPoint(k,eld,del);
	k++;
      }
    } 
  }

  can->cd(3);
  g3->Draw("AP");
  g3->GetXaxis()->SetTitle("elevation drive (deg)");
  g3->GetYaxis()->SetTitle("#Delta elevation CCD (deg)");
  can->cd(4);
  g4->Sort();
  g4->Draw("AP");
  g4->GetXaxis()->SetTitle("elevation drive (deg)");
  g4->GetYaxis()->SetTitle("#Delta azimuth CCD (deg)");
can->SaveAs(nam);
}


