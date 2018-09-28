//global constants
const Float_t f = TMath::Pi()/180.0;
Double_t AZ = 0;
Double_t DAZ = 2;
Double_t cutEl = 35;

//functions for el and az
Double_t funcEl(Double_t el,Double_t az,Double_t az0,Double_t el0,Double_t phi, Double_t psi){
  az = az*f;
  el = el*f;
  az0 = az0*f;
  el0 = el0*f;
  phi = phi*f;
  psi = psi*f;
  //Versuch mit 4 Parametern
  double X = +cos(az)*(cos(el)*cos(az0)*cos(el0)-sin(el)*sin(el0))+sin(az)*sin(az0)*cos(el0);
  double Y = sin(az)*(cos(el)*cos(az0)*cos(el0)-sin(el)*sin(el0))-cos(az)*sin(az0)*cos(el0);
  double Z = sin(el)*cos(az0)*cos(el0)+cos(el)*sin(el0);
  //return asin(Z+psi*X+phi*Y)/f;  
  return (asin(Z)+phi*sin(az+psi))/f;
}

Double_t funcAz(Double_t el,Double_t az,Double_t az0,Double_t el0,Double_t phi, Double_t psi){
  az = az*f;
  el = el*f;
  az0 = az0*f;
  el0 = el0*f;
  phi = phi*f;
  psi = psi*f;
  double X = cos(az)*(cos(el)*cos(az0)*cos(el0)-sin(el)*sin(el0))+sin(az)*sin(az0)*cos(el0);
  double Y = sin(az)*(cos(el)*cos(az0)*cos(el0)-sin(el)*sin(el0))-cos(az)*sin(az0)*cos(el0);
  double Z = sin(el)*cos(az0)*cos(el0)+cos(el)*sin(el0);
  //return atan2(Y-phi*psi*X-phi*Z,X+psi*X+phi*Y)/f;
  return atan2(Y,X)/f;
}

//funtion to minimize (chisq)
Double_t chi(Double_t az0, Double_t el0,Double_t phi, Double_t psi){
  Double_t result=0;
  TFile *file = TFile::Open("ntuple2nt_v12.root");
  TNtuple*nt = (TNtuple*)file->Get("run341_ccd3_tpoint_0_00_nt_ntuple");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t ps = a[5];
    if( fabs(ps-11.03) > 0.05 ) continue; 
    Double_t azd = a[0];
    Double_t eld = a[1];
    Double_t az = a[2];
    Double_t el = a[3];
    Double_t del = funcEl(eld,azd,az0,el0,phi,psi)-el;
    if (del < -180){
      del+=360;
    }
    if (del > 180){
      del-=360;
    }
    Double_t daz = funcAz(eld,azd,az0,el0,phi,psi)-az;
    if (daz < -180){
      daz+=360;
    }
    if (daz > 180){
      daz-=360;
    }
    result+=pow(del,2)+pow(daz,2);
  }
  file->Close();
  return result;
}

//for minimization
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  result = chi(par[0], par[1],par[2],par[3]);
}

void fitD2C4(){
  TFile *file= TFile::Open("ntuple2nt_v12.root");
  TNtuple*nt = (TNtuple*)file->Get("run341_ccd3_tpoint_0_00_nt_ntuple");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  int k = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t ps = a[5];
    if( fabs(ps-11.03) > 0.05 ) continue;
    k++;
  }
  const int kk=k;
  Double_t az_vec[kk],el_vec[kk],azd_vec[kk],eld_vec[kk],daz_vec[kk],del_vec[kk],azerr[kk],elerr[kk];
  
  //Einlesen der Daten
  k=0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t ps = a[5];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Double_t azd = a[0];
    Double_t eld = a[1];
    Double_t az = a[2];
    if (az < -180){
      az+=360;
    }
    Double_t el = a[3];
    azd_vec[k]=azd;
    eld_vec[k]=eld;
    az_vec[k]=az;
    el_vec[k]=el;
    k++;
  }

  //minimize
  TFitter* minimizer = new TFitter(4);
  { 
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
  }
  minimizer->SetFCN(minuitFunction);
  minimizer->SetParameter(0,"az0",0,10,0,0);
  minimizer->SetParameter(1,"el0",0,5,0,0);
  minimizer->SetParameter(2,"phi",0,0.1,0,0);
  minimizer->SetParameter(3,"psi",0,0.1,0,0);
  minimizer->ExecuteCommand("SIMPLEX",0,0);
  minimizer->ExecuteCommand("MIGRAD",0,0);
  double az0 = minimizer->GetParameter(0);
  double el0 = minimizer->GetParameter(1);
  double phi = minimizer->GetParameter(2);
  double psi = minimizer->GetParameter(3);
  double minimum = chi(az0,el0,phi,psi);
  //double az0 = -12.048;
  //double el0 = -1.18437;
  //double phi = 0.125;
  //double psi = -0.07;
  //phi = -0.15;
  //psi = 0;
  //double minimum = chi(az0,el0,phi,psi);
  std::cout<<"az0 = "<<az0<<std::endl;
  std::cout<<"el0 = "<<el0<<std::endl;
  std::cout<<"phi = "<<phi<<std::endl;
  std::cout<<"psi = "<<psi<<std::endl;
  std::cout<<"Minimum chi^2 = "<<minimum<<std::endl;

  //calculate differences
  for(int i=0; i<kk; i++){
    del_vec[i] = funcEl(eld_vec[i],azd_vec[i],az0,el0,phi,psi)-el_vec[i];//,el0,phi,psi)-el_vec[i];
    if (del_vec[i] < -180){
      del_vec[i]+=360;
    }
    if (del_vec[i] > 180){
      del_vec[i]-=360;
    }
  }
				
  for(int i=0; i<kk; i++){
    daz_vec[i] = funcAz(eld_vec[i],azd_vec[i],az0,el0,phi,psi)-az_vec[i];
    //daz_vec[i] = azd_vec[i]-az_vec[i];
    if (daz_vec[i] < -180){
      daz_vec[i]+=360;
    }
    if (daz_vec[i] > 180){
      daz_vec[i]-=360;
    }
  }

  Double_t chisq=0;
  for(int i=0; i<kk; i++){
    chisq+=pow(del_vec[i],2);
  }
  std::cout<<"chisqtest = "<<chisq<<std::endl;

  //plot
  TCanvas* can = new TCanvas("plots","Plots",0,0,800,600);
  TString nam("run341D2C.png");
  TString tit1("fit drive to CCD");
  TString tit2("el0 = ");
  //tit2 += Int_t(AZ);
  tit2 += el0;
  tit2 += ", az0 = ";
  tit2 += az0;
  can->Divide(2,2);
  TGraph* g_delel=new TGraph(kk,eld_vec,del_vec);
  can->cd(1);
  g_delel->SetMarkerStyle(20);
  g_delel->SetMarkerSize(0.80);
  g_delel->GetXaxis()->SetTitle("elevation drive (deg)");
  g_delel->GetYaxis()->SetTitle("#Delta elevation CCD (deg)");
  g_delel->SetTitle(tit1);
  g_delel->Draw("AP");
  TGraph* g_delaz=new TGraph(kk,azd_vec,del_vec);
  can->cd(2);
  g_delaz->SetMarkerStyle(20);
  g_delaz->SetMarkerSize(0.80);
  g_delaz->GetXaxis()->SetTitle("azimuth drive (deg)");
  g_delaz->GetYaxis()->SetTitle("#Delta elevation CCD (deg)");
  g_delaz->SetTitle(tit2);
  g_delaz->Draw("AP");
  TGraph* g_dazel=new TGraph(kk,eld_vec,daz_vec);
  can->cd(3);
  g_dazel->SetMarkerStyle(20);
  g_dazel->SetMarkerSize(0.80);
  g_dazel->GetXaxis()->SetTitle("elevation drive (deg)");
  g_dazel->GetYaxis()->SetTitle("#Delta azimuth CCD (deg)");
  g_dazel->Draw("AP");
  TGraph* g_dazaz=new TGraph(kk,azd_vec,daz_vec);
  can->cd(4);
  g_dazaz->SetMarkerStyle(20);
  g_dazaz->SetMarkerSize(0.80);
  g_dazaz->GetXaxis()->SetTitle("azimuth drive (deg)");
  g_dazaz->GetYaxis()->SetTitle("#Delta azimuth CCD (deg)");
  g_dazaz->Draw("AP");
  can->SaveAs(nam);
  file->Close();
  }


