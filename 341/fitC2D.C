//global constants
const Float_t f = TMath::Pi()/180.0;
Double_t AZ =0;
Double_t DAZ =2;
Double_t cutEl=90;

//functions for el and az
Double_t funcEl(Double_t el,Double_t az,Double_t az0,Double_t el0){
  az = az*f;
  el = el*f;
  az0 = az0*f;
  el0 = el0*f;
  double temp1 = cos(az0)*cos(el0)-sqrt(pow(cos(az0)*cos(el0),2)+pow(sin(el0),2)-pow(sin(el),2));
  double temp2 = sin(el0)+sin(el);
  return 2*atan(temp1/temp2)/f;
}

Double_t funcAz(Double_t el,Double_t az,Double_t az0,Double_t el0){
  az = az*f;
  el = el*f;
  az0 = az0*f;
  el0 = el0*f;
  double temp1 = cos(az0)*cos(el0)-sqrt(pow(cos(az0)*cos(el0),2)+pow(sin(el0),2)-pow(sin(el),2));
  double temp2 = sin(el0)+sin(el);
  el = 2*atan(temp1/temp2);
  double X = cos(az)*(cos(el)*cos(az0)*cos(el0)-sin(el)*sin(el0))+sin(az)*sin(az0)*cos(el0);
  double Y = -sin(az)*(cos(el)*cos(az0)*cos(el0)-sin(el)*sin(el0))+cos(az)*sin(az0)*cos(el0);
  return atan2(-Y,X)/f;
}

//funtion to minimize (chisq)
Double_t chi(Double_t az0, Double_t el0){
  Double_t result=0;
  Double_t temp3=0;
  Double_t temp4=0; 
  TFile *file= TFile::Open("ntuple2nt_v12.root");
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
    if( eld < cutEl ){
      Double_t del = funcEl(el,az,az0,el0);
      if (del < -180){
	del+=360;
      }
      if (del > 180){
	del-=360;
      }
      //std::cout<<del;
      Double_t daz = funcAz(el,az,az0,el0);
      if (daz < -180){
	daz+=360;
      }
      if (daz > 180){
	daz-=360;
      }
      temp3=sin(eld*f)*sin(del*f)+cos(eld*f)*cos(del*f)*cos((azd-daz)*f);
      if (fabs(temp3)>1){
	temp3=0;
      }
      else{
	temp3=acos(temp3);
      }
      result+=pow(temp3,2);
    }
  }
  file->Close();
  return result;
}

//for minimization
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  result = chi(par[0], par[1]);
}

void fitC2D(){
  TFile *file= TFile::Open("ntuple2nt_v12.root");
  TNtuple*nt = (TNtuple*)file->Get("run341_ccd3_tpoint_0_00_nt_ntuple");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  int k = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t ps = a[5];
    if( fabs(ps-11.03) > 0.05 ) continue;  
    Double_t azd = a[0];
    Double_t eld = a[1];
    if( eld < cutEl ){
      k++;
    }
  }
  const int kk=k;
  Double_t  az_vec[kk], el_vec[kk],azd_vec[kk],eld_vec[kk],daz_vec[kk],del_vec[kk],azerr[kk],elerr[kk];
  
  //Einlesen der Daten
  k=0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t ps = a[5];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Double_t azd = a[0];
    Double_t eld = a[1];
    if( eld < cutEl ){
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
  }
  //test
  std::cout<<sizeof(eld_vec)/sizeof(eld_vec[0])<<"vs"<<kk<<std::endl;

  //minimize
  TFitter* minimizer = new TFitter(2);
  { 
   double p1 = -1;
   minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
  }
  minimizer->SetFCN(minuitFunction);
  minimizer->SetParameter(0,"az0",12,1,0,0);//0,10
  minimizer->SetParameter(1,"el0",0,1,0,0);//0,5
  minimizer->ExecuteCommand("SIMPLEX",0,0);
  minimizer->ExecuteCommand("MIGRAD",0,0);
  double az0 = minimizer->GetParameter(0);
  double el0 = minimizer->GetParameter(1);
  double minimum = chi(az0,el0);
  std::cout<<"az0 = "<<az0<<std::endl;
  std::cout<<"el0 = "<<el0<<std::endl;
  std::cout<<"Minimum chi^2 = "<<minimum<<std::endl;

  //Matrizen berechnen
  const int npar = 2;
  double matrix[npar][npar];
  gMinuit->mnemat(&matrix[0][0],npar);
  TMatrixD* fCovar = new TMatrixD(npar,npar,&matrix[0][0]);
  std::cout << "Error matrix" << std::endl;
  fCovar->Print();
  double sigma[npar];
  for(int i=0;i<npar;i++){
    sigma[i]=sqrt((*fCovar)[i][i]);
  }
  for(int i=0;i<npar;i++){
    for(int k=0;k<npar;k++){
      double s = sigma[i]*sigma[k];
      (*fCovar)[i][k] = (*fCovar)[i][k]/s;
    }
  }
  std::cout << "Correlation matrix" << std::endl;
  fCovar->Print();

  //calculate differences
  for(int i=0; i<kk; i++){
    del_vec[i] = eld_vec[i]-funcEl(el_vec[i],az_vec[i],az0,el0);
    //del_vec[i] = eld_vec[i]-el_vec[i];
    if (del_vec[i] < -180){
      del_vec[i]+=360;
    }
    if (del_vec[i] > 180){
      del_vec[i]-=360;
    }
  }				
  for(int i=0; i<kk; i++){
    daz_vec[i] = azd_vec[i]-funcAz(el_vec[i],az_vec[i],az0,el0);
    //daz_vec[i] = azd_vec[i]-az_vec[i];
    if (daz_vec[i] < -180){
      daz_vec[i]+=360;
    }
    if (daz_vec[i] > 180){
      daz_vec[i]-=360;
    }
  }

  //chitest
  Double_t chisq=0;
  for(int i=0; i<kk; i++){
    chisq+=pow(del_vec[i],2)+pow(daz_vec[i],2);
  }
  std::cout<<"errorbars = "<<sqrt(chisq/(N-2))<<std::endl;

  //plot
  TCanvas* can = new TCanvas("plotsC","PlotsC",0,0,800,600);
  can->Divide(2,2);
  TString nam("run341C2D.png");
  TString tit1("fit CCD to drive");
  TString tit2("el0 = ");
  tit2 += el0;
  tit2 += ", az0 = ";
  tit2 += az0;
  TGraph* g_delel=new TGraph(kk,el_vec,del_vec);
  can->cd(1);
  g_delel->SetMarkerStyle(20);
  g_delel->SetMarkerSize(0.80);
  g_delel->GetXaxis()->SetTitle("elevation center (deg)");
  g_delel->GetYaxis()->SetTitle("#Delta elevation CCD (deg)");
  g_delel->SetTitle(tit1);
  g_delel->Draw("AP");
  TGraph* g_delaz=new TGraph(kk,az_vec,del_vec);
  can->cd(2);
  g_delaz->SetMarkerStyle(20);
  g_delaz->SetMarkerSize(0.80);
  g_delaz->GetXaxis()->SetTitle("azimuth center (deg)");
  g_delaz->GetYaxis()->SetTitle("#Delta elevation CCD (deg)");
  g_delaz->SetTitle(tit2);
  g_delaz->Draw("AP");
  TGraph* g_dazel=new TGraph(kk,el_vec,daz_vec);
  can->cd(3);
  g_dazel->SetMarkerStyle(20);
  g_dazel->SetMarkerSize(0.80);
  g_dazel->GetXaxis()->SetTitle("elevation center (deg)");
  g_dazel->GetYaxis()->SetTitle("#Delta azimuth CCD (deg)");
  g_dazel->Draw("AP");
  TGraph* g_dazaz=new TGraph(kk,az_vec,daz_vec);
  can->cd(4);
  g_dazaz->SetMarkerStyle(20);
  g_dazaz->SetMarkerSize(0.80);
  g_dazaz->GetXaxis()->SetTitle("azimuth center (deg)");
  g_dazaz->GetYaxis()->SetTitle("#Delta azimuth CCD (deg)");
  g_dazaz->Draw("AP");
  can->SaveAs(nam);
  file->Close();
  }

