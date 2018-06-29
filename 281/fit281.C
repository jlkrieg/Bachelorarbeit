const Float_t f = TMath::Pi()/180.0;
double AZ = 90;

Double_t fitEl(Double_t* x,Double_t* par){
  Double_t el = x[0];
  Double_t az = x[1];
  Double_t az0 = par[0];
  Double_t el0 = par[1];
  Double_t Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  return asin(Z)/f;
}

Double_t funcEl(Double_t el,Double_t az,Double_t az0,Double_t el0){
  //double az = AZ;
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  return asin(Z)/f;
}

Double_t fitAz(Double_t* arg,Double_t* par){
  double el = arg[0];
  double az = arg[1];
  double az0 = par[0];
  double el0 = par[1];
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  return atan2(-Y,X)/f;
}

Double_t funcAz(Double_t el,Double_t az,Double_t az0,Double_t el0){
  //double az = AZ;
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  return atan2(-Y,X)/f; 
}

Double_t chi(Double_t az0, Double_t el0){
  Double_t cutEl=35; //grosse Ausreisser unter 10deg
  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  int k = 0;
  Double_t result=0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t u = a[0];
    //if( u == U ) continue; //evtl Hysterese
    Float_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue; 
    Float_t azd = a[1]; 
    Float_t eld = a[2];
    Float_t az = a[3];
    Float_t el = a[4];
    if( eld > cutEl ){ //grosse Ausreisser unter 10deg
      result+=pow((el-funcEl(eld,azd,az0,el0)),2)+pow((az-funcAz(eld,azd,az0,el0)),2);
    }
  }
  return result;
}

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  result = chi(par[0], par[1]);
}

void fit281(){
  double U =0;
  double DAZ =2;
  double cutEl=35;
  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();

  TCanvas* can = new TCanvas("can","Plots",0,0,800,600);
  TString nam("run281_az");
  nam += Int_t(AZ);
  nam += ".png";
  can->Divide(2,2);

  int k = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t u = a[0];
    //if( u == U ) continue; //evtl Hysterese
    Float_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;  
    Float_t eld = a[2];
    Float_t azd = a[1];
    //if ( azd = AZ ) continue;
    if( eld > cutEl ){ //grosse Ausreisser unter 10deg
      k++;
    }
  }
  const int kk=k;
  //std::cout<<kk<<std::endl;
  Double_t az_vec[kk],el_vec[kk],azd_vec[kk],eld_vec[kk],daz_vec[kk],del_vec[kk],azerr[kk],elerr[kk];
  //Einlesen der Daten
  k=0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t u = a[0];
    //if( u == U ) continue; //evtl Hysterese
    Float_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Float_t azd = a[1];    
    Float_t eld = a[2];
    //if ( azd = AZ ) continue;
    if( eld > cutEl ){ //grosse Ausreisser unter 10deg
      Float_t az = a[3];
      Float_t el = a[4];
      Float_t flag = a[5];
      azd_vec[k]=azd;
      eld_vec[k]=eld;
      az_vec[k]=az;
      el_vec[k]=el;
      k++;
    }
  }
  //minimize
  //function to minimize chi(el_vec,az_vec,eld_vec,azd_vec,kk,n,par)
  TFitter* minimizer = new TFitter(2);
  { 
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
  }
  minimizer->SetFCN(minuitFunction);
  // Define the parameters
         //   arg1 – parameter number
  //   arg2 – parameter name
  //   arg3 – first guess at parameter value
  //   arg4 – estimated distance to minimum
  //   arg5, arg6 – ignore for now
  minimizer->SetParameter(0,"az0",12,1,0,0);
  minimizer->SetParameter(1,"el0",-1.4,0.5,0,0);
  // Run the simplex minimizer to get close to the minimum
  minimizer->ExecuteCommand("SIMPLEX",0,0);
  // Run the migrad minimizer (an extended Powell's method) to improve the
  // fit.
  minimizer->ExecuteCommand("MIGRAD",0,0);
  // Get the best fit values
  double bestX = minimizer->GetParameter(0);
  double bestY = minimizer->GetParameter(1);
  // Get the function value at the best fit.
  double minimum = chi(bestX, bestY);
  std::cout<<"az0 = "<<bestX<<std::endl;
  std::cout<<"el0 = "<<bestY<<std::endl;
  std::cout<<minimum<<std::endl;
  //testparams
  //bestX=12;
  //bestY=-1.4;
  //calculate differences
  for(int i=0; i<kk; i++){
    del_vec[i] = funcEl(eld_vec[i],azd_vec[i],bestX,bestY)-el_vec[i];
  }
  for(int i=0; i<kk; i++){
    daz_vec[i] = funcAz(eld_vec[i],azd_vec[i],bestX,bestY)-az_vec[i];
  }
  //chitest
  Double_t chisq=0;
  for(int i=0; i<kk; i++){
    chisq+=pow(del_vec[i],2)+pow(daz_vec[i],2);
  }
  std::cout<<"chisqtest = "<<chisq<<std::endl;
  //plot
  TGraph* g_delel=new TGraph(kk,el_vec,del_vec);
  can->cd(1);
  g_delel->SetMarkerStyle(20);
  g_delel->SetMarkerSize(0.80);
  g_delel->Draw("AP");
  TGraph* g_delaz=new TGraph(kk,az_vec,del_vec);
  can->cd(2);
  g_delaz->SetMarkerStyle(20);
  g_delaz->SetMarkerSize(0.80);
  g_delaz->Draw("AP");
  TGraph* g_dazel=new TGraph(kk,el_vec,daz_vec);
  can->cd(3);
  g_dazel->SetMarkerStyle(20);
  g_dazel->SetMarkerSize(0.80);
  g_dazel->Draw("AP");
  TGraph* g_dazaz=new TGraph(kk,az_vec,daz_vec);
  can->cd(4);
  g_dazaz->SetMarkerStyle(20);
  g_dazaz->SetMarkerSize(0.80);
  g_dazaz->Draw("AP");
}


