//global constants
const Float_t f = TMath::Pi()/180.0;
Double_t AZ =0;
Double_t U =0;
Double_t DAZ =2;
Double_t cutEl=35;

//functions for el and az
Double_t funcEl(Double_t el,Double_t az,Double_t az0,Double_t el0){
  double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  return asin(Z)/f;
  // az=az*f;
  // el=el*f;
  // az0=az0*f;
  // el0=el0*f;
  // Double_t r = pow(cos(az0)*cos(el0),2) + pow(cos(el0)*sin(az0),2) + pow(cos(az0)*sin(el0),2) + pow(sin(az0)*sin(el0),2);
  // Double_t z = cos(az0)*cos(el0)*sin(el)-cos(az0)*cos(el)*sin(el0);
  // return asin(z/r)/f;
}

Double_t funcAz(Double_t el,Double_t az,Double_t az0,Double_t el0){
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  return atan2(-Y,X)/f;
  // az=az*f;
  // el=el*f;
  // az0=az0*f;
  // el0=el0*f;
  // Double_t x=cos(az)*cos(az0)*cos(el)*cos(el0)+cos(az)*cos(az0)*sin(el)*sin(el0)+sin(az)*(pow(cos(el0),2)*sin(az0)+sin(az0)*pow(sin(el0),2));
  // Double_t y=cos(az0)*cos(el)*cos(el0)*sin(az)+cos(az0)*sin(az)*sin(el)*sin(el0)+cos(az)*(pow(cos(el0),2)*sin(az0)+sin(az0)*pow(sin(el0),2));
  // return atan2(-y,x)/f;
}

//funtion to minimize (chisq)
Double_t chi(Double_t az0, Double_t el0){
  Double_t result=0;
  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t u = a[0];
    Double_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue; 
    Double_t azd = a[1];
    //if (azd == AZ){
      Double_t eld = a[2];
      Double_t az = a[3];
      Double_t el = a[4];
      if( eld > cutEl ){
	Double_t del = funcEl(el,az,az0,el0)-eld;
	if (del < -180){
	  del+=360;
	}
	Double_t daz = funcAz(el,az,az0,el0)-eld;
	if (daz < -180){
	  daz+=360;
	}
	result+=pow((del-eld),2);//+pow((daz-azd),2);
	//result+=pow((el-funcEl(eld,azd,az0,el0)),2)+pow((az-funcAz(eld,azd,az0,el0)),2); //eld abh
      }
      // }
  }
  return result;
}

//for minimization
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  result = chi(par[0], par[1]);
}

void fit281(){

  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  int k = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;  
    Double_t azd = a[1];
    //if ( azd == AZ ){
      Double_t eld = a[2];
      if( eld > cutEl ){
	k++;
      }
      //}
  }
  const int kk=k;
  Double_t az_vec[kk],el_vec[kk],azd_vec[kk],eld_vec[kk],daz_vec[kk],del_vec[kk],azerr[kk],elerr[kk];
  
  //Einlesen der Daten
  k=0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t u = a[0];
    Double_t ps = a[6];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Double_t azd = a[1];
    //if ( azd == AZ ){
      Double_t eld = a[2];
      if( eld > cutEl ){
	Double_t az = a[3];
	if (az < -180){
	  az+=360;
	}
	Double_t el = a[4];
	Float_t flag = a[5];
	azd_vec[k]=azd;
	eld_vec[k]=eld;
	az_vec[k]=az;
	el_vec[k]=el;
	k++;
      }
      //}
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
  minimizer->SetParameter(0,"az0",0,10,0,0);
  minimizer->SetParameter(1,"el0",0,5,0,0);
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
  std::cout<<"Minimum chi^2 = "<<minimum<<std::endl;

  //testparams
  // for drive dependency az0 = 12.1029 el0 = -1.23133
  //bestX=12.1029;
  //bestY=-1.23133;

  //calculate differences
  for(int i=0; i<kk; i++){
    //del_vec[i] = funcEl(el_vec[i],az_vec[i],bestX,bestY)-eld_vec[i];
    //del_vec[i] = funcEl(eld_vec[i],azd_vec[i],bestX,bestY)-el_vec[i]; //eld abh
    del_vec[i] = el_vec[i]-eld_vec[i];
    if (del_vec[i] < -180){
      del_vec[i]+=360;
    }
  }				
  for(int i=0; i<kk; i++){
    //daz_vec[i] = funcAz(el_vec[i],az_vec[i],bestX,bestY)-azd_vec[i];
    //daz_vec[i] = funcAz(eld_vec[i],azd_vec[i],bestX,bestY)-az_vec[i]; //eld abh
    daz_vec[i] = az_vec[i]-azd_vec[i];
    if (daz_vec[i] < -180){
      daz_vec[i]+=360;
    }
  }

  //chitest
  Double_t chisq=0;
  for(int i=0; i<kk; i++){
    chisq+=pow(del_vec[i],2)+pow(daz_vec[i],2);
  }
  std::cout<<"chisqtest = "<<chisq<<std::endl;

  //plot
  TCanvas* can = new TCanvas("plots","Plots",0,0,800,600);
  TString nam("run281_az");
  nam += Int_t(AZ);
  nam += ".png";
  can->Divide(2,2);
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


