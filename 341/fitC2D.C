//global constants
const Float_t f = TMath::Pi()/180.0;
Double_t AZ =0;
Double_t DAZ =2;
Double_t cutEl=90;

//functions for el and az
Double_t funcEl(Double_t el,Double_t az,Double_t az0,Double_t el0){
  //double temp1 = el;
  //double temp2 = az;
  //el = el0;
  //az = az0;
  //el0 = temp1;
  //az0 = temp2;
  //double Z = sin(el*f)*cos(az0*f)*cos(el0*f)+cos(el*f)*sin(el0*f);
  double Z = sin(el*f)*az0+cos(el*f)*el0;  
return asin(Z)/f;
}

Double_t funcAz(Double_t el,Double_t az,Double_t az0,Double_t el0){
  double X = +cos(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-sin(az*f)*sin(az0*f)*cos(el0*f);
  double Y = -sin(az*f)*(cos(el*f)*cos(az0*f)*cos(el0*f)-sin(el*f)*sin(el0*f))-cos(az*f)*sin(az0*f)*cos(el0*f);
  return atan2(-Y,X)/f;
}

//funtion to minimize (chisq)
Double_t chi(Double_t az0, Double_t el0){
  Double_t result=0;
  TFile *file= TFile::Open("ntuple2nt_v12.root");
  TNtuple*nt = (TNtuple*)file->Get("run341_ccd3_tpoint_0_00_nt_ntuple");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Double_t ps = a[5];
    if( fabs(ps-11.03) > 0.05 ) continue; 
    Double_t azd = a[0];
    //if (azd == AZ){
      Double_t eld = a[1];
      Double_t az = a[2];
      Double_t el = a[3];
      if( eld < cutEl ){
	Double_t del = funcEl(el,az,az0,el0)-eld;
	if (del < -180){
	  del+=360;
	}
	if (del > 180){
	  del-=360;
	}
	Double_t daz = funcAz(el,az,az0,el0)-azd;
	if (daz < -180){
	  daz+=360;
	}
	if (daz > 180){
	  daz-=360;
	}
	result+=pow(del,2);//+pow(daz,2);
	}
      // }
  }
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
    //if ( azd == AZ ){
      Double_t eld = a[1];
      if( eld < cutEl ){
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
    Double_t ps = a[5];
    if( fabs(ps-11.03) > 0.05 ) continue;
    Double_t azd = a[0];
    //if ( azd == AZ ){
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
      //}
  }

  //minimize
  //function to minimize chi(el_vec,az_vec,eld_vec,azd_vec,kk,n,par)
  //TFitter* minimizer = new TFitter(2);
  //{ 
  //  double p1 = -1;
  //  minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
  //}
  //minimizer->SetFCN(minuitFunction);
  // Define the parameters
  //   arg1 – parameter number
  //   arg2 – parameter name
  //   arg3 – first guess at parameter value
  //   arg4 – estimated distance to minimum
  //   arg5, arg6 – ignore for now
  //minimizer->SetParameter(0,"az0",1,0.1,0,0);//0,10
  //minimizer->SetParameter(1,"el0",0,0.1,0,0);//0,5
  // Run the simplex minimizer to get close to the minimum
  //minimizer->ExecuteCommand("SIMPLEX",0,0);
  // Run the migrad minimizer (an extended Powell's method) to improve the fit.
  //minimizer->ExecuteCommand("MIGRAD",0,0);
  // Get the best fit values
  //double el0 = minimizer->GetParameter(0);
  //double az0 = minimizer->GetParameter(1);
  // Get the function value at the best fit.
  //double minimum = chi(el0, az0);
  //std::cout<<"az0 = "<<el0<<std::endl;
  //std::cout<<"el0 = "<<az0<<std::endl;
  //std::cout<<"Minimum chi^2 = "<<minimum<<std::endl;

  //testparams
  // for drive dependency az0 = 12.1029 el0 = -1.23133
  double bestX=1.0198;//1.02 1.02252 ;
  double bestY=0.0225;//0.023 0.0206697;
  double az0 = 6.6;
  double el0 = -10.48;

  //calculate differences
  for(int i=0; i<kk; i++){
    del_vec[i] = eld_vec[i]-funcEl(el_vec[i],az_vec[i],bestX,bestY);
    //del_vec[i] = eld_vec[i]-el_vec[i];
    //del_vec[i] = eld_vec[i];
    if (del_vec[i] < -180){
      del_vec[i]+=360;
    }
    if (del_vec[i] > 180){
      del_vec[i]-=360;
    }
  }				
  for(int i=0; i<kk; i++){
    daz_vec[i] = funcAz(el_vec[i],az_vec[i],el0,az0)-azd_vec[i];
    //daz_vec[i] = azd_vec[i]-az_vec[i];
    if (daz_vec[i] < -180){
      daz_vec[i]+=360;
    }
    if (daz_vec[i] > 180){
      daz_vec[i]-=360;
    }
  }
  std::cout<<"test"<<funcEl(65,0,1,0)<<std::endl;
  //chitest
  Double_t chisq=0;
  for(int i=0; i<kk; i++){
    chisq+=pow(del_vec[i],2)+pow(daz_vec[i],2);
  }
  std::cout<<"chisqtest = "<<chisq<<std::endl;

  //plot
  TCanvas* can = new TCanvas("plots","Plots",0,0,800,600);
  can->Divide(2,2);
  TString nam("run341C2D.png");
  TString tit1("fit CCD to drive");
  TString tit2("el0 = ");
  //tit2 += Int_t(AZ);
  tit2 += el0;
  tit2 += ", az0 = ";
  tit2 += az0;
  TGraph* g_delel=new TGraph(kk,el_vec,del_vec);
  can->cd(1);
  g_delel->SetMarkerStyle(20);
  g_delel->SetMarkerSize(0.80);
  g_delel->GetXaxis()->SetTitle("elevation drive (deg)");
  g_delel->GetYaxis()->SetTitle("#Delta elevation CCD (deg)");
  g_delel->SetTitle(tit1);
  g_delel->Draw("AP");
  TGraph* g_delaz=new TGraph(kk,az_vec,del_vec);
  can->cd(2);
  g_delaz->SetMarkerStyle(20);
  g_delaz->SetMarkerSize(0.80);
  g_delaz->GetXaxis()->SetTitle("azimuth drive (deg)");
  g_delaz->GetYaxis()->SetTitle("#Delta elevation CCD (deg)");
  g_delaz->SetTitle(tit2);
  g_delaz->Draw("AP");
  TGraph* g_dazel=new TGraph(kk,el_vec,daz_vec);
  can->cd(3);
  g_dazel->SetMarkerStyle(20);
  g_dazel->SetMarkerSize(0.80);
  g_dazel->GetXaxis()->SetTitle("elevation drive (deg)");
  g_dazel->GetYaxis()->SetTitle("#Delta azimuth CCD (deg)");
  g_dazel->Draw("AP");
  TGraph* g_dazaz=new TGraph(kk,az_vec,daz_vec);
  can->cd(4);
  g_dazaz->SetMarkerStyle(20);
  g_dazaz->SetMarkerSize(0.80);
  g_dazaz->GetXaxis()->SetTitle("azimuth drive (deg)");
  g_dazaz->GetYaxis()->SetTitle("#Delta azimuth CCD (deg)");
  g_dazaz->Draw("AP");
  can->SaveAs(nam);
  file->Close();
  }


