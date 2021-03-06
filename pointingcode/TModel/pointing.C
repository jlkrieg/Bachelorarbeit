#include <iostream>
#include <fstream>
#include <vector>
#include <TROOT.h>
#include <TMinuit.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1F.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSpline.h>
#include <TLine.h>
#include <TMatrixD.h>

vector <Double_t> elc_vec, eld_vec, azc_vec, azd_vec; // Data from the measurements
vector <Double_t> psi_vec; // Angular distance
vector <Double_t> err_psi; // Errors for the measurements

// Conversion between degrees and radians, and the opposite:
const Double_t deg2rad = TMath::Pi()/180;
const Double_t rad2deg = 180/TMath::Pi();
const Double_t mSiz = .75;

class PointingModel {  
public:
  PointingModel(int np,Double_t* p,const char* modelname) : npar(np){
    par = new Double_t[np];
    SetParams(p);
    name = std::string(modelname);
  }
  virtual ~PointingModel(){ if(par) delete[] par; }
  virtual Double_t deltaAz(Double_t azc,Double_t elc){ return 0; }
  virtual Double_t deltaEl(Double_t azc,Double_t elc){ return 0; }
  int GetNpar(){ return npar; }
  virtual Double_t* getInitialValues(){ return 0; }
  void SetParams(Double_t* p){
    for(int i=0;i<npar;i++){
      par[i] = (p!=0) ? p[i] : 0.0;
    }
  }
  std::string GetName(){ return name; }
  void read(const char *file){
    std::string s;
    ifstream f(file);
    f >> s;
    if( s!=name ){
      std::cout << "This must not happen: " << s << " vs " << name << std::endl;
      exit(0);
    }
    int n = 0;
    f >> n;
    if( n!=npar ){
      std::cout << "This must not happen: " << n << " vs " << npar << std::endl;
      exit(0);
    }
    for(int i=0;i<npar;i++){
      f >> par[i];
    }
    f.close();
  }
  void write(const char* file){
    TString fname(file);
    fname += name;
    fname += ".pm";
    ofstream f(fname.Data());
    f << name << std::endl;
    f << npar << std::endl;
    for(int i=0;i<npar;i++)
      f << par[i] << std::endl;
    f.close();
  }
  //private:
  Double_t* par;
  int npar;
  std::string name;
};

class BasicModel4 : public PointingModel {
public:
  BasicModel4(Double_t* p=0,const char* nam="BasicModel4") : PointingModel(4,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){
    return tan(elc)*(par[1] + par[2]*cos(azc) + par[3]*sin(azc));
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc);
  }
};

class BasicModel5 : public PointingModel {
public:
  BasicModel5(Double_t* p=0,const char* nam="BasicModel5") : PointingModel(5,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){
    return tan(elc)*(par[1] + par[2]*cos(azc) + par[3]*sin(azc));
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};

class BasicModel6 : public PointingModel {
public:
  BasicModel6(Double_t* p=0,const char* nam="BasicModel6") : PointingModel(6,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){
    return tan(elc)*(par[1] + par[2]*cos(azc) + par[3]*sin(azc))+par[5]*pow(sin(azc),2);
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};

class BasicModel8 : public PointingModel {
public:
  BasicModel8(Double_t* p=0,const char* nam="BasicModel8") : PointingModel(8,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){
    return tan(elc)*(par[1] + par[2]*cos(azc) + par[3]*sin(azc))+par[5]*pow(sin(azc),2)+par[6]*sin(azc-par[7]);
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};

class BasicModel9 : public PointingModel {
public:
  BasicModel9(Double_t* p=0,const char* nam="BasicModel9") : PointingModel(9,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){
    return tan(elc)*(par[1] + par[2]*cos(azc) + par[3]*sin(azc))+par[5]*pow(sin(azc),2)+par[6]*cos(azc-par[7])+par[8];
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};

class BasicModel10 : public PointingModel {
public:
  BasicModel10(Double_t* p=0,const char* nam="BasicModel10") : PointingModel(10,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){ 
    Double_t ret = tan(elc)*(par[1]+par[2]*cos(azc) + par[3]*sin(azc));
    return ret + par[5]*cos(azc-par[6])+par[7]*elc+par[8]*pow(elc,2)+par[9]*pow(sin(azc),2);

  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};

class BasicModel11 : public PointingModel {
public:
  BasicModel11(Double_t* p=0,const char* nam="BasicModel11") : PointingModel(11,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){
    Double_t ret = tan(elc)*(par[1] + par[2]*cos(azc) + par[3]*sin(azc))+par[5]*pow(sin(azc),2)+par[6]*cos(azc-par[7])+par[8];
    return ret + par[9]*elc+par[10]*elc*elc;
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};

class BasicModel11V2 : public PointingModel {
public:
  BasicModel11V2(Double_t* p=0,const char* nam="BasicModel11V2") : PointingModel(11,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){ 
    Double_t ret = tan(elc)*(par[1]+par[2]*cos(azc) + par[3]*sin(azc));
    return ret + par[5]*cos(azc-par[6])+par[7]*azc+par[8]*cos(elc)+par[9]*pow(sin(azc),2)-par[10]*sin(azc);
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
    return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};


class BasicModel11V3 : public PointingModel {
public:
  BasicModel11V3(Double_t* p=0,const char* nam="BasicModel11V3") : PointingModel(11,p,nam) {}
  Double_t deltaAz(Double_t azc,Double_t elc){ 
//     Double_t ret = tan(elc)*(par[1]+ (par[2]*cos(azc)+par[5]*cos(azc-par[6])+ par[3]*sin(azc)));
//     return ret + par[7]*(exp(elc)-sin(azc))+par[8]*cos(elc)+par[9]*pow(sin(azc),2);
    Double_t ret = tan(elc)*(par[1]+par[2]*cos(azc) + par[3]*sin(azc));
    return ret + par[5]*cos(azc-par[6])*tan(elc)+par[7]*sin(elc)+par[8]*cos(elc)+par[9]*pow(sin(azc),2)-par[10]*sin(azc);
  }
  Double_t deltaEl(Double_t azc,Double_t elc){
//     return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
     return par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*elc;
  }
};


class TModel : public PointingModel {
public:
  TModel(Double_t* p=0,const char* nam="TModel") : PointingModel(11,p,nam) {}
  
  Double_t deltaEl(Double_t az,Double_t el){
    Double_t ret = par[0] - par[2]*sin(az) + par[3]*cos(az) + par[4]*el;
    return ret;
  }
  
  Double_t deltaAz(Double_t az,Double_t el){
    Double_t ret = tan(el)*(par[1]+par[2]*cos(az) + par[3]*sin(az));
    ret += par[5]*cos(el);
    ret += par[6]*sin(el);
    ret += par[7]*sin(az);
    ret += par[8]*cos(az);
    ret += par[9]*sin(2*az);
    ret += par[10]*cos(2*az);
    return ret;
  }
  
};

class UModel : public PointingModel {
public:
  UModel(Double_t* p=0,const char* nam="UModel") : PointingModel(16,p,nam) {}
  
  Double_t deltaEl(Double_t az,Double_t el){
    Double_t ret = par[0]+(el-asin(cos(par[1])*sin(el)));
    ret += par[4]*sin(az)+par[5]*cos(az);
    ret += par[8]*sin(2*az)+par[9]*cos(2*az);
    ret += par[10]*sin(el)+par[11]*cos(el);
    ret += par[14]*sin(2*el)+par[15]*cos(2*el);
    return ret;
  }
  
  Double_t deltaAz(Double_t az,Double_t el){
    Double_t ret = par[2]-(el-asin(cos(par[3])*sin(el)));
    ret += par[6]*sin(az)+par[7]*cos(az);
    ret += par[12]*sin(el)+par[13]*cos(el);
    return ret;
  }

  Double_t* getInitialValues(){
    static Double_t arr[] = { 0.8*deg2rad, 14./180*TMath::Pi(), -12*deg2rad, 8./180*TMath::Pi(),
			      0.1*deg2rad,-0.1*deg2rad,0.04*deg2rad, 0.01, 0.01,-0.01,0.01,0.01,0.01,0.01,0.01,0.01  };
    return arr;
  }
};


double delta(double az1,double az2){ //expects deg!
  if(az1<0) az1 +=  360;
  if(az2<0) az2 +=  360;
  double d = az1-az2;
  if( d<-180 )
    d = 360 + d;
  else if( d>180 )
    d = d - 360;
  return d; //returns deg
}

Double_t func1(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Double_t function1 = (-par[2]*sin(xx))+(par[3]*cos(xx))+par[0];
  return function1;
}

Double_t func2(Double_t *x, Double_t *par){
  Double_t xx = x[0]; 
  Double_t function2 = (par[2]*cos(xx))+(par[3]*sin(xx))+par[1];
  return function2;
}

// Compare actual and predicted position for the drive system (version 1):
Double_t fit_func(float elc, float eld, float azc, float azd, Double_t *par) {
  //Double_t del = (par[0] - par[2]*sin(azc) + par[3]*cos(azc));
  Double_t del = (par[0] - par[2]*sin(azc) + par[3]*cos(azc) + par[4]*(elc-1.2));
  Double_t daz = tan(elc)*(par[1] + par[2]*cos(azc) + par[3]*sin(azc));  
  // Create the new coordinates eld_p(eld') and azd_p(azd'):
  Double_t eld_p = elc + del;
  Double_t azd_p = azc + daz;
  // Calculate the angular distance between the original and the new point:
  Double_t value = cos(eld_p-eld) - cos(eld_p)*cos(eld)*(1-cos(azd_p-azd));
  value = acos(value);
  return value;
}

Double_t GetDiff(PointingModel* point,Double_t azc,Double_t elc,Double_t azd,Double_t eld,Double_t* diffEl=0,Double_t* diffAz=0){
  Double_t del = point->deltaEl(azc,elc);
  Double_t daz = point->deltaAz(azc,elc); 
  // Create the new coordinates eld_p(eld') and azd_p(azd'):
  Double_t eld_p = elc + del;
  Double_t azd_p = azc + daz;
  // Calculate the angular distance between the original and the new point:
  Double_t value = cos(eld_p-eld) - cos(eld_p)*cos(eld)*(1-cos(azd_p-azd));
  if(diffEl) *diffEl = (eld-elc)-del;
  if(diffAz) *diffAz = ((azd-azc)-daz);
  return acos(value);
}

PointingModel* point = 0;
// Calculate chi square:
void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  Double_t chisq = 0;
  point->SetParams(par);
  const int n = azc_vec.size();
  for (int i=0;i<n; i++) {
    Double_t delta = GetDiff(point,azc_vec[i],elc_vec[i],azd_vec[i],eld_vec[i])/err_psi[i];
    chisq += delta*delta;
  }
  f = chisq;
}

//Set-up for Simulated Data
const double phix = 0.03*deg2rad;
const double phiy = 0.1*deg2rad;
const double offset_del = 0.5*deg2rad;
const double offset_daz = 0.25*deg2rad;
const double sigma = 0.01*deg2rad;
const Double_t error_psi = 0.02*deg2rad; // Assumed to be the same for all points

void FitModel(PointingModel* point,double meanDeltaEl,double meanDeltaAz){
  // Start of minimization:
  // Initialize TMinuit with a maximum of 6 parameters:
  const int npar = point->GetNpar();
  TMinuit *ptMinuit = new TMinuit(npar);
 
  ptMinuit->SetPrintLevel(3);
  ptMinuit->SetFCN(calc_chi_square);

  
  Double_t arglist[10]; 
  Int_t ierflg = 0;
  arglist[0] = 1;
  ptMinuit->mnexcm("SET ERRORS",arglist,1,ierflg);
  Double_t vstart[] = {meanDeltaEl*deg2rad, meanDeltaAz*deg2rad, 0.1*deg2rad, 0.1*deg2rad,0.07,-0.02,0.015,0.014,-0.9,0.2,-0.3};
  Double_t step[] = {0.01*deg2rad, 0.01*deg2rad, 0.01*deg2rad, 0.01*deg2rad,0.05,0.001,0.001,0.001,0.001,0.001,0.001};

  int NP = sizeof(vstart)/sizeof(vstart[0]);
  Double_t *ini = point->getInitialValues();
  for(int i=0;i<npar;i++){
    TString aname("a");
    aname += i;
    if( i<NP ){
      // Set starting values and step sizes for parameters.
      if(ini==0) 
	ptMinuit->mnparm(i, aname.Data(), vstart[i], step[i], -2,2,ierflg);
      else {
	double vary = 5*fabs(ini[i]);
	double min = ini[i]-vary;
	double max = ini[i]+vary;
	ptMinuit->mnparm(i, aname.Data(), ini[i], 0.1*vary,min,max,ierflg);
	std::cout << "US " << i << " " << aname.Data() << " " << ini[i] << " " << 0.1*vary << " " << min << ".." << max << std::endl;
      }
    } else
      ptMinuit->mnparm(i, aname.Data(), 1e-5*gRandom->Gaus(),1e-6, -.1,.1,ierflg);
  }
#if 0  
  if(npar>4) ptMinuit->mnparm(4, "a4", vstart[4], step[4], -0.5,0.5,ierflg);
  if(npar>5) ptMinuit->mnparm(5, "a5", vstart[5], step[5], -0.5,0.5,ierflg);
  if(npar>6) ptMinuit->mnparm(6, "a6", vstart[6], step[6], -0.5,0.5,ierflg);
  if(npar>7) ptMinuit->mnparm(7, "a7", vstart[7], step[7], -0.5,0.5,ierflg);
  if(npar>8) ptMinuit->mnparm(8, "a8", vstart[8], step[8], -0.5,0.5,ierflg);
  if(npar>9) ptMinuit->mnparm(9, "a9", vstart[9], step[9], -0.5,0.5,ierflg);
  if(npar>10) ptMinuit->mnparm(10, "a10", vstart[10], step[10], -0.5,0.5,ierflg);
#endif

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  ptMinuit->mnexcm("MIGRAD",arglist,2,ierflg);


  // Results of minimization:
  // Get the parameter values:
  Double_t par_a[npar];
  Double_t par_a_err[npar];
  for(int i=0;i<npar;i++){
    ptMinuit->GetParameter(i,par_a[i],par_a_err[i]);
    printf("%02i) %7.4f pm %7.4f deg, sigma = %5.2f\n",i,par_a[i]*rad2deg,par_a_err[i]*rad2deg,fabs(par_a[i])/par_a_err[i]);
  }
  for(int i=0;i<npar;i++){
    printf("  g->SetPoint(%i,%i,%.10f); g->SetPointError(%i,0,%.10f);\n",i,i,par_a[i]*rad2deg,i,par_a_err[i]*rad2deg);
  }
  point->SetParams(par_a);
  
  // For access to these parameters, use:
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  std:: cout << std::endl << "Print results from minuit" << std::endl << std::endl;
  std::cout << "Minimum chi square = " << amin << std::endl;
  std::cout << "Estimated vert. distance to min. = " << edm << std::endl;
  std::cout << "Number of variable parameters = " << nvpar << std::endl;
  std::cout << "Highest number of parameters defined by user = " << nparx << std::endl;
  std::cout << "Status of covariance matrix = " << icstat << std::endl;

  ptMinuit->mnprin(3,amin);

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
  
}


Double_t polyfunc(Double_t* x,Double_t* par){
  const int Poly = int(par[0]+0.5);
  if( Poly<3 ){
    std::cout << "BAD BAD IDEA" << std::endl;
    exit(0);
  }
  const int N = (Poly & 1) ? Poly : Poly-1;
  Double_t para[Poly+1];
  Double_t ret = 0;
  int k = 1;
  para[N]=0;
  for(int i=0;i<=Poly;i++){
    if( i!=N ){
      para[i] = par[k];
      //std::cout << "para[" << i << "] = par[" << k << "]" << std::endl;
      if( i & 1 ){
	//std::cout << "b para[" << N << "]=" << para[N] <<std::endl;
	para[N] -= para[i]*pow(TMath::Pi(),i-1);
	//std::cout << "a para[" << N << "]=" << para[N] <<std::endl;
      }
    }
    k++;
  }
  para[N] /= pow(TMath::Pi(),N-1);
  for(int i=0;i<=Poly;i++){
    ret += para[i]*pow(x[0],i);
    //std::cout << i << ") " << para[i];
    //if( i==N ) std::cout << " * ";
    //std::cout << std::endl;
  }
  return ret;
}

TF1* makePoly(const int N,const char* name){
  TF1* f = new TF1(name,polyfunc,-TMath::Pi(),+TMath::Pi(),N+1);
  f->FixParameter(0,N);
  const int Npar = N;
  for(int i=0;i<Npar;i++)
    f->SetParameter(1+i,1e-3*gRandom->Gaus());
  f->SetLineColor(kGreen);
  return f;
}

double mround(Double_t r){
  double rr(r);
  int sign = 1;
  const double fak = 1.0;
  if( r<0 ){
    r = -r;
    sign = -1;
  }
  double s = 1+int(log10(r*fak)+0.5);
  double z = pow(10.0,s);
  double ret = sign*round(r*fak*z)/z;
  //std::cout << "MROUND " << rr << " -> " << ret << std::endl;
  return ret;
}

void pointing_tel(const char* name="run038_0102_tpoint_0_00_nt_ntuple",bool fit=true){

//pointing_in.root
//nt

// ntuple2nt.root
// run038_0102_tpoint_0_00_nt_ntuple
// run039_0102_tpoint_0_00_nt_ntuple
// run091_tpoint_0_00_nt_ntuple
 // run093_01_tpoint_0_00_nt_ntuple
// run095_01_tpoint_0_00_nt_ntuple

  //CREATE
  //point = new UModel(); 
  point = new TModel();
  TGraphErrors* g_del = new TGraphErrors();
  TGraphErrors* g_daz = new TGraphErrors();
  TGraphErrors* g_del_el = new TGraphErrors();
  TGraphErrors* g_daz_el = new TGraphErrors();
  
  TNtuple* nt = (TNtuple*)gDirectory->Get(name);
  if(nt==0){
    std::cout << "No ntuple " << name << " found" << std::endl;
    return;
  }
  const Int_t N = nt->GetEntries();   //Check number of entries

 //Create a pointer to the arguments in nt:
  const Float_t* a = nt->GetArgs(); // Has to be Float_t 
  int ID = 1;

  int gindex = 0;
  Double_t meanDeltaAz = 0;
  Double_t meanDeltaEl = 0;
  for (int i=0; i<N; i++){
    nt->GetEntry(i);   //Get the correct entry
    int id = int(a[4]+0.5);   //Check image ID, 1 or 2
    if(id>0 && id!=ID) continue;    //Only use image 1 for each star
    Double_t azd = a[0];
    Double_t eld = a[1];
    Double_t azc = a[2];
    Double_t elc = a[3];
    Double_t ps = a[5];
    double Delta = delta(azd,azc);    
    double dist = a[6];    
    std::cout << i << ") " << azd << " " << eld << " " << azc << " " << elc << " " << Delta << std::endl;
    if( !(fabs(ps-11.03)<0.03 || fabs(ps-21.9)<0.015 || fabs(ps)<1e-5 ) ){
      std::cout << "Skip, pixel scale is " << ps << std::endl;
      continue;
    }
    if( false && dist>2 ){
      std::cout << "Skip, dist is " << dist << std::endl;
      continue;
    }
    if( fabs(dist-257.5)<2.5 || dist>480 ){
      std::cout << "Skip, t is " << dist << std::endl;
      continue;
    }
    if( eld<11 ){
      std::cout << "Skip, elevation is " << eld << std::endl;
      continue;
    }
    
    
    //Store the data in global variables:
    azd_vec.push_back(deg2rad*azd);
    eld_vec.push_back(deg2rad*eld);
    azc_vec.push_back(deg2rad*azc);
    elc_vec.push_back(deg2rad*elc);
    err_psi.push_back(error_psi);

    g_del->SetPoint(gindex,azc*deg2rad,eld-elc);
    g_del->SetPointError(gindex,0,sigma*rad2deg);
    g_daz->SetPoint(gindex,azc*deg2rad,Delta);    
    g_daz->SetPointError(gindex,0,sigma*rad2deg);
    g_del_el->SetPoint(gindex,elc*deg2rad,eld-elc);
    g_del_el->SetPointError(gindex,0,sigma*rad2deg);
    g_daz_el->SetPoint(gindex,elc*deg2rad,Delta);
    g_daz_el->SetPointError(gindex,0,sigma*rad2deg);

    meanDeltaAz += Delta/tan(elc*deg2rad);
    meanDeltaEl += eld-elc; 
    gindex++;
  }
  if( gindex>0 ){
    meanDeltaAz /= Double_t(gindex);
    meanDeltaEl /= Double_t(gindex);
  }
  std::cout << "meanDeltaAz = " << meanDeltaAz << " deg" << std::endl;
  std::cout << "meanDeltaEl = " << meanDeltaEl << " deg" << std::endl;
  g_daz->SetMarkerStyle(20);
  g_daz->SetMarkerSize(mSiz);
  g_daz_el->SetMarkerStyle(20);
  g_daz_el->SetMarkerSize(mSiz);
  g_del->SetMarkerStyle(20);
  g_del->SetMarkerSize(mSiz);
  g_del_el->SetMarkerStyle(20);
  g_del_el->SetMarkerSize(mSiz);

  if( fit )
    FitModel(point,meanDeltaEl,meanDeltaAz);
  else
    point->read("run095_01_tpoint_0_00_nt_ntuple_TModel.pm");
  
  std::vector< Double_t > vreso;
  std::vector< Double_t > vreso_el;
  std::vector< Double_t > vreso_az;
  Double_t diffEl(0),diffAz(0);
  for(int i=0;i<gindex;i++){
    Double_t delta = GetDiff(point,azc_vec[i],elc_vec[i],azd_vec[i],eld_vec[i],&diffEl,&diffAz);
    Double_t v = delta*rad2deg*3600;
    vreso.push_back(v*v);
    vreso_el.push_back(diffEl*rad2deg*3600);
    vreso_az.push_back(diffAz*rad2deg*3600);
  }
  sort(vreso.begin(),vreso.end());
  sort(vreso_el.begin(),vreso_el.end());
  sort(vreso_az.begin(),vreso_az.end());
  sort(vreso.begin(),vreso.end());
  TH1F* hreso = new TH1F("hreso","hreso",gindex/2,0,1.1*vreso.back());
  hreso->GetXaxis()->SetTitle("dN/d#Theta^2 (arcsec^2)");
  hreso->GetYaxis()->SetTitle("entries");
  double dis = 20.0*(vreso_el.back()-vreso_el.front())/gindex;
  TH1F* hreso_el = new TH1F("hreso_el","hreso_el",gindex/2,vreso_el.front()-dis,vreso_el.back()+dis);
  hreso_el->GetXaxis()->SetTitle("#Delta el (arcsec)");
  hreso_el->GetYaxis()->SetTitle("entries");
  dis = 20.0*(vreso_az.back()-vreso_az.front())/gindex;
  TH1F* hreso_az = new TH1F("hreso_az","hreso_az",gindex/2,vreso_az.front()-dis,vreso_az.back()+dis);
  hreso_az->GetXaxis()->SetTitle("#Delta az (arcsec)");
  hreso_az->GetYaxis()->SetTitle("entries");
  const double thresh = vreso.size()*0.683;
  double resolution = -1;
  int n=0;
  for(std::vector< Double_t >::iterator i=vreso.begin();i!=vreso.end();i++){
    hreso->Fill(*i);
    if( resolution<0 && n>thresh ){
      resolution = sqrt(*i);
    }
    n++;
  }
  for(std::vector< Double_t >::iterator i=vreso_el.begin();i!=vreso_el.end();i++) hreso_el->Fill(*i);
  for(std::vector< Double_t >::iterator i=vreso_az.begin();i!=vreso_az.end();i++) hreso_az->Fill(*i);
  gStyle->SetOptFit(1111111111);
  hreso_el->Fit("gaus");
  hreso_az->Fit("gaus");
  Double_t minDel = 1e5;
  Double_t maxDel = -1e5;
  Double_t minDaz = 1e5;
  Double_t maxDaz = -1e5;
  TGraphErrors* k_del = new TGraphErrors(gindex);
  TGraphErrors* k_del_el = new TGraphErrors(gindex);
  TGraphErrors* k_daz = new TGraphErrors(gindex);
  TGraphErrors* k_daz_el = new TGraphErrors(gindex);
  k_daz->SetTitle("");
  k_daz->SetMarkerStyle(20);
  k_daz->SetMarkerSize(mSiz);
  k_daz_el->SetTitle("");
  k_daz_el->SetMarkerStyle(20);
  k_daz_el->SetMarkerSize(mSiz);
  k_del->SetTitle("");
  k_del->SetMarkerStyle(20);
  k_del->SetMarkerSize(mSiz);
  k_del_el->SetTitle("");
  k_del_el->SetMarkerStyle(20);
  k_del_el->SetMarkerSize(mSiz);

  TFile f6("delme.root","recreate");
  g_del_el->Write("g_del_el");
  g_daz_el->Write("g_daz_el");
  f6.Close();

  for(int i=0;i<gindex;i++){
    Double_t azc,elc,del;
    g_del->GetPoint(i,azc,del);
    g_del_el->GetPoint(i,elc,del);
    if( del<minDel ) minDel = del;
    if( del>maxDel ) maxDel = del;
    k_del->SetPoint(i,azc,del);
    k_del->SetPointError(i,g_del->GetErrorX(i),g_del->GetErrorY(i));
    k_del_el->SetPoint(i,elc,del);
    k_del_el->SetPointError(i,g_del_el->GetErrorX(i),g_del_el->GetErrorY(i));
    Double_t model = point->deltaEl(azc,elc)*rad2deg;
    g_del->SetPoint(i,azc,del-model);
    g_del_el->SetPoint(i,elc,del-model);
  }
  for(int i=0;i<gindex;i++){
    Double_t azc,elc,daz;
    g_daz->GetPoint(i,azc,daz);
    g_daz_el->GetPoint(i,elc,daz);
    Double_t fak = 1.0;
    Double_t pdaz = daz*fak;
    if( pdaz<minDaz ) minDaz = pdaz;
    if( pdaz>maxDaz ) maxDaz = pdaz;
    k_daz->SetPoint(i,azc,pdaz);
    k_daz->SetPointError(i,g_daz->GetErrorX(i),g_daz->GetErrorY(i)*fak);
    k_daz_el->SetPoint(i,elc,pdaz);
    k_daz_el->SetPointError(i,g_daz_el->GetErrorX(i),g_daz_el->GetErrorY(i)*fak);
    Double_t model = point->deltaAz(azc,elc)*rad2deg;
    g_daz->SetPoint(i,azc,daz-model);
    g_daz_el->SetPoint(i,elc,daz-model);
  }
  TF1* myf = new TF1("myf","[0]*sin(x-[1])",-TMath::Pi(),TMath::Pi());
  myf->SetParameter(0,-0.01);
  myf->SetParameter(1,0.01);
  myf->SetLineColor(kBlue);
  TF1* myg = new TF1("myg","[0]*sin([1]*x)",-TMath::Pi(),TMath::Pi());
  myg->SetParameter(0,-0.01);
  myg->SetParameter(1,0.01);
  myg->SetLineColor(kGreen);

  const Double_t max = 1.5; //.15
  TCanvas * c5 = new TCanvas("c5","deltaEl and deltaAz",1200,800); 
  c5->Divide(2,2);
  c5->cd(1);
  g_del->SetMinimum(-max);
  g_del->SetMaximum(max);
  g_del->GetXaxis()->SetTitle("azimuth (rad)");
  g_del->GetYaxis()->SetTitle("#Delta el (deg)");
  g_del->Draw("AP");
  g_del->Fit("pol1");
  TF1* tst3 = makePoly(7,"tst3");
  //g_del->Fit("tst3");
  c5->cd(2);
  g_del_el->SetMinimum(-max);
  g_del_el->SetMaximum(max);
  g_del_el->GetXaxis()->SetTitle("elevation (rad)");
  g_del_el->GetYaxis()->SetTitle("#Delta el (deg)");
  g_del_el->Draw("AP");
  g_del_el->Fit("pol1");
  //g_del_el->Fit("tst3");
  	
  TFile f6c("delme_el.root","recreate");
  g_del_el->Write("g_del_el");
  g_daz_el->Write("g_daz_el");
  f6c.Close();

  
  c5->cd(3);
  g_daz->SetMinimum(-max);
  g_daz->SetMaximum(max);
  g_daz->GetXaxis()->SetTitle("azimuth (rad)");
  g_daz->GetYaxis()->SetTitle("#Delta az cos(el) (deg)");
  g_daz->Draw("AP");  
  g_daz->Fit("pol1");
  //g_daz->Fit("tst3");    
  c5->cd(4);
  g_daz_el->SetMinimum(-max);
  g_daz_el->SetMaximum(max);
  g_daz_el->GetXaxis()->SetTitle("elevation (rad)");
  g_daz_el->GetYaxis()->SetTitle("#Delta az cos(el) (deg)");
  g_daz_el->Draw("AP");
  std::cout << "** FIT" << std::endl;    
  g_daz_el->Fit("pol1");
  //g_daz_el->Fit("tst3");  
  std::cout << "** FIT done" << std::endl;  
  c5->cd(0);
  std::cout << "Number of Entries: " << N <<  std::endl << std::endl;
  TString cname(name);
  cname += "_";
  cname += TString(point->GetName());
  if( fit )
    cname += "_1";
  else
    cname += "_0";
  c5->SaveAs(cname+"_del_daz.png");

  TCanvas* c6 = new TCanvas("c6","c6",0,0,800,600);
  TString htit("68.3% below ");
  char help[80];
  sprintf(help,"%.1f",resolution);
  htit += TString(help);
  htit += " arcsec";
  hreso->SetTitle(htit);
  const double hmax = hreso->GetMaximum();
  hreso->SetMaximum(hmax+1);
  hreso->Draw();
  TLine* l6 = new TLine(resolution*resolution,0,resolution*resolution,hmax);
  l6->SetLineColor(kRed);
  l6->SetLineWidth(2);
  l6->Draw();
  TString c6name(cname);
  c6->SaveAs(c6name+"_hreso.png");
  if( fit ){
    TString pm_name(name);
    pm_name += "_";
    point->write(pm_name.Data());
  }

  TCanvas* c8 = new TCanvas("c8","c8",0,0,1000,500);
  c8->Divide(2,1);
  {
    const double hmax = hreso_az->GetMaximum()>hreso_el->GetMaximum() ? hreso_az->GetMaximum() : hreso_el->GetMaximum();
    hreso_el->SetMaximum(hmax+2);
    hreso_az->SetMaximum(hmax+2);
  }
  c8->cd(1); hreso_el->Draw();
  c8->cd(2); hreso_az->Draw();
  c8->cd(0);
  TString c8name(cname);
  c8->SaveAs(c8name+"_hreso_elaz.png");
  
  //********** show data ****************
  std::cout << "Daz " << minDaz << ".." << maxDaz << std::endl;
  std::cout << "Del " << minDel << ".." << maxDel << std::endl;
  double marginDaz = (maxDaz-minDaz)/4.0;
  minDaz = mround(minDaz-marginDaz);
  maxDaz = mround(maxDaz+marginDaz);
  double marginDel = (maxDel-minDel)/4.0;
  minDel = mround(minDel-marginDel);
  maxDel = mround(maxDel+marginDel);
  TCanvas * c7 = new TCanvas("c7","Data",1200,800); 
  c7->Divide(2,2);
  c7->cd(1);
  k_del->SetMinimum(minDel);
  k_del->SetMaximum(maxDel);
  k_del->GetXaxis()->SetTitle("azimuth (rad)");
  k_del->GetYaxis()->SetTitle("#Delta el (deg)");
  k_del->Draw("AP");
  c7->cd(2);
  k_del_el->SetMinimum(minDel);
  k_del_el->SetMaximum(maxDel);
  k_del_el->GetXaxis()->SetTitle("elevation (rad)");
  k_del_el->GetYaxis()->SetTitle("#Delta el (deg)");
  k_del_el->Draw("AP");
  c7->cd(3);
  k_daz->SetMinimum(minDaz);
  k_daz->SetMaximum(maxDaz);
  k_daz->GetXaxis()->SetTitle("azimuth (rad)");
  k_daz->GetYaxis()->SetTitle("#Delta az cos(el) (deg)");
  k_daz->Draw("AP");  
  c7->cd(4);
  k_daz_el->SetMinimum(minDaz);
  k_daz_el->SetMaximum(maxDaz);
  k_daz_el->GetXaxis()->SetTitle("elevation (rad)");
  k_daz_el->GetYaxis()->SetTitle("#Delta az cos(el) (deg)");
  k_daz_el->Draw("AP");
  c7->cd(0);
  TString cname2(cname);
  c7->SaveAs(cname2+"_data.png");
}

void pointing(const char* nt = "run038_0102_tpoint_0_00_nt_ntuple",bool fit=true){
  gStyle->SetOptStat(1);
  gStyle->SetOptTitle(1);
  pointing_tel(nt,fit);
}
