#include <iostream>
#include <math.h>
#include <TGraph.h>
#include <TFile.h>

Double_t rdx(Double_t el,Double_t az){
  el=el*M_PI/180;
  az=az*M_PI/180;
  return cos(el)*sin(az);
}
Double_t rdy(Double_t el,Double_t az){
  el=el*M_PI/180;
  az=az*M_PI/180;
  return cos(el)*cos(az);
}
Double_t rdz(Double_t el,Double_t az){
  el=el*M_PI/180;
  return sin(el);
}
Double_t rcx(Double_t el,Double_t az, Double_t phi0){
  el=el*M_PI/180;
  az=az*M_PI/180;
  phi0=phi0*M_PI/180;
  return sin(phi0)*cos(az)+cos(phi0)*cos(el)*sin(az);
}
Double_t rcy(Double_t el,Double_t az, Double_t phi0){
  el=el*M_PI/180;
  az=az*M_PI/180;
  phi0=phi0*M_PI/180;
  return -sin(phi0)*sin(az)+cos(phi0)*cos(el)*cos(az);
}
Double_t rcz(Double_t el,Double_t az, Double_t phi0){
  el=el*M_PI/180;
  az=az*M_PI/180;
  phi0=phi0*M_PI/180;
  return cos(phi0)*sin(el);
}
Double_t deltaEl(Double_t el,Double_t az, Double_t phi0){
  Double_t temp=(acos(rdz(el,az))-acos(rcz(el,az,phi0)))*180/M_PI;
  if (temp>180) temp=temp-360;
  if (temp<-180) temp=temp+360;
  return temp;
}
Double_t deltaAz(Double_t el,Double_t az, Double_t phi0){
  Double_t temp = (atan2(rdy(el,az),rdx(el,az))-atan2(rcy(el,az,phi0),rcx(el,az,phi0)))*180/M_PI;
  if (temp>180) temp=temp-360;
  if (temp<-180) temp=temp+360;
  return temp;
}

void pointing(){
  Double_t phi0=5;
  int n=91;
  int m=180;
  int m0=-180;
  Double_t el[n],az[m-m0],a[n],b[n],c[n],d[n],e[m-m0],f[m-m0],g[m-m0],h[m-m0],j[n],k[n],l[n],o[n],p[m-m0],q[m-m0],r[m-m0],s[m-m0];
  for (int i=0;i<=n;i++){
    el[i]=i;
    a[i]=deltaAz(i,0,phi0);
    b[i]=deltaAz(i,90,phi0);
    c[i]=deltaAz(i,180,phi0);
    d[i]=deltaAz(i,270,phi0);
  }
  for (int i=0;i<(m-m0);i++){
    az[i]=i+m0;
    e[i]=deltaAz(0,i+m0,phi0);
    f[i]=deltaAz(22.5,i+m0,phi0);
    g[i]=deltaAz(45,i+m0,phi0);
    h[i]=deltaAz(67.5,i+m0,phi0);
    //std::cout<<i<<std::endl;
    //std::cout<<rcy(0,i,phi0)/rcx(0,i,phi0)<<std::endl;
    //std::cout<<atan(rcy(0,i,phi0)/rcx(0,i,phi0))<<std::endl;
  }
  for (int i=0;i<=n;i++){
    el[i]=i;
    j[i]=deltaEl(i,0,phi0);
    k[i]=deltaEl(i,90,phi0);
    l[i]=deltaEl(i,180,phi0);
    o[i]=deltaEl(i,270,phi0);
  }
  for (int i=0;i<(m-m0);i++){
    az[i]=i+m0;
    p[i]=deltaEl(0,i+m0,phi0);
    q[i]=deltaEl(22.5,i+m0,phi0);
    r[i]=deltaEl(45,i+m0,phi0);
    s[i]=deltaEl(67.5,i+m0,phi0);
  }
  TGraph* dazVel=new TGraph(n,el,b);
  dazVel->SetName("dazVel");
  //dazVel->Draw("A*");
  TGraph* dazVaz=new TGraph(m-m0,az,g);
  dazVaz->SetName("dazVaz");
  //dazVaz->Draw("A*");
  TGraph* delVel=new TGraph(n,el,k);
  delVel->SetName("delVel");
  //delVel->Draw("A*");
  TGraph* delVaz=new TGraph(m-m0,az,r);
  delVaz->SetName("delVaz");
  delVaz->Draw("A*");
  TFile myfile("pointing.root","recreate");
  dazVel->Write();
  dazVaz->Write();
  delVel->Write();
  delVaz->Write();
  myfile.Close();
}
