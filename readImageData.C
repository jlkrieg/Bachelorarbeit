#include <iostream>
#include <TStyle.h>
#include <TH2.h>
#include <TFITS.h>
#include <TImage.h>
#include <TH1.h>
#include <TH1F.h>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <TNtuple.h>
#include <TFile.h>

void readImageData()
{
  //Dateinamen aus Verzeichnis einlesen
  fstream f;
  f.open("run199datanames.txt",fstream::out | fstream::trunc);
  char pfad[] ="/datr1/HESS/schwanke/cta/run199/";
  DIR *verzeichnis;
  struct dirent *datei;
  verzeichnis = opendir(pfad);
  while (datei = readdir(verzeichnis)) {
    f << datei->d_name << std::endl;
  }
  f.close();


  int i=0;
  TNtuple *ntuple = new TNtuple("ntuple", "Auswertung","gain:exposure:median:width", 32000);
  std::ifstream infile("run199datanames.txt");


  // Bilder einlesen
  for (std::string line; getline(infile,line);){
    i++;
    std::string str= "/datr1/HESS/schwanke/cta/run199/"+line;
    const char *fname = str.c_str();
    // Gain und Belichtungszeit einlesen
    int gain;
    int exposure;
    if ((line.find(".fit")!=std::string::npos) && not (line.find("_02_")!=std::string::npos)){
      if (line.find("gain0")!=std::string::npos){
	gain=0;
      }
      else if (line.find("gain7")!=std::string::npos){
	gain=7;
      }
      else if (line.find("gain14")!=std::string::npos){
	gain=14;
      }
      else if (line.find("gain21")!=std::string::npos){
	gain=21;
      }
      else {
	std::cout<<"Errory"<<std::endl;
	return;
      }
      if (line.find("__6.0_")!=std::string::npos){
	exposure=1;
      }
      else if (line.find("__7.0_")!=std::string::npos){
	exposure=10;
      }
      else if (line.find("__7.3_")!=std::string::npos){
	exposure=20;
      }
      else {
	std::cout<<"ERROR"<<std::endl;
	return;
      }
      //std::cout<<"gain"<<gain<<std::endl;
      //std::cout<<"expousre"<<exposure<std::endl;
      //std::cout<<fname<<std::endl;
      
      //Bild als Histogramm einlesen
      TFITSHDU *hdu = new TFITSHDU(fname);
      TH2* h = dynamic_cast<TH2*>(hdu->ReadAsHistogram());
      TH1F* g=new TH1F("g","",256,-0.5,255.5);
      for(int x=1;x<=h->GetNbinsX();x++){
	for(int y=1;y<=h->GetNbinsY();y++){
	  g->Fill((h->GetBinContent(x,y)-1));
	}
      }
      delete hdu;
      delete h;
      //Median berechnen
      int n=0;
      double med=0;
      int p=1360*1024;
      for(int x=0;x<=g->GetNbinsX()-1;x++){
	n=n+(int)(g->GetBinContent(x));
	med++;
	if( n>=p/2 ){
	  med=med+(p/2-n)/g->GetBinContent(x)-1;
	  //std::cout << "Median: " << med << std::endl;
	  break;
	}
      }
      //Breite berechnen
      n=0;
      double wl=0;
      double wr=0;
      double w;
      double pl=p*0.1585;
      double pr=p*0.8415;
      for(int x=0;x<=g->GetNbinsX()-1;x++){
	n=n+(int)g->GetBinContent(x);
	wl=wl++;
	if(n>pl){
	  wl=wl+(pl-n)/g->GetBinContent(x)-1; 
	  break;
	}
      }
      n=0;
      for(int x=0;x<=g->GetNbinsX()-1;x++){
	n=n+(int)g->GetBinContent(x);
	wr=wr++;
	if(n>pr){
	  wr=wr+(pr-n)/g->GetBinContent(x)-1;
	  break;
	}
      }
      delete g;
      w=wr-wl;
      //std::cout<<w<<std::endl;

      //Daten in nTuple schreiben
      ntuple->TNtuple::Fill(gain,exposure,med,w);
     }
    //     if (i>90){
    //       break;
    //     }
     std::cout<<i<<std::endl;
  }
  TFile myfile("imgdata.root","recreate");
  ntuple->Write();
  myfile.Close();
  ntuple->Print();
  ntuple->Draw("median:exposure");
  return;
}
