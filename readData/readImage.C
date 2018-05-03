#include <iostream>
#include <TStyle.h>
#include <TH2.h>
#include <TFITS.h>
#include <TImage.h>
#include <TH1.h>
#include <TH1F.h>


void readImage(const char* fname="/datr1/HESS/schwanke/cta/run103/field_109_-70.2_44.0__7.3_000000_01_2016-09-07T19:32:14.962000_CCD_3.fits"){
  TFITSHDU *hdu = new TFITSHDU(fname);
  TH2* h = dynamic_cast<TH2*>(hdu->ReadAsHistogram());
  gStyle->SetOptLogz(1);
  h->Draw("colz");
//   int count = 0;
//   std::cout << "Bins in x: " << h->GetNbinsX() << std::endl;
//   std::cout << "Bins in y: " << h->GetNbinsY() << std::endl;
// //   for(int x=1;x<=h->GetNbinsX();x++){
// //     for(int y=1;y<=h->GetNbinsY();y++){
// //       std::cout << h->GetBinContent(x,y) << std::endl;
// //       if( ++count>11 ) return;
// //     }
// //   }
//   TH1F* h2=new TH1F("h2","h2",256,-0.5,255.5);
//    for(int x=0;x<=h->GetNbinsX()-1;x++){
//     for(int y=0;y<=h->GetNbinsY()-1;y++){
//       h2->Fill((h->GetBinContent(x,y)-1));
//       if (h->GetBinContent(x,y)>250){
//       std::cout<<h->GetBinContent(x,y)<<std::endl;
//       }
//     }
//   }
//   gStyle->SetOptLogy(1);
//   h2->Draw();
//   int n=0;
//   int p=1360*1024;
//   double med=0;
//   for(int x=0;x<=h2->GetNbinsX()-1;x++){
//     n=n+(h2->GetBinContent(x));
//     med++;
//     if( n>=p/2 ){
//       med  = med+(p/2-n)/h2->GetBinContent(x)-1;
//       std::cout << "Median: " << med << std::endl;
//       break;
//     }
//   }
//   n=0;
//   double wl=0;
//   double wr=0;
//   double pl=p*0.1585;
//   double pr=p*0.8415;
//   for(int x=0;h2->GetNbinsX();x++){
//     n=n+h2->GetBinContent(x);
//     wl=wl++;
//     if(n>pl){
//       wl=wl+(pl-n)/h2->GetBinContent(x)-1; 
//       break;
//     }
//   }
//   n=0;
//   for(int x=0;x<=h2->GetNbinsX()-1;x++){
//     n=n+h2->GetBinContent(x);
//     wr=wr++;
//     if(n>pr){
//       wr=wr+(pr-n)/h2->GetBinContent(x)-1;
//       break;
//     }
//   }
//   std::cout<<wl<<std::endl;
//   std::cout<<wr<<std::endl;
//   std::cout<<wr-wl<<std::endl;
//   //std::cout<<h2->GetNbinsX()<<std::endl;
}
