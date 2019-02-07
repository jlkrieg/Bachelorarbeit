void dataPlot(){
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
  Double_t azc_vec[kk],elc_vec[kk],azd_vec[kk],eld_vec[kk],daz_vec[kk],del_vec[kk];
  
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
      azc_vec[k]=az;
      elc_vec[k]=el;
      k++;
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
  TCanvas* can = new TCanvas("plots1","Plots1",0,0,800,600);
  can->Divide(2,2);
  TString nam("data1.png");
  TGraph* g_delel1=new TGraph(kk,elc_vec,eld_vec);
  can->cd(1);
  g_delel1->SetMarkerStyle(20);
  g_delel1->SetMarkerSize(0.80);
  g_delel1->GetXaxis()->SetTitle("elevation center (deg)");
  g_delel1->GetYaxis()->SetTitle("elevation drive (deg)");
  g_delel1->Draw("AP");
  TGraph* g_delaz1=new TGraph(kk,azc_vec,eld_vec);
  can->cd(2);
  g_delaz1->SetMarkerStyle(20);
  g_delaz1->SetMarkerSize(0.80);
  g_delaz1->GetXaxis()->SetTitle("azimuth center (deg)");
  g_delaz1->GetYaxis()->SetTitle("elevation drive (deg)");
  g_delaz1->Draw("AP");
  TGraph* g_dazel1=new TGraph(kk,elc_vec,azd_vec);
  can->cd(3);
  g_dazel1->SetMarkerStyle(20);
  g_dazel1->SetMarkerSize(0.80);
  g_dazel1->GetXaxis()->SetTitle("elevation center (deg)");
  g_dazel1->GetYaxis()->SetTitle("azimuth drive (deg)");
  g_dazel1->Draw("AP");
  TGraph* g_dazaz1=new TGraph(kk,azc_vec,azd_vec);
  can->cd(4);
  g_dazaz1->SetMarkerStyle(20);
  g_dazaz1->SetMarkerSize(0.80);
  g_dazaz1->GetXaxis()->SetTitle("azimuth center (deg)");
  g_dazaz1->GetYaxis()->SetTitle("azimuth drive (deg)");
  g_dazaz1->Draw("AP");
  can->SaveAs(nam);

  TCanvas* can = new TCanvas("plots2","Plots2",0,0,800,600);
  can->Divide(2,2);
  TString nam2("data2.png");
  TGraph* g_delel2=new TGraph(kk,elc_vec,del_vec);
  can->cd(1);
  g_delel2->SetMarkerStyle(20);
  g_delel2->SetMarkerSize(0.80);
  g_delel2->GetXaxis()->SetTitle("elevation center (deg)");
  g_delel2->GetYaxis()->SetTitle("#Delta elevation (deg)");
  g_delel2->Draw("AP");
  TGraph* g_delaz2=new TGraph(kk,azc_vec,del_vec);
  can->cd(2);
  g_delaz2->SetMarkerStyle(20);
  g_delaz2->SetMarkerSize(0.80);
  g_delaz2->GetXaxis()->SetTitle("azimuth center (deg)");
  g_delaz2->GetYaxis()->SetTitle("#Delta elevation (deg)");
  g_delaz2->Draw("AP");
  TGraph* g_dazel2=new TGraph(kk,elc_vec,daz_vec);
  can->cd(3);
  g_dazel2->SetMarkerStyle(20);
  g_dazel2->SetMarkerSize(0.80);
  g_dazel2->GetXaxis()->SetTitle("elevation center (deg)");
  g_dazel2->GetYaxis()->SetTitle("#Delta azimuth (deg)");
  g_dazel2->Draw("AP");
  TGraph* g_dazaz2=new TGraph(kk,azc_vec,daz_vec);
  can->cd(4);
  g_dazaz2->SetMarkerStyle(20);
  g_dazaz2->SetMarkerSize(0.80);
  g_dazaz2->GetXaxis()->SetTitle("azimuth center (deg)");
  g_dazaz2->GetYaxis()->SetTitle("#Delta azimuth (deg)");
  g_dazaz2->Draw("AP");
  can->SaveAs(nam2);

  TCanvas* can = new TCanvas("plots3","Plots3",0,0,800,600);
  can->Divide(2,2);
  TString nam3("data3.png");
  TGraph* g_el=new TGraph(kk,elc_vec,azc_vec,del_vec);
  can->cd(1);
  g_el->SetMarkerStyle(20);
  g_el->SetMarkerSize(0.80);
  g_el>GetXaxis()->SetTitle("elevation center (deg)");
  g_el>GetYaxis()->SetTitle("azimuth center (deg)");
  g_el->GetZaxis()->SetTitle("#Delta elevation (deg)");
  g_el->Draw("AP");
  TGraph* g_az=new TGraph2D(kk,elc_vec,azc_vec,daz_vec);
  can->cd(2);
  g_az->SetMarkerStyle(20);
  g_az->SetMarkerSize(0.80);
  g_az>GetXaxis()->SetTitle("elevation center (deg)");
  g_az>GetYaxis()->SetTitle("azimuth center (deg)");
  g_az->GetZaxis()->SetTitle("#Delta elevation (deg)");
  g_az->Draw("AP");
  can->SaveAs(nam3);

  file->Close();
}
