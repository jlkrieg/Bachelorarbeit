void readtest(){ 
  double U =0;
  double DAZ =2;
  double cutEl=35;
  double AZ =0;
  TNtuple*nt = (TNtuple*)gDirectory->Get("nt");
  Int_t N = nt->GetEntries();
  Float_t* a = nt->GetArgs();
  int k = 0;
  for(int i=0;i<N;i++){
    nt->GetEntry(i);
    Float_t u = a[0];
    //if( u == U ) continue; //evtl Hysterese
    Float_t ps = a[6];
    //std::cout<<ps<<std::endl;
    if( fabs(ps-11.03) > 0.05 ) continue;  
    Double_t eld = a[2];
    Double_t azd = a[1];
    Double_t az = a[3];
    if( azd == AZ ) {// continue;
      std::cout<<azd<<std::endl;
      if( eld > cutEl ){ //grosse Ausreisser unter 10deg
	k++;
	}
    }
  }
  std::cout<<k<<std::endl;
}
