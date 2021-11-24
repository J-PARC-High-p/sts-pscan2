void show_hist(const char* filename="pscan/pscan_20211102_RedFEB8.root"){
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  TFile* file = new TFile(filename);

  int ch = 0;
  int ch_max = 127;

  TString outputname;
  {
    TString tmp = filename;
    int ipos = tmp.Last('.');
    tmp.Remove(ipos);
    outputname = tmp + ".pdf";
  }

  // SHOW h_quality_%d.
  c1->Print(outputname + "[","pdf");
  while(ch <= ch_max){
    for( int itmp = 0; itmp < 16;itmp ++){
      if (itmp == 0 ) {
	c1->Clear();
	c1->Divide(4,4);
      }
      c1->cd(itmp+1);
      TString name = TString::Format("h_quality_%d",ch);
      TH1D* hist = (TH1D*)file->Get(name);
      hist->SetTitle(name);
      hist->Draw();
      ch++;
      if ( ch > ch_max ) break;
    }
    c1->Draw();
    c1->Update();
    c1->Print(outputname,"PDF");
    //getchar();
  }

  // SHOW scurve and dcnt

  int d_max = 31;
  for( int ch = 0;ch<=ch_max;ch++){
    int d = 0;
    while(d <= d_max){
      for( int itmp = 0;itmp<8;itmp++){
	if ( itmp==0) {
	  c1->Clear();
	  c1->Divide(2,4);
	}
	c1->cd(itmp+1);
	TString name = TString::Format("h_scurve_%d_%d",ch,d);
	TH1D* hist = (TH1D*)file->Get(name);
	hist->SetTitle(name);
	hist->Draw("HIST");
	
	name = TString::Format("h_d_%d_%d",ch,d);
	TH1D* histd = (TH1D*)file->Get(name);
	histd->SetLineColor(kRed);
	histd->Draw("SAME");
	
	d++;
	if ( d > d_max ) break;
      }
      c1->Draw();
      c1->Update();
      c1->Print(outputname,"PDF");
      //getchar();
    }
  }
  c1->Clear();
  c1->Print(outputname + "]","PDF");
}
