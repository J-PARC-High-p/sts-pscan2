#include <fstream>
#include <string>

using namespace std;

void show_hist(const char* filename="pscan/pscan_20211102_RedFEB8.root"){
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  TFile* file = new TFile(filename);

  TString filename_trim = filename;
  int ipos = filename_trim.Last('.');
  if ( ipos >= 0 ) {
    filename_trim.Remove(ipos);
    filename_trim = filename_trim + "_trim.txt";
  }else{
    cout << "unexpected filename" << endl;
    return;
  }

  ofstream out(filename_trim);

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

  TF1* func_calib = new TF1("func_calib","pol1");
  // 20-240. 
  //func_calib->SetParameters(200,-6);
  func_calib->SetParameters(240,-7.33333333);
  const double trim_amp = 18/12.5;
  
  while(ch <= ch_max){
    for( int itmp = 0; itmp < 16;itmp ++){
      if (itmp == 0 ) {
	c1->Clear();
	c1->Divide(4,4);
      }
      c1->cd(itmp+1);

      /// FastTrim
      TString name = TString::Format("h_d_%d_31",ch);
      TH1D* hist_fast = (TH1D*)file->Get(name);
      hist_fast->Fit("gaus");
      TF1* func = hist_fast->GetFunction("gaus");
      int fast_trim = 36;
      if ( func != NULL ){
	double var = func->GetParameter(0);
	fast_trim = -(var-20.)*0.675 + 36;
      }


      string str;
      str += string(TString::Format("ch:%4d",ch));
      name = TString::Format("h_quality_%d",ch);
      TH1D* hist = (TH1D*)file->Get(name);
      hist->SetTitle(name);
      hist->Draw();
      
      for(int i = 0;i<31; i++){
	int ibin = hist->FindBin(i+0.001);
	double var = hist->GetBinContent(ibin);
	int adj = 128.-(func_calib->Eval(i) - var)*trim_amp;
	str += string(TString::Format("%5d",adj));
      }

      str += string(TString::Format("%5d",fast_trim)); // dummy entry for FASTth.
      out << str << endl;
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
	//c1->GetPad(itmp+1)->SetLogy();
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
  out.close();
}
