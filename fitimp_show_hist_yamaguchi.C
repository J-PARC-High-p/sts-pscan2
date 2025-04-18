// edited by S.Ochiai 2023.02.03

#include <fstream>
#include <string>
#include <TMath.h>

using namespace std;

double binomial_err(double entry){
  double tmp = 0;
  if(entry == 0){
    tmp = 1;
  }
  else if(entry == 200){
    tmp = 199;
  }
  else{
    tmp = entry;
  }
  double p = tmp/200.;
  int n = 200;

  return 1/sqrt(n*p*(1-p));
}

double ENC(double amp_cal){
  double qin = 0.05216*amp_cal; //sigma[amp_cal]->sigma[fC]
  double enc = qin / 1.6e-4; //ENC = e??? noise coulomb
  return enc;

}

void show_hist(const char* filename="pscan/pscan_20211102_RedFEB8.root"){
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  TFile* file = new TFile(filename);
  
  TH2F *sig_ch[128];
  for(int i=0;i<128;i++){
    sig_ch[i] = new TH2F("sig_ch","sig_ch",32,0,32,60,0,60);
  }

  TH1D *sig_hist[128];
  for(int i=0;i<128;i++){
    sig_hist[i] = new TH1D("sig_hist","sig_hist",60,0,60);
  }

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
      hist_fast->Fit("gaus","0");
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
	if ( itmp==0 ) {
	  c1->Clear();
	  c1->Divide(2,4);
	}
	c1->cd(itmp+1);
	TString name = TString::Format("h_scurve_%d_%d",ch,d);
	TH1D* hist_scurve = (TH1D*)file->Get(name);
	hist_scurve->SetTitle(name);
	
	//give the binomial error
	//TGraphErrors *graph = new TGraphErrors;
	//int nbins = hist_scurve->GetXaxis()->GetNbins();
	//for(int i=1;i<nbins;i+=4){
	//double y = hist_scurve->GetBinContent(i);
	//double x = i;
	//double ex = 0;
	///double ey = binomial_err(y); 
	//graph->SetPoint(i, x, y);
	//graph->SetPointError(i, ex, ey);
	//cout << "x;" << x << " y;" << y << endl;
	//}
	//graph->SetMarkerStyle(20);
	//graph->SetMarkerSize(1);
	//graph->Draw("AP");

	TF1* erf = new TF1("erf", "200*ROOT::Math::gaussian_cdf(x,[0],[1])",0,300);
	erf->SetParameters(20,150);
	//c1->GetPad(itmp+1)->SetLogy();
	hist_scurve->Rebin(4);
	hist_scurve->Draw("HIST");
	hist_scurve->Fit("erf","L","",0,255);//Fitting #1
	
	//----------------------added by S.Ochiai 2023.02.03------------------
	TF1* erf0 = new TF1("erf0", "erf", 0, 300);
	erf0->SetParameters(erf->GetParameters());
	erf0->Draw("same");

	Double_t sigma, x0;
	Double_t xlow;
	Double_t chi2_previous, chi2_current;
	
	//+++++++++Fitting Iteration+++++++++
	for ( Int_t i = 2; i < 11; i++ ){

		chi2_previous = erf->GetChisquare();
		sigma = erf->GetParameter(0);
		x0 = erf->GetParameter(1);
		xlow = x0 - 2.*sigma;

		erf->SetParameters(sigma, x0);
		hist_scurve->Fit("erf", "L", "", xlow, 255);
		chi2_current = erf->GetChisquare();
		
		if( chi2_current > chi2_previous || TMath::Abs( chi2_previous - chi2_current ) < chi2_previous * 0.01){
			break;
		}
	} //+++++++++Fitting Iteration end+++++
	
	//----------------------added by S.Ochiai 2023.02,03 end0-------------
	

	//graph->Fit("erf");
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(1);
	//double chi2 = erf->GetChisquare();
	double chi2 = chi2_previous;
	//double sigma = erf->GetParameter(0);
	int bin1 = hist_scurve->GetBinContent(1);
	sig_ch[ch]->Fill(d,sigma);
	if(chi2 <= 20 && bin1 < 200){
	  sig_hist[ch]->Fill(sigma);
	}

	TF1 *erf2 = new TF1("erf2","erf",0,300);
	//erf2->SetParameters(erf->GetParameters());
	erf2->SetParameters(sigma, x0);
	erf2->SetLineWidth(1);
	erf2->SetLineColor(kBlack);
	erf2->SetLineStyle(7);
	erf2->Draw("same");
	
	TF1* erf3 = new TF1("erf3", "erf", xlow, 255);
	erf3->SetParameters(sigma, x0);
	erf3->SetLineColor(kViolet);
	erf3->Draw("same");
	/*
	name = TString::Format("h_d_%d_%d",ch,d);
	TH1D* histd = (TH1D*)file->Get(name);
	histd->SetLineColor(kRed);
	histd->Draw("SAME");
	*/
	
	d++;
	if ( d > d_max ) break;
      }
      c1->Draw();
      c1->Update();
      c1->Print(outputname,"PDF");
      //getchar();
    }
  }

  //Draw the distribution of sigma of every channel
  ch = 0;
  while(ch <= ch_max){
    for(int itemp=0;itemp<8;itemp++){
      if(itemp == 0){
	c1->Clear();
	c1->Divide(2,4);
      }
      c1->cd(itemp+1);
      TString hist_name = TString::Format("sigma_ch_%d",ch);
      sig_ch[ch]->Draw();
      sig_ch[ch]->SetTitle(hist_name);
      gStyle->SetMarkerStyle(20);
      gStyle->SetMarkerColor(1);
      gStyle->SetMarkerSize(1);
      ch++;
      if(ch > ch_max) break;
    }
    c1->Draw();
    c1->Update();
    c1->Print(outputname,"PDF");
}

  //Show the histgram of sigma
ch = 0;
 TH2D *enc_ch = new TH2D("enc_ch","enc_ch",128,0,128,5000,0,5000);
  while(ch <= ch_max){
    for(int itemp=0;itemp<8;itemp++){
      if(itemp == 0){
	c1->Clear();
	c1->Divide(2,4);
      }
      c1->cd(itemp+1);
      TString hist_name = TString::Format("sigma_ch_%d",ch);
      sig_hist[ch]->Draw();
      sig_hist[ch]->SetTitle(hist_name);
      sig_hist[ch]->Fit("gaus","","",1,25);
      //if(ch!=6){
      TF1 *func3 = sig_hist[ch]->GetFunction("gaus");
      if ( func3 != NULL ) {
	double mean = func3->GetParameter(1);
	double enc = ENC(mean);
	//cout << enc << endl;
	enc_ch->Fill(ch,enc);
      }
	//}
      gStyle->SetOptStat(1);
      gStyle->SetOptFit(1);
      ch++;
      if(ch > ch_max) break;
    }
    c1->Draw();
    c1->Update();
    c1->Print(outputname,"PDF");
}

  c1->Clear();
  c1->cd(1);
  enc_ch->Draw();
  c1->Draw();  
  c1->Update();
  c1->Print(outputname, "PDF");
  c1->Print(outputname + "]","PDF");
  out.close();
}
