#define _GLIBCXX_USE_C99 1
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdint.h>
#include <unistd.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TProfile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TF1.h>
#include <TBrowser.h>

using std::ostringstream;
using namespace std;

void test_sts_disp_trim_dev(){
  
  
  int read_flag;
  
  TFile *file1;
  TTree *t_data;
  
  //std::string filename = "ftrim_171011_feb_b_module_build_asic_addr_0_5228188_vp_40_226.0_sts_elect";
  //std::string filename = "pscan_171012_setup2_feb-b_2_asic_addr_0_5326186_holes_temp__temp_24";
  std::string filename = "pscan_171018_setup1_asic_addr_0_5324186_holes";
  //std::string filename = "pscan_171009_setup2_feb-b_2_asic_addr_0_5326186_holes_config_normal_capboard_2_erni_1";
  std::string extension = ".root"; 
  std::string dir = "pscan_files/";
  std::string dir2 = "root_files/";
  
  const char *root_filename = (dir2+filename+extension).c_str();
  const char *datafile_1 =(dir+filename+".txt").c_str();
  const char *cap_file = (dir2+filename+".txto").c_str();
 // const char *datafile_1 ="pscan_files/pscan_170912_setup2_feb-b_8_asic_addr_0_5632186_sts_electrons.txt"
  //read_flag =1; 
  //file1 = new TFile(root_filename,"recreate"); 
  // -------------------------------------------------------------
  // 				Variables
  // -------------------------------------------------------------
 
  int ch;
  int ch_min = 0;
  int ch_max = 128;
  int ch_step =1;
  
  int d;
  int d_counter;
  int d_counter1;
  int d_min = 0;
  int d_max = 31;
  
  int grp;
  int grp_min = 0;
  int grp_max = 4;
  
  int vp;
  int ivp = 0;
  // input values for the vp scan. Please modify this values acording to the pscan vp range; values should be correct or s-curves will be shifted. *Needs improvements
  int vp_min = 0;
  int vp_max = 240;
  int vp_step = 1;
  uint32_t vcnt[130][32][250];   
  uint32_t cnt_val;
  
  // variables for data processing
  double mean;
  double sigma;
  double sigma_e;
  
  double sum_mean;
  double sum_delta;
  double sum_sig;
  double sum_sige;
  
  int d_cnt;
  int a;
  double vpe;
  

  
  //-------------------------------------------------------
  //  Channel and group for testing
  //-------------------------------------------------------
  int ch_sel = 54;						// Channel under test-> Display s-curves, ADC linearity and fitting, Disc. fitting
  grp = 2;							// Channel group to display. 
  int d_alt =20;						// Not used for the moment
  int ch_comp[10] = {6, 33, 79, 72, 75, 91, 94, 101, 110, 121};  	// 3 different channel connection type to compare. TYpe here the channel where the noise should be evaluated. 
  //int ch_comp[9] = {6, 20, 28, 34, 60, 76, 82, 92, 114};
  
  
  // ------------------------------------------------------------------
  // 		checking up if rootfile existspscan_171009_setup2_feb-b_2_asic_addr_0_5326186_holes_config_off_dis_capboard_1_erni_1
  // ------------------------------------------------------------------
 
  ifstream scanfile;
  scanfile.open(datafile_1);		// Input here the file name you want to read. 
  ifstream rootfile1;
  ofstream myfile;
  myfile.open(cap_file);
  
  //const char *filename1 = "pscan_170411_asic0_sensor_on_120_2.root";
  rootfile1.open(root_filename);
  
   if (!rootfile1.good()){
     rootfile1.close();
     file1 = new TFile(root_filename,"recreate"); 
     cout<<"Creating root_file with name: "<<root_filename<<endl;
//     t_data = new TTree("t_data","pscan_sensor_on_100");
     read_flag = 1;
   }
   else {
     //file1 = new TFile(root_filename,"update");; 
     cout<<"Reading root_file with name: "<<root_filename<<endl;
//     t_data = (TTree*)file1->Get("t_data");
     read_flag = 0;
//     rootfile1.close();
   }
  
  // --------------------------------------------------
  // 			Histograms
  // --------------------------------------------------
  
  TH1F *hdcnt[128][32];   // Differential count 
  TH1F *hscurve[128][32]; // S-Curves
  TH1F *hdisc_sig [128][32];  // Discriminators thr values for a specific channel 
  TH1F *hdisc_mean [128][32];  // Discriminators thr values for a specific channel 
 
  
  TH1F*h_mean  = new TH1F("hfit_m","hfit_m",128,0,128);
  TH1F*h_adc_mean = new TH1F("h_adc_mean","h_adc_mean",100,-3,0);
  TH1F*h_adc_q = new TH1F("h_adc_q","Chi-square",50,0,1);
  TH1F*h_sigma = new TH1F("hfit_s","hfit_s",128,0,128);
  TH1F *h_aux1 = new TH1F("h_aux1","h_aux1",31,0,31);
  TH1F *h_aux1_c = new TH1F("h_aux1_c","h_aux1_c",31,0,31);
  TH1F *h_aux2 = new TH1F("h_aux2","h_aux2",31,0,31);
  
  TH2F*hmeane  = new TH2F("hmeane", "hmean_calc", 128,0,128, 31,0,31);
  TH2F*hsige   = new TH2F("hsige", "hsig_calc", 128,0,128, 31,0,31);
    
  TH1F*hcalc_s  = new TH1F("hcalc_s","hcalc_s",128,0,128);
  TH1F*hfit_s  = new TH1F("hfit_s","hfit_s",128,0,128);
  TH2F*hmeanf  = new TH2F("hmeanf", "hmean_fit", 128,0,128, 31,0,31);
  TH2F*hsigef  = new TH2F("hsig_fit", "hsig_fit", 128,0,128, 31,0,31);
  
  TH1F*hfit_s2 = new TH1F("hfit_s2","hfit_s2",31,0,31); 
  TH1F*hcalc_s2 = new TH1F("hcalc_s2","hcalc_s2",31,0,31); 
  
  //-----------------------------------------
  //		fitting window
  //-----------------------------------------
  
  float amp_cal_min = 40.;   					// feedback from the get trim calibration process. Important to have these two values to make the proper cut around the peak (To avoid double peaks in the fitting)
  float amp_cal_max = 226.;  					// feedback from the get trim calibration process. These values can be taken from the amp_cal that was used in the get_trim process
  int thr_min =0;
  int thr_max =0;
  int thr_val = 0;						
  int thr_val1 =12;						// thr_value to select the range around the fitting peak. For ASIC alone (thr_val = 10, sigma ~2; so 5*sigma ~10); For ASIC + sensor (thr_val ~20,30 or 40 to cover whole S-curve)
  float vp_d = (amp_cal_max-amp_cal_min)/(d_max-d_min);
  
  int vp_set[31]= {0,};
  for (d = d_min; d<d_max;d++){
    vp_set[30-d] =(amp_cal_min + vp_d*(d-d_min));
  }
  
  int thr_min1 =0;
  int thr_max1 =0;

  
  Double_t f_s1mean = 0;
  Double_t f_mean = 0; 
  Double_t f_s1sigma = 0;
  Double_t f_sigma = 0;
  Double_t f_sigma2 = 0;
  sum_sige =0;
  // ----------------------------------------------------------------------
  //				Reading file
  // ----------------------------------------------------------------------

  if (read_flag == 1){
    
  std::cout<<"------------- Clearing counting array -----------"<<std::endl;
  
  for (ch=ch_min;ch<ch_max; ch++) {
    for (d =d_min; d<d_max; d++) {
      for (vp = vp_min;vp<vp_max; vp+=vp_step) {
	vcnt[ch][d][vp] = 0;
      }
    } 
  }
  
  int fvp;
  int fch;
  std::string line;
  std::string ele;
  
  std::cout<<"------------- Reading data file: -----------"<<std::endl;
  
  for (vp =vp_min; vp<vp_max; vp+=vp_step) {
    for (ch =ch_min;ch<ch_max; ch++){
      std::getline (scanfile,line);
      
      std::istringstream iss(line);
      iss>>ele;
      iss>>fvp;
      iss>>ele;
      iss>>ele;
      cout << "vp " << fvp << "   ch " << ch <<  ":  ";
      for (d=d_min;d<d_max;d++) {
	iss>>vcnt[ch][d][ivp];
	if (vcnt[ch][d][ivp] > 500){vcnt[ch][d][ivp] = 500;} 		// Comment this line to see the full s_curves; this is just for cutting double pulses after understand it. the value (i.e. 400) depends on the number of injected pulses. This is defined in the pscan file
	printf("%4d ", vcnt[ch][d][ivp]);
      }
      cout<<endl;  
    }
    ivp++;
  }
  scanfile.close();
 
  
  for (ch=ch_min;ch<ch_max;ch+=ch_step){
    //writefile<<ch<<'\t';
    thr_min,thr_max = (0,0);
    d_counter = 0;
    d_counter1 = 0;
    
    
    for (d=d_min;d <d_max;d++) {
      //cout <<"  ch: "<<ch<<"  disc: "<<d<< " ";
      char *histoname = new char[16];
      sprintf(histoname, "h_d_%d_%d",ch,d);
      char *histoscurve = new char[16];
      sprintf(histoscurve, "h_scurve_%d_%d",ch,d);
      char *histodisc = new char[16];
      sprintf(histodisc, "h_disc_%d",d);
      hdcnt[ch][d]=new TH1F(histoname,"",300,0,300);
      hscurve[ch][d]=new TH1F(histoscurve,"",300,0,300);
      hdisc_sig[ch][d] = new TH1F(histodisc,"",31, 0,31);
      hdisc_mean[ch][d] = new TH1F(histodisc,"",31,0,31);
      ivp = 1;
      a = 0;
      sum_mean = 0;
      sum_delta = 0.000001;
      
      for ( vp = vp_min +vp_step; vp <vp_max; vp += vp_step ) {
	d_cnt = vcnt[ch][d][ivp]-vcnt[ch][d][ivp-1];
	sum_delta += d_cnt;
	//sum_mean += vp*d_cnt;
	//printf("%5d", vcnt[ch][d][ivp]);
	vpe = (vp-0.5*vp_step) * 15. / 256 * 6258;
	hdcnt[ch][d]->Fill(vp,d_cnt);
	hscurve[ch][d]->Fill(vp,vcnt[ch][d][ivp]);
	ivp++;
      }

    

      //cout<<endl;
     
      // ------------ Fitting. Sigma 1. Fitting curves ---------------- //
      // ------------Fitting ranges for disc. -------------------------//
      thr_val = thr_val1;
      thr_min = (int)(vp_set[d]-thr_val); 
      thr_max = (int)(vp_set[d]+thr_val); 

      
      fit_state0:
      TFitResultPtr f_s1; 
      f_s1= hdcnt[ch][d]->Fit("gaus","SR","",thr_min,thr_max); // fitting disc in the range thr_min-thr_max
      //f_s1= hdcnt[ch][d]->Fit("gaus","SW"); // fitting disc in the range thr_min-thr_max
      Int_t fitstatus = f_s1;
      if ( fitstatus == 0 && hdcnt[ch][d]->GetFunction("gaus")->GetChisquare()!=0) {
	
	f_s1mean+= f_s1->Parameter(1);
	f_sigma2 = f_s1->Parameter(2);
	f_mean = f_s1->Parameter(1);
	f_sigma = f_s1->Parameter(2);
	//writefile<<f_mean<<'\t';
	
	if (f_mean !=0 ){
	  hdisc_sig[ch][d]->Fill(d,f_sigma2);
	  hdisc_mean[ch][d]->Fill(d,f_mean);
	  if (ch == ch_sel) {
	    hfit_s2->Fill(d,f_sigma2);
	    h_aux1->Fill(d,f_mean);
	    //h_aux2->Fill(d,mean);
	  }
	}
	else {

	  while (thr_val > 5){
	    thr_val = thr_val-5;
	    thr_min = (int)(vp_set[d]-thr_val); 
	    thr_max = (int)(vp_set[d]+thr_val); 
	    goto fit_state0;
	  }
	}  
	
	if (f_s1->Parameter(2) >0){
	//if (ch<64 && ch%2 ==0 && d>6 && f_s1->Parameter(2) >0){
	  f_s1sigma+= f_s1->Parameter(2);
	  d_counter++; 
	}
	/*else ()
	{
	  f_s1sigma+= f_s1->Parameter(2);
	  d_counter++;
	} */

      }
      
      else {
	  cout << "fit not ok" << endl;
	 // writefile<<0<<'\t';
      }     
      cout<<endl; 
      
    
    hmeanf->Fill(ch,d,f_mean);
    hsigef->Fill(ch,d,f_sigma);
    f_s1mean = f_s1mean/d_max;	
 
    
    // --------------------- Calculating sigma --------------------------//
    thr_val = thr_val1;
    sum_delta = 0.0000001;
    sum_sig = 0;
    d_cnt = 0;
    
   // range where I will look for the Scurve. Expanded range compare to the fit
    
   thr_min1 = (int)(vp_set[d]-vp_min-thr_val); 
   if (thr_min1 <=0){thr_min1 = 1;}
   thr_max1 = (int)(vp_set[d]-vp_min+thr_val);
   if (thr_max1 >vp_max-vp_min){thr_max1 = vp_max-vp_min;} 
   
//    if (ch == 10){
//      myfile<<"disc: "<<d<<"\t"<<d<<"\t";
//      myfile<<thr_min1<<"\t"<<vp_set[d]<<thr_max1<<endl;
//   }
//   myfile.close();

  ivp = thr_min1;
   
   for ( vp = thr_min1; vp <thr_max1; vp += vp_step ) {
	d_cnt = vcnt[ch][d][ivp]-vcnt[ch][d][ivp-1]; 
	sum_delta += d_cnt;
	sum_mean += (vp+vp_min)*d_cnt;	
	ivp+=1;
	
   }
      mean = sum_mean/sum_delta;
      if ((mean<0.)||(mean>250.)) {
	hmeane->Fill(ch,d,1);
      }
      else {
	hmeane->Fill(ch,d,mean);
      }
   sum_delta = 0.0000001;
   
   ivp = thr_min1;
   
   for ( vp = thr_min1; vp <thr_max1; vp += vp_step ) {
      d_cnt = vcnt[ch][d][ivp]-vcnt[ch][d][ivp-1];
      sum_delta += d_cnt;
      sum_sig +=((vp+vp_min)-mean)*((vp+vp_min)-mean)*d_cnt;
      ivp++;
   }
   sigma = sqrt(sum_sig/sum_delta);
   
   
   if (ch == ch_sel) {
    h_aux1_c->Fill(d,mean);
    hcalc_s2->Fill(d,sigma);
   }
   
  // if (mean >0 && sigma >0 && d>6){
   if (sigma >0){
    sum_sige +=sigma;
    d_counter1++;
    }
   hsige->Fill(ch,d,sigma); 
  } 

  //-----------------------------------------------------------------------//
   
   if (d_counter !=0){f_s1sigma = f_s1sigma/d_counter;}
   else {f_s1sigma = f_s1sigma/1;}
   hfit_s->Fill(ch,f_s1sigma);
   
   // writing file for analyzing capacitance-noise effect
    if (ch%2 ==0){
      myfile<<ch<<"\t"<<f_s1sigma<<endl;
    }
   
      
   if (d_counter1 !=0){sum_sige = sum_sige/d_counter1;}
   else {sum_sige = sum_sige/1;}
   hcalc_s->Fill(ch,sum_sige);
 
   
  // writefile<<endl;
  }
  
  myfile.close();
  
  
// -------------------------------------------------
//		Printing out histograms 
// -------------------------------------------------
  TCanvas *c1 =new TCanvas("c1","c1");
  c1->cd();
  thr_min = (int)(vp_set[0]-thr_val); 
  thr_max = (int)(vp_set[0]+thr_val);
  hdcnt[ch_sel][0]->Draw("");
  hdcnt[ch_sel][0]->Fit("gaus","SR","", thr_min,thr_max);
  
 
  int n =2;
  TCanvas *c2 =new TCanvas("c2","c2");
  c2->Divide(6,6);
  c2->cd(1);
  hdcnt[ch_sel][30]->Draw("");
  hdcnt[ch_sel][30]->Fit("gaus","SWR","",0,100);
  for (d=29; d>=d_min;d--){
    c2->cd(n);
    hdcnt[ch_sel][d]->Draw("");
    thr_min = (int)(vp_set[d]-thr_val); 
    thr_max = (int)(vp_set[d]+thr_val); 
    n++;
  }
  
 
 n=0;
  gStyle->SetOptStat(0);
  TCanvas *c3 =new TCanvas("c3","S-Curves_selected channel",1400,600);
  c3->Divide(2,1);
  c3->cd(1)->SetGrid();
  hscurve[ch_sel][30]->Draw("");
  hscurve[ch_sel][30]->GetXaxis()->SetTitle("Pulse amplitude [au]");
  hscurve[ch_sel][30]->GetYaxis()->SetTitle("Number of counts");
  hscurve[ch_sel][30]->SetLineWidth(2);
   for (d=d_min+1;d <31;d++) {
	n = 30-d;
        hscurve[ch_sel][n]->Draw("SAME");
	//hscurve[ch_sel][d]->SetLineColor(n);
	hscurve[ch_sel][n]->SetLineWidth(2);
    } 
 
  c3->cd(2)->SetGrid();

  TPad *pad1 = new TPad("pad1","pad1",0,0.30,1,0.95);
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.30);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(4.0);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd()->SetGrid();
  h_aux1->Draw("P");
  gStyle->SetErrorX(0);
  h_aux1->Fit("pol1");
  h_aux1->SetMarkerStyle(kFullCircle);
  h_aux1->GetXaxis()->SetTitle("Discriminator Number");
  h_aux1->GetYaxis()->SetTitle("Discriminator Thr [au]"); 
  
  float par0 = h_aux1->GetFunction("pol1")->GetParameter(0);
  float par1 = h_aux1->GetFunction("pol1")->GetParameter(1);
  float residuals = 0;
  
  pad2->cd()->SetGrid();
  for(d=d_min; d<d_max; d++){
    residuals = par0+par1*(d+0.5) - h_aux1->GetBinContent(d+1);
    h_aux2->Fill(d,residuals);
  } 
  h_aux2->Draw("PE");
  h_aux2->SetMarkerStyle(kFullCircle);
  h_aux2->GetXaxis()->SetLabelSize(0.06);
  h_aux2->GetYaxis()->SetLabelSize(0.06);
  h_aux2->GetYaxis()->SetTitle(""); 
  h_aux2->GetXaxis()->SetTitle("Discriminator number"); 
  h_aux2->GetXaxis()->SetTitleSize(0.06);
  h_aux2->GetYaxis()->SetTitle("Residuals");
  h_aux2->GetYaxis()->SetTitleSize(0.06);
  
  c3->Write();


   int i=1;
   n = 0;
   TCanvas *c4 =new TCanvas("c4","S-Curves_selected group");
   c4->Divide(6,5); 
   for (ch=grp;ch<ch_max;ch+=4){
     c4->cd(i);
     hscurve[ch][0]->Draw("HIST");
     hscurve[ch][0]->GetXaxis()->SetTitle("Pulse amplitude [au]");
     hscurve[ch][0]->GetYaxis()->SetTitle("Number of counts");
     hscurve[ch][0]->GetYaxis()->SetRangeUser(0,3000);
     for (d=d_min+1;d <31;d++) {
	n = 30-d;
	hscurve[ch][d]->Draw("SAME");
	//hscurve[ch][d]->SetLineColor(d);
	}
     i++;
   }
   
   
  TCanvas *c5 =new TCanvas("c5","Mean_and sigma_channels");
  c5->Divide(3,2);
  c5->cd(1);
  hmeanf->Draw("COLZ");
  hmeanf->GetZaxis()->SetRangeUser(0,255);

  c5->cd(2);
  hmeane->Draw("COLZ");
  hmeane->GetZaxis()->SetRangeUser(0,255);
  
  c5->cd(3)->SetGrid();
  ch =0;
  TH1D *py[128];
  for (int j = 1; j<=128;j++){
    ch =j-1;
    char *h_adc_linearity = new char[16];
    sprintf(h_adc_linearity, "h_quality_%d",ch);
    py[ch] = hmeanf->ProjectionY(h_adc_linearity,j,j);
  }
  py[0]->Draw("");
  py[0]->GetXaxis()->SetTitle("ADC value");
  py[0]->GetYaxis()->SetTitle("Signal Vp amplitude [au]");
  for (ch = ch_min+2; ch<ch_max;ch+=ch_step){
    py[ch]->Draw("same");
    //py[ch]->SetLineColor(ch);
    
  }
  
  c5->cd(4);
  hsigef->Draw("COLZ");
  hsigef->GetZaxis()->SetRangeUser(0,5);
  
  c5->cd(5);
  hsige->Draw("COLZ");
  hsige->GetZaxis()->SetRangeUser(0,5);
  
  c5->cd(6)->SetGrid();
  ch =0;
  for (int ch =ch_min; ch<ch_max;ch+=ch_step){
    py[ch]->Fit("pol1","R","",0,31);
    float par1_adc = py[ch]->GetFunction("pol1")->GetParameter(1);
    float par2_adc = py[ch]->GetFunction("pol1")->GetChisquare();
    h_adc_q->Fill(par2_adc);
    h_adc_mean->Fill(par1_adc);  
   
  }
  h_adc_mean->Draw("");
  h_adc_q->Draw("");
  h_adc_q->SetFillColor(kGreen-6);
  h_adc_q->SetLineColor(kGreen-6);
  h_adc_q->Fit("gaus");
  
  c5->Write();
  

  

  n=0;
  TCanvas *c6 =new TCanvas("c6","Noise level/discriminator");
  c6->Divide(2,1);
  c6->cd(1)->SetGrid();
  h_aux1->Draw("PTEXT");
  h_aux1->Fit("pol1");
  h_aux1->SetMarkerColor(kRed);
  h_aux1->GetXaxis()->SetTitle("Discriminator Number");
  h_aux1->GetYaxis()->SetTitle("Discriminator Thr [au]");
  h_aux1_c->Draw("PSAMETEXT");
  h_aux1_c->Fit("pol1");
  h_aux1_c->SetMarkerColor(kBlue);
  h_aux1_c->SetMarkerStyle(kFullCircle);
    
  c6->cd(2)->SetGrid();
   hfit_s2->Draw("HISTTEXT");
   hfit_s2->SetLineColor(kRed);
   hfit_s2->GetXaxis()->SetTitle("Discriminator Number");
   hfit_s2->GetYaxis()->SetTitle("Sigma [au]");
   hcalc_s2->Draw("HISTTEXTsame");
   hcalc_s2->SetLineColor(kBlue);
   
   
   TCanvas*c_comp = new TCanvas("c_comp","c_comp");
   c_comp->Divide(3,3);
   
   for (int j=0; j<3; j++){
    c_comp->cd(j+1)->SetGrid();
    hscurve[ch_comp[j]][30]->Draw("");
    hscurve[ch_comp[j]][30]->GetXaxis()->SetTitle("Pulse amplitude [au]");
    hscurve[ch_comp[j]][30]->GetYaxis()->SetTitle("Number of counts");
    hscurve[ch_comp[j]][30]->SetLineWidth(2);
    for (d=d_min+1;d <31;d++) {
      n = 30-d;
      hscurve[ch_comp[j]][n]->Draw("SAME");
      //hscurve[ch_sel][d]->SetLineColor(n);
      hscurve[ch_comp[j]][n]->SetLineWidth(2);
    } 
    int k = j+4;
    c_comp->cd(k)->SetGrid();
    hdisc_mean[ch_comp[j]][0]->Draw("P");
    hdisc_mean[ch_comp[j]][0]->SetMarkerStyle(8);
    hdisc_mean[ch_comp[j]][0]->SetMarkerSize(0.7);
    hdisc_mean[ch_comp[j]][0]->SetMarkerColor(kRed);
    hdisc_mean[ch_comp[j]][0]->GetXaxis()->SetTitle("Discriminator number");
    hdisc_mean[ch_comp[j]][0]->GetYaxis()->SetTitle("Disc. thr [au]");
    hdisc_mean[ch_comp[j]][0]->GetYaxis()->SetRangeUser(0,260);
    hdisc_mean[ch_comp[j]][0]->SetLineWidth(2);
    for (d=d_min+1;d <31;d++) {
      n = 30-d;
      hdisc_mean[ch_comp[j]][d]->Draw("PSAME");
      //hdisc_sig[ch_sel][d]->SetLineColor(n);
      hdisc_mean[ch_comp[j]][d]->SetLineWidth(2);
      hdisc_mean[ch_comp[j]][d]->SetMarkerStyle(8);
      hdisc_mean[ch_comp[j]][d]->SetMarkerSize(0.7);
      hdisc_mean[ch_comp[j]][d]->SetMarkerColor(kRed);
    }  
    
    int z = j+7;
    c_comp->cd(z)->SetGrid();
    hdisc_sig[ch_comp[j]][0]->Draw("P");
    hdisc_sig[ch_comp[j]][0]->SetMarkerStyle(8);
    hdisc_sig[ch_comp[j]][0]->SetMarkerSize(0.7);
    hdisc_sig[ch_comp[j]][0]->GetXaxis()->SetTitle("Discriminator number");
    hdisc_sig[ch_comp[j]][0]->GetYaxis()->SetTitle("ENC [au]");
    hdisc_sig[ch_comp[j]][0]->GetYaxis()->SetRangeUser(0,10);
    hdisc_sig[ch_comp[j]][0]->SetLineWidth(2);
    for (d=d_min+1;d <31;d++) {
      n = 30-d;
      hdisc_sig[ch_comp[j]][d]->Draw("PSAME");
      //hdisc_sig[ch_sel][d]->SetLineColor(n);
      hdisc_sig[ch_comp[j]][d]->SetLineWidth(2);
      hdisc_sig[ch_comp[j]][d]->SetMarkerStyle(8);
      hdisc_sig[ch_comp[j]][d]->SetMarkerSize(0.7);
    }  
  }

  TCanvas *c_noise = new TCanvas("c_noise","c_noise");
  c_noise->cd()->SetGrid();
  hfit_s->Draw("");
  hfit_s->GetXaxis()->SetTitle("Channel number");
  hfit_s->GetYaxis()->SetTitle("ENC [au]");
  hfit_s->SetLineWidth(2);
  hfit_s->SetLineColor(kRed);
  //hfit_s->GetYaxis()->SetRangeUser(0,20);
  hcalc_s->Draw("same");
  hcalc_s->SetLineWidth(2);
  hcalc_s->SetLineColor(kBlue);
  //hfit_s->SaveAs("noise.pdf");
  hfit_s->Write();
  hcalc_s->Write();
  
  TLegend *l_noise = new TLegend(0.5,0.7,0.88,0.88);
  l_noise->SetHeader(""); // option "C" allows to center the header
  l_noise->SetLineColor(kWhite);
  l_noise->AddEntry(hfit_s,"Noise_fitting","l");
  l_noise->AddEntry(hcalc_s,"Noise_calculated_RMS","l");
  l_noise->Draw();
  
  
  
  
  cout<<"Channels noise comparison"<<endl;
  cout<<"-----------------------------------"<<endl;
    for (int j =0; j<10; j++){
    cout<<"ch:"<<"\t"<<ch_comp[j]<<"\t"<<hfit_s->GetBinContent(ch_comp[j]+1)<<endl;
  }
  cout<<"-----------------------------------"<<endl;

  cout<<" "<<endl;
  cout<<"-------------------------------"<<endl;
  cout<< "disc"<<"\t"<<"thr_min"<<"\t"<<"vp_set[d]"<<"\t"<<"thr_max"<<endl;
  for (d=d_min; d<d_max; d++){
  thr_min = (int)(vp_set[d]-thr_val); 
  thr_max = (int)(vp_set[d]+thr_val);
   
  cout<<d<<"\t"<< thr_min<<"\t"<<vp_set[d]<<"\t"<<thr_max<<endl;
  
  }
  
  //file1->Close("");
  }
  
  
  if (read_flag == 0){
    TBrowser*my_browser = new TBrowser("my_browser","");
    //file1->Draw("");   
  }
 
}
