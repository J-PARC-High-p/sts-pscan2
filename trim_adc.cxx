#include "trim_adc.hxx"
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
#include "TObject.h"

#define KAZ_DEBUG 1

ClassImp(trim_adc);

//------------------------------------------+-----------------------------------
//! Default constructor.
trim_adc::trim_adc(TString filenameData)
: ch(0),
  ch_min(0),
  ch_max(127),
  ch_step(1),
  d(0),
  d_counter(0),
  d_counter1(0),
  d_min(0),
  d_max(31),
  d_step(0),
  grp(0),
  grp_min(0),
  grp_max(4),
  grp_step(0),
  vp(0),
  ivp(0),
  vp_min(0),
  vp_max(250),
  cut_db_pulses(100),
  vp_step(1),   
  cnt_val(0),
  mean(0),
  sigma(0),
  sigma_e(0), 
  sum_mean(0),
  sum_delta(0.0001),
  sum_sig(0),
  sum_sige(0), 
  d_cnt(0),
  a(0),
  vpe(0),
  thr_val1(0),
  amp_cal_max(0),
  amp_cal_min(0),
  filename_data(filenameData),
  read_flag(-1),
  soft_flag(0),
  f_s1mean(0),
  f_mean(0),
  f_mean_d1(0),
  f_mean_d20(0),
  f_mean_d30(0),
  f_s1sigma(0),
  f_sigma(0),
  f_sigma2(0),
  f_figma_fast(0),
  f_mean_fast(0),
  // ----- extras ----
  f_s1mean_erfc(0),
  f_mean_erfc(0),
  f_mean_d1_erfc(0),
  f_mean_fast_erfc(0),
  f_s1sigma_erfc(0),
  f_sigma_erfc(0),
  f_sigma2_erfc(0),
  
  thr_min(0),
  thr_max(0),
  ch_sel(0)
{
 for (ch=ch_min;ch<=ch_max; ch++) {
    for (d =d_min; d<d_max+1; d++) {
      for (vp = vp_min;vp<vp_max; vp+=vp_step) {
	vcnt[ch][d][vp] = 0;
	vcnt_soft[ch][d][vp] = 0;
      }
    } 
  }
  for (d=d_min;d<=d_max;d++){
     vp_set[d]=0;
  }
}

//------------------------------------------+-----------------------------------
//! Destructor.
trim_adc::~trim_adc()
{

}
//------------------------------------------+-----------------------------------
//! Open and check root file
bool trim_adc::Check_root_file()
{
   
  filename_root = (filename_data+ ".root");
  filename_val = (filename_data+ ".txto");
  filename_data = (filename_data+".txt");
  scanfile.open(filename_data);		// Input here the file name you want to read. 
  if (! scanfile.is_open() ) {
    cout << " cannot open file " << filename_data << endl;
    return false;
  }
  ifstream rootfile;
  
  
  rootfile.open(filename_root);
  
   if (!rootfile.good()){
     rootfile.close();
     file1 = new TFile(filename_root,"recreate"); 
     cout<<"Creating root_file with name: "<<filename_root<<endl;
//     t_data = new TTree("t_data","pscan_sensor_on_100");
     read_flag = 1;
   }
   else {
     file1 = new TFile(filename_root,"update");; 
     cout<<"Reading root_file with name: "<<filename_root<<endl;
//     t_data = (TTree*)file1->Get("t_data");
     read_flag = 0;
   }
   
   if (read_flag ==-1) return false;
   return true;
}

//------------------------------------------+-----------------------------------
//! Close root_file
void trim_adc::Close_root_file()
{
   file1->Close();
   scanfile.close();
}
//------------------------------------------+-----------------------------------
//! Create histograms
void trim_adc::Create_histo()
{
  h_mean  = new TH1F("hfit_m","hfit_m",128,0,128);
  h_adc_mean = new TH1F("h_adc_mean","h_adc_mean",1000,-10,10);
  h_adc_slope = new TH1F("h_adc_gain_spread","h_adc_gain_spread",300,0,0.6);
  h_adc_slope_ch = new TH1F("h_adc_gain_ch","h_adc_gain_ch",128,0,128);
  h_adc_int = new TH1F("h_adc_int","Intercept",500,150,250);
  h_sigma = new TH1F("hfit_s","hfit_s",128,0,128);
  h_aux1 = new TH1F("h_aux1","h_aux1",31,0,31);
  h_aux1_erfc = new TH1F("h_aux1_erfc","h_aux1_erfc",31,0,31);
  h_aux1_c = new TH1F("h_aux1_c","h_aux1_c",31,0,31);
  h_aux2 = new TH1F("h_aux2","h_aux2",31,0,31);
  
  
  hmeane  = new TH2F("hmeane", "hmean_calc", 128,0,128, 31,0,31);
  hsige   = new TH2F("hsige", "hsig_calc", 128,0,128, 31,0,31);
  
  hmeane_erfc  = new TH2F("hmeane_erfc", "hmeane_erfc", 128,0,128, 31,0,31);
  hsige_erfc   = new TH2F("hsige_erfc", "hsige_erfc", 128,0,128, 31,0,31);
    
  hcalc_s  = new TH1F("hcalc_s","hcalc_s",128,0,128);
  hfit_s  = new TH1F("hfit_s","hfit_s",128,0,128);
  hfit_s_erfc = new TH1F("hfit_s_erfc","hfit_s_erfc",128,0,128);
  
  hmean_f_disc_30  = new TH1F("hmean_f_disc_30","hmean_f_disc_30",128,0,128);
  hfit_s_disc_30  = new TH1F("hfit_s_disc_30","hfit_s_disc_30",128,0,128);
  hmeanf  = new TH2F("hmeanf", "hmean_fit", 128,0,128, 31,0,31);
  hsigef  = new TH2F("hsig_fit", "hsig_fit", 128,0,128, 31,0,31);
  
  hfit_s2 = new TH1F("hfit_s2","hfit_s2",31,0,31); 
  hfit_s2_erfc = new TH1F("hfit_s2_erfc","hfit_s2_erfc",31,0,31); 
  hcalc_s2 = new TH1F("hcalc_s2","hcalc_s2",31,0,31); 
  
  hnoise_fast = new TH1F("hnoise_fast", "hnoise_fast", 128,0,128);
  hmean_fast = new TH1F("hmean_fast","hmean_fast",128,0,128);
  h_diff_fast_adc = new TH1F("h_diff_fast_adc","h_diff_fast_adc",200, -20,20);
  hmean_fast_disc = new TH1F("hmean_fast_disc","hmean_fast_disc",100, 0,50);

  
}

//------------------------------------------+-----------------------------------
//! Initilaizing ch

bool trim_adc::Init_ch(int chMin, int chMax, int chStep)
{
  ch_min = chMin;
  ch_max = chMax;
  if (chStep !=1) ch_step = chStep;
  if (chMin<0 | chMin>chMax | chMax>=128) return false;
  return true;
}

//------------------------------------------+-----------------------------------
//! Initilizing grp

bool trim_adc::Init_grp(int grpMin, int grpMax, int grpStep)
{
  grp_min = grpMin;
  grp_max = grpMax;
  if (grpStep !=1) grp_step = grpStep;
  if (grpMin<0 | grpMin>grpMax | grpMax>4) return false;
  return true;
}

//------------------------------------------+-----------------------------------
//! Initializing disc

bool trim_adc::Init_d(int dMin, int dMax, int dStep)
{
  d_min = dMin;
  d_max = dMax;
  if (dStep !=1) d_step = dStep;
  if (dMin<0 | dMin>dMax | dMax>31) return false;
  return true;
}

//------------------------------------------+-----------------------------------
//! Initializing vp

bool trim_adc::Init_vp(int vpMin, int vpMax, int vpStep)
{

  vp_min = vpMin;
  vp_max = vpMax;
  if (vpStep !=1) vp_step = vpStep;
  if (vpMin<0 | vpMin>vpMax | vpMax>255) return false;
  return true;

}
//------------------------------------------+-----------------------------------
//! Initializing fitting windows

bool trim_adc::Fitting_windows(int ampcalMin, int ampcalMax, int width)
{

  thr_val1=width;
  amp_cal_max= ampcalMax;
  amp_cal_min=ampcalMin;
 					// feedback from the get trim calibration process. These values can be taken from the amp_cal that was used in the get_trim process
  int thr_min =0;
  int thr_max =0;
  int thr_val = 0;							// thr_value to select the range around the fitting peak. For ASIC alone (thr_val = 10, sigma ~2; so 5*sigma ~10); For ASIC + sensor (thr_val ~20,30 or 40 to cover whole S-curve)
  float vp_d = (amp_cal_max-amp_cal_min)/(d_max-d_min);
  
  for (d = d_min; d<=d_max;d++){
    vp_set[30-d] =(amp_cal_min + vp_d*(d-d_min));
  }
  return true;
}

//------------------------------------------+-----------------------------------
//! Reading data file
bool trim_adc::Reading_file(int cut_db_pulses_user)
{
  cut_db_pulses = cut_db_pulses_user;
  if (read_flag != 1) return false;

  int fvp;
  int fch;
  std::string line;
  std::string ele;
  
  std::cout<<"------------- Reading data file: -----------"<<std::endl;
  
  for (vp =vp_min; vp<vp_max; vp+=vp_step) {
    for (ch =ch_min;ch<=ch_max; ch++){
      if ( ! std::getline (scanfile,line) ) {
	cout << "cannot read lines anymore " << endl;
	return false;
      }
      
      std::istringstream iss(line);
      //      cout << line << endl;
      iss>>ele;
      iss>>fvp;
      iss>>ele;
      iss>>ele;
      cout << "vp " << fvp << "   ch " << ch <<  ":  ";
      for (d=d_min;d<=d_max;d++) {
	//for (d=d_min;d<d_max;d++){
	iss>>vcnt[ch][d][ivp];
	if (vcnt[ch][d][ivp] > cut_db_pulses) vcnt[ch][d][ivp] = cut_db_pulses;		// Comment this line to see the full s_curves; this is just for cutting double pulses after understand it. the value (i.e. 400) depends on the number of injected pulses. This is defined in the pscan file
	printf("%4d ", vcnt[ch][d][ivp]);
      }
      cout<<endl;  
    }
    ivp++;
  }
  std::cout<<"------------- Finished reading data file: -----------"<<std::endl;
  return true;
}


//------------------------------------------+-----------------------------------
//! Analyzing the data
void trim_adc::Analysis(int cut_db_pulses_user, int width_user, int dcut_min_user,int dcut_max_user, int ch_sel_user, bool fit, bool calc, bool fast_fit)
{
#ifdef KAZ_DEBUG
  cout << "The begining of trim_adc::Analysis" << endl;
#endif
  ch_sel = ch_sel_user;
//   val_file.open(filename_val);
  
  Soft_val(true);
  
  cut_db_pulses = cut_db_pulses_user;
  for (ch=ch_min;ch<=ch_max;ch+=ch_step){
    //val_file << ch << "\t\t";
    d_counter = 0;
    d_counter1 = 0;
    d_counter_erfc = 0;
    
    
    for (d=d_min;d <=d_max;d++) {
      TString name;
      name = TString::Format("h_d_%d_%d",ch,d);
      hdcnt[ch][d]=new TH1F(name,"",300,0,300);
      
      name = TString::Format("h_scurve_%d_%d",ch,d);
      hscurve[ch][d]=new TH1F(name,"",300,0,300);
      
      name = TString::Format("h_disc_sig_%d_%d",ch,d);
      hdisc_sig[ch][d] = new TH1F(name,"",31, 0,31);

      name = TString::Format("h_disc_mean_%d_%d",ch,d);
      hdisc_mean[ch][d] = new TH1F(name,"",31,0,31);
      ivp = 1;
      a = 0;
      sum_mean = 0;
      sum_delta = 0.000001;
      
      for ( vp = vp_min +vp_step; vp <vp_max; vp += vp_step ) {
#ifdef KAZ_DEBUG
	//	cout << "DEBUG vp = " << vp << endl;
#endif
	d_cnt = vcnt_soft[ch][d][ivp]-vcnt_soft[ch][d][ivp-1];
	//if (d_cnt <0 | d_cnt > 35) d_cnt =0;
	if (d_cnt <0) d_cnt =0;
	//sum_delta += d_cnt;
	//sum_mean += vp*d_cnt;
	//printf("%5d", vcnt[ch][d][ivp]);
	vpe = (vp-0.5*vp_step) * 15. / 256 * 6258;
	//if (abs(d_cnt) > 30) d_cnt = 5;
	hdcnt[ch][d]->Fill(vp,d_cnt);
	hscurve[ch][d]->Fill(vp,vcnt_soft[ch][d][ivp]);
	ivp++;
      }
      
      
      // for debbuging...
      
      if (ch==ch_sel){
	hscurve[ch][d]->Write();
	hdcnt[ch][d]->Write();
	
      }
      
      //if (fit == true && d<31) Fit_values(width_user);
      //if (calc == true && d<31) Calc_values(width_user);
      if (fit == true && d<=30) {
	Fit_values(width_user, dcut_min_user, dcut_max_user); 
	Fit_values_erfc(width_user, dcut_min_user, dcut_max_user, cut_db_pulses);}
      if (calc == true && d<=30) Calc_values(width_user, dcut_min_user, dcut_max_user);
      
      if (fast_fit == true && d==31) Fitting_Fast(width_user);
      
    }
    
   if (d_counter !=0){f_s1sigma = f_s1sigma/d_counter;}
   else {f_s1sigma = f_s1sigma/1;}
   hfit_s->Fill(ch,f_s1sigma*350);
   
   if (d_counter_erfc !=0){f_s1sigma_erfc = f_s1sigma_erfc/d_counter_erfc;}
   else {f_s1sigma_erfc = f_s1sigma_erfc/1;}
   hfit_s_erfc->Fill(ch,f_s1sigma_erfc*350);
   
   
   if (d_counter1 !=0){sum_sige = sum_sige/d_counter1;}
   else {sum_sige = sum_sige/1;}
   hcalc_s->Fill(ch,sum_sige*350);
   
   // writing file for analyzing capacitance-noise effect
//    if (ch%2 ==0 && ch<64 ){
//      val_file << ch << "\t\t"<<f_s1sigma<<"\t\t"<<f_mean_d1<<"\t\t"<<f_mean_d20<<"\t\t"<<f_mean_fast<<"\t\t"<<f_figma_fast/350.<<endl;
//    }
   
     //val_file << ch <<"\t\t"<<f_mean_d1<<"\t\t"<<f_mean_d30<< "\t\t"<<f_s1sigma*350<<endl;
      //val_file << ch <<"\t"<<f_s1sigma*350<<endl;
      //val_file << ch <<"\t"<<f_s1sigma*350<<endl;
     
     
 }
 
//  val_file.close();
  
}



//------------------------------------------+-----------------------------------
//! Analyzing the data. Fitting S-curves
bool trim_adc::Soft_val(bool soft_flag)
{
  if (soft_flag==true){
    for (ch=ch_min;ch<=ch_max;ch+=ch_step){
      //for (d=d_min;d <d_max;d++){
      for (d=d_min;d <d_max+1;d++){
	 ivp = 2;
	 int soft_count =0;
	 int d_cnt0 = 0;
	 int d_cnt1 = 0;
	 
	 vcnt_soft[ch][d][0]= vcnt[ch][d][0];
	 vcnt_soft[ch][d][1]= vcnt[ch][d][1];
	 vcnt_soft[ch][d][vp_max-vp_min-1]= vcnt[ch][d][vp_max-vp_min-1];
	 vcnt_soft[ch][d][vp_max-vp_min-2]= vcnt[ch][d][vp_max-vp_min-2];
	 for ( vp = vp_min +2*vp_step; vp <vp_max-2*vp_step; vp += vp_step ){
	   d_cnt0 = vcnt[ch][d][ivp]-vcnt[ch][d][ivp-1];
	   d_cnt1 = vcnt[ch][d][ivp+1]-vcnt[ch][d][ivp];
	   if (d_cnt0 < -10 &&  d_cnt1 > 10) {vcnt_soft[ch][d][ivp] = int((vcnt[ch][d][ivp-1]+vcnt[ch][d][ivp+1])/2.); }
	   else if (d_cnt0 >10 &&  d_cnt1 <-10) {vcnt_soft[ch][d][ivp] = int((vcnt[ch][d][ivp-1]+vcnt[ch][d][ivp+1])/2.);}
	   else {
	     vcnt_soft[ch][d][ivp] = vcnt[ch][d][ivp];
	     
	   //vcnt_soft[ch][d][ivp] = int((vcnt[ch][d][ivp-1]+vcnt[ch][d][ivp]+vcnt[ch][d][ivp+1])/3);
	   //vcnt_soft[ch][d][ivp] = int((vcnt[ch][d][ivp-2]+vcnt[ch][d][ivp-1]+vcnt[ch][d][ivp]+vcnt[ch][d][ivp+1]+vcnt[ch][d][ivp+2])/5);
	   //vcnt_soft[ch][d][ivp] = int((vcnt[ch][d][ivp-1]+vcnt[ch][d][ivp])/2);
	  }
	  ivp++;
	}
      }
    }
  }
  
  if (soft_flag == false){
    for (ch=ch_min;ch<=ch_max;ch+=ch_step){
       //for (d=d_min;d <d_max;d++){
	for (d=d_min;d <d_max+1;d++){
	 ivp = 0;
	 for ( vp = vp_min; vp <vp_max; vp += vp_step ){
	   vcnt_soft[ch][d][ivp] = vcnt[ch][d][ivp];
	  ivp++;
	}
      }
    }
    
  }
  
  return true;
}




//------------------------------------------+-----------------------------------
//! Analyzing the data. Fitting S-curves
bool trim_adc::Fit_values(int width, int d_cut_min, int d_cut_max)
{
  d_cut_min = 30-d_cut_min;
  d_cut_max = 30-d_cut_max;
  
  if (d_cut_min < 0 | d_cut_max > 30) return false;
  else { 
    if (d>=d_cut_max && d<=d_cut_min) {
      int binmax = hdcnt[ch][d]->GetMaximumBin();
      float x =hdcnt[ch][d]->GetXaxis()->GetBinCenter(binmax); 
      
      /*
	thr_min = (int)(vp_set[d]-width); 
	thr_max = (int)(vp_set[d]+width); 
      */
      
      thr_min = (int)(x-width); 
      thr_max = (int)(x+width); 
      
      //if (thr_min< (int)(vp_set[d]-width)){ thr_min = vp_set[d]-width;}
      //if (thr_max> (int)(vp_set[d]+width)){ thr_max = vp_set[d]+width;}
      
      if (d >28 ){thr_min = 0.; thr_max = 80.;}
      
      
    fit_state0:
      TFitResultPtr f_s1; 
      f_s1= hdcnt[ch][d]->Fit("gaus","SWLQR","",thr_min,thr_max); // fitting disc in the range thr_min-thr_max
      //f_s1= hdcnt[ch][d]->Fit("gaus","SWL"); // fitting disc in the range thr_min-thr_max
      Int_t fitstatus = f_s1;
      float mean_tmp =0;
      float sigma_tmp =0;
      
      if ( fitstatus == 0 && hdcnt[ch][d]->GetFunction("gaus")->GetChisquare()!=0) {
	
	/*mean_tmp = f_s1->Parameter(1);
	  sigma_tmp = f_s1->Parameter(2);
	  TF1 *f1 = new TF1("f1", "gaus", mean_tmp-4*sigma_tmp, mean_tmp+4*sigma_tmp);
	  f1->SetParameter(1,mean_tmp);
	  
	  hdcnt[ch][d]->Fit("f1","SWLQR","",mean_tmp-4*sigma_tmp, mean_tmp+4*sigma_tmp );
	  
	  
	  f_s1mean+= hdcnt[ch][d]->GetFunction("f1")->GetParameter(1);
	  f_sigma2 = hdcnt[ch][d]->GetFunction("f1")->GetParameter(2);
	  f_mean = hdcnt[ch][d]->GetFunction("f1")->GetParameter(1);
	  f_sigma = hdcnt[ch][d]->GetFunction("f1")->GetParameter(2);
	*/
	
	
	f_s1mean+= f_s1->Parameter(1);
	f_sigma2 = f_s1->Parameter(2);
	f_mean = f_s1->Parameter(1);
	f_sigma = f_s1->Parameter(2);
	
	
	if (d==30) f_mean_d1=f_s1->Parameter(1);
	if (d==10) f_mean_d20=f_s1->Parameter(1);
	if (d==0) f_mean_d30=f_s1->Parameter(1);
	
	//writefile<<f_mean<<'\t';
	
	if (f_mean!=0 && f_sigma2 >0.5){
	  hdisc_sig[ch][d]->Fill(d,f_sigma2);
	  hdisc_mean[ch][d]->Fill(d,f_mean);
	  if (ch == ch_sel) {
	    hfit_s2->Fill(d,f_sigma2);
	    h_aux1->SetBinContent(d+1,f_mean);
	    h_aux1->SetBinError(d+1,f_sigma2); 
	    //h_aux2->Fill(d,mean);
	  }
	}
	else {
	  while (width > 5){
	    width = width-5;
	    thr_min = (int)(vp_set[d]-width); 
	    thr_max = (int)(vp_set[d]+width); 
	    goto fit_state0;
	  }
	}  
	
	if (f_sigma2 >0.5){
	  f_s1sigma+= f_sigma2;
	  d_counter++; 
	}
	//delete f1;
      }
      
      else {
	cout << "fit not ok" << endl;
      }     
      cout<<endl; 
      
      hmeanf->Fill(ch,d,f_mean);
      hsigef->Fill(ch,d,f_sigma);
      f_s1mean = f_s1mean/d_max;
    }
  }
#ifdef KAZ_DEBUG
  cout << "Fit_values returns true for now. check if it is ok." << endl;
#endif
  return true;
}


//------------------------------------------+-----------------------------------
//! Analyzing the data. Fitting S-curves
bool trim_adc::Fit_values_erfc(int width, int d_cut_min, int d_cut_max, int cut_db_pulses_user)
{
  d_cut_min = 30-d_cut_min;
  d_cut_max = 30-d_cut_max;
  
  int plateau = cut_db_pulses_user;
  
  if (d_cut_min < 0 | d_cut_max > 30) return false;
  
  else { 
    if (d>=d_cut_max && d<=d_cut_min) {
      int binmax = hdcnt[ch][d]->GetMaximumBin();
      float x =hdcnt[ch][d]->GetXaxis()->GetBinCenter(binmax); 
      
      /*
      thr_min = (int)(vp_set[d]-width); 
      thr_max = (int)(vp_set[d]+width); 
      */
      
      width+=15;
      thr_min = (int)(x-width); 
      thr_max = (int)(x+width); 
      
      
      
      //if (thr_min< (int)(vp_set[d]-2*width)) thr_min = vp_set[d]-3*width;
      //if (thr_max> (int)(vp_set[d]+2*width)) thr_max = vp_set[d]+3*width;
      if (d >28 ){thr_min = 0.; thr_max = 80.;}
      
      TF1 *f_erfc = new TF1("f_erfc", "[0]-[1]*TMath::Erfc((x-[2])/(sqrt(2)*[3]))",0,250);
      
      fit_state0:
      TFitResultPtr f_s1; 
      //f_s1= hdcnt[ch][d]->Fit("f_erfc","SWLQR","",thr_min,thr_max); // fitting disc in the range thr_min-thr_max
      f_erfc->SetParameters(plateau,plateau/2,x,5);
      f_s1= hscurve[ch][d]->Fit("f_erfc", "SWR", "",thr_min, thr_max); // fitting disc in the range thr_min-thr_max
      Int_t fitstatus = f_s1;
      float mean_tmp =0;
      float sigma_tmp =0;
      
      if ( fitstatus == 0 ) {
	
	
	f_s1mean_erfc+= f_s1->Parameter(2);
	f_sigma2_erfc = f_s1->Parameter(3);
	f_mean_erfc = f_s1->Parameter(2);
	f_sigma_erfc = f_s1->Parameter(3);
	
	/*
	if (d==30) f_mean_d1=f_s1->Parameter(3);
	if (d==10) f_mean_d20=f_s1->Parameter(3);
	if (d==0) f_mean_d30=f_s1->Parameter(3);
	*/
	if (ch == ch_sel) {
	    hfit_s2_erfc->Fill(d,f_sigma2_erfc);
	    h_aux1_erfc->SetBinContent(d+1,f_mean_erfc);
	    h_aux1_erfc->SetBinError(d+1,f_sigma2_erfc); 
	    //h_aux2->Fill(d,mean);
	  }
	
	if (f_sigma2_erfc >0.5){
	  f_s1sigma_erfc+= f_sigma2_erfc;
	  d_counter_erfc++; 
	}
	//delete f1;
      }
           
      else {
	  cout << "fit not ok" << endl;
      }     
      cout<<endl; 
      
    hmeane_erfc->Fill(ch,d,f_mean_erfc);
    hsige_erfc->Fill(ch,d,f_sigma_erfc);
    f_s1mean_erfc = f_s1mean_erfc/d_max;
    
    delete f_erfc;
    }
  }
#ifdef KAZ_DEBUG
  cout << "Fit_values_erfc returns true for now. check if it is ok." << endl;
#endif
  return true;
}

//------------------------------------------+-----------------------------------
//! Analyzing the data. Calculating values of a Normal distribution
bool trim_adc::Calc_values(int width, int d_cut_min, int d_cut_max)
{
  d_cut_min = 30-d_cut_min;
  d_cut_max = 30-d_cut_max;
  
  if (d_cut_min < 0 | d_cut_max > 30) return false;
  
  else { 
    if (d>=d_cut_max && d<=d_cut_min) {
    sum_delta = 0.0000001;
    sum_sig = 0;
    d_cnt = 0;
    
   // range where I will look for the Scurve. Expanded range compare to the fit
    
   thr_min = (int)(vp_set[d]-vp_min-width); 
   if (thr_min <=0){thr_min = 1;}
   thr_max = (int)(vp_set[d]-vp_min+width);
   if (thr_max >vp_max-vp_min){thr_max = vp_max-vp_min;} 
  
  thr_min = 0;
  thr_max = 250;
  ivp = thr_min;
   
   for ( vp = thr_min; vp <thr_max; vp += vp_step ) {
	d_cnt = vcnt_soft[ch][d][ivp]-vcnt_soft[ch][d][ivp-1]; 
	//if (d_cnt <0 | d_cnt >35) d_cnt =0;
	if (d_cnt <0) d_cnt =0;
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
   
   ivp = thr_min;
   
   for ( vp = thr_min; vp <thr_max; vp += vp_step ) {
      d_cnt = vcnt_soft[ch][d][ivp]-vcnt_soft[ch][d][ivp-1];
      //if (d_cnt <0 | d_cnt >35) d_cnt =0;
      if (d_cnt <0) d_cnt =0;
      sum_delta += d_cnt;
      sum_sig +=((vp+vp_min)-mean)*((vp+vp_min)-mean)*d_cnt;
      ivp++;
   }
   sigma = sqrt(sum_sig/sum_delta);
   
   
   if (ch == ch_sel) {
    h_aux1_c->Fill(d,mean);
    hcalc_s2->Fill(d,sigma);
   }

   
   if (sigma >0){
  // if (sigma >0){
    sum_sige +=sigma;
    d_counter1++;
    }
    
   hsige->Fill(ch,d,sigma); 
    }
  }

#ifdef KAZ_DEBUG
  cout << "Calc_values returns true for now. check if it is ok." << endl;
#endif
  return true;
}


//------------------------------------------+-----------------------------------
//! Analyzing data. Fast discriminator
void trim_adc::Fitting_Fast(int width){
  width = width - 10;
  thr_min = (int)(vp_set[30]-width); 
  thr_max = (int)(vp_set[30]+width);
  
  TFitResultPtr f_s2; 
  f_s2= hdcnt[ch][31]->Fit("gaus","SWLQR","",20,70);
  Int_t fitstatus = f_s2;
      if ( fitstatus == 0 && hdcnt[ch][31]->GetFunction("gaus")->GetChisquare()!=0){
	f_mean = f_s2->Parameter(1);
	f_sigma = f_s2->Parameter(2)*350;
	f_figma_fast = f_s2->Parameter(2)*350;
	f_mean_fast = f_s2->Parameter(1);
	hmean_fast->Fill(ch, f_mean);
	hnoise_fast->Fill(ch,f_sigma);
      }   
}

//------------------------------------------+-----------------------------------
//! Displaying histograms_ADC
void trim_adc::Display_histo_adc(int width_user, int d_cut_min, int d_cut_max,  int* ch_comp, int grp_sel)
{
  
  val_file.open(filename_val);
  
  d_cut_min = 30-d_cut_min;
  d_cut_max = 30-d_cut_max;
  
  TCanvas *c1 =new TCanvas("c1","c1");
  c1->cd();
  thr_min = (int)(vp_set[0]-width_user); 
  thr_max = (int)(vp_set[0]+width_user);
  hdcnt[ch_sel][0]->Draw("");
  hdcnt[ch_sel][0]->Fit("gaus","SWQR","", thr_min,thr_max);
  //hdcnt[ch_sel][0]->Fit("gaus","SWL");
  c1->Close();
  
 
  int n =2;
  TCanvas *c2 =new TCanvas("c2","c2");
  c2->Divide(6,6);
  c2->cd(1);
  hdcnt[ch_sel][30]->Draw("");
  hdcnt[ch_sel][30]->Fit("gaus","SWQR","",10,80);
  //hdcnt[ch_sel][30]->Fit("gaus","SW");
  for (d=29; d>=d_cut_max;d--){
    c2->cd(n);
    hdcnt[ch_sel][d]->Draw("");
    thr_min = (int)(vp_set[d]-width_user); 
    thr_max = (int)(vp_set[d]+width_user); 
    n++;
  }
  c2->Write();
  c2->Close();
 
  
  n =2;
  TCanvas *c2_erfc =new TCanvas("c2_erfc","c2_erfc");
  //TF1 *f_erfc = new TF1("f_erfc", "[0]-[1]*TMath::Erfc((x-[2])/(sqrt(2)*[3]))",0,250);
  c2_erfc->Divide(6,6);
  c2_erfc->cd(1);
  hscurve[ch_sel][30]->Draw("");
  hscurve[ch_sel][30]->Fit("f_erfc","SWQRN","",10,80);
  //hdcnt[ch_sel][30]->Fit("gaus","SW");
  for (d=29; d>=d_cut_max;d--){
    c2_erfc->cd(n);
    hscurve[ch_sel][d]->Draw("");
    thr_min = (int)(vp_set[d]-width_user); 
    thr_max = (int)(vp_set[d]+width_user); 
    
    n++;
  }
  c2_erfc->Write();
  c2_erfc->Close();
  
  n=0;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c_scurves =new TCanvas("c_scurves","S-Curves_selected channel",1200,900);
  c_scurves->cd(1)->SetGrid();
  hscurve[ch_sel][30]->Draw("");
  hscurve[ch_sel][30]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
  hscurve[ch_sel][30]->GetYaxis()->SetTitle("Number of counts");
  hscurve[ch_sel][30]->SetLineWidth(2);
   for (d=d_min+1;d <31;d++) {
	n = 30-d;
        hscurve[ch_sel][n]->Draw("SAME");
	//hscurve[ch_sel][d]->SetLineColor(n);
	hscurve[ch_sel][n]->SetLineWidth(2);
    } 
  c_scurves->Write();
  c_scurves->Close();
  
  
 n=0;
  gStyle->SetOptStat(0);
  TCanvas *c3 =new TCanvas("c3","S-Curves_selected channel",1400,600);
  c3->Divide(2,1);
  c3->cd(1)->SetGrid();
  hscurve[ch_sel][30]->Draw("");
  hscurve[ch_sel][30]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
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
  h_aux1->SetMarkerColor(kRed);
  h_aux1->Draw("PE");
  h_aux1_erfc->Draw("PESAME");
  h_aux1_erfc->SetMarkerColor(kBlue);
  gStyle->SetErrorX(0);
  h_aux1->Fit("pol1", "R", "",d_cut_max, d_cut_min+1);
  h_aux1->SetMarkerStyle(kFullCircle);
  h_aux1_erfc->SetMarkerStyle(kFullCircle);
  h_aux1->GetXaxis()->SetTitle("Discriminator Number");
  h_aux1->GetYaxis()->SetTitle("Discriminator Thr [amp_cal_units]"); 
  h_aux1->Write();
  h_aux1_erfc->Write();
  
  float par0 = h_aux1->GetFunction("pol1")->GetParameter(0);
  float par1 = h_aux1->GetFunction("pol1")->GetParameter(1);
  float residuals = 0;
  
  pad2->cd()->SetGrid();
  for(d=d_cut_max; d<=d_cut_min; d++){
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
  h_aux2->Write();
  c3->Write();
  c3->Close();


   int i=1;
   n = 0;
   TCanvas *c4 =new TCanvas("c4","S-Curves_selected group");
   c4->Divide(6,5); 
   for (ch=grp_sel;ch<=ch_max;ch+=4){
     c4->cd(i);
     hscurve[ch][0]->Draw("HIST");
     hscurve[ch][0]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
     hscurve[ch][0]->GetYaxis()->SetTitle("Number of counts");
     hscurve[ch][0]->GetYaxis()->SetRangeUser(0,500);
     for (d=d_min+1;d <31;d++) {
	n = 30-d;
	hscurve[ch][d]->Draw("SAME");
	//hscurve[ch][d]->SetLineColor(d);
	}
     i++;
   }
   c4->Write();
   c4->Close();
   
  TCanvas *c5 =new TCanvas("c5","Mean_and sigma_channels");
  c5->Divide(3,2);
  c5->cd(1);
  hmeanf->Draw("COLZ");
  hmeanf->GetXaxis()->SetTitle("Channel number");
  hmeanf->GetYaxis()->SetTitle("ADC discriminator");
  hmeanf->GetZaxis()->SetRangeUser(0,255);

  c5->cd(2);
  hmeane->Draw("COLZ");
  hmeane->GetXaxis()->SetTitle("Channel number");
  hmeane->GetYaxis()->SetTitle("ADC discriminator");
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
  py[0]->GetXaxis()->SetTitle("ADC value LSB");
  py[0]->GetYaxis()->SetTitle("Signal Vp amplitude [amp_cal_units]");
  for (ch = ch_min+2; ch<=ch_max;ch+=ch_step){
    py[ch]->Draw("same");
    //py[ch]->SetLineColor(ch);
    
  }
  
  c5->cd(4);
  hsigef->Draw("COLZ");
  hsigef->GetZaxis()->SetRangeUser(0,8);
  hsigef->GetXaxis()->SetTitle("Channel number");
  hsigef->GetYaxis()->SetTitle("ADC discriminator");
  
  c5->cd(5);
  hsige->Draw("COLZ");
  hsige->GetZaxis()->SetRangeUser(0,8);
  hsige->GetXaxis()->SetTitle("Channel number");
  hsige->GetYaxis()->SetTitle("ADC discriminator");
  
  c5->cd(6)->SetGrid();
  ch =0;
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);
  for (ch =ch_min; ch<=ch_max;ch+=ch_step){
    //if (ch>0 && ch< 127){
    py[ch]->Fit("pol1","R","",d_cut_max,d_cut_min+1);
    float par1_adc = py[ch]->GetFunction("pol1")->GetParameter(1);
    float par0_adc = py[ch]->GetFunction("pol1")->GetParameter(0);
    float par2_adc = py[ch]->GetFunction("pol1")->GetChisquare();
    h_adc_slope->Fill(-1*par1_adc*0.056);
    h_adc_slope_ch->Fill(ch,-1*par1_adc*0.056);
    h_adc_int->Fill(par0_adc);
   
   
    
    float dnl = 0;
    for (d=d_min-1;d<d_max+1; d++){
      dnl =  py[ch]->GetBinContent(d+1)-(par1_adc*d+par0_adc);
      h_adc_mean->Fill(dnl);
    }
    val_file<<ch<<"\t"<<hfit_s_erfc->GetBinContent(ch+1)<<"\t\t"<<-1*par1_adc<<"\t\t"<<hnoise_fast->GetBinContent(ch+1)<<endl;
  //}
  }
  
  h_adc_mean->Draw("");
  h_adc_mean->SetFillColor(kGreen-6);
  h_adc_mean->SetLineColor(kGreen-6);
  h_adc_mean->Fit("gaus", "SWL");
  h_adc_mean->GetXaxis()->SetTitle("Residuals [amp_cal_units]");
  h_adc_mean->GetYaxis()->SetTitle("Entries");
  h_adc_mean->Write();
  h_adc_slope->SetTitle("ADC gain uniformity; ADC gain [fC/LSB]; Entries");
  h_adc_slope->Fit("gaus");
  h_adc_slope->Write();
  h_adc_slope_ch->SetMarkerStyle(kFullCircle);
  h_adc_slope_ch->Draw("P");
  h_adc_slope_ch->SetTitle("ADC gain uniformity;Channel; ADC gain [fC/LSB]");
  h_adc_slope_ch->Write();
  
  h_adc_int->Draw("");
  h_adc_int->SetFillColor(kGreen-6);
  h_adc_int->SetLineColor(kGreen-6);
  h_adc_int->Fit("gaus");
  h_adc_int->Write();
  
  c5->Write();
  c5->Close();
  
  TCanvas *c_adc_lin =new TCanvas("c_adc_lin","ASIC_ADC_linearity");
  c_adc_lin->cd()->SetGrid();
  py[0]->Draw("");
  py[0]->GetXaxis()->SetTitle("ADC value");
  py[0]->GetYaxis()->SetTitle("Threshold level [amp_cal units]");
  cout << "KAZ KAZ  watch out, decision based on ch_max maybe wrong" << endl;
  for (ch = ch_min+2; ch<=ch_max;ch+=ch_step){
    py[ch]->Draw("same");    
  }
  
  c_adc_lin->Write();
  c_adc_lin->Close();

  n=0;
  TCanvas *c6 =new TCanvas("c6","Noise level/discriminator");
  c6->Divide(2,1);
  c6->cd(1)->SetGrid();
  h_aux1->Draw("PTEXT");
  h_aux1->Fit("pol1");
  h_aux1->SetMarkerColor(kRed);
  h_aux1->GetXaxis()->SetTitle("Discriminator Number");
  h_aux1->GetYaxis()->SetTitle("Discriminator Thr [amp_cal_units]");
  h_aux1_c->Draw("PSAMETEXT");
  h_aux1_c->Fit("pol1");
  h_aux1_c->SetMarkerColor(kBlue);
  h_aux1_c->SetMarkerStyle(kFullCircle);
    
  c6->cd(2)->SetGrid();
   hfit_s2->Draw("HISTTEXT");
   hfit_s2->SetLineColor(kRed);
   hfit_s2->GetXaxis()->SetTitle("Discriminator Number");
   hfit_s2->GetYaxis()->SetTitle("Sigma [electrons]");
   hcalc_s2->Draw("HISTTEXTsame");
   hcalc_s2->SetLineColor(kBlue);
   
   c6->Close();
   
   
   TCanvas*c_comp = new TCanvas("c_comp","c_comp");
   c_comp->Divide(4,3);
   
   ofstream test_file;
   test_file.open("vref_t_184.txt");
   
   for (int j=0; j<4; j++){
    c_comp->cd(j+1)->SetGrid();
    hscurve[ch_comp[j]][30]->Draw("");
    hscurve[ch_comp[j]][30]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
    hscurve[ch_comp[j]][30]->GetYaxis()->SetTitle("Number of counts");
    hscurve[ch_comp[j]][30]->SetLineWidth(2);
    for (d=d_min+1;d <31;d++) {
      n = 30-d;
      hscurve[ch_comp[j]][n]->Draw("SAME");
      //hscurve[ch_sel][d]->SetLineColor(n);
      hscurve[ch_comp[j]][n]->SetLineWidth(2);
    } 
    int k = j+5;
    c_comp->cd(k)->SetGrid();
    hdisc_mean[ch_comp[j]][0]->Draw("P");
    hdisc_mean[ch_comp[j]][0]->SetMarkerStyle(8);
    hdisc_mean[ch_comp[j]][0]->SetMarkerSize(0.7);
    hdisc_mean[ch_comp[j]][0]->SetMarkerColor(kRed);
    hdisc_mean[ch_comp[j]][0]->GetXaxis()->SetTitle("Discriminator number");
    hdisc_mean[ch_comp[j]][0]->GetYaxis()->SetTitle("Disc. thr [amp_cal_units]");
    hdisc_mean[ch_comp[j]][0]->GetYaxis()->SetRangeUser(0,260);
    hdisc_mean[ch_comp[j]][0]->SetLineWidth(2);
    
    test_file<<"ch"<<"\t"<<ch_comp[j]<<"\n";
    test_file<<hdisc_mean[ch_comp[j]][0]->GetBinContent(1)<<endl;
    for (d=d_min+1;d <31;d++) {
      n = 30-d;
      hdisc_mean[ch_comp[j]][d]->Draw("PSAME");
      //hdisc_sig[ch_sel][d]->SetLineColor(n);
      hdisc_mean[ch_comp[j]][d]->SetLineWidth(2);
      hdisc_mean[ch_comp[j]][d]->SetMarkerStyle(8);
      hdisc_mean[ch_comp[j]][d]->SetMarkerSize(0.7);
      hdisc_mean[ch_comp[j]][d]->SetMarkerColor(kRed);
      test_file<<hdisc_mean[ch_comp[j]][d]->GetBinContent(d+1)<<endl;
    }  
    
    int z = j+9;
    c_comp->cd(z)->SetGrid();
    hdisc_sig[ch_comp[j]][0]->Draw("P");
    hdisc_sig[ch_comp[j]][0]->SetMarkerStyle(8);
    hdisc_sig[ch_comp[j]][0]->SetMarkerSize(0.7);
    hdisc_sig[ch_comp[j]][0]->GetXaxis()->SetTitle("Discriminator number");
    hdisc_sig[ch_comp[j]][0]->GetYaxis()->SetTitle("ENC [amp_cal_units]");
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
  test_file.close();
  c_comp->Write();
  c_comp->Close();

  TCanvas *c_noise = new TCanvas("c_noise","c_noise");
  c_noise->cd()->SetGrid();
  
  hfit_s->Draw("");
  hfit_s->GetXaxis()->SetTitle("Channel number");
  hfit_s->GetYaxis()->SetTitle("ENC [electrons]");
  hfit_s->SetLineWidth(2);
  hfit_s->SetLineColor(kRed);
  //hfit_s->GetYaxis()->SetRangeUser(0,20);
  hcalc_s->Draw("same");
  hcalc_s->SetLineWidth(2);
  hcalc_s->SetLineColor(kBlue);
  //hfit_s->SaveAs("noise.pdf");
  hfit_s->Write();
  hcalc_s->Write();
  hfit_s_erfc->Write();
  
  TLegend *l_noise = new TLegend(0.5,0.7,0.88,0.88);
  l_noise->SetHeader(""); // option "C" allows to center the header
  l_noise->SetLineColor(kWhite);
  l_noise->AddEntry(hfit_s,"Noise_fitting","l");
  l_noise->AddEntry(hcalc_s,"Noise_calculated_RMS","l");
  l_noise->Draw();
  
  c_noise->Close();
  
  TCanvas *c_all_scurves = new TCanvas("c_all_scurves","c_all_scurves", 1600, 800);
  c_all_scurves->cd()->SetGrid();
  hscurve[0][0]->Draw("");
  hscurve[0][0]->GetXaxis()->SetTitle("Vp [amp_cal_units]");
  hscurve[0][0]->GetYaxis()->SetTitle("Entries");
  cout << "WATCH OUT . DECISION regarding ch_max maybe wrong" << endl;
  for (ch = ch_min; ch<=ch_max; ch++){
    for (d=d_min; d<d_max+1; d++){
      hscurve[ch][d]->Draw("same");
      //hscurve[ch][d]->SetLineColor(d);
    }
  }
  c_all_scurves->Write();
  c_all_scurves->Close();
  
  val_file.close();
  
}

//------------------------------------------+-----------------------------------
//! Displaying histograms_FAST_discriminator
void trim_adc::Display_histo_fast(int width_user, int* ch_comp){
  
  thr_min = (int)(vp_set[30]-30); 
  thr_max = (int)(vp_set[30]+50);  
  
  TCanvas *c_fast =new TCanvas("c_fast","ADC vs FAST comparison", 1600,1000);
  c_fast->Divide(2,2);
  
  c_fast->cd(1)->SetGrid();
  hscurve[ch_sel][30]->Draw("");
  hscurve[ch_sel][30]->SetTitle("S-curves comparison");
  hscurve[ch_sel][30]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
  hscurve[ch_sel][30]->GetYaxis()->SetTitle("Number of counts");
  hscurve[ch_sel][30]->GetXaxis()->SetRangeUser(vp_min, vp_min+80);
  hscurve[ch_sel][30]->SetLineWidth(2);
  hscurve[ch_sel][30]->SetLineColor(kRed);
  hscurve[ch_sel][31]->Draw("SAME");
  hscurve[ch_sel][31]->SetLineColor(kBlue);
  hscurve[ch_sel][31]->SetLineWidth(2);
  
  TLegend *l_scurve = new TLegend(0.65,0.65,0.85,0.85);
  l_scurve->SetHeader("S-curves for test_channel"); // option "C" allows to center the header
  l_scurve->SetLineColor(kWhite);
  l_scurve->AddEntry(hscurve[ch_sel][30],"S-curve_ADC_disc_30","l");
  l_scurve->AddEntry(hscurve[ch_sel][31],"S-curve_FAST_disc","l");
  l_scurve->Draw();
  
  c_fast->cd(2)->SetGrid();
  TPad *pad1_h = new TPad("pad1","ADC",0,0,0.45,0.95,0);
  TPad *pad2_h = new TPad("pad2","FAST",0.5,0,0.95,0.95,0);
  pad1_h->SetBorderMode(0);
  pad1_h->SetTopMargin(0.05);
  pad1_h->SetBottomMargin(4.0);
  pad2_h->SetBorderMode(0);
  pad2_h->SetTopMargin(0.05);
  pad2_h->SetBottomMargin(4.0);
  pad1_h->Draw();
  pad2_h->Draw();
  pad1_h->cd()->SetGrid();  
  hdcnt[ch_sel][30]->Draw("");
  hdcnt[ch_sel][30]->SetLineColor(kRed);
  hdcnt[ch_sel][30]->SetFillColor(kRed-10);
  hdcnt[ch_sel][30]->Fit("gaus", "SWQR", "",thr_min, thr_max);
//hdcnt[ch_sel][30]->Fit("gaus","SWLR","",vp_min+10, vp_max-20);
  hdcnt[ch_sel][30]->SetTitle("ADC disc_30");
  hdcnt[ch_sel][30]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
  hdcnt[ch_sel][30]->GetXaxis()->SetRangeUser(vp_min, vp_min+80);
  pad2_h->cd()->SetGrid();  
  hdcnt[ch_sel][31]->Draw("");
  hdcnt[ch_sel][31]->SetLineColor(kBlue);
  hdcnt[ch_sel][31]->SetFillColor(kBlue-9);
  hdcnt[ch_sel][31]->Fit("gaus", "SWQR", "",thr_min, thr_max);
  //hdcnt[ch_sel][31]->Fit("gaus","SWLR","",vp_min+10, vp_max-20);
  hdcnt[ch_sel][31]->SetTitle("FAST disc");
  hdcnt[ch_sel][31]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
  hdcnt[ch_sel][31]->GetXaxis()->SetRangeUser(vp_min, vp_min+80);
  
  
  c_fast->cd(3)->SetGrid();  
  TPad *pad1_m = new TPad("pad1_m","ADC",0,0,0.75,0.95,0);
  TPad *pad2_m = new TPad("pad2_m","FAST",0.76,0,1,0.95,0);
  pad1_m->SetBorderMode(0);
  pad1_m->SetTopMargin(0.05);
  pad1_m->SetBottomMargin(4.0);
  pad2_m->SetBorderMode(0);
  pad2_m->SetTopMargin(0.05);
  pad2_m->SetBottomMargin(4.0);
  pad1_m->Draw();
  pad2_m->Draw();
  pad1_m->cd()->SetGrid(); 
  for(ch=ch_min; ch<=ch_max; ch++){
   hmean_f_disc_30->Fill(ch,hmeanf->GetBinContent(ch+1,31)); 
  }
  hmean_f_disc_30->Draw("");
  hmean_f_disc_30->GetXaxis()->SetTitle("Channel number");
  hmean_f_disc_30->GetYaxis()->SetTitle("Discriminator_thr [amp_cal_units]");
  hmean_f_disc_30->GetYaxis()->SetRangeUser(35,65);
  hmean_f_disc_30->SetLineColor(kRed);
  hmean_f_disc_30->SetLineWidth(2);
  hmean_fast->Draw("SAME");
  hmean_fast->SetLineColor(kBlue);
  hmean_fast->SetLineWidth(2);
  
  TLegend *l_mean = new TLegend(0.2,0.75,0.5,0.93);
  l_mean->SetHeader(""); // option "C" allows to center the header
  l_mean->SetLineColor(kWhite);
  l_mean->AddEntry(hmean_f_disc_30,"Threshold_ADC_disc_30","l");
  l_mean->AddEntry(hmean_fast,"Threshold_FAST_disc","l");
  l_mean->Draw();
  
  Double_t diff = 0;
  Double_t mean_adc = 0;
  Double_t mean_fast = 0;
 
  pad2_m->cd()->SetGrid(); 
  
  for (int j=1; j<=ch_max; j++){
    mean_adc = hmean_f_disc_30->GetBinContent(j);
    mean_fast = hmean_fast->GetBinContent(j);
    diff = mean_adc- mean_fast;
    h_diff_fast_adc->Fill(diff);
  }
  h_diff_fast_adc->Draw("");
  //h_diff_fast_adc->GetXaxis()->SetRangeUser(-5,5);
  h_diff_fast_adc->SetLineWidth(2);
  h_diff_fast_adc->SetLineColor(kBlack);
  h_diff_fast_adc->SetFillColor(kGray+1);
  //h_diff_fast_adc->SetFillStyle();
  //h_diff_fast_adc->Fit("gaus");
  
  c_fast->cd(4)->SetGrid();
  
  for(ch=ch_min; ch<=ch_max; ch++){
   hfit_s_disc_30->Fill(ch,hsigef->GetBinContent(ch+1,31)*350); 
  }
  hfit_s_disc_30->Draw("");
  hfit_s_disc_30->GetXaxis()->SetTitle("Channel number");
  hfit_s_disc_30->GetYaxis()->SetTitle("ENC [electrons]");
  hfit_s_disc_30->SetLineWidth(2);
  hfit_s_disc_30->SetLineColor(kRed);
  hnoise_fast->Draw("SAME");
  hnoise_fast->SetLineColor(kBlue);
  hnoise_fast->SetLineWidth(2);
  hnoise_fast->Write();
  
  TLegend *l_noise_fast = new TLegend(0.2,0.7,0.58,0.88);
  l_noise_fast->SetHeader(""); // option "C" allows to center the header
  l_noise_fast->SetLineColor(kWhite);
  l_noise_fast->AddEntry(hmean_f_disc_30,"Noise_ADC_disc_30","l");
  l_noise_fast->AddEntry(hnoise_fast,"Noise_FAST_disc","l");
  l_noise_fast->Draw();
  
  c_fast->Write();
  c_fast->Close();
  
  
  
  TCanvas *c_fast_disc =new TCanvas("c_fast_disc","Fast_discriminator_channels", 1600,1000);
  c_fast_disc->Divide(2,1);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);
  c_fast_disc->cd(1)->SetGrid();
  hscurve[0][31]->Draw("L");
  hscurve[0][31]->GetXaxis()->SetRangeUser(vp_min, vp_min+80);
  hscurve[0][31]->SetTitle("S-curves fast_disc");
  hscurve[0][31]->GetXaxis()->SetTitle("Pulse amplitude [amp_cal_units]");
  hscurve[0][31]->GetYaxis()->SetTitle("Number of counts");
  hscurve[0][31]->SetLineColor(kRed);
  hscurve[0][31]->SetLineWidth(2);
  for(ch =ch_min+1; ch<=ch_max; ch++){
    hscurve[ch][31]->Draw("LSAME");
    hscurve[ch][31]->SetLineWidth(2);
    if (ch%4 == 0) hscurve[ch][31]->SetLineColor(kRed);
    else if (ch%4 == 1) hscurve[ch][31]->SetLineColor(kBlue);
    else if (ch%4 == 2) hscurve[ch][31]->SetLineColor(kGreen);
    else hscurve[ch][31]->SetLineColor(kMagenta);
  }
  
 c_fast_disc->cd(2)->SetGrid();
 for (int j = 1; j<=ch_max; j++){
    hmean_fast_disc->Fill(hmean_fast->GetBinContent(j));
  }
 hmean_fast_disc->Draw("");
 hmean_fast_disc->Fit("gaus", "SWLQR", "",thr_min, thr_max);
 hmean_fast_disc->GetXaxis()->SetTitle("Threshold fast disc [amp_cal_units]");
 hmean_fast_disc->GetYaxis()->SetTitle("Entries");
 hmean_fast_disc->SetLineColor(kBlack);
 hmean_fast_disc->SetLineWidth(2);
 hmean_fast_disc->SetFillColor(kGray+1);
 
 
 c_fast_disc->Write();
 c_fast_disc->Close();
 

 
}

//------------------------------------------+-----------------------------------
//! Displaying selected values
void trim_adc::Display_values(int* ch_comp)
{
  
  cout<<"Channels noise comparison"<<endl;
  cout<<"-----------------------------------"<<endl;
    for (int j =0; j<10; j++){
    cout<<"ch:"<<"\t"<<ch_comp[j]<<"\t"<< hfit_s->GetBinContent(ch_comp[j]+1) <<endl;
  }
  cout<<"-----------------------------------"<<endl;

  
}
  
