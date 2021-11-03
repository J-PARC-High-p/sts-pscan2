#ifndef _trim_adc_HXX
#define _trim_adc_HXX

#include "TFile.h"
#include <string>
#include <TH1F.h>
#include <TH2F.h>
#include <unistd.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include "TObject.h"

class trim_adc
{
  public:
    explicit trim_adc(TString filenameData = "z");
    ~trim_adc();
  
  private:

  int ch;
  int ch_min;
  int ch_max;
  int ch_step;
  
  int d;
  int d_counter;
  int d_counter1;
  int d_counter_erfc;
  int d_min;
  int d_max;
  int d_step;
  int cut_db_pulses;
  
  int grp;
  int grp_min;
  int grp_max;
  int grp_step;
  
  int vp;
  int ivp;
  // input values for the vp scan. Please modify this values acording to the pscan vp range; values should be correct or s-curves will be shifted. *Needs improvements
  int vp_min;
  int vp_max;
  int vp_step;
  int vcnt[130][32][260];   
  int vcnt_soft[130][32][260];
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
  
  float amp_cal_min;   		
  float amp_cal_max;  						
  int thr_val1;					
  int vp_set[31];
  
  int read_flag;
  int soft_flag;
  
  
  double f_s1mean;
  double f_mean;
  double f_mean_d1;
  double f_mean_d20;
  double f_mean_d30;
  double f_mean_fast;
  double f_s1sigma;
  double f_sigma;
  double f_sigma2;
  double f_figma_fast;
  
  // ----- extras ----
  double f_s1mean_erfc;
  double f_mean_erfc;
  double f_mean_d1_erfc;
  double f_mean_fast_erfc;
  double f_s1sigma_erfc;
  double f_sigma_erfc;
  double f_sigma2_erfc;
  
  
  int thr_min; 
  int thr_max;
  
  int ch_sel;
  
  TFile *file1;
  TString filename_data;
  TString filename_root;
  TString filename_val;
  ifstream scanfile;
  ofstream val_file;
  
  TH1F *hdcnt[128][32];   // Differential count 
  TH1F *hscurve[128][32]; // S-Curves
  TH1F *hdisc_sig [128][32];  // Discriminators thr values for a specific channel 
  TH1F *hdisc_mean [128][32];
  
  TH1F *h_mean;
  TH1F *h_adc_mean;
  TH1F *h_adc_slope;
  TH1F *h_adc_slope_ch;
  TH1F *h_adc_int;
  TH1F *h_sigma;
  TH1F *h_aux1;
  TH1F *h_aux1_erfc;
  TH1F *h_aux1_c;
  TH1F *h_aux2;
  
  
  TH2F *hmeane;
  TH2F *hsige;
  
  TH2F *hmeane_erfc;
  TH2F *hsige_erfc;
  
  TH1F *hcalc_s;
  TH1F *hfit_s;
  TH1F *hfit_s_erfc;
  TH2F *hmeanf;
  TH2F *hsigef;
  
  TH1F *hfit_s2; 
  TH1F *hfit_s2_erfc;
  TH1F *hcalc_s2; 
  TH1F *hnoise_fast;
  TH1F *hmean_fast;
  TH1F *h_diff_fast_adc;
  TH1F *hmean_f_disc_30;
  TH1F *hfit_s_disc_30;
  TH1F *hmean_fast_disc;
  
  


public:

  bool Init_ch(int chMin, int chMax, int chStep=1);

  bool Init_grp(int grphMin, int grpMax, int grpStep=1);

  bool Init_d(int dMin, int dMax, int dStep=1);

  bool Init_vp(int vpMin, int vpMax, int vpStep=1);
  bool Fitting_windows(int ampcalMin, int ampcalMax, int width=10);
  
  bool Check_root_file();
  
  void Close_root_file();
  void Create_histo();
  bool Reading_file(int cut_db_pulses_user);
  void Analysis(int cut_db_pulses_user,int width_user, int dcut_min_user,int dcut_max_user, int ch_sel_user, bool fit=true, bool calc=true, bool fast_fit =false);
  
  int  Get_fitting_values(int i){ return vp_set[i];};
  
  void Display_histo_adc(int width_user, int d_cut_min, int d_cut_max, int* ch_comp,  int grp_sel);
  void Display_histo_fast(int width_user, int* ch_comp);
  void Display_values(int* ch_comp);
  
private:
  bool Fit_values(int width, int d_cut_min, int d_cut_max);
  bool Fit_values_erfc(int width, int d_cut_min, int d_cut_max, int cut_db_pulses_user);
  bool Calc_values(int width, int d_cut_min, int d_cut_max); 
  void Fitting_Fast(int width);
  bool Soft_val(bool soft_flag=true);
  
public:
  ClassDef(trim_adc,0)

};
#endif
