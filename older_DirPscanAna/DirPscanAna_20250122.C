// edited by S.Ochiai 2023.10.07
/**
 !!! allenc.root should be overwritten, but for some reason,
     it doesnt get updated unless you delete it first.
 When you provide the directory containing the root files of pscan as argumnents,
 this code create a hist of all ASIC's ENC values sorted by ASIC number,
 and PDFs of the pscan results, trim files for each ASIC individually.

 root allenc.root
 c_all->Draw();
**/
#include <fstream>
#include <string>
#include <TMath.h>
#include <dirent.h>
#include <picojson/picojson.h>

using Vector4D = std::vector<std::vector<std::vector<std::vector<double>>>>;

double binomial_error(double entry){
    double tmp = 0;
    if(entry == 0){
        tmp = 1;
    }else if(entry == 200){
        tmp = 199;
    }else{
        tmp = entry;
    }
    double p = tmp/200.;
    int n = 200;

    return sqrt(n*p*(1-p));
}

double ENC(double amp_cal){
    double qin =0.056*amp_cal; //sigma[amp_cal]->sigma[fC]
    double enc = qin / 1.6e-4; //ENC = e??? noise coulomb

    return enc;
}

int ConvertModuleID(string str_module_id) {

    std::vector<string> k_module_ids = {"101", "102", "103", "104", "106", "107", "108", "109", "206", "207"};

    for (int i = 0; i < k_module_ids.size(); ++i) {
        if (k_module_ids[i] == str_module_id) {
            return i;
        }
    }
    std::cerr << "Invalid module_id: " << str_module_id << std::endl;
    return -1;
}

void DirPscanAna(const char* dirpath){
    std::vector<string> K_MODULE_IDs = {"101", "102", "103", "104", "106", "107", "108", "109", "206", "207"};
    const int NumOfModules = 10;
    const int NumOfStrips = 128;
    const int NumOfAllStrips = 1024;
    const int NumOfADCch = 32;
    const int NumOfPol = 2;
    int RangeY = 5000;


    //adc_calib
    //std::vector<std::vector<std::vector<std::vector<float>>>> calib_data;//mod,p/n,ch,adc,fc
    //Vector4D calib_data(NumOfModules, std::vector<std::vector<std::vector<std::vector<double>>>>(2,std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(2, std::vector<double>(2, 0.0)))));

    Vector4D calib_data(NumOfModules, std::vector<std::vector<std::vector<double>>>(NumOfPol,std::vector<std::vector<double>>(NumOfAllStrips, std::vector<double>(NumOfADCch, -1))));

    //adc_thr/////////////////////////////////////////////////////////
    TH2D* h_all_N_adc_thr = new TH2D("h_all_N_adc_thr","h_all_N_adc_thr", NumOfAllStrips,0,NumOfAllStrips,100,0,100);
    h_all_N_adc_thr->GetXaxis()->SetTitle("CH");
    h_all_N_adc_thr->GetYaxis()->SetTitle("amp_cal");
    h_all_N_adc_thr->SetMarkerStyle(8);
    h_all_N_adc_thr->SetMarkerSize(0.8);
    h_all_N_adc_thr->SetMarkerColor(kRed);

    TH2D* h_all_P_adc_thr = new TH2D("h_all_P_adc_thr","h_all_P_adc_thr", NumOfAllStrips,0,NumOfAllStrips,100,0,100);
    h_all_P_adc_thr->GetXaxis()->SetTitle("CH");
    h_all_P_adc_thr->GetYaxis()->SetTitle("amp_cal");
    h_all_P_adc_thr->SetMarkerStyle(8);
    h_all_P_adc_thr->SetMarkerSize(0.8);
    h_all_P_adc_thr->SetMarkerColor(kRed);

    //fast_thr/////////////////////////////////////////////////////////////
    TH2D* h_all_N_fast_thr = new TH2D("h_all_N_fast_thr","h_all_N_fast_thr", NumOfAllStrips,0,NumOfAllStrips,100,0,100);
    h_all_N_fast_thr->GetXaxis()->SetTitle("CH");
    h_all_N_fast_thr->GetYaxis()->SetTitle("amp_cal");
    h_all_N_fast_thr->SetMarkerStyle(8);
    h_all_N_fast_thr->SetMarkerSize(0.8);
    h_all_N_fast_thr->SetMarkerColor(kBlue);

    TH2D* h_all_P_fast_thr = new TH2D("h_all_P_fast_thr","h_all_P_fast_thr", NumOfAllStrips,0,NumOfAllStrips,100,0,100);
    h_all_P_fast_thr->GetXaxis()->SetTitle("CH");
    h_all_P_fast_thr->GetYaxis()->SetTitle("amp_cal");
    h_all_P_fast_thr->SetMarkerStyle(8);
    h_all_P_fast_thr->SetMarkerSize(0.8);
    h_all_P_fast_thr->SetMarkerColor(kBlue);

    TH2D* h_all_N_fast_thr_dcnt = new TH2D("h_all_N_fast_thr_dcnt","h_all_N_fast_thr_dcnt", NumOfAllStrips,0,NumOfAllStrips,100,0,100);
    h_all_N_fast_thr_dcnt->GetXaxis()->SetTitle("CH");
    h_all_N_fast_thr_dcnt->GetYaxis()->SetTitle("amp_cal");
    h_all_N_fast_thr_dcnt->SetMarkerStyle(8);
    h_all_N_fast_thr_dcnt->SetMarkerSize(0.8);
    h_all_N_fast_thr_dcnt->SetMarkerColor(kGreen);

    TH2D* h_all_P_fast_thr_dcnt = new TH2D("h_all_P_fast_thr_dcnt","h_all_P_fast_thr_dcnt", NumOfAllStrips,0,NumOfAllStrips,100,0,100);
    h_all_P_fast_thr_dcnt->GetXaxis()->SetTitle("CH");
    h_all_P_fast_thr_dcnt->GetYaxis()->SetTitle("amp_cal");
    h_all_P_fast_thr_dcnt->SetMarkerStyle(8);
    h_all_P_fast_thr_dcnt->SetMarkerSize(0.8);
    h_all_P_fast_thr_dcnt->SetMarkerColor(kGreen);

    //adc_noise/////////////////////////////////////////////////////////
    TH2D* h_all_N_adc_noise = new TH2D("h_all_N_adc_noise","h_all_N_adc_noise", NumOfAllStrips,0,NumOfAllStrips,RangeY,0,RangeY);
    h_all_N_adc_noise->GetXaxis()->SetTitle("CH");
    h_all_N_adc_noise->GetYaxis()->SetTitle("enc");
    h_all_N_adc_noise->SetMarkerStyle(8);
    h_all_N_adc_noise->SetMarkerSize(0.8);
    h_all_N_adc_noise->SetMarkerColor(kRed);

    TH2D* h_all_P_adc_noise = new TH2D("h_all_P_adc_noise","h_all_P_adc_noise", NumOfAllStrips,0,NumOfAllStrips,RangeY,0,RangeY);
    h_all_P_adc_noise->GetXaxis()->SetTitle("CH");
    h_all_P_adc_noise->GetYaxis()->SetTitle("enc");
    h_all_P_adc_noise->SetMarkerStyle(8);
    h_all_P_adc_noise->SetMarkerSize(0.8);
    h_all_P_adc_noise->SetMarkerColor(kRed);

    //fast_noise/////////////////////////////////////////////////////////////
    TH2D* h_all_N_fast_noise = new TH2D("h_all_N_fast_noise","h_all_N_fast_noise", NumOfAllStrips,0,NumOfAllStrips,RangeY,0,RangeY);
    h_all_N_fast_noise->GetXaxis()->SetTitle("CH");
    h_all_N_fast_noise->GetYaxis()->SetTitle("enc");
    h_all_N_fast_noise->SetMarkerStyle(8);
    h_all_N_fast_noise->SetMarkerSize(0.8);
    h_all_N_fast_noise->SetMarkerColor(kBlue);

    TH2D* h_all_P_fast_noise = new TH2D("h_all_P_fast_noise","h_all_P_fast_noise", NumOfAllStrips,0,NumOfAllStrips,RangeY,0,RangeY);
    h_all_P_fast_noise->GetXaxis()->SetTitle("CH");
    h_all_P_fast_noise->GetYaxis()->SetTitle("enc");
    h_all_P_fast_noise->SetMarkerStyle(8);
    h_all_P_fast_noise->SetMarkerSize(0.8);
    h_all_P_fast_noise->SetMarkerColor(kBlue);

    TH2D* h_all_N_fast_noise_dcnt = new TH2D("h_all_N_fast_noise_dcnt","h_all_N_fast_noise_dcnt", NumOfAllStrips,0,NumOfAllStrips,RangeY,0,RangeY);
    h_all_N_fast_noise_dcnt->GetXaxis()->SetTitle("CH");
    h_all_N_fast_noise_dcnt->GetYaxis()->SetTitle("enc");
    h_all_N_fast_noise_dcnt->SetMarkerStyle(8);
    h_all_N_fast_noise_dcnt->SetMarkerSize(0.8);
    h_all_N_fast_noise_dcnt->SetMarkerColor(kGreen);

    TH2D* h_all_P_fast_noise_dcnt = new TH2D("h_all_P_fast_noise_dcnt","h_all_P_fast_noise_dcnt", NumOfAllStrips,0,NumOfAllStrips,RangeY,0,RangeY);
    h_all_P_fast_noise_dcnt->GetXaxis()->SetTitle("CH");
    h_all_P_fast_noise_dcnt->GetYaxis()->SetTitle("enc");
    h_all_P_fast_noise_dcnt->SetMarkerStyle(8);
    h_all_P_fast_noise_dcnt->SetMarkerSize(0.8);
    h_all_P_fast_noise_dcnt->SetMarkerColor(kGreen);

    TString out_thr_noise_pdf = string(dirpath) + "thr_noise.pdf";
    TString N_out_root = string(dirpath) + "rlt_N_out.root";
    TString P_out_root = string(dirpath) + "rlt_P_out.root";
    TString adc_calib = string(dirpath) + "ADCcalib.dat";

    DIR *dir;
    struct dirent *dp;
    dir = opendir(dirpath);
    if ( dir == NULL ) { return 1; };
    dp = readdir(dir);

    int Addr2AsicNo_PA_Nside[8] = {1, 0, 3, 2, 5, 4, 7, 6};
    int Addr2AsicNo_PA_Pside[8] = {7, 6, 5, 4, 3, 2, 1, 0};
    int Addr2AsicNo_PB_Nside[8] = {7, 6, 5, 4, 3, 2, 1, 0};
    int Addr2AsicNo_PB_Pside[8] = {1, 0, 3, 2, 5, 4, 7, 6};

    string RootFilePath;
    string GetAddr;
    int Addr;
    int AsicNo;
    const char* filename;
    string EfuseStr;
    string Pol;
    string Module_ID;
    string Module_type;

	//Open EfuseMap.json
    //ifstream ifs("./EfuseMap.json", ios::in);
    // ifstream ifs("./EfuseMap_run0e.json", ios::in);
    ifstream ifs("./EfuseMap_20250122_2021.json", ios::in);
    if(ifs.fail()){
        //cerr << "Failed to read EfuseMap.json" << endl;
        cerr << "Failed to read EfuseMap_20250122_2021.json" << endl;
        return 1;
    }
    const string json((istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
    ifs.close();
    picojson::value v;
    const string err = picojson::parse(v, json);
    if(err.empty() == false){
        cerr << err << endl;
        return 2;
    }
    picojson::object& obj = v.get<picojson::object>();


    while( dp != NULL ){//--1111111-- Dir Scan Loop Layer1 --1111111--
	    filename = dp->d_name;
        string RootFilePath = strrchr(filename, '.');
        if( RootFilePath == ".root"){//--222222-- if(root file) Layer2 --2222222--

            //Get info from file name
            GetAddr = (string)filename;
            GetAddr = GetAddr[GetAddr.find("addr")+4];
            Addr = stoi(GetAddr);
            EfuseStr = (string)filename;
            size_t found = EfuseStr.find("XA");
            if (found != std::string::npos) {
               EfuseStr = EfuseStr.substr(found, 28);
            } else {
            // "XA" が見つからない場合の処理
             EfuseStr = EfuseStr.substr(EfuseStr.find("na"), 28);
            }
            RootFilePath = (string)dirpath + (string)filename;
            filename = RootFilePath.c_str();

            //Get info from JSON file
            Pol = obj[EfuseStr].get<picojson::object>()["Pol"].get<string>();
            Module_ID = obj[EfuseStr].get<picojson::object>()["Module_ID"].get<string>();
            Module_type = Module_ID.substr(Module_ID.find("P"), 2);

            //convert Addr to Asic No
            if(Module_type == "PA"){
                 if(Pol=="N-side"){
                     AsicNo = Addr2AsicNo_PA_Nside[Addr];
                 }
                 else if(Pol=="P-side"){
                     AsicNo = Addr2AsicNo_PA_Pside[Addr];
                 }
            }
            else if(Module_type == "PB"){
                 if(Pol=="N-side"){
                     AsicNo = Addr2AsicNo_PB_Nside[Addr];
                 }
                 else if(Pol=="P-side"){
                     AsicNo = Addr2AsicNo_PB_Pside[Addr];
                 }
            }

            TCanvas* c1 = new TCanvas("c1","c1",800,800);
            TFile* file = new TFile(filename);
            cout << file << endl;

            TH2F *sig_ch[NumOfStrips];
            for(int i=0;i<NumOfStrips;i++){
                sig_ch[i] = new TH2F("sig_ch","sig_ch",32,0,32,60,0,60);
            }

            TH1D *sig_adc_noise_hist[NumOfStrips];
            for(int i=0;i<NumOfStrips;i++){
                sig_adc_noise_hist[i] = new TH1D("sig_adc_noise_hist","sig_adc_noise_hist",60,0,60);
            }

            TString filename_trim = filename;
            int ipos = filename_trim.Last('.');
            if ( ipos >= 0 ) {
                filename_trim.Remove(ipos);
                filename_trim = filename_trim + "_trim.txt";
            }else{
                cout << "unexpected filename" << endl;
                return 1;
            }

            ofstream out_trim(filename_trim);
			TString out_calib_adc_name = filename;
            out_calib_adc_name.Remove(ipos);
            out_calib_adc_name = out_calib_adc_name + "_ADCcalib.txt";
            ofstream out_calib_adc(out_calib_adc_name);

            int ch = 0;
            int ch_max = NumOfStrips-1;
            //std::cout << "==========================ch_max: " << ch_max << std::endl;
            int excepted_sigma_cnt = 0; //sum of cutted sigma

            TString out_pdf=filename;
            out_pdf.Remove(ipos);
            out_pdf = out_pdf + ".pdf";

            // SHOW h_quality_%d.
            c1->Print(out_pdf + "[","pdf");

            TF1* func_calib = new TF1("func_calib","pol1");
            // 20-240.
            //func_calib->SetParameters(200,-6);
            func_calib->SetParameters(240,-7.33333333);
            const double trim_amp = 18/12.5;

            TH2D* hist_ch_fast_thr_dcnt = new TH2D("hist_ch_fast_thr_dcnt","hist_ch_fast_thr_dcnt", 128,0,128,100,0,100);
            hist_ch_fast_thr_dcnt->GetXaxis()->SetTitle("CH");
            hist_ch_fast_thr_dcnt->GetYaxis()->SetTitle("amp_cal");
            hist_ch_fast_thr_dcnt->SetMarkerStyle(8);
            hist_ch_fast_thr_dcnt->SetMarkerSize(0.8);
            hist_ch_fast_thr_dcnt->SetMarkerColor(kGreen);

            TH2D* hist_ch_fast_noise_dcnt = new TH2D("hist_ch_fast_noise_dcnt","hist_ch_fast_noise_dcnt", 128,0,128,5000,0,5000);//
            hist_ch_fast_noise_dcnt->GetXaxis()->SetTitle("CH");
            hist_ch_fast_noise_dcnt->GetYaxis()->SetTitle("enc");
            hist_ch_fast_noise_dcnt->SetMarkerStyle(8);
            hist_ch_fast_noise_dcnt->SetMarkerSize(0.8);
            hist_ch_fast_noise_dcnt->SetMarkerColor(kGreen);

            //trim & calib///////////////////////////////////////////////////////////////////
            while(ch <= ch_max){//--3333333-- Layer3-1 --3333333--
                for( int itmp = 0; itmp < 16;itmp ++){
                    if (itmp == 0 ) {
                        c1->Clear();
                        c1->Divide(4,4);
                    }
                    c1->cd(itmp+1);

                    /// FastTrim
                    TString name = TString::Format("h_d_%d_31",ch);
                    TH1D* hist_fast = (TH1D*)file->Get(name);
                    hist_fast->Fit("gaus","LQ0");
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
					//------------added by R.Yamada 2023.11.11
                    TF1 *fitFunc = new TF1("fitFunc", "pol1",5,25);
                    hist->Fit("fitFunc");
                    fitFunc->Draw("same");
					double intere =0.056*fitFunc->GetParameter(0);
                    double slope =0.056*fitFunc->GetParameter(1);
                    string str_ch;
                    str_ch = str;
                    for(int i=0;i<31; i++){
                        double charge_fc = (30-i)*slope+intere;
                        str_ch += string(TString::Format(" %9.6f",charge_fc));
                        string str_charge_fc = string(TString::Format(" %4.2f",charge_fc));//小数点以下2桁
                        float charge_fc_2 = std::stof(str_charge_fc);
                        string Module_name = Module_ID.substr(0, 3);
                        int k_mod = ConvertModuleID(Module_name);
                        if(Pol=="N-side"){
                            calib_data[k_mod][1][ch+NumOfStrips*AsicNo][i+1] = charge_fc;
                        }else if(Pol=="P-side"){
                            calib_data[k_mod][0][ch+NumOfStrips*AsicNo][i+1] = charge_fc;
                        }
                    }
                    out_calib_adc << str_ch << endl;

                    for(int i = 0;i<31; i++){
                        int ibin = hist->FindBin(i+0.001);
                        double var = hist->GetBinContent(ibin);
                        int adj = 128.-(func_calib->Eval(i) - var)*trim_amp;
                        str += string(TString::Format(" %5d",adj));
                    }

                    str += string(TString::Format(" %5d",fast_trim)); // dummy entry for FASTth
                    out_trim << str << endl;
                    ch++;
                    if ( ch > ch_max ) break;
                }
                c1->Draw();
                c1->Update();
                c1->Print(out_pdf,"PDF");
            }//--3333333 Layer3-1 END --3333333--
            /*c1->Clear();
            c1->cd(1);
            hist_ch_fast_thr_dcnt->Draw();
            c1->Update();
            c1->Print(out_pdf,"PDF");

            c1->Clear();
            c1->cd(1);
            hist_ch_fast_noise_dcnt->Draw();
            c1->Update();
            c1->Print(out_pdf,"PDF");*/

            // SHOW scurve and dcnt
            TH2D* hist_ch_adc_thr = new TH2D("hist_ch_adc_thr","hist_ch_adc_thr", 128,0,128,100,0,100);
            hist_ch_adc_thr->GetXaxis()->SetTitle("CH");
            hist_ch_adc_thr->GetYaxis()->SetTitle("amp_cal");
            hist_ch_adc_thr->SetMarkerStyle(8);
            hist_ch_adc_thr->SetMarkerSize(0.8);
            hist_ch_adc_thr->SetMarkerColor(kRed);

            TH2D* hist_ch_fast_thr = new TH2D("hist_ch_fast_thr","hist_ch_fast_thr", 128,0,128,100,0,100);
            hist_ch_fast_thr->GetXaxis()->SetTitle("CH");
            hist_ch_fast_thr->GetYaxis()->SetTitle("amp_cal");
            hist_ch_fast_thr->SetMarkerStyle(8);
            hist_ch_fast_thr->SetMarkerSize(0.8);
            hist_ch_fast_thr->SetMarkerColor(kBlue);

            TH2D *hist_ch_adc_noise = new TH2D("hist_ch_adc_noise","hist_ch_adc_noise",128,0,128,5000,0,5000);//
            hist_ch_adc_noise->GetXaxis()->SetTitle("CH");
            hist_ch_adc_noise->GetYaxis()->SetTitle("enc");
            hist_ch_adc_noise->SetMarkerStyle(8);
            hist_ch_adc_noise->SetMarkerSize(0.8);
            hist_ch_adc_noise->SetMarkerColor(kRed);

            TH2D* hist_ch_fast_noise = new TH2D("hist_ch_fast_noise","hist_ch_fast_noise", 128,0,128,5000,0,5000);//
            hist_ch_fast_noise->GetXaxis()->SetTitle("CH");
            hist_ch_fast_noise->GetYaxis()->SetTitle("enc");
            hist_ch_fast_noise->SetMarkerStyle(8);
            hist_ch_fast_noise->SetMarkerSize(0.8);
            hist_ch_fast_noise->SetMarkerColor(kBlue);

            int d_max = 31;
            for( int ch = 0;ch<=ch_max;ch++){//--3333333-- Layer3-2 --3333333--
                int d = 0;
                while(d <= d_max){//--4444444-- Layer4-1 --4444444--
                    for( int itmp = 0;itmp<8;itmp++){//--5555555-- Layer5-1 --5555555--
                    if ( itmp==0 ) {
                            c1->Clear();
                            c1->Divide(2,4);
                        }
                        c1->cd(itmp+1);

                        TString name = TString::Format("h_scurve_%d_%d",ch,d);
                        TH1D* hist_scurve = (TH1D*)file->Get(name);
                        hist_scurve->SetTitle(name);

                        //give the binomial error : eddited by S.Ochiai 2023.03.01
                        Double_t error[300];
                        int nbins = hist_scurve->GetXaxis()->GetNbins();
                        for(int i=1; i<nbins; i++){
                            double y = hist_scurve->GetBinContent(i);
                            error[i-1] = binomial_error(y);
                        }
                        hist_scurve->SetError(error); //Set Binomial Error to BinError

                        TF1* erf = new TF1("erf", "200*ROOT::Math::gaussian_cdf(x,[0],[1])",0,300);
                        erf->SetParameters(20,160);
                        hist_scurve->Rebin(4);
                        hist_scurve->Draw("HIST");
                        hist_scurve->Fit("erf","LQ","",0,255);//Fitting #1

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
                            hist_scurve->Fit("erf", "LQ", "", xlow, 255);
                            chi2_current = erf->GetChisquare();
                            if( chi2_current > chi2_previous || TMath::Abs( chi2_previous - chi2_current ) < chi2_previous * 0.01){
                                break;
                            }
                         } //+++++++++Fitting Iteration end+++++

                        //-----------------added by S.Ochiai 2023.02,03 end-------------
                        gStyle->SetOptStat(1);
                        gStyle->SetOptFit(1);

                        double chi2 = chi2_previous;
                        int NDF = erf->GetNDF();//自由度
                        int bin1 = hist_scurve->GetBinContent(1);
                        sig_ch[ch]->Fill(d,sigma);

                        if(chi2/NDF <= 20 && bin1 < 200 && sigma<15 && d!=31){//chi2/NDF <= 20 ->モデルの適応具合
                            sig_adc_noise_hist[ch]->Fill(sigma);//adcのsigmaのみ。d=0-30を射影
                            if(d==30){
                                hist_ch_adc_thr->Fill(ch,x0);
                                if(Pol=="N-side"){
                                    h_all_N_adc_thr->Fill(ch+NumOfStrips*AsicNo, x0);
                                }else if(Pol=="P-side"){
                                    h_all_P_adc_thr->Fill(ch+NumOfStrips*AsicNo, x0);
                                }
                            }
                        }else{excepted_sigma_cnt++;}

                        if(d==d_max){
                            double sigma_enc = ENC(sigma);
                            hist_ch_fast_thr->Fill(ch,x0);
                            hist_ch_fast_noise->Fill(ch,sigma_enc);
                            if(Pol=="N-side"){
                                h_all_N_fast_thr->Fill(ch+NumOfStrips*AsicNo, x0);
                                h_all_N_fast_noise->Fill(ch+NumOfStrips*AsicNo, sigma_enc);
                            }else if(Pol=="P-side"){
                                h_all_P_fast_thr->Fill(ch+NumOfStrips*AsicNo, x0);
                                h_all_P_fast_noise->Fill(ch+NumOfStrips*AsicNo, sigma_enc);
                            }

                            TString name = TString::Format("h_d_%d_31",ch);
                            TH1D* hist_fast = (TH1D*)file->Get(name);
                            hist_fast->Rebin(4);
                            hist_fast->Fit("gaus","LQR0","",25,100);
                            TF1* func = hist_fast->GetFunction("gaus");
                            int fast_trim = 36;
                            if ( func != NULL ){
                                double var = func->GetParameter(0);
                                double fast_thr_dcnt = func->GetParameter(1);
                                double fast_noise_dcnt_amp = func->GetParameter(2);
                                double fast_noise_dcnt = ENC(fast_noise_dcnt_amp);
                                hist_ch_fast_thr_dcnt->Fill(ch,fast_thr_dcnt);
                                hist_ch_fast_noise_dcnt->Fill(ch,fast_noise_dcnt);
                                if(Pol=="N-side"){
                                    h_all_N_fast_thr_dcnt->Fill(ch+NumOfStrips*AsicNo, fast_thr_dcnt);
                                    h_all_N_fast_noise_dcnt->Fill(ch+NumOfStrips*AsicNo, fast_noise_dcnt);
                                }else if(Pol=="P-side"){
                                    h_all_P_fast_thr_dcnt->Fill(ch+NumOfStrips*AsicNo, fast_thr_dcnt);
                                    h_all_P_fast_noise_dcnt->Fill(ch+NumOfStrips*AsicNo, fast_noise_dcnt);
                                }

                                TF1* dcnt_gaus = new TF1("dcnt_gaus", "gaus", 25, 100);
                                dcnt_gaus->SetParameters(var,fast_thr_dcnt,fast_noise_dcnt_amp);
                                dcnt_gaus->SetLineColor(kGreen);
                                hist_fast->SetLineColor(kBlack);
                                hist_fast->Draw("same HIST");
                                dcnt_gaus->Draw("same");
                            }
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

                        d++;
                        if ( d > d_max ) break;
		            }//--5555555-- Layer5-1 END --5555555--
                    c1->Draw();
                    c1->Update();
                    c1->Print(out_pdf,"PDF");
                }//--4444444-- Layer4-1 END --4444444--
            }//--3333333-- Layer3-2 END -- 3333333--

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
                c1->Print(out_pdf,"PDF");
            }

            //Noise
            ch = 0;

            while(ch <= ch_max){
                for(int itemp=0;itemp<8;itemp++){
                    if(itemp == 0){
                        c1->Clear();
                        c1->Divide(2,4);
                    }
                    c1->cd(itemp+1);
                    TString hist_name = TString::Format("sigma_adc_noise_ch_%d",ch);
                    sig_adc_noise_hist[ch]->Draw();
                    sig_adc_noise_hist[ch]->SetTitle(hist_name);
                    sig_adc_noise_hist[ch]->Fit("gaus","SLQR","",1,25);
                    TF1 *func3 = sig_adc_noise_hist[ch]->GetFunction("gaus");

                    if ( func3 != NULL ) {
                        double mean = func3->GetParameter(1);
                        double enc = ENC(mean);
                        hist_ch_adc_noise->Fill(ch,enc);
                        if(Pol=="N-side"){
                            h_all_N_adc_noise->Fill(ch+NumOfStrips*AsicNo, enc);
                        }else if(Pol=="P-side"){
                            h_all_P_adc_noise->Fill(ch+NumOfStrips*AsicNo, enc);
                        }
                    }
                    gStyle->SetOptStat(1);
                    gStyle->SetOptFit(1);
                    ch++;
                    if(ch > ch_max) break;
                }
                c1->Draw();
                c1->Update();
                c1->Print(out_pdf,"PDF");
            }

            //thresold///////////////////////////////////////////////////////////////
            c1->Clear();
            c1->cd(1);

            TH1F *bg_thr=new TH1F("bg_thr","adc&fast_thr",128,0,128);
            bg_thr->SetXTitle("strip_ch");
            bg_thr->SetYTitle("amp_cal");
            bg_thr->SetMaximum(100);
	        bg_thr->SetMinimum(0);
            bg_thr->SetStats(0);
            bg_thr->Draw();
            hist_ch_adc_thr ->Draw("same");
            hist_ch_fast_thr ->Draw("same");

            TLegend *lg_thr=new TLegend(0.5,0.15,0.75,0.3);
            lg_thr->AddEntry(hist_ch_adc_thr,"adc_thr","p");
            lg_thr->AddEntry(hist_ch_fast_thr,"fast_thr","p");
            lg_thr->SetBorderSize(0);
            lg_thr->SetFillColor(0);
            lg_thr->SetTextSize(0.04);
            lg_thr->SetBorderSize(1);
            lg_thr->Draw("same");
            c1->Update();
            c1->Print(out_pdf, "PDF");

            //noise///////////////////////////////////////////////////////////////////
            c1->Clear();
            c1->cd(1);

            TH1F *bg_noise=new TH1F("bg_thr","adc&fast_noise",128,0,128);
            bg_noise->SetXTitle("strip_ch");
            bg_noise->SetYTitle("ENC");
            bg_noise->SetMaximum(5000);
	        bg_noise->SetMinimum(0);
            bg_noise->SetStats(0);
            bg_noise->Draw();
            hist_ch_adc_noise->Draw("same");
            hist_ch_fast_noise->Draw("same");

            TLegend *lg_noise=new TLegend(0.5,0.55,0.75,0.75);
            lg_noise->AddEntry(hist_ch_adc_thr,"adc_noise","p");
            lg_noise->AddEntry(hist_ch_fast_thr,"fast_noise","p");
            lg_noise->SetBorderSize(0);
            lg_noise->SetFillColor(0);
            lg_noise->SetTextSize(0.04);
            lg_noise->SetBorderSize(1);
            lg_noise->Draw("same");

            c1->Update();
            c1->Print(out_pdf, "PDF");
            c1->Print(out_pdf +"]","PDF");

            out_trim.close();
            for (int ch=0; ch<NumOfStrips; ch++){
                 delete sig_ch[ch];
                 delete sig_adc_noise_hist[ch];
            }
            delete func_calib;
            delete c1;
            delete hist_ch_adc_thr;
            delete hist_ch_fast_thr;
            delete hist_ch_adc_noise;
            delete hist_ch_fast_noise;
            delete file;
        }//--2222222-- if(root file) Layer2 END --2222222--
	  dp = readdir(dir);
    }//--1111111-- Dir Scan Loop Layer1 END--1111111--
    if ( dir != NULL ){ closedir(dir); }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //TFile* AllADC_CALIB_N = new TFile(N_adc_calib, "RECREATE");
    ofstream adc_calib_outFile(adc_calib);
    if (!adc_calib_outFile) {
        std::cerr << "could not open file!" << std::endl;
        return 1;
    }

    // データを書き込む
    for (int module = 0; module < NumOfModules; ++module) {
        for (int pol = 0; pol < NumOfPol; ++pol) {
            for (int strip_ch = 0; strip_ch < NumOfAllStrips; ++strip_ch) {
                for (int adc_ch = 0; adc_ch < NumOfADCch; ++adc_ch) {
                    if (calib_data[module][pol][strip_ch][adc_ch] != -1) {
                        adc_calib_outFile << K_MODULE_IDs[module] << " " << pol << " " << strip_ch << " " << adc_ch << " " << calib_data[module][pol][strip_ch][adc_ch] << std::endl;
                    }
                }
            }
        }
    }

    adc_calib_outFile.close();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    TFile* AllRLT_N = new TFile(N_out_root, "RECREATE");

    TCanvas* c_all_N_th = new TCanvas("c_all_N_th", "c_all_N_th", 800, 600);
    c_all_N_th->cd(1);
    h_all_N_adc_thr->Draw("");
    h_all_N_fast_thr->Draw("same");
    h_all_N_fast_thr_dcnt->Draw("same");
    const int NumOfLines = 7;
    TLine* Line_N[NumOfLines];
    for(int i=0; i<NumOfLines; i++){
        Line_N[i] = new TLine(NumOfStrips*(i+1), 0, NumOfStrips*(i+1), 100);
        Line_N[i]->SetLineWidth(3);
        Line_N[i]->SetLineColor(kViolet);
        Line_N[i]->Draw("same");
    }
    c_all_N_th->Draw();
    c_all_N_th->Update();

    TCanvas* c_all_N_noise = new TCanvas("c_all_N_noise", "c_all_N_noise", 800, 600);
    c_all_N_noise->cd(1);
    h_all_N_adc_noise->Draw("");
    h_all_N_fast_noise->Draw("same");
    h_all_N_fast_noise_dcnt->Draw("same");
    TLine* Line_N_n[NumOfLines];
    for(int i=0; i<NumOfLines; i++){
        Line_N_n[i] = new TLine(NumOfStrips*(i+1), 0, NumOfStrips*(i+1), 5000);
        Line_N_n[i]->SetLineWidth(3);
        Line_N_n[i]->SetLineColor(kViolet);
        Line_N_n[i]->Draw("same");
    }
    c_all_N_noise->Draw();
    c_all_N_noise->Update();

    AllRLT_N->cd();
    h_all_N_adc_thr->Write();
    h_all_N_fast_thr->Write();
    h_all_N_adc_noise->Write();
    h_all_N_fast_noise->Write();
    c_all_N_th->Write();
    c_all_N_noise->Write();
    AllRLT_N->Close();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    TFile* AllRLT_P = new TFile(P_out_root, "RECREATE");

    TCanvas* c_all_P_th = new TCanvas("c_all_P_th", "c_all_P_th", 800, 600);
    c_all_P_th->cd();
    h_all_P_adc_thr->Draw("");
    h_all_P_fast_thr->Draw("same");
    h_all_P_fast_thr_dcnt->Draw("same");
    TLine* Line_P[NumOfLines];
    for(int i=0; i<NumOfLines; i++){
        Line_P[i] = new TLine(NumOfStrips*(i+1), 0, NumOfStrips*(i+1), 100);
        Line_P[i]->SetLineWidth(3);
        Line_P[i]->SetLineColor(kViolet);
        Line_P[i]->Draw("same");
    }
    c_all_P_th->Draw();
    c_all_P_th->Update();

    TCanvas* c_all_P_noise = new TCanvas("c_all_P_noise", "c_all_P_noise", 800, 600);
    c_all_P_noise->cd(1);
    h_all_P_adc_noise->Draw("");
    h_all_P_fast_noise->Draw("same");
    h_all_P_fast_noise_dcnt->Draw("same");
    TLine* Line_P_n[NumOfLines];
    for(int i=0; i<NumOfLines; i++){
        Line_P_n[i] = new TLine(NumOfStrips*(i+1), 0, NumOfStrips*(i+1), 5000);
        Line_P_n[i]->SetLineWidth(3);
        Line_P_n[i]->SetLineColor(kViolet);
        Line_P_n[i]->Draw("same");
    }
    c_all_P_noise->Draw();
    c_all_P_noise->Update();

    AllRLT_P->cd();
    h_all_P_adc_thr->Write();
    h_all_P_fast_thr->Write();
    h_all_P_adc_noise->Write();
    h_all_P_fast_noise->Write();
    c_all_P_th->Write();
    c_all_P_noise->Write();
    AllRLT_P->Close();

}
