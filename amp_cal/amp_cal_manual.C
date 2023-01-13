TGraph* gr;
void amp_cal_manual()
{
  gr = new TGraph("amp_cal_manual.txt");
  gr->GetXaxis()->SetTitle("amp_cal");
  gr->GetYaxis()->SetTitle("Qin[fC]");
  gr->SetMarkerStyle(20);
  gr->Draw("ALP");
}
