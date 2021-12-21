#include <fstream>
#include <string>

using namespace std;

void make_trim_raise_LowADC(const char* filename){
  ifstream in(filename);
  char buf[1024];
  
  int threshold[128];
  while(true){
    in.getline(buf,1024);
    if ( in.eof() ) break;
    int ch,var;
    int results = sscanf(buf,"%d %d",&ch,&var);
    if ( results != 2 ){
      cout << "strange format : " << buf << endl;
      continue;
    }
    threshold[ch] = var;
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
  
  for( int ch = 0;ch < 128; ch++){
    string str;
    str += string(TString::Format("ch:%4d",ch));
    for(int i = 0;i<31; i++){
      int adj = 128;
      if ( i >= threshold[ch] ) adj = 1;
      str += string(TString::Format("%5d",adj));
    }
    int fast_trim = 36;
    str += string(TString::Format("%5d",fast_trim)); // dummy entry for FASTth.
    out << str << endl;
  }
  out.close();
}
