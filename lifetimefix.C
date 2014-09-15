#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <iterator>
#include "TString.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"

using namespace std;

TRandom randomForLifetime;

int numdiff=0;
vector<TString> fixlines(vector<string> lines, float lifetime);
string fixformat(Double_t val);
bool linecheck(TString before, TString after);
float exponentialLifetime(float properLifetime);

void lifetimefix(std::string fileloc, float lifetime){
  gStyle->SetOptStat(0);
  unsigned sz = fileloc.size();
  string outbase=fileloc;
  outbase.resize(sz-4);

  char lifetimestr[4];
  sprintf(lifetimestr,"%i",int(lifetime));
  lifetime/=1000;

  std::string line;
  ifstream colorfile;
  colorfile.open(fileloc.c_str());
  ofstream fixedfile;
  char oname[250];
  sprintf(oname,"%s_lifetime%sum_thrown.lhe",outbase.c_str(),lifetimestr);
  fixedfile.open(oname);
  
  //string test;
  //string prevline;
  int prevnum=0;
  vector<string> linestofix;
  vector<TString> fixedlines;
  int a=0;
  while (true){
    //if(a>10) assert(0);
    getline(colorfile,line);
    if(colorfile.eof()) break;
    istringstream linestr(line);
    string testval;
    linestr>>testval;
    int npart;
    if(!(istringstream(testval) >> npart)) npart=-99;
    if(npart==-99 || npart==2212 || prevnum==2212){
      fixedfile<<line<<endl; // Output to file
      prevnum=npart;
    }
    else{
      if(prevnum==-99){
        linestofix.push_back(line);
        for(int i=0;i<npart;i++){
          getline(colorfile,line);
          linestofix.push_back(line);
        }
        fixedlines = fixlines(linestofix,lifetime);
        for(int i=0;i<npart+1;i++){
          fixedfile<<fixedlines[i]<<endl;
        }
        for(int i=0;i<npart+1;i++){
          linestofix.pop_back();
        }
        if(a%1000==0) cout<<"Event number: "<<a<<endl;
        prevnum=npart;
        a++;
      }
      else{
	      cout<<"PROBLEM"<<endl;
      }
    }
  }
  cout<<a<<" events processed"<<endl;
}

vector<TString> fixlines(vector<string> lines, float lifetime){
  vector<TString> Tlines;
  for(unsigned int k=0;k<lines.size();k++){
    Tlines.push_back((TString)lines[k]);
  }
  
  int idup[15], istup[15], mothup[15][2], icolup[15][2];
  float pup[15][5],vtimup[15],spinup[15];
  int npart,procid;
  float weight, scale, alphaqed, alphaqcd;
  istringstream linevals;
  linevals.str(lines[0]);
  linevals >> npart >> procid >> weight >> scale >> alphaqed >> alphaqcd; 
  for(int k=0;k<15;k++){
    if(k<npart){
      istringstream onelinevals(lines[k+1]);
      onelinevals >> idup[k] >> istup[k] >> mothup[k][0] >> mothup[k][1] >> icolup[k][0] >> icolup[k][1];
      for(int i=0;i<5;i++){
        onelinevals >> pup[k][i];
      }
      onelinevals >> vtimup[k] >> spinup[k];
    }
    else{
      idup[k]=0;
      istup[k]=0;
      mothup[k][0]=0;
      mothup[k][1]=0;
      icolup[k][0]=0;
      icolup[k][1]=0;
      for(int i=0;i<5;i++) pup[k][i]=0.;
      vtimup[k]=0.;
      spinup[k]=0.;
    }
  }

  int Higgs=0;
  for(int k=0;k<npart;++k){
    if(abs(idup[k])==25) Higgs=k;
  }

  TString temp;
  //Using proper lifetime
  //string fix1 = fixformat(lifetime);
  //Using thrown lifetime
  string fix1 = fixformat(exponentialLifetime(lifetime));

  //if(fix1.size()!=13){ //POWHEG Base
  if(fix1.size()!=15){ //POWHEG+JHUGen
    cout<<fix1<<endl;
    cout<<"Length of lifetime string is wrong. Exiting..."<<endl;
    assert(0);
  }
  temp=Tlines[Higgs+1];
  //POWHEG Base
  //Tlines[Higgs+1].Replace(123,13,fix1.c_str());
  //JHUGen w/ or w/o POWHEG
  Tlines[Higgs+1].Replace(95,15,fix1.c_str());  

  return Tlines;
}

string fixformat(Double_t val){
  std::stringstream fixstream;

  char buffer[10]="";
  char expstr[1]="";
  int expval = int(floor(log10(fabs(val))));
  if(fabs(int(log10(fabs(val)))-log10(fabs(val)))<0.00001){
    if(expval<=0) expval++;
    //else if (expval>0) expval--;
  }
  float dec = fabs(val)/pow(10,expval);
  //POWHEG Base
  //sprintf(buffer,"%.5f",dec);
  //JHUGen w/ or w/o POWHEG
  sprintf(buffer,"%.7f",dec);
  sprintf(expstr,"%i",abs(expval));
  if(val>=0.) fixstream<<"  ";
  if(val<0.) fixstream<<" -";
  fixstream<<buffer;
  if(expval>=0) fixstream<<"E+0";
  if(expval<0) fixstream<<"E-0";
  fixstream<<expstr;

  return fixstream.str();
}

bool linecheck(TString before, TString after){
  bool okay=true;
  TString pxbefore = TString(before(40,3));
  TString pxafter = TString(after(40,3));
  TString pybefore = TString(before(54,3));
  TString pyafter = TString(after(54,3));
  TString pzbefore = TString(before(68,3));
  TString pzafter = TString(after(68,3));
  if(pxbefore!=pxafter || pybefore!=pyafter || pzbefore!=pzafter) okay=false;

  if(!okay){
    cout<<before<<endl;
    cout<<after<<endl;
  }

  return okay;
}

float exponentialLifetime(float properLifetime){
  float expLifetime;
  expLifetime = randomForLifetime.Exp(properLifetime);
  return expLifetime;
}
