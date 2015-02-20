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
#include "TSystem.h"

using namespace std;

vector<TString> fixlines(vector<string> lines);

void scalefix(string fileloc){
  unsigned sz = fileloc.size();
  string outbase=fileloc;
  outbase.resize(sz-4);

  std::string line;
  ifstream colorfile;
  colorfile.open(fileloc.c_str());
  ofstream fixedfile;
  char oname[250];
  sprintf(oname,"%s_new.lhe",outbase.c_str());
  fixedfile.open(oname);
  
  string test;
  string prevline;
  int prevnum;
  vector<string> linestofix;
  vector<TString> fixedlines;
  int a=0;
  while (true){
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
	fixedlines = fixlines(linestofix);
	for(int i=0;i<npart+1;i++){
	  fixedfile<<fixedlines[i]<<endl;
	}
	//linestofix.erase(linestofix.begin(),linestofix.begin()+npart);
	for(int i=0;i<npart+1;i++){
	  linestofix.pop_back();
	}
	prevnum=npart;
	a++;
      }
      else{
	cout<<"PROBLEM"<<endl;
      }
    }
  }
  cout<<a<<" events processed"<<endl;

  char command[250];
  sprintf(command,"rm %s",fileloc.c_str());
  gSystem->Exec(command);
  sprintf(command,"mv %s %s",oname,fileloc.c_str());
  gSystem->Exec(command);  
}

vector<TString> fixlines(vector<string> lines){
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

  string scalefix;
  char buffer[7];
  char expstr[1];
  int exp = int(floor(log10(scale*2)));
  float dec = 2*scale/pow(10,exp);
  sprintf(buffer,"%.5f",dec);
  sprintf(expstr,"%i",exp);

  scalefix+=buffer;
  scalefix+="E+0";
  scalefix+=expstr;
  
  if(scalefix.size()!=11) cout<<"WARNING"<<endl;
  Tlines[0].Replace(29,11,scalefix.c_str());

  return Tlines;
}
