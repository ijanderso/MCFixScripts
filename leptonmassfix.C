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

using namespace std;

const float Zmass = 91.19;
int numdiff=0;
TH1F* electrondiff1 = new TH1F("electrondiff1","electrondiff1",100,-0.02,0.02);
TH1F* electrondiff2 = new TH1F("electrondiff2","electrondiff2",100,-0.02,0.02);
TH1F* electrondiffinit1 = new TH1F("electrondiffinit1","electrondiffinit1",100,-0.02,0.02);
TH1F* electrondiffinit2 = new TH1F("electrondiffinit2","electrondiffinit2",100,-0.02,0.02);
TH2F* electrondiff12D = new TH2F("electrondiff12D","electrondiff12D",100,0.,200.,100,-0.02,0.02);
TH2F* electrondiff22D = new TH2F("electrondiff22D","electrondiff22D",100,0.,200.,100,-0.02,0.02);
TH2F* electrondiffinit12D = new TH2F("electrondiffinit12D","electrondiffinit12D",100,0.,200.,100,-0.02,0.02);
TH2F* electrondiffinit22D = new TH2F("electrondiffinit22D","electrondiffinit22D",100,0.,200.,100,-0.02,0.02);

TH1F* muondiff1 = new TH1F("muondiff1","muondiff1",100,-0.02,0.02);
TH1F* muondiff2 = new TH1F("muondiff2","muondiff2",100,-0.02,0.02);
TH1F* muondiffinit1 = new TH1F("muondiffinit1","muondiffinit1",100,-0.02,0.02);
TH1F* muondiffinit2 = new TH1F("muondiffinit2","muondiffinit2",100,-0.02,0.02);
TH2F* muondiff12D = new TH2F("muondiff12D","muondiff12D",100,0.,200.,100,-0.02,0.02);
TH2F* muondiff22D = new TH2F("muondiff22D","muondiff22D",100,0.,200.,100,-0.02,0.02);
TH2F* muondiffinit12D = new TH2F("muondiffinit12D","muondiffinit12D",100,0.,200.,100,-0.02,0.02);
TH2F* muondiffinit22D = new TH2F("muondiffinit22D","muondiffinit22D",100,0.,200.,100,-0.02,0.02);
vector<TString> fixlines(vector<string> lines);
string makefix(TLorentzVector lep, Double_t mass);
string fixformat(Double_t val);
bool linecheck(TString before, TString after);
TLorentzVector makemassless(TLorentzVector base);

void leptonmassfix(std::string fileloc){
  gStyle->SetOptStat(0);
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
        fixedlines = fixlines(linestofix);
        for(int i=0;i<npart+1;i++){
          fixedfile<<fixedlines[i]<<endl;
        }
        //linestofix.erase(linestofix.begin(),linestofix.begin()+npart);
        for(int i=0;i<npart+1;i++){
          linestofix.pop_back();
        }
        //cout<<a<<endl;
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
  cout<<numdiff<<" events with initial |m_cal-m_listed|/E > 1%"<<endl;

  /*muondiff1->Scale(1./muondiff1->Integral());
  muondiff2->Scale(1./muondiff2->Integral());
  muondiffinit1->Scale(1./muondiffinit1->Integral());
  muondiffinit2->Scale(1./muondiffinit2->Integral());
  muondiff1->SetTitle("");
  muondiff1->GetXaxis()->SetTitle("#Deltam/E");
  muondiff1->GetXaxis()->CenterTitle();

  TCanvas* c = new TCanvas("c","c",800,800);
  c->cd();
  muondiff1->SetLineWidth(2);
  muondiff1->SetLineColor(kRed);
  muondiff1->Draw();
  muondiff2->SetLineWidth(2);
  muondiff2->SetLineColor(kBlue);
  muondiff2->Draw("same");
  muondiffinit1->SetLineWidth(2);
  muondiffinit1->SetLineColor(kRed);
  muondiffinit1->SetLineStyle(2);
  muondiffinit1->Draw("same");
  muondiffinit2->SetLineWidth(2);
  muondiffinit2->SetLineColor(kBlue);
  muondiffinit2->SetLineStyle(2);
  muondiffinit2->Draw("same");
  c->SaveAs("muon_massdiscrepancies_frac_1D.png");
  c->SaveAs("muon_massdiscrepancies_frac_1D.eps");
  c->SaveAs("muon_massdiscrepancies_frac_1D.pdf");
  c->SaveAs("muon_massdiscrepancies_frac_1D.root");
  c->SaveAs("muon_massdiscrepancies_frac_1D.C");

  muondiff12D->Scale(1./muondiff12D->Integral());
  muondiff22D->Scale(1./muondiff22D->Integral());
  muondiffinit12D->Scale(1./muondiffinit12D->Integral());
  muondiffinit22D->Scale(1./muondiffinit22D->Integral());
  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->cd();
  muondiff12D->SetTitle("");
  muondiff12D->GetXaxis()->SetTitle("E_{#mu} [GeV]");
  muondiff12D->GetXaxis()->CenterTitle();
  muondiff12D->GetYaxis()->SetTitle("#Deltam/E");
  muondiff12D->GetYaxis()->CenterTitle();
  muondiff12D->Draw("COLZ");
  c2->SaveAs("muon1_massdiscrepancy_frac_2D.png");
  c2->SaveAs("muon1_massdiscrepancy_frac_2D.eps");
  c2->SaveAs("muon1_massdiscrepancy_frac_2D.pdf");
  c2->SaveAs("muon1_massdiscrepancy_frac_2D.root");
  c2->SaveAs("muon1_massdiscrepancy_frac_2D.C");
  muondiff22D->SetTitle("");
  muondiff22D->GetXaxis()->SetTitle("E_{#mu} [GeV]");
  muondiff22D->GetXaxis()->CenterTitle();
  muondiff22D->GetYaxis()->SetTitle("#Deltam/E");
  muondiff22D->GetYaxis()->CenterTitle();
  muondiff22D->Draw("COLZ");
  c2->SaveAs("muon2_massdiscrepancy_frac_2D.png");
  c2->SaveAs("muon2_massdiscrepancy_frac_2D.eps");
  c2->SaveAs("muon2_massdiscrepancy_frac_2D.pdf");
  c2->SaveAs("muon2_massdiscrepancy_frac_2D.root");
  c2->SaveAs("muon2_massdiscrepancy_frac_2D.C");
  muondiffinit12D->SetTitle("");
  muondiffinit12D->GetXaxis()->SetTitle("E_{#mu} [GeV]");
  muondiffinit12D->GetXaxis()->CenterTitle();
  muondiffinit12D->GetYaxis()->SetTitle("#Deltam/E");
  muondiffinit12D->GetYaxis()->CenterTitle();
  muondiffinit12D->Draw("COLZ");
  c2->SaveAs("muoninit1_massdiscrepancy_frac_2D.png");
  c2->SaveAs("muoninit1_massdiscrepancy_frac_2D.eps");
  c2->SaveAs("muoninit1_massdiscrepancy_frac_2D.pdf");
  c2->SaveAs("muoninit1_massdiscrepancy_frac_2D.root");
  c2->SaveAs("muoninit1_massdiscrepancy_frac_2D.C");
  muondiffinit22D->SetTitle("");
  muondiffinit22D->GetXaxis()->SetTitle("E_{#mu} [GeV]");
  muondiffinit22D->GetXaxis()->CenterTitle();
  muondiffinit22D->GetYaxis()->SetTitle("#Deltam/E");
  muondiffinit22D->GetYaxis()->CenterTitle();
  muondiffinit22D->Draw("COLZ");
  c2->SaveAs("muoninit2_massdiscrepancy_frac_2D.png");
  c2->SaveAs("muoninit2_massdiscrepancy_frac_2D.eps");
  c2->SaveAs("muoninit2_massdiscrepancy_frac_2D.pdf");
  c2->SaveAs("muoninit2_massdiscrepancy_frac_2D.root");
  c2->SaveAs("muoninit2_massdiscrepancy_frac_2D.C");*/

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
    //cout<<lines[k]<<endl;
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

  vector<int> electrons;
  vector<int> muons;
  vector<int> taus;
  vector<int> leptons;
  for(int k=0;k<npart;++k){
    if(abs(idup[k])==11) electrons.push_back(k);
    if(abs(idup[k])==13) muons.push_back(k);
    if(abs(idup[k])==15) taus.push_back(k);  
    if(abs(idup[k])==11 || abs(idup[k])==13 || abs(idup[k])==15) leptons.push_back(k);
  }

  Double_t m_el = 0.000510998;
  double m_mu = 0.105658;
  double a;

  bool problem=false;
  //float dm1,dm2,dmfix1,dmfix2;
  TString temp1,temp2;
  if(electrons.size()==2 && muons.size()==2 && taus.size()==0){
    TLorentzVector lep1,lep2,lep1pr,lep2pr;
    lep1.SetPxPyPzE(pup[electrons[0]][0],pup[electrons[0]][1],pup[electrons[0]][2],pup[electrons[0]][3]);
    lep2.SetPxPyPzE(pup[electrons[1]][0],pup[electrons[1]][1],pup[electrons[1]][2],pup[electrons[1]][3]);
    lep1=makemassless(lep1);
    lep2=makemassless(lep2);
    a=sqrt(1.-2*pow(m_el,2)/lep1.Dot(lep2));
    lep1pr=(1+a)/2.*lep1+(1-a)/2.*lep2;
    lep2pr=(1+a)/2.*lep2+(1-a)/2.*lep1;
    if(fabs(lep1pr.M()-m_el)/lep1pr.E()>0.01 || fabs(lep2pr.M()-m_el)/lep2pr.E()>0.01) problem=true;
    electrondiff1->Fill((lep1pr.M()-m_el)/lep1pr.E());
    electrondiff2->Fill((lep2pr.M()-m_el)/lep2pr.E());
    electrondiffinit1->Fill(lep1.M()/lep1.E());
    electrondiffinit2->Fill(lep2.M()/lep2.E());
    electrondiff12D->Fill(lep1pr.E(),(lep1pr.M()-m_el)/lep1pr.E());
    electrondiff22D->Fill(lep2pr.E(),(lep2pr.M()-m_el)/lep2pr.E());
    electrondiffinit12D->Fill(lep1.E(),lep1.M()/lep1.E());
    electrondiffinit22D->Fill(lep2.E(),lep2.M()/lep2.E());
    string fix1 = makefix(lep1pr,m_el);
    string fix2 = makefix(lep2pr,m_el);
    if(fix1.size()!=70 || fix2.size()!=70) cout<<"WARNING"<<endl;
    temp1=Tlines[electrons[0]+1];
    temp2=Tlines[electrons[1]+1];
    Tlines[electrons[0]+1].Replace(29,70,fix1.c_str());
    Tlines[electrons[1]+1].Replace(29,70,fix2.c_str());
    electrons.pop_back();
    electrons.pop_back();

    if(lep1.E()<1. || lep2.E()<1) cout<<"WARNING!"<<endl;

    lep1.SetPxPyPzE(pup[muons[0]][0],pup[muons[0]][1],pup[muons[0]][2],pup[muons[0]][3]);
    lep2.SetPxPyPzE(pup[muons[1]][0],pup[muons[1]][1],pup[muons[1]][2],pup[muons[1]][3]);
    lep1=makemassless(lep1);
    lep2=makemassless(lep2);
    a=sqrt(1.-2*pow(m_mu,2)/lep1.Dot(lep2));
    lep1pr=(1+a)/2.*lep1+(1-a)/2.*lep2;
    lep2pr=(1+a)/2.*lep2+(1-a)/2.*lep1;
    fix1 = makefix(lep1pr,m_mu);
    fix2 = makefix(lep2pr,m_mu);
    if(fix1.size()!=70 || fix2.size()!=70) cout<<"WARNING"<<endl;
    temp1=Tlines[muons[0]+1];
    temp2=Tlines[muons[1]+1];
    Tlines[muons[0]+1].Replace(29,70,fix1.c_str());
    Tlines[muons[1]+1].Replace(29,70,fix2.c_str()); 
    if(fabs(lep1pr.M()-m_mu)/lep1pr.E()>0.01){
      problem=true;
    }
    if(fabs(lep2pr.M()-m_mu)/lep2pr.E()>0.01){
      problem=true;
    }
    muondiff1->Fill((lep1pr.M()-m_mu)/lep1pr.E());
    muondiff2->Fill((lep2pr.M()-m_mu)/lep2pr.E());
    muondiffinit1->Fill(lep1.M()/lep1.E());
    muondiffinit2->Fill(lep2.M()/lep2.E());
    muondiff12D->Fill(lep1pr.E(),(lep1pr.M()-m_mu)/lep1pr.E());
    muondiff22D->Fill(lep2pr.E(),(lep2pr.M()-m_mu)/lep2pr.E());
    muondiffinit12D->Fill(lep1.E(),lep1.M()/lep1.E());
    muondiffinit22D->Fill(lep2.E(),lep2.M()/lep2.E());
    muons.pop_back();
    muons.pop_back();

    if(lep1.E()<1. || lep2.E()<1) cout<<"WARNING!"<<endl;

    if(problem) numdiff++;

  } else if(electrons.size()==4 || muons.size()==4){
    double m_l;
    if(electrons.size()==4) m_l=m_el;
    if(muons.size()==4) m_l=m_mu;
    TLorentzVector l1_minus, l1_plus, l2_minus, l2_plus;
    int l1p,l1m,l2p,l2m;

    if(idup[leptons[0]]>0){
      l1_minus.SetPxPyPzE(pup[leptons[0]][0], pup[leptons[0]][1], pup[leptons[0]][2], pup[leptons[0]][3]);
      if(idup[leptons[1]]>0){
        l2_minus.SetPxPyPzE(pup[leptons[1]][0], pup[leptons[1]][1], pup[leptons[1]][2], pup[leptons[1]][3]);
        l1_plus.SetPxPyPzE(pup[leptons[2]][0], pup[leptons[2]][1], pup[leptons[2]][2], pup[leptons[2]][3]);
        l2_plus.SetPxPyPzE(pup[leptons[3]][0], pup[leptons[3]][1], pup[leptons[3]][2], pup[leptons[3]][3]);
        l1p=2; l1m=0; l2p=3; l2m=1;
      }
      else if(idup[leptons[2]]>0){
        l1_plus.SetPxPyPzE(pup[leptons[1]][0], pup[leptons[1]][1], pup[leptons[1]][2], pup[leptons[1]][3]);
        l2_minus.SetPxPyPzE(pup[leptons[2]][0], pup[leptons[2]][1], pup[leptons[2]][2], pup[leptons[2]][3]);
        l2_plus.SetPxPyPzE(pup[leptons[3]][0], pup[leptons[3]][1], pup[leptons[3]][2], pup[leptons[3]][3]);
        l1p=1; l1m=0; l2p=3; l2m=2;
      }
      else{
        l1_plus.SetPxPyPzE(pup[leptons[1]][0], pup[leptons[1]][1], pup[leptons[1]][2], pup[leptons[1]][3]);
        l2_plus.SetPxPyPzE(pup[leptons[2]][0], pup[leptons[2]][1], pup[leptons[2]][2], pup[leptons[2]][3]);
        l2_minus.SetPxPyPzE(pup[leptons[3]][0], pup[leptons[3]][1], pup[leptons[3]][2], pup[leptons[3]][3]);
        l1p=1; l1m=0; l2p=2; l2m=3;
      }
    }
    else{
      l1_plus.SetPxPyPzE(pup[leptons[0]][0], pup[leptons[0]][1], pup[leptons[0]][2], pup[leptons[0]][3]);
      if(idup[leptons[1]]<0){
        l2_plus.SetPxPyPzE(pup[leptons[1]][0], pup[leptons[1]][1], pup[leptons[1]][2], pup[leptons[1]][3]);
        l1_minus.SetPxPyPzE(pup[leptons[2]][0], pup[leptons[2]][1], pup[leptons[2]][2], pup[leptons[2]][3]);
        l2_minus.SetPxPyPzE(pup[leptons[3]][0], pup[leptons[3]][1], pup[leptons[3]][2], pup[leptons[3]][3]);
        l1p=0; l1m=2; l2p=1; l2m=3;
      }
      else if(idup[leptons[2]]<0){
        l1_minus.SetPxPyPzE(pup[leptons[1]][0], pup[leptons[1]][1], pup[leptons[1]][2], pup[leptons[1]][3]);
        l2_plus.SetPxPyPzE(pup[leptons[2]][0], pup[leptons[2]][1], pup[leptons[2]][2], pup[leptons[2]][3]);
        l2_minus.SetPxPyPzE(pup[leptons[3]][0], pup[leptons[3]][1], pup[leptons[3]][2], pup[leptons[3]][3]);
        l1p=0; l1m=1; l2p=2; l2m=3;
      }
      else{
        l1_minus.SetPxPyPzE(pup[leptons[1]][0], pup[leptons[1]][1], pup[leptons[1]][2], pup[leptons[1]][3]);
        l2_minus.SetPxPyPzE(pup[leptons[2]][0], pup[leptons[2]][1], pup[leptons[2]][2], pup[leptons[2]][3]);
        l2_plus.SetPxPyPzE(pup[leptons[3]][0], pup[leptons[3]][1], pup[leptons[3]][2], pup[leptons[3]][3]);           
        l1p=0; l1m=1; l2p=3; l2m=2;
      }
    }

    TLorentzVector Z1 = l1_minus + l1_plus;
    TLorentzVector Z2 = l2_minus + l2_plus;
    TLorentzVector Z1alt = l1_minus + l2_plus;
    TLorentzVector Z2alt = l2_minus + l1_plus;

    if(fabs(Z1.M()-Zmass)>fabs(Z2.M()-Zmass)) swap(Z1,Z2);
    if(fabs(Z1alt.M()-Zmass)>fabs(Z2alt.M()-Zmass)) swap(Z1alt,Z2alt);
    if(fabs(Z1alt.M()-Zmass)<fabs(Z1.M()-Zmass)){
      swap(l1m,l2m);
      swap(l1_minus,l2_minus);
    }

    TLorentzVector lep1,lep2,lep1pr,lep2pr;
    lep1.SetPxPyPzE(pup[leptons[l1p]][0],pup[leptons[l1p]][1],pup[leptons[l1p]][2],pup[leptons[l1p]][3]);
    lep2.SetPxPyPzE(pup[leptons[l1m]][0],pup[leptons[l1m]][1],pup[leptons[l1m]][2],pup[leptons[l1m]][3]);
    lep1=makemassless(lep1);
    lep2=makemassless(lep2);
    a=sqrt(1.-2*pow(m_l,2)/lep1.Dot(lep2));
    lep1pr=(1+a)/2.*lep1+(1-a)/2.*lep2;
    lep2pr=(1+a)/2.*lep2+(1-a)/2.*lep1;
    string fix1 = makefix(lep1pr,m_l);
    string fix2 = makefix(lep2pr,m_l);
    if(fix1.size()!=70 || fix2.size()!=70) cout<<"WARNING"<<endl;
    temp1=Tlines[leptons[l1p]+1];
    temp2=Tlines[leptons[l1m]+1];
    Tlines[leptons[l1p]+1].Replace(29,70,fix1.c_str());
    Tlines[leptons[l1m]+1].Replace(29,70,fix2.c_str());
    if(fabs(lep1pr.M()-m_l)/lep1pr.E()>0.01 || fabs(lep2pr.M()-m_l)/lep2pr.E()>0.01) problem=true;

    if(lep1.E()<1. || lep2.E()<1) cout<<"WARNING!"<<endl;

    lep1.SetPxPyPzE(pup[leptons[l2p]][0],pup[leptons[l2p]][1],pup[leptons[l2p]][2],pup[leptons[l2p]][3]);
    lep2.SetPxPyPzE(pup[leptons[l2m]][0],pup[leptons[l2m]][1],pup[leptons[l2m]][2],pup[leptons[l2m]][3]);
    lep1=makemassless(lep1);
    lep2=makemassless(lep2);
    a=sqrt(1.-2*pow(m_l,2)/lep1.Dot(lep2));
    lep1pr=(1+a)/2.*lep1+(1-a)/2.*lep2;
    lep2pr=(1+a)/2.*lep2+(1-a)/2.*lep1;
    fix1 = makefix(lep1pr,m_l);
    fix2 = makefix(lep2pr,m_l);
    if(fix1.size()!=70 || fix2.size()!=70) cout<<"WARNING"<<endl;
    temp1=Tlines[leptons[l2p]+1];
    temp2=Tlines[leptons[l2m]+1];
    Tlines[leptons[l2p]+1].Replace(29,70,fix1.c_str());
    Tlines[leptons[l2m]+1].Replace(29,70,fix2.c_str());
    if(fabs(lep1pr.M()-m_l)/lep1pr.E()>0.01 || fabs(lep2pr.M()-m_l)/lep2pr.E()>0.01) problem=true;

    if(lep1.E()<1. || lep2.E()<1) cout<<"WARNING!"<<endl;

    if(problem) numdiff++;

    if(electrons.size()==4){
      electrons.pop_back();
      electrons.pop_back();
      electrons.pop_back();
      electrons.pop_back();
    }
    if(muons.size()==4){
      muons.pop_back();
      muons.pop_back();
      muons.pop_back();
      muons.pop_back();
    }
    leptons.pop_back();
    leptons.pop_back();
    leptons.pop_back();
    leptons.pop_back();

  } else{
    cout<<"Other lepton settings not yet functional."<<endl;
    return Tlines;
  }
  
  return Tlines;
}

string makefix(TLorentzVector lep, Double_t mass){
  std::stringstream fix;

  fix<<fixformat(lep.Px());
  fix<<fixformat(lep.Py());
  fix<<fixformat(lep.Pz());
  fix<<fixformat(lep.E());
  fix<<fixformat(mass);

  return fix.str();
}

string fixformat(Double_t val){
  std::stringstream fixstream;

  char buffer[10]="";
  char expstr[1]="";
  int exp = int(floor(log10(fabs(val))));
  float dec = fabs(val)/pow(10,exp+1);
  sprintf(buffer,"%.6f",dec);
  if(exp>-1) sprintf(expstr,"%i",abs(exp)+1);
  if(exp<=-1) sprintf(expstr,"%i",abs(exp)-1);
  if(val>=0.) fixstream<<"  ";
  if(val<0.) fixstream<<" -";
  fixstream<<buffer;
  if(exp>=-1) fixstream<<"E+0";
  if(exp<-1) fixstream<<"E-0";
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

TLorentzVector makemassless(TLorentzVector base){
  TLorentzVector output;
  double efactor,pfactor;
  efactor=sqrt(1.-(pow(base.E(),2)-pow(base.P(),2))/(2*pow(base.E(),2)));
  pfactor=sqrt(1.+(pow(base.E(),2)-pow(base.P(),2))/(2*pow(base.P(),2)));
  output.SetPxPyPzE(pfactor*base.Px(),pfactor*base.Py(),pfactor*base.Pz(),efactor*base.E());  
  return output;
}
