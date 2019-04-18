void get_Kstar2K_pp502(TGraphErrors* g,TBox* b){
  double x=4.30,xl=0.11,xu=0.32;//https://aliceinfo.cern.ch/Notes/node/603, version: 13 January 2017
  xl=xl/x/3*pow(x,1./3.);
  xu=xu/x/3*pow(x,1./3.);
  x=pow(x,1./3.);

  double v=0.3291,t=0.0024,y=0.0314;//Ranbir, 25 September 2017

  g->SetPoint(0,x,v);
  g->SetPointError(0,0.,t);
  b->SetX1(x-xl);
  b->SetY1(v-y);
  b->SetX2(x+xu);
  b->SetY2(v+y);

  return;
}

void get_phi2K_pp502(TGraphErrors* g,TBox* b){
  double x=4.30,xl=0.11,xu=0.32;//https://aliceinfo.cern.ch/Notes/node/603, version: 13 January 2017
  xl=xl/x/3*pow(x,1./3.);
  xu=xu/x/3*pow(x,1./3.);
  x=pow(x,1./3.);

  //0.11290 +/- 0.000803 (stat.) +/- 0.010415 (sys.); Sushanta, 4 July 2017
  double v=0.11290,t=0.000803,y=0.010415;

  g->SetPoint(0,x,v);
  g->SetPointError(0,0.,t);
  b->SetX1(x-xl);
  b->SetY1(v-y);
  b->SetX2(x+xu);
  b->SetY2(v+y);

  return;
}


void get_KStar2K_pp13Mult(TGraphErrors* gt,TBox** by,TGraphErrors* gu){
  double X[9]={26.18, 20.16, 16.40, (14.00+12.28)/2, 10.31, 8.24, 6.62, 4.77, 2.76};
  double UX[9]={0.55, 0.41, 0.34, (0.29+0.25)/2, 0.21, 0.17, 0.13, 0.09, 0.05};
  double V[9]={0.29711, 0.292132, 0.300112, 0.300922, 0.295984, 0.308614, 0.309597, 0.314755, 0.356309};
  double T[9]={0.0152817, 0.00775211, 0.00846854, 0.00650299, 0.0069524, 0.00775581, 0.00794056, 0.00650175, 0.00669031};
  double Y[9]={0.0302605, 0.0289287, 0.0300567, 0.0301463, 0.0293947, 0.0308925, 0.0315497, 0.0327565, 0.0385581};
  double U[9]={0.0163327, 0.01566, 0.0160537, 0.0134308, 0.0165804, 0.0142741, 0.0156264, 0.0157139, 0.0183867};

  for(int j=0;j<9;j++){
    UX[j]=UX[j]/X[j]/3.;
    X[j]=pow(X[j],1./3.);
    UX[j]*=X[j];

    gt->SetPoint(j,X[j],V[j]);
    gt->SetPointError(j,0.,T[j]);
    by[j]=new TBox(X[j]-UX[j],V[j]-Y[j],X[j]+UX[j],V[j]+Y[j]);
    gu->SetPoint(j,X[j],V[j]);
    gu->SetPointError(j,UX[j],U[j]);
  }

  return;
}


void get_phi2K_pp13Mult(TGraphErrors* gt,TBox** by,TGraphErrors* gu){
  double X[10]={26.18, 20.16, 16.40, 14.00, 12.28, 10.31, 8.24, 6.62, 4.77, 2.76};
  double UX[10]={0.55, 0.41, 0.34, 0.29, 0.25, 0.21, 0.17, 0.13, 0.09, 0.05};
  double V[10]={1.298808214e-01, 1.239350180e-01, 1.245949487e-01, 1.183999937e-01, 1.191009856e-01, 1.178998652e-01, 1.206026122e-01, 1.100308117e-01, 1.117142718e-01, 1.112157540e-01};
  double T[10]={4.884498452e-03, 2.245981927e-03, 2.232748912e-03, 2.202173612e-03, 2.555975587e-03, 1.964021921e-03, 2.741911936e-03, 2.776588136e-03, 2.463968924e-03, 3.839918927e-03};
  double Y[10]={1.055515048e-02, 1.001183125e-02, 1.008259878e-02, 9.907249429e-03, 9.673362131e-03, 9.668879000e-03, 1.006892976e-02, 9.540144244e-03, 1.037493420e-02, 1.177919284e-02};
  double U[10]={3.458048837e-03, 3.116890345e-03, 3.194780784e-03, 3.839738880e-03, 3.163422355e-03, 3.408822820e-03, 3.175601271e-03, 2.925871538e-03, 3.911558479e-03, 5.737222706e-03};

  for(int j=0;j<10;j++){
    UX[j]=UX[j]/X[j]/3.;
    X[j]=pow(X[j],1./3.);
    UX[j]*=X[j];

    gt->SetPoint(j,X[j],V[j]);
    gt->SetPointError(j,0.,T[j]);
    by[j]=new TBox(X[j]-UX[j],V[j]-Y[j],X[j]+UX[j],V[j]+Y[j]);
    gu->SetPoint(j,X[j],V[j]);
    gu->SetPointError(j,UX[j],U[j]);
  }

  return;
}


void get_Kstar2K_pPb502(TGraphErrors* gt,TBox** by,TGraphErrors* gu){
  // dN/deta
  // double x[5]={ 35.6, 23.2, 16.1, 9.8, 4.4};
  // double ux[5]={ 0.8,  0.5,  0.4, 0.2, 0.1};
  // (dNch/deta)^1/3
  double x[5]={TMath::Power(35.6,1./3.),TMath::Power(23.2,1./3.),TMath::Power(16.1,1./3.),TMath::Power(9.8,1./3.),TMath::Power(4.4,1./3.)};
  double udndeta[5]={0.8, 0.5, 0.4, 0.2, 0.1};

  int j;
  
  double ux[5];
  for(j=0;j<5;j++) ux[j]=TMath::Power(x[j],-2.)*udndeta[j]/3.;

  Double_t y[5]  = {0.270104, 0.288596, 0.297597, 0.304737, 0.324889};
  Double_t ty[5] = {0.003762, 0.003818, 0.003986, 0.004150, 0.005191};
  Double_t yy[5] = {0.026760, 0.027251, 0.027792, 0.027855, 0.027892};
  Double_t uncorr[5] = {0.025504, 0.025839, 0.026319, 0.026312, 0.026134};

  double corr[5];
  double uncorr[5];
  for(j=0;j<5;j++){ 
    uncorr[j] = y[j] * TMath::Sqrt(TMath::Power((yy[j]/y[j]),2.0) - 0.0391*0.0391);
    corr[j] = 0.0391 * y[j];
  }

  for(j=0;j<5;j++){
    gt->SetPoint(j,x[j],y[j]);
    gt->SetPointError(j,ux[j],ty[j]);
    by[j]=new TBox(x[j]-ux[j],y[j]-yy[j],x[j]+ux[j],y[j]+yy[j]);
    gu->SetPoint(j,x[j],y[j]);
    gu->SetPointError(j,ux[j],uncorr[j]);
  }
  return;
}


void get_phi2K_pPb502(TGraphErrors* gt,TBox** by,TGraphErrors* gu){
  double x[7]={ TMath::Power(45.15,1./3.), TMath::Power(36.22,1./3.), TMath::Power(30.46,1./3.),TMath::Power(23.24,1./3.), TMath::Power(16.08,1./3.), TMath::Power(9.82,1./3.), TMath::Power(4.4,1./3.)};
  double udndeta[7]={1., 0.8, 0.67, 0.51, 0.35, 0.21, 0.1};
  int j;

  double ux[7];
  for(j=0;j<7;j++) ux[j]=TMath::Power(x[j],-2.)*udndeta[j]/3.;
                                                                                                           
  Double_t y[7] = {0.1290, 0.1241, 0.1254, 0.1250, 0.1213, 0.1143, 0.1160};//value
  Double_t ty[7] = {0.0013, 0.0013, 0.0010, 0.0008, 0.0009, 0.0012, 0.0018};//stat. err
  Double_t yy[7] = {0.0126, 0.0112, 0.0110, 0.0107, 0.0102, 0.0089, 0.0110};//syst. err
  Double_t uncorr[7] = {0.0076, 0.0057, 0.0053, 0.0053, 0.0052, 0.0046, 0.0078};//uncorr

  for(int j=0;j<7;j++){
    gt->SetPoint(j,x[j],y[j]);
    gt->SetPointError(j,ux[j],ty[j]);
    by[j]=new TBox(x[j]-ux[j],y[j]-yy[j],x[j]+ux[j],y[j]+yy[j]);
    gu->SetPoint(j,x[j],y[j]);
    gu->SetPointError(j,ux[j],uncorr[j]);
  }

  return;
}


void get_Kstar2K_PbPb502(TGraphErrors* gt,TBox** by,TGraphErrors* gu){
  double X[8],UX[8],V[8],T[8],Y[8],U[8];
  int j=0;
  X[j]=(2035.43+1850.15+1666.32+1505.11)/4; UX[j]=(52.43+55.23+47.89+43.71)/4; j++;
  X[j]=1180.48; UX[j]=31.18; j++;
  X[j]=786.45; UX[j]=20.44; j++;
  X[j]=511.66; UX[j]=15.2; j++;
  X[j]=318.15; UX[j]=12.27; j++;
  X[j]=183.33; UX[j]=8.07; j++;
  X[j]=96.27; UX[j]=5.81; j++;
  X[j]=44.89; UX[j]=3.44; j++;

  j=0;
  V[j]=0.159; T[j]=0.0111; Y[j]=0.0208; U[j]=0.; j++;
  V[j]=0.193; T[j]=0.0109; Y[j]=0.0243; U[j]=0.; j++;
  V[j]=0.210; T[j]=0.0094; Y[j]=0.0250; U[j]=0.; j++;
  V[j]=0.217; T[j]=0.0094; Y[j]=0.0299; U[j]=0.; j++;
  V[j]=0.225; T[j]=0.0111; Y[j]=0.0256; U[j]=0.; j++;
  V[j]=0.226; T[j]=0.0105; Y[j]=0.0239; U[j]=0.; j++;
  V[j]=0.247; T[j]=0.0119; Y[j]=0.0352; U[j]=0.; j++;
  V[j]=0.260; T[j]=0.0118; Y[j]=0.0278; U[j]=0.; j++;

  for(int j=0;j<8;j++){
    UX[j]=UX[j]/X[j]/3.;
    X[j]=pow(X[j],1./3.);
    UX[j]*=X[j];

    gt->SetPoint(j,X[j],V[j]);
    gt->SetPointError(j,UX[j],T[j]);
    by[j]=new TBox(X[j]-UX[j],V[j]-Y[j],X[j]+UX[j],V[j]+Y[j]);
    gu->SetPoint(j,X[j],V[j]);
    gu->SetPointError(j,UX[j],U[j]);
  }

  return;
}


void get_phi2K_PbPb502(TGraphErrors* gt,TBox** by,TGraphErrors* gu){
  double X[8],UX[8],V[8],T[8],Y[8],U[8];
  int j=0;
  X[j]=(2035.43+1850.15+1666.32+1505.11)/4; UX[j]=(52.43+55.23+47.89+43.71)/4; j++;
  X[j]=1180.48; UX[j]=31.18; j++;
  X[j]=786.45; UX[j]=20.44; j++;
  X[j]=511.66; UX[j]=15.2; j++;
  X[j]=318.15; UX[j]=12.27; j++;
  X[j]=183.33; UX[j]=8.07; j++;
  X[j]=96.27; UX[j]=5.81; j++;
  X[j]=44.89; UX[j]=3.44; j++;

  j=0;
  V[j]=0.120; T[j]=0.0012; Y[j]=0.0110; U[j]=0.; j++;
  V[j]=0.121; T[j]=0.0013; Y[j]=0.0106; U[j]=0.; j++;
  V[j]=0.128; T[j]=0.0014; Y[j]=0.0115; U[j]=0.; j++;
  V[j]=0.134; T[j]=0.0014; Y[j]=0.0116; U[j]=0.; j++;
  V[j]=0.130; T[j]=0.0015; Y[j]=0.0117; U[j]=0.; j++;
  V[j]=0.134; T[j]=0.0017; Y[j]=0.0118; U[j]=0.; j++;
  V[j]=0.128; T[j]=0.0018; Y[j]=0.0126; U[j]=0.; j++;
  V[j]=0.121; T[j]=0.0023; Y[j]=0.0115; U[j]=0.; j++;

  for(int j=0;j<8;j++){
    UX[j]=UX[j]/X[j]/3.;
    X[j]=pow(X[j],1./3.);
    UX[j]*=X[j];

    gt->SetPoint(j,X[j],V[j]);
    gt->SetPointError(j,UX[j],T[j]);
    by[j]=new TBox(X[j]-UX[j],V[j]-Y[j],X[j]+UX[j],V[j]+Y[j]);
    gu->SetPoint(j,X[j],V[j]);
    gu->SetPointError(j,UX[j],U[j]);
  }

  return;
}

/*************************************************************/
/*************************************************************/

void Kstar2K_phi2K_Run2(){
  myOptions();
  gROOT->ForceStyle();

  //array of colors for pp vs mult,  pPb,    PbPb
  Color_t color[20] = {kMagenta+2, kRed, kOrange-3, //K*/K value + sys
		       kMagenta-9, kRed-9, kOrange-9, //K*/K uncorr
		       kGreen+2, kBlue+2, kTeal-6,   //phi/K value+sys 
		       kSpring+2, kAzure-9, kCyan-8,
		       kBlack, kMagenta-8, kGreen, kGreen-9, 
		       kRed+1,  kBlue, kAzure-9, kMagenta}; //pp min bias INEL is black

  int j;

  TCanvas* c=new TCanvas("c","",10,10,800,800);
  c->SetFillColor(0);
  c->SetTickx(1);
  c->SetTicky(1);
  c->Draw();
  c->cd();

  TH1F* hf=new TH1F("frame","",1, 0.,13.);
  hf->SetLineColor(1);

  hf->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT^{1/3}_{|#it{#eta_{lab}}| < 0.5}");
  hf->GetXaxis()->SetTitleSize(0.06);
  hf->GetXaxis()->SetLabelSize(0.045);
  hf->GetXaxis()->SetTitleOffset(1.25);
  hf->SetNdivisions(510,"x");

  hf->SetYTitle("Ratio to K");
  hf->SetMinimum(0.0);
  hf->SetMaximum(0.5);
  hf->GetYaxis()->SetTitleSize(0.06);
  hf->GetYaxis()->SetLabelSize(0.045);
  hf->GetYaxis()->SetTitleOffset(1.0);
  hf->SetNdivisions(508,"y");

  double ml=0.13,mt=0.02,mr=0.01,mb=0.18;

  TPad* p2=new TPad("p2","",0,0,1,1);
  myPadSetUp(p2,ml,mt,mr,mb);
  p2->SetFillColor(0);
  p2->SetTickx(1);
  p2->SetTicky(1);

  //-----

  TGraphErrors* pkt=new TGraphErrors(1);
  TBox* pky=new TBox();
  get_Kstar2K_pp502(pkt,pky);

  pkt->SetLineColor(TColor::GetColor("#aa00ff"));//color[12]); 
  pkt->SetMarkerColor(TColor::GetColor("#aa00ff"));//color[12]);
  pkt->SetMarkerStyle(21); 
  pkt->SetMarkerSize(1.8);

  pky->SetLineColor(TColor::GetColor("#aa00ff"));//color[12]); 
  pky->SetFillStyle(0);

  //-----

  TGraphErrors* ppt=new TGraphErrors(1);
  TBox* ppy=new TBox();
  get_phi2K_pp502(ppt,ppy);

  ppt->SetLineColor(TColor::GetColor("#aa00ff")); 
  ppt->SetMarkerColor(TColor::GetColor("#aa00ff"));
  ppt->SetMarkerStyle(20);
  ppt->SetMarkerSize(1.8);

  ppy->SetLineColor(TColor::GetColor("#aa00ff")); 
  ppy->SetFillStyle(0);

  //-----
  
  TGraphErrors* xkt=new TGraphErrors(9);
  TBox* xky[9];
  TGraphErrors *xku=new TGraphErrors(9);
  get_KStar2K_pp13Mult(xkt,xky,xku);
  
  xkt->SetLineColor(TColor::GetColor("#ffaa00"));
  xkt->SetMarkerColor(TColor::GetColor("#ffaa00")); 
  xkt->SetMarkerStyle(21);
  xkt->SetMarkerSize(1.8);

  for(j=0;j<xkt->GetN();j++){
    xky[j]->SetFillStyle(0);
    xky[j]->SetLineColor(TColor::GetColor("#ffaa00"));
  }
  
  xku->SetMarkerStyle(0);
  ////####xku->SetMarkerColor(kGreen+2);
  xku->SetLineColor(0);
  xku->SetFillColor(kGreen-8);
  xku->SetFillStyle(1001);

  //-----
  
  TGraphErrors* xpt=new TGraphErrors(10);
  TBox* xpy[10];
  TGraphErrors* xpu=new TGraphErrors(10);
  get_phi2K_pp13Mult(xpt,xpy,xpu);

  xpt->SetLineColor(TColor::GetColor("#ffaa00"));
  xpt->SetMarkerColor(TColor::GetColor("#ffaa00")); 
  xpt->SetMarkerStyle(20);
  xpt->SetMarkerSize(1.8);

  for(j=0;j<xpt->GetN();j++){
    xpy[j]->SetFillStyle(0);
    xpy[j]->SetLineColor(TColor::GetColor("#ffaa00"));
  }
  
  xpu->SetMarkerStyle(0);
  xpu->SetLineColor(0);
  xpu->SetFillColor(kGreen-8);
  xpu->SetFillStyle(1001);

  //-----

  TGraphErrors* akt=new TGraphErrors(5);
  TBox* aky[5];
  TGraphErrors* aku=new TGraphErrors(5);
  get_Kstar2K_pPb502(akt,aky,aku);

  akt->SetLineColor(color[17]); 
  akt->SetMarkerColor(color[17]); 
  akt->SetMarkerStyle(25); 
  akt->SetMarkerSize(1.8);

  for(j=0;j<akt->GetN();j++){
    aky[j]->SetFillStyle(0);
    aky[j]->SetFillColor(color[0]);
    aky[j]->SetLineColor(color[17]);
  }

  aku->SetLineColor(0); 
  aku->SetMarkerColor(0); 
  aku->SetMarkerStyle(0);
  aku->SetFillColor(color[18]);
  aku->SetFillStyle(1001);

  //-----
  
  TGraphErrors* apt=new TGraphErrors(7);
  TBox* apy[7];
  TGraphErrors* apu=new TGraphErrors(7);
  get_phi2K_pPb502(apt,apy,apu);

  apt->SetLineColor(color[17]);
  apt->SetMarkerColor(color[17]); 
  apt->SetMarkerStyle(4);
  apt->SetMarkerSize(1.8);

  for(j=0;j<apt->GetN();j++){
    apy[j]->SetFillStyle(0);
    apy[j]->SetFillColor(color[17]);
    apy[j]->SetLineColor(color[17]);
  }

  apu->SetLineColor(0); 
  apu->SetMarkerColor(0); 
  apu->SetMarkerStyle(0);
  apu->SetFillColor(color[18]);
  apu->SetFillStyle(1001);

  //-----

  TGraphErrors* lkt=new TGraphErrors(8);
  TBox* lky[8];
  TGraphErrors* lku=new TGraphErrors(8);
  get_Kstar2K_PbPb502(lkt,lky,lku);

  lkt->SetLineColor(1);
  lkt->SetMarkerColor(1); 
  lkt->SetMarkerStyle(21); 
  lkt->SetMarkerSize(1.8);

  lku->SetLineColor(color[1]); 
  lku->SetFillColor(color[4]);

  for(j=0;j<lkt->GetN();j++){
    lky[j]->SetLineColor(1); 
    lky[j]->SetFillColor(0); 
    lky[j]->SetFillStyle(0);
  }

  //-----

  TGraphErrors* lpt=new TGraphErrors(8);
  TBox* lpy[8];
  TGraphErrors* lpu=new TGraphErrors(8);
  get_phi2K_PbPb502(lpt,lpy,lpu);

  lpt->SetLineColor(1); 
  lpt->SetMarkerColor(1); 
  lpt->SetMarkerStyle(20); 
  lpt->SetMarkerSize(1.8);

  lpu->SetLineColor(color[1]); 
  lpu->SetFillColor(color[4]);

  for(j=0;j<lpt->GetN();j++){
    lpy[j]->SetLineColor(1); 
    lpy[j]->SetFillStyle(0);
  }

  //-----

  c->cd();
  SetStyle();
  p2->Draw();
  p2->cd();
  hf->Draw();

  //lku->Draw("e2same");
  for(j=0;j<lkt->GetN();j++) if(lky[j]) lky[j]->Draw();
  lkt->Draw("pzsame");

  //lpu->Draw("e2same");
  for(j=0;j<lpt->GetN();j++) if(lpy[j]) lpy[j]->Draw();
  lpt->Draw("pzsame");

  //aku->Draw("e2same");
  for(j=0;j<akt->GetN();j++) if(aky[j]) aky[j]->Draw();
  akt->Draw("pzsame");

  //apu->Draw("e2same");
  for(j=0;j<apt->GetN();j++) if(apy[j]) apy[j]->Draw();
  apt->Draw("pzsame");

  //xku->Draw("e2same");
  for(j=0;j<xkt->GetN();j++) if(xky[j]) xky[j]->Draw();
  xkt->Draw("pzsame");

  //xpu->Draw("e2same");
  for(j=0;j<xpt->GetN();j++) if(xpy[j]) xpy[j]->Draw();
  xpt->Draw("pzsame");

  pky->Draw();
  pkt->Draw("pzsame");

  ppy->Draw();
  ppt->Draw("pzsame");

  /*
  TLine* smk=new TLine(9.7, 29.56/100.5, 11.9, 29.56/100.5);
  smk->SetLineColor(kGray+2);
  smk->SetLineWidth(4);  
  smk->SetLineStyle(1);

  TLine* smp=new TLine(9.7,12.16/100.5,11.9,12.16/100.5);
  smp->SetLineColor(kGray+2);
  smp->SetLineWidth(4);
  smp->SetLineStyle(7);
  smp->Draw();
  smk->Draw();

  TLegend *lth=new TLegend(0.3, 0.25, 0.8, 0.30);//, | y | < 0.5");
  myLegendSetUp(lth,0.04);
  lth->SetNColumns(2);
  lth->SetColumnSeparation(0.02);
  lth->AddEntry(smk," ","l");//thermal
  lth->AddEntry(smp," GC Thermal model (T=156 MeV)","l"); //thermal
  lth->Draw();
  */
  
  //Legends
  TMarker *mks[4]; 
  TMarker *mphi[4]; 
  TMarker *mphi_kstar1[4]; 
  TMarker *mphi_kstar2[4];
  TMarker *mks_pPb[4]; //kn
  TMarker *mphi_pPb[4];
  
  mks[0]  = new TMarker(11.0, 0.44, 21);   
  mphi[0] = new TMarker(12.0, 0.44, 20);
  mks[2]  = new TMarker(11.0, 0.41, 21);
  mphi[2] = new TMarker(12.0, 0.41, 20);
  mks[1]  = new TMarker(11.0, 0.38, 25);
  mphi[1] = new TMarker(12.0, 0.38, 4);
  mks[3]  = new TMarker(11.0, 0.35, 21);
  mphi[3] = new TMarker(12.0, 0.35, 20);
  
  for (int j=0; j<4; j++) {
    mks[j]->SetMarkerSize(1.8);
    mphi[j]->SetMarkerSize(1.8);
  }

  mks[0]->SetMarkerColor(TColor::GetColor("#aa00ff"));//color[12]);//pp INEL
  mphi[0]->SetMarkerColor(TColor::GetColor("#aa00ff"));//pp INEL
  mks[1]->SetMarkerColor(color[17]);//p-Pb
  mphi[1]->SetMarkerColor(color[17]);//p-Pb
  mks[2]->SetMarkerColor(TColor::GetColor("#ffaa00"));//pp 13 TeV mult
  mphi[2]->SetMarkerColor(TColor::GetColor("#ffaa00"));//pp 13 TeV mult
  mks[3]->SetMarkerColor(1); //Pb-Pb
  mphi[3]->SetMarkerColor(1); //Pb-Pb
  
  TPaveText * alice = new TPaveText(0.83, 0.91, 0.3, 0.93, "NDC");
  alice->SetFillColor(kWhite);
  alice->SetFillStyle(0);
  alice->SetBorderSize(0);
  alice->SetTextAlign(31);
  alice->SetTextSize(0.04);
  alice->SetLineWidth(2);
  alice->InsertText("#bf{ALICE Preliminary}");

  TPaveText * particles = new TPaveText(0.81, 0.90, 0.98, 0.94, "NDC");
  particles->SetFillColor(kWhite);
  particles->SetFillStyle(0);
  particles->SetBorderSize(0);
  particles->SetTextAlign(11);
  particles->SetTextSize(0.04);
  particles->InsertText("K*^{0}/K #phi/K");

  TPaveText * systems = new TPaveText(0.5, 0.715, 0.83, 0.907, "NDC");
  systems->SetFillColor(kWhite);
  systems->SetFillStyle(0);
  systems->SetBorderSize(0);
  systems->SetTextAlign(32);
  systems->SetTextSize(0.04);

  systems->InsertText("pp 5.02 TeV INEL");
  systems->InsertText("pp 13 TeV");
  systems->InsertText("p-Pb 5.02 TeV [1]");
  systems->InsertText("Pb-Pb 5.02 TeV");

  TPaveText* paverr = new TPaveText(0.32, 0.2, 0.8, 0.25,"NDC");
  paverr->InsertText("Uncertainties: stat. (bars), syst. (box)");
  paverr->SetBorderSize(0);
  paverr->SetFillColor(kWhite);
  paverr->SetTextFont(42);
  paverr->SetTextAlign(22);
  paverr->SetTextSize(0.035);
  paverr->SetFillStyle(0);  

  for (int j=0; j<4; j++) {
    mks[j]->Draw();
    mphi[j]->Draw();
  }
  
  particles->Draw();
  alice->Draw();
  systems->Draw();
  paverr->Draw();
  
  Double_t lx = 0.6;
  Double_t ly = 0.8;

  ///epos
  /*
  TLegend* ll4=new TLegend(lx-0.21,ly-0.51,lx+0.04,ly-0.43);
  myLegendSetUp(ll4,0.04);
  ll4->SetFillColor(0);
  ll4->SetBorderSize(0);
  ll4->AddEntry(lnt_epos,"EPOS3","l");
  ll4->Draw();
  */

  TLatex* t1=new TLatex(0.52,0.68,"[1]: #it{EPJC} #bf{76} 245 (2016)");
  t1->SetNDC();
  t1->SetTextSize(0.03);
  t1->Draw();

  c->SaveAs("Kstar2K_phi2K_Run2.pdf");
  c->SaveAs("Kstar2K_phi2K_Run2.eps");

  return;
}

void SetStyle(){
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X"); 
  //gStyle->SetTickLength(0.04,"Y"); 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
		  int currentMarkerStyle=20, int currentMarkerColor=0,
		  int currentLineStyle=1, int currentLineColor=0){
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  return;
}

void myOptions(Int_t lStat=0){
  // Set gStyle
  int font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.0,"xyz");  
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }

}
