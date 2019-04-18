void phi2K_energy(){
  myOptions();
  gROOT->ForceStyle();

  TCanvas* c=new TCanvas("c","",10,10,450,350);
  c->SetFillColor(0);
  c->Draw();
  c->cd();

  double ml=0.13,mr=0.01,mb=0.16,mt=0.01;
  int j;

  TH1F* hf=new TH1F("frame","",1,6.,20000.);
  hf->SetLineColor(1);

  hf->SetXTitle("#sqrt{#it{s}_{NN}} (GeV)");
  hf->GetXaxis()->SetTitleSize(0.07);
  hf->GetXaxis()->SetLabelSize(0.05);
  hf->GetXaxis()->SetTitleOffset(1.);
  hf->SetNdivisions(505,"x");

  hf->SetYTitle("#phi/K");
  hf->GetYaxis()->SetTitleSize(0.07);
  hf->GetYaxis()->SetLabelSize(0.05);
  hf->GetYaxis()->SetTitleOffset(0.9*0.93);
  hf->SetNdivisions(505,"y");

  double lt=0.06;

  TPad* p4=new TPad("p4","p4",0,0,1,1);
  myPadSetUp(p4,ml,mt,mr,mb);
  p4->SetFillColor(0);
  p4->SetTicky(1);
  p4->SetLogx();

  c->cd();
  p4->Draw();
  p4->cd();

  hf->GetXaxis()->SetLabelOffset(0.);
  hf->SetMinimum(0.);
  hf->SetMaximum(0.23);
  hf->SetNdivisions(509,"y");
  hf->Draw();

  TLatex* cent4=new TLatex(0.25,0.94,"d-Au, p-Pb, and A-A Data: High Mult.");
  cent4->SetNDC();
  cent4->SetTextSize(lt);
  cent4->Draw();

  double euw=0.05;
  double euw1=0.08;

  TGraphErrors* g41=new TGraphErrors(1);
  g41->SetMarkerStyle(4); g41->SetMarkerColor(4); g41->SetMarkerSize(1.3); g41->SetLineColor(4);
  TGraphErrors* b41=new TGraphErrors(1);
  b41->SetLineColor(4);
  get_phi_NA49_pp(g41,b41,euw);

  TGraphErrors* g42=new TGraphErrors(3);
  g42->SetMarkerStyle(21); g42->SetMarkerColor(4); g42->SetMarkerSize(1.3); g42->SetLineColor(4);
  TGraphErrors* b42=new TGraphErrors(3);
  b42->SetFillColor(TColor::GetColor("#aaaadd"));
  get_phi_NA49_PbPb(g42,b42,euw);

  TGraphErrors* g43=new TGraphErrors(1);
  g43->SetMarkerStyle(4); g43->SetMarkerColor(1); g43->SetMarkerSize(1.3); g43->SetLineColor(1);
  TGraphErrors* b43=new TGraphErrors(1);
  b43->SetLineColor(1);
  get_phi_STAR_pp(g43,b43,euw);

  TGraphErrors* g44=new TGraphErrors(1);
  g44->SetMarkerStyle(27); g44->SetMarkerColor(1); g44->SetMarkerSize(1.5); g44->SetLineColor(4);
  TGraphErrors* b44=new TGraphErrors(1);
  b44->SetLineColor(1);
  get_phi_STAR_dAu(g44,b44,euw);

  TGraphErrors* g45=new TGraphErrors(2);
  g45->SetMarkerStyle(33); g45->SetMarkerColor(1); g45->SetMarkerSize(1.5); g45->SetLineColor(1);
  TGraphErrors* b45=new TGraphErrors(2);
  b45->SetFillColor(17);
  get_phi_STAR_CuCu(g45,b45,euw);

  TGraphErrors* g46=new TGraphErrors(3);
  g46->SetMarkerStyle(8); g46->SetMarkerColor(1); g46->SetMarkerSize(1.3); g46->SetLineColor(1);
  TGraphErrors* b46=new TGraphErrors(3);
  b46->SetLineColor(1);
  get_phi_STAR_AuAu(g46,b46,euw);

  TGraphErrors* g47=new TGraphErrors(1);
  g47->SetMarkerStyle(4); g47->SetMarkerColor(TColor::GetColor("#238e23")); g47->SetMarkerSize(1.3); g47->SetLineColor(1);
  TGraphErrors* b47=new TGraphErrors(1);
  b47->SetLineColor(TColor::GetColor("#238e23"));
  get_phi_PHENIX_pp(g47,b47,euw);

  TGraphErrors* g48=new TGraphErrors(1);
  g48->SetMarkerStyle(8); g48->SetMarkerColor(TColor::GetColor("#238e23")); g48->SetMarkerSize(1.3); g48->SetLineColor(TColor::GetColor("#238e23"));
  TGraphErrors* b48=new TGraphErrors(1);
  b48->SetFillColor(TColor::GetColor("#aaddaa"));
  get_phi_PHENIX_AuAu(g48,b48,euw);

  TGraphErrors* g49=new TGraphErrors(4);
  g49->SetMarkerStyle(4); g49->SetMarkerColor(2); g49->SetMarkerSize(1.3); g49->SetLineColor(2);
  TGraphAsymmErrors* b49=new TGraphAsymmErrors(2);
  b49->SetFillColor(TColor::GetColor("#ddaaaa"));
  get_phi_ALICE_pp(g49,b49,euw);

  TGraphErrors* g4a=new TGraphErrors(1);
  g4a->SetMarkerStyle(21); g4a->SetMarkerColor(2); g4a->SetMarkerSize(1.3); g4a->SetLineColor(2);
  TGraphErrors* b4a=new TGraphErrors(1);
  b4a->SetFillColor(TColor::GetColor("#ddaaaa"));
  get_phi_ALICE_PbPb(g4a,b4a,euw);

  TGraphErrors* g4b=new TGraphErrors(1);
  g4b->SetMarkerStyle(25); g4b->SetMarkerColor(2); g4b->SetMarkerSize(1.3); g4b->SetLineColor(2);
  TGraphAsymmErrors* b4b=new TGraphAsymmErrors(1);
  b4b->SetFillColor(TColor::GetColor("#ddaaaa"));
  get_phi_ALICE_pPb(g4b,b4b,euw);

  TLine* sm4=new TLine(1600.,12.16/100.5,4000.,12.16/100.5);
  sm4->SetLineWidth(4);
  sm4->SetLineStyle(2);

  b48->Draw("e2same");
  b42->Draw("e2same");
  b45->Draw("e2same");
  plot_box(b41);
  plot_box(b43);
  plot_box(b44);
  plot_box(b46);
  plot_box(b47);
  b49->Draw("e2same");
  b4a->Draw("e2same");
  b4b->Draw("e2same");

  sm4->Draw();

  g41->Draw("pzsame");
  g42->Draw("pzsame");
  g47->Draw("pzsame");
  g48->Draw("pzsame");
  g43->Draw("pzsame");
  g44->Draw("pzsame");
  g45->Draw("pzsame");
  g46->Draw("pzsame");
  g49->Draw("pzsame");
  g4a->Draw("pzsame");
  g4b->Draw("pzsame");

  TGraph* q=new TGraph(1);
  q->SetMarkerColor(0);

  TGraph* g4c=new TGraph(1);
  g4c->SetMarkerStyle(21); g4c->SetMarkerSize(1.3);

  double x=0.57,y=0.2,dx=0.3,dy=0.06,dl=0.005;

  TLegend* l1=new TLegend(x,y,x+dx,y+3*dy);
  myLegendSetUp(l1,lt);
  l1->AddEntry(q,"pp","p");
  l1->AddEntry(q,"d-Au","p");
  l1->AddEntry(q,"p-Pb","p");
  l1->Draw();

  TLegend* l2=new TLegend(x,y+dl,x+dx,y+3*dy+dl);
  myLegendSetUp(l2,lt);
  l2->AddEntry(g43," ","p");
  l2->AddEntry(g44," ","p");
  l2->AddEntry(g4b," ","p");
  l2->Draw();

  x=0.76;

  TLegend* l3=new TLegend(x,y,x+dx,y+3*dy);
  myLegendSetUp(l3,lt);
  l3->AddEntry(q,"Cu-Cu","p");
  l3->AddEntry(q,"Au-Au","p");
  l3->AddEntry(q,"Pb-Pb","p");
  l3->Draw();

  TLegend* l4=new TLegend(x,y+dl,x+dx,y+3*dy+dl);
  myLegendSetUp(l4,lt);
  l4->AddEntry(g45," ","p");
  l4->AddEntry(g46," ","p");
  l4->AddEntry(g4c," ","p");
  l4->Draw();

  x=0.58; y=0.9;
  TLegend* l5=new TLegend(x,y-3*dy,x+dx,y);
  myLegendSetUp(l5,lt);
  l5->AddEntry(q,"Grand Canonical","p");
  l5->AddEntry(sm4,"Thermal Model","l");
  l5->AddEntry(q,"#it{T}_{ch} = 156 MeV","p");
  l5->Draw();

  TLatex* t1=new TLatex(24.,0.13,"NA49");
  t1->SetTextAlign(22);
  t1->SetTextSize(lt);
  t1->SetTextColor(4);
  t1->Draw();

  TLatex* t2=new TLatex(65.,0.105,"STAR");
  t2->SetTextAlign(22);
  t2->SetTextSize(lt);
  t2->Draw();

  TLatex* t3=new TLatex(200.,0.07,"PHENIX");
  t3->SetTextAlign(22);
  t3->SetTextSize(lt);
  t3->SetTextColor(TColor::GetColor("#238e23"));
  t3->Draw();

  TLatex* t4=new TLatex(2760.,0.08,"ALICE");
  t4->SetTextAlign(22);
  t4->SetTextSize(lt);
  t4->SetTextColor(2);
  t4->Draw();

  TLatex* t5=new TLatex(10.,0.044,"#bf{ALICE Preliminary}");
  t5->SetTextAlign(12);
  t5->SetTextSize(0.05);
  t5->Draw();

  TLatex* t6=new TLatex(10.,0.031,"pp #sqrt{#it{s}} = 5.02, 8 and 13 TeV");
  t6->SetTextAlign(12);
  t6->SetTextSize(0.05);
  t6->Draw();

   TLatex* t7=new TLatex(10.,0.014,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  t7->SetTextAlign(12);
  t7->SetTextSize(0.05);
  t7->Draw();

  p4->RedrawAxis();

  c->Print("phi2K_energy.eps");
  c->Print("phi2K_energy.gif");

  return;
}


void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetTextAlign(11);
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


void get_phi_NA49_pp(TGraphErrors* gv,TGraphErrors* gu,double euw){
  double x=17.3,y=0.07878,u=0.0227;
  gv->SetPoint(0,x,y);
  gu->SetPoint(0,x,y); gu->SetPointError(0,euw*x,u);
  return;
}


void get_phi_NA49_PbPb(TGraphErrors* gv,TGraphErrors* gu,double euw){
  //phi: NA49, Phys. Rev. C 78, 044907 (2008), Table VI
  //K-: NA49, Phys. Rev. C 66, 054902 (2002), Table II
  double x,phi,phiT,phiY,k,kT,kY,r,rT,rY;
  for(int j=0;j<3;j++){
    if(!j){
      x=8.8;
      phi=2.55; phiT=0.17; phiY=0.19;
      k=19.2; kT=0.5; kY=1.;
    }else if(j==1){
      x=12.3;
      phi=4.04; phiT=0.19; phiY=0.31;
      k=32.4; kT=0.6; kY=1.6;
    }else if(j==2){
      x=17.3;
      phi=8.46; phiT=0.38; phiY=0.33;
      k=51.9; kT=1.9; kY=3.;
    }
    r=phi/k; rT=standard_err(phi,k,phiT,kT); rY=standard_err(phi,k,phiY,kY);
    gv->SetPoint(j,x,r); gv->SetPointError(j,0.,rT);
    gu->SetPoint(j,x,r); gu->SetPointError(j,euw*x,rY);
  }
  return;
}


void get_phi_STAR_pp(TGraphErrors* gv,TGraphErrors* gu,double euw){
  //point from STAR, Phys. Rev. C 79, 064903 (2009) Table III, http://drupal.star.bnl.gov/STAR/files/starpublications/127/data.html
  //referring to Phys. Lett. B 612 (2005) 181 (for phi in p+p) and Phys. Rev. Lett. 92, 112301 (2004) (for K)
  double x=200.,y=0.12414,u=0.025;
  gv->SetPoint(0,x,y);
  gu->SetPoint(0,x,y); gu->SetPointError(0,euw*x,u);
  return;
}


void get_phi_STAR_dAu(TGraphErrors* gv,TGraphErrors* gu,double euw){
  //phi/K-: STAR, Phys. Rev. C 79, 064903 (2009): https://drupal.star.bnl.gov/STAR/files/starpublications/127/data.html
  double x=200./(1+3*euw),y=0.13993,uy=0.01742;
  gv->SetPoint(0,x,y); gv->SetPointError(0,0.,0.);
  gu->SetPoint(0,x,y); gu->SetPointError(0,euw*x,uy);
  return;
}


void get_phi_STAR_CuCu(TGraphErrors* gt,TGraphErrors* gy,double euw){
  //STAR, Phys. Lett. B 673, 183 (2009): Table 3
  double x,phi,phiT,phiY,k,kT,kY,r,rT,rY;
  for(int j=0;j<2;j++){
    if(!j){
      x=62.4/(1+1.5*euw);
      phi=1.3; phiT=0.03; phiY=0.2;
      k=8.536; kT=0.129; kY=0.982;
    }else if(j==1){
      x=200.;
      phi=2.3; phiT=0.03; phiY=0.3;
      k=13.239; kT=0.068; kY=1.417;
    }
    r=phi/k;
    rT=r*sqrt(pow(phiT/phi,2)+pow(kT/k,2));
    rY=r*sqrt(pow(phiY/phi,2)+pow(kY/k,2));
    gt->SetPoint(j,x,r); gt->SetPointError(j,euw*x,rT);
    gy->SetPoint(j,x,r); gy->SetPointError(j,euw*x,rY);
  }
  return;
}


void get_phi_STAR_AuAu(TGraphErrors* gv,TGraphErrors* gu,double euw){
  //points from STAR, Phys. Rev. C 79, 064903 (2009) Table III, http://drupal.star.bnl.gov/STAR/files/starpublications/127/data.html
  //STAR, Phys. Rev. C 65 041901 (2002) (for 130 GeV)
  double x,y,u;
  for(int j=0;j<3;j++){
    if(!j){x=62.4*(1+1.5*euw); y=0.13799; u=0.01906;}
    else if(j==1){x=130.; y=0.14; u=0.0399;}
    else if(j==2){x=200.*(1+3*euw); y=0.16061; u=0.03078;}
    gv->SetPoint(j,x,y);
    gu->SetPoint(j,x,y); gu->SetPointError(j,euw*x,u);
  }
  return;
}


void get_phi_PHENIX_pp(TGraphErrors* gv,TGraphErrors* gu,double euw){
  //phi from Phys. Rev. D 83, 052004 (2011), Table IX
  //Email form Viktor Riabov on 24-7-2012: "For the ratio systematic error of 9.7% should be also subtracted from corresponding error for phi."
  //K from Phys. Rev. C 83, 064903 (2011), Table VI
  double x=200./(1+3*euw),y=0.0998613037,u=0.0109547419;
  gv->SetPoint(0,x,y);
  gu->SetPoint(0,x,y); gu->SetPointError(0,euw*x,u);
  return;
}


void get_phi_PHENIX_AuAu(TGraphErrors* gv,TGraphErrors* gu,double euw){
  //ratio from Phys. Rev. C 72, 014903 (2005), http://www.phenix.bnl.gov/phenix/WWW/info/data/ppg016_data.html
  //double x=200.,y=0.0911,ut=0.0141,uy=0.0150;//phi/<K+->
  double x=200.*(1+3*euw),y=0.094,ut=0.014,uy=0.016;//phi/K-
  gv->SetPoint(0,x,y); gv->SetPointError(0,0.,ut);
  gu->SetPoint(0,x,y); gu->SetPointError(0,euw*x,uy);
  return;
}


void get_phi_ALICE_pp(TGraphErrors* gv,TGraphAsymmErrors* gu,double euw){
  double x,phi,phiT,phiL,phiU,k,kT,kY,r,rT,rL,rU,rNL,rNU;

  x=900.;
  phi=0.021; phiT=0.004; phiL=phiU=0.00297;
  k=0.182; kT=0.004; kY=0.015;

  r=phi/k; rT=standard_err(phi,k,phiT,kT); rL=standard_err(phi,k,phiL,kY); rU=standard_err(phi,k,phiU,kY);
  gv->SetPoint(0,x,r); gv->SetPointError(0,0.,rT);
  gu->SetPoint(0,x,r); gu->SetPointError(0,euw*x,euw*x,rL,rU);

  //x=3100.;
  
  x=2760.;
  r=0.113; rT=0.001; rL=rU=0.013;
  gv->SetPoint(1,x/1.075,r); gv->SetPointError(1,0.,rT);
  gu->SetPoint(1,x/1.075,r); gu->SetPointError(1,euw*x,euw*x,rL,rU);

   x=5020.;
 r=1.1290014e-01; rT=8.03423138993941298e-04; rL=rU=1.04153386996841731e-02;
  gv->SetPoint(2,x/1.05,r); gv->SetPointError(2,0.,rT);
  gu->SetPoint(2,x/1.05,r); gu->SetPointError(2,euw*x,euw*x,rL,rU);


  x=7000.;
  phi=0.032; phiT=0.0004; phiL=sqrt(pow(0.0035,2)-pow(0.03*phi,2)); phiU=sqrt(pow(0.004,2)-pow(0.062*phi,2));
  k=0.286; kT=0.0002; kY=0.016;

  r=phi/k; rT=standard_err(phi,k,phiT,kT); rL=standard_err(phi,k,phiL,kY); rU=standard_err(phi,k,phiU,kY);
  rNL=r*(1.-(1.+0.062)/(1.+0.07));
  rNU=r*((1.-0.03)/(1-0.04)-1.);
  rL=sqrt(rL*rL+rNL*rNL);
  rU=sqrt(rU*rU+rNU*rNU);

  gv->SetPoint(3,x,r); gv->SetPointError(3,0.,rT);
  gu->SetPoint(3,x,r); gu->SetPointError(3,euw*x,euw*x,rL,rU);

  x=8000.;
  phi=2*0.0335; phiT=2*0.00025; phiL=phiU=2*0.0030;
  k=0.585904; kT=0.0; kY=0.026626;

  r=phi/k; rT=standard_err(phi,k,phiT,kT); rL=standard_err(phi,k,phiL,kY); rU=standard_err(phi,k,phiU,kY);
  gv->SetPoint(4,x,r); gv->SetPointError(4,0.,rT);
  gu->SetPoint(4,x,r); gu->SetPointError(4,euw*x,euw*x,rL,rU);

  x=13000.;
  r=1.1193380e-01; rT=9.5289277e-04; rL=rU=9.6845125e-03;

  gv->SetPoint(5,x,r); gv->SetPointError(5,0.,rT);
  gu->SetPoint(5,x,r); gu->SetPointError(5,euw*x,euw*x,rL,rU);

  return;
}


void get_phi_ALICE_pPb(TGraphErrors* gv,TGraphAsymmErrors* gu,double euw){
  double x,r,rT,rL,rU;

  x=5020.;
  r=0.1290; rT=0.0013; rL=rU=0.0126;

  gv->SetPoint(0,x*1.05,r); gv->SetPointError(0,0.,rT);
  gu->SetPoint(0,x*1.05,r); gu->SetPointError(0,euw*x,euw*x,rL,rU);

  return;
}


void get_phi_ALICE_PbPb(TGraphErrors* gv,TGraphErrors* gu,double euw){
  double x=2760.,r=1.266071640e-01,rT=4.277173479e-03,rY=1.4486718e-02;
  gv->SetPoint(0,x*1.075,r);
  gv->SetPointError(0,0.,rT);
  gu->SetPoint(0,x*1.075,r);
  gu->SetPointError(0,euw*x,rY);

  //from Ajay Kumar Dash
  double x=5020.,r=0.120,rT=0.0012,rY=0.0110;
  gv->SetPoint(1,x/1.2,r);
  gv->SetPointError(1,0.,rT);
  gu->SetPoint(1,x/1.2,r);
  gu->SetPointError(1,euw*x,rY);

  return;
}


void plot_box(const TGraphErrors* g,int style=0){
  int N=g->GetN();
  const int n=N;
  TBox* b[1000];
  double x1,x2,y1,y2;
  for(int j=0;j<N;j++){
    g->GetPoint(j,x1,y1);
    x2=x1+g->GetErrorX(j);
    x1-=g->GetErrorX(j);
    y2=y1+g->GetErrorY(j);
    y1-=g->GetErrorY(j);
        
    b[j]=new TBox(x1,y1,x2,y2);
    if(!style){
      b[j]->SetLineColor(g->GetLineColor());
      b[j]->SetLineWidth(g->GetLineWidth());
      b[j]->SetLineStyle(g->GetLineStyle());
      b[j]->SetFillStyle(0);
    }else{
      b[j]->SetLineColor(g->GetFillColor());
      b[j]->SetFillColor(g->GetFillColor());
      b[j]->SetFillStyle(1);
    }
    b[j]->Draw();
  }
  return;
}


double standard_err(double A,double B,double EA,double EB){
  double ER=0.;
  if(B>0.) ER=sqrt(EA*EA/(B*B)+A*A*EB*EB/(B*B*B*B));
  return ER;
}
