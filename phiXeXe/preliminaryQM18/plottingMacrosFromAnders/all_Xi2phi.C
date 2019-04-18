TBox *by[20],*bu[20],*bl[20],*bp[20];

void all_Xi2phi(){
  gROOT->LoadMacro("/Users/aknospe/Documents/Austin/resonance/3xphi/overhead.C");
  gROOT->LoadMacro("/Users/aknospe/akroot/AKToolsDraw.C");
  AKstyle();
  myOptions();
  gROOT->ForceStyle();

  ncb--;

  TCanvas* c=new TCanvas("c","",10,10,420,400);
  c->SetFillColor(0);
  c->Draw();
  c->cd();

  double ml=0.16,mr=0.01;
  TPad* p=new TPad("p","",0,0,1,1);
  myPadSetUp(p,ml,0.01,mr,0.16);
  p->SetFillColor(0);
  p->Draw();
  p->cd();
  p->SetTicky(1);
  p->SetLogx();

  int j;

  TH1F* hf=new TH1F("frame","",1,0.8,3000.);
  hf->SetMinimum(0.);
  hf->SetMaximum(0.85);
  hf->SetXTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5}");
  hf->GetXaxis()->SetTitleSize(0.07);
  //hf->GetXaxis()->SetTitleOffset(1.);
  hf->GetXaxis()->SetLabelSize(0.05);
  hf->GetXaxis()->SetLabelOffset(0.);
  hf->SetNdivisions(505,"x");
  hf->SetYTitle("#Xi/#phi");
  hf->GetYaxis()->SetTitleSize(0.07);
  hf->GetYaxis()->SetLabelSize(0.05);
  hf->GetYaxis()->SetTitleOffset(1.1);
  hf->SetNdivisions(505,"y");
  hf->Draw();

  //-----

  TGraphErrors* g=new TGraphErrors(1);
  get_Xi2phi(g);

  AKMarker(g,8,"orange");
  for(cb=0;cb<ncb;cb++){by[cb]->SetLineColor(TColor::GetColor("#ffaa00")); by[cb]->SetFillStyle(0);}
  for(cb=0;cb<ncb;cb++){bu[cb]->SetFillColor(TColor::GetColor("#00dd00")); bu[cb]->SetFillStyle(1001);}

  //-----

  TGraphErrors* glt=new TGraphErrors(5);
  //TGraphErrors* glu=new TGraphErrors(5);
  get_PbPb(glt);

  AKMarker(glt,21,1);
  //glu->SetFillColor(TColor::GetColor("#ffc0cb"));
  for(cb=0;cb<8;cb++){bl[cb]->SetLineColor(1); bl[cb]->SetFillStyle(0);}

  //-----

  TGraphErrors* gpt=new TGraphErrors(7);
  TGraphErrors* gpu=new TGraphErrors(7);
  get_pPb(gpt,gpu);

  AKMarker(gpt,27,4);
  gpu->SetFillColor(TColor::GetColor("#aaaaff"));
  for(cb=0;cb<7;cb++){bp[cb]->SetLineColor(4); bp[cb]->SetFillStyle(0);}

  //-----

  c->cd();
  p->Draw();
  p->cd();

  hf->Draw();

  //for(cb=0;cb<ncb;cb++) bu[cb]->Draw();
  for(cb=0;cb<ncb;cb++) by[cb]->Draw();
  //glu->Draw("e2same");
  //gpu->Draw("e2same");
  for(cb=0;cb<8;cb++) bl[cb]->Draw();
  for(cb=0;cb<7;cb++) bp[cb]->Draw();
  g->Draw("pzsame");
  glt->Draw("pzsame");
  gpt->Draw("pzsame");

  //giy->Draw("e2same");
  //git->Draw("pzsame");

  double lx=0.38,ldx=0.25,ly=0.18,ldy=0.08,lt=0.06,lo=0.008;
  TLegend *l1=new TLegend(lx,ly,lx+ldx,ly+3*ldy); myLegendSetUp(l1,lt);//////3
  TLegend *l2=new TLegend(lx,ly+lo,lx+ldx,ly+3*ldy+lo); myLegendSetUp(l2,lt);//////3
  TGraph* dummy=new TGraph(1); dummy->SetMarkerColor(0);

  l1->AddEntry(dummy,"pp #sqrt{#it{s}} = 13 TeV","p"); l2->AddEntry(g," ","p");
  l1->AddEntry(dummy,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","p"); l2->AddEntry(gpt," ","p");
  l1->AddEntry(dummy,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","p"); l2->AddEntry(glt," ","p");
  l1->Draw();
  l2->Draw();

  TGraphErrors* g1=new TGraphErrors(1); AKMarker(g1,1,1);
  TBox* b1=new TBox(0,0,1,1); b1->SetLineColor(1); b1->SetFillColor(0);
  TBox* b2=new TBox(0,0,1,1); b2->SetFillColor(16);

  lx=0.6; ly=0.42; ldx=0.25; lo=0.009; ldy=0.05; lt=0.04;
  TLegend* l3=new TLegend(lx,ly,lx+ldx,ly+3*ldy); myLegendSetUp(l3,lt);
  TLegend* l4=new TLegend(lx,ly+lo,lx+ldx,ly+3*ldy+lo); myLegendSetUp(l4,lt);
  l3->AddEntry(dummy,"stat.","p"); l4->AddEntry(g1," ","le");
  l3->AddEntry(dummy,"sys.","p"); l4->AddEntry(b1," ","f");
  l3->AddEntry(dummy,"uncorr. sys.","p"); l4->AddEntry(b2," ","f");
  l3->Draw();
  l4->Draw();

  lx=0.19; ly=0.98; ldy=0.06; lt=0.04;
  TLatex* t1=new TLatex(lx,ly,"#bf{ALICE Preliminary}");
  t1->SetTextAlign(13);
  t1->SetNDC();
  t1->SetTextSize(lt);
  t1->Draw();

  TLatex* t2=new TLatex(lx,ly-ldy,"|#it{y}| < 0.5 ");
  t2->SetTextAlign(13);
  t2->SetNDC();
  t2->SetTextSize(lt);
  //t2->Draw();

  c->cd();
  draw_sqrt(0.611,0.1935,0.0825,0.0673,0.11,2);
  draw_sqrt(0.577,0.2715,0.0825,0.0673,0.11,2);
  draw_sqrt(0.5295,0.3683,0.035,0.05125,0.095,2);

  c->SaveAs("all_Xi2phi.eps");
  c->SaveAs("all_Xi2phi.gif");
  c->SaveAs("all_Xi2phi.pdf");
  c->SaveAs("all_Xi2phi.root");

  return;
}


void draw_sqrt(double x,double y,double dx,double dy,double s,double w){
  double f=0.0011;
  TLine* a=new TLine(x,y,x+s*dy,y+dy);
  a->SetLineWidth(w);
  TLine* b=new TLine(x+s*dy-f,y+dy,x+s*dy+dx-f,y+dy);
  b->SetLineWidth(w);
  a->Draw();
  b->Draw();
  return;
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


void get_Xi2phi(TGraphErrors* g){
  TFile* f=TFile::Open("/Users/aknospe/Documents/Austin/resonance/3xphi/COMBINE/particle_ratios.root");
  TFile* fm=TFile::Open("/Users/aknospe/Documents/Austin/resonance/3xphi/mult.root");
  if(!f || !fm) return;

  TH3D* h=(TH3D*) f->Get("phi_dNdyH_over_Xix");
  if(!h) return;

  TH2D* m=(TH2D*) fm->Get("mult");
  if(!m) return;

  double x,dx,v,t,y,u;
  int j;
  for(j=0;j<h->GetNbinsZ()-1;j++){
    if(j!=2){
      gbb(m,j+1,1,x,dx);
      gbb(h,1148,1,j+1,v,t);
      y=gbc(h,1148,2,j+1);
      u=AKquad_sum_sg(y,-gbc(h,1148,"corr",j+1));
    }

    t/=v*v;
    y/=v*v;
    u/=v*v;
    v=1/v;

    g->SetPoint(j,x,v);
    g->SetPointError(j,0.,t);
    by[j]=new TBox(x-dx,v-y,x+dx,v+y);
    bu[j]=new TBox(x-dx,v-u,x+dx,v+u);
  }

  f->Close();
  fm->Close();

  return;
}


void get_PbPb(TGraphErrors* gt){
  TFile* f=TFile::Open("XiOverPhi.root");
  if(!f) return;

  TGraphErrors* T=(TGraphErrors*) f->Get("gRatioStat");
  TGraphErrors* Y=(TGraphErrors*) f->Get("gRatioSyst");
  if(!T || !Y) return;

  double x,ux,v,t,y;
  for(int j=0;j<8;j++){
    T->GetPoint(j,x,v);
    t=T->GetErrorY(j);
    y=Y->GetErrorY(j);
    ux=Y->GetErrorX(j);

    gt->SetPoint(j,x,v);
    gt->SetPointError(j,0.,t);
    bl[j]=new TBox(x-ux,v-y,x+ux,v+y);
  }

  f->Close();
  return;
}


void get_PbPb276(TGraphErrors* gt,TGraphErrors* gu){
  int j;
  double x[5],ux[5],v[5],t[5],u[5],b[5],X,UX,y;

  j=0;
  x[j]=0.; ux[j]=0.; v[j]=3.697881; t[j]=0.1148481; u[j]=0.4812985; b[j]=0.04467803; j++;
  x[j]=9.885357396; ux[j]=1.336099003e-01; v[j]=3.520333; t[j]=0.09534332; u[j]=0.448539; b[j]=0.06063911; j++;
  x[j]=8.13; ux[j]=0.0958; v[j]=3.720423; t[j]=0.07080089; u[j]=0.4959548; b[j]=0.08654105; j++;
  x[j]=5.896; ux[j]=0.0719; v[j]=3.842173; t[j]=0.07047916; u[j]=0.4928767; b[j]=0.132317; j++;
  x[j]=3.814; ux[j]=0.06873; v[j]=4.216453; t[j]=0.1401695; u[j]=0.569145; b[j]=0.1702702; j++;

  for(j=1;j<5;j++){
    X=pow(x[j],3);
    UX=X*3*ux[j]/x[j];
    x[j]=X;
    ux[j]=UX;
  }

  X=(pow(1.169850713e+01,3)+pow(1.089711227e+01,3))/2;
  UX=(pow(1.169850713e+01,3)*3*1.266071640e-01/1.169850713e+01 + pow(1.089711227e+01,3)*3*1.298899748e-01/1.089711227e+01)/2;
  x[0]=X;
  ux[0]=UX;

  //-----

  for(j=0;j<5;j++){
    v[j]/=2.; t[j]/=2.; u[j]/=2.; b[j]/=2.;
    t[j]/=v[j]*v[j];
    u[j]=AKquad_sum_sg(u[j],-b[j]);
    y=AKquad_sum(u[j],0.090009*v[j],0.03004808*v[j]);
    u[j]/=v[j]*v[j];
    y/=v[j]*v[j];
    v[j]=1/v[j];
    
    gt->SetPoint(j,x[j],v[j]);
    gt->SetPointError(j,0.,t[j]);
    gu->SetPoint(j,x[j],v[j]);
    gu->SetPointError(j,ux[j],u[j]);
    bl[j]=new TBox(x[j]-ux[j],v[j]-y,x[j]+ux[j],v[j]+y);
  }
}


void get_pPb(TGraphErrors* gt,TGraphErrors* gu){
  double x[7],ux[7];
  double a[7],at[7],au[7],ac[7],ay;
  double b[7],bt[7],bu[7],by[7];
  double c[7],ct[7],cu[7],cy[7];
  double d,dt,du,dy,r,rt,ru,ry;

  int j=0;
  a[j]=0.3767; at[j]=0.0037; au[j]=0.0196; ac[j]=0.0226; j++;
  a[j]=0.2878; at[j]=0.0030; au[j]=0.0138; ac[j]=0.0173; j++;
  a[j]=0.2436; at[j]=0.0020; au[j]=0.0118; ac[j]=0.0146; j++;
  a[j]=0.1850; at[j]=0.0012; au[j]=0.0095; ac[j]=0.0111; j++;
  a[j]=0.12285; at[j]=0.00085; au[j]=0.00640; ac[j]=0.00737; j++;
  a[j]=0.06948; at[j]=0.00069; au[j]=0.00369; ac[j]=0.00417; j++;
  a[j]=0.02973; at[j]=0.00045; au[j]=0.00228; ac[j]=0.00178; j++;

  j=0;
  b[j]=0.1180; bt[j]=0.0013; bu[j]=0.0025; by[j]=0.0082; j++;
  b[j]=0.0923; bt[j]=0.0011; bu[j]=0.0019; by[j]=0.0063; j++;
  b[j]=0.0746; bt[j]=0.0007; bu[j]=0.0015; by[j]=0.0052; j++;
  b[j]=0.0545; bt[j]=0.0004; bu[j]=0.0010; by[j]=0.0038; j++;
  b[j]=0.0356; bt[j]=0.0005; bu[j]=0.0007; by[j]=0.0026; j++;
  b[j]=0.0198; bt[j]=0.0003; bu[j]=0.0005; by[j]=0.0015; j++;
  b[j]=0.0072; bt[j]=0.0003; bu[j]=0.0003; by[j]=0.0009; j++;

  j=0;
  c[j]=0.1174; ct[j]=0.0013; cu[j]=0.0025; cy[j]=0.0082; j++;
  c[j]=0.0938; ct[j]=0.0011; cu[j]=0.0019; cy[j]=0.0075; j++;
  c[j]=0.0754; ct[j]=0.0007; cu[j]=0.0015; cy[j]=0.0060; j++;
  c[j]=0.0555; ct[j]=0.0005; cu[j]=0.0010; cy[j]=0.0047; j++;
  c[j]=0.0370; ct[j]=0.0004; cu[j]=0.0007; cy[j]=0.0039; j++;
  c[j]=0.0201; ct[j]=0.0003; cu[j]=0.0005; cy[j]=0.0017; j++;
  c[j]=0.0076; ct[j]=0.0002; cu[j]=0.0003; cy[j]=0.0010; j++;

  j=6;
  x[j]=(4.2+4.4)/2; ux[j]=(4.4-4.2)/2; j--;
  x[j]=(9.6+10.0)/2; ux[j]=(10.0-9.6)/2; j--;
  x[j]=(15.7+16.5)/2; ux[j]=(16.5-15.7)/2; j--;
  x[j]=(22.7+23.7)/2; ux[j]=(23.7-22.7)/2; j--;
  x[j]=(29.8+31.2)/2; ux[j]=(31.2-29.8)/2; j--;
  x[j]=(35.4+37.0)/2; ux[j]=(37.0-35.4)/2; j--;
  x[j]=(44.0+46.0)/2; ux[j]=(46.0-44.0)/2; j--;

  for(j=0;j<7;j++){
    ay=AKquad_sum(au[j],ac[j]);
    d=b[j]+c[j];
    dt=AKquad_sum(bt[j],ct[j]);
    du=bu[j]+cu[j];
    dy=by[j]+cy[j];

    r=d/a[j];
    rt=standard_err(d,a[j],dt,at[j]);
    ru=standard_err(d,a[j],du,au[j]);
    ry=standard_err(d,a[j],dy,ay[j]);

    gt->SetPoint(j,x[j],r);
    gt->SetPointError(j,0.,rt);
    gu->SetPoint(j,x[j],r);
    gu->SetPointError(j,ux[j],ru);
    bp[j]=new TBox(x[j]-ux[j],r-ry,x[j]+ux[j],r+ry);
  }

  return;
}
