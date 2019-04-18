void PlotPhiKStarToKRatio()
{
  
  
  Double_t dNdYKaon1[] ={245.821,168.392,114.858, 75.5912,46.6784,26.3643,13.2937,6.01793,2.25409};
  Double_t dNdYKaonstat1[] ={0.4809 ,0.326334 ,0.275685 ,0.200172 ,0.149861 ,0.107479 ,0.0670287, 0.0431206, 0.0206383};
  Double_t dNdYKaonsys1[] ={ 11.5, 6.86674, 4.62415, 3.09558, 1.99945, 1.1433, 0.757718, 0.329975, 0.138256};

  Double_t dNdYKStar[] = {19.5218, 16.2772, 12.0617, 8.21309, 5.26205, 2.97302, 1.6436, 0.783306, };
  Double_t dNdYKStarStat[] = {1.35978, 0.918699, 0.537506, 0.354556, 0.258265, 0.138427, 0.0788077, 0.0350609, };
  Double_t dNdYKStarSys[] = {2.38185, 1.93861, 1.35253, 1.07999, 0.553478, 0.287313, 0.214377, 0.0718446, };
  
  const Int_t kCenks = 8;
  Double_t dNchdY[kCenks] ={1764,1180 ,786 ,512  ,318  ,183  ,96.3 ,44.9};
  Double_t dNchdY1[kCenks] ={1764,1180 ,786 ,512  ,318  ,183  ,96.3 ,44.9};
  Double_t ErdNchdY1[kCenks] ={49.75,31  ,20  ,15  ,12  , 8  , 5.8, 3.4};

  Double_t dNdYPhi[] = {14.7221, 10.1724, 7.32674, 5.07992, 3.03874, 1.7717, 0.854083, 0.362676, };
  Double_t dNdYPhiStat[] = {0.145732, 0.105384, 0.0755311, 0.0516598, 0.0326515, 0.0209748, 0.0113248, 0.00646689, };
  Double_t dNdYPhiSys[] = {1.15703, 0.788032, 0.588263, 0.385745, 0.239898, 0.135048, 0.0682854, 0.0281975, };
 
  Double_t PhiByK[kCenks];
  Double_t ErPhiByKStat[kCenks];
  Double_t ErPhiByKSys[kCenks];

  Double_t KStarByK[kCenks];
  Double_t ErKStarByKStat[kCenks];
  Double_t ErKStarByKSys[kCenks];
  
  for(Int_t i = 0; i < kCenks; i++)
    {
      //if(i < kCen-1){
      dNdYKaon1[i]      = 0.5*dNdYKaon1[i];
      dNdYKaonstat1[i]   = 0.5*dNdYKaonstat1[i];
      dNdYKaonsys1[i]   = 0.5*dNdYKaonsys1[i];
      
      ErKStarByKStat[i] = RatioError(dNdYKStar[i],dNdYKStarStat[i],dNdYKaon1[i],dNdYKaonstat1[i]);
      ErKStarByKSys[i] = RatioError(dNdYKStar[i],dNdYKStarSys[i],dNdYKaon1[i],dNdYKaonsys1[i]);
      KStarByK[i] = dNdYKStar[i]/dNdYKaon1[i];
      
      ErPhiByKStat[i]  = RatioError(dNdYPhi[i],dNdYPhiStat[i],dNdYKaon1[i],dNdYKaonstat1[i]);
      ErPhiByKSys[i]  = RatioError(dNdYPhi[i],dNdYPhiSys[i],dNdYKaon1[i],dNdYKaonsys1[i]);
      PhiByK[i]  = dNdYPhi[i]/dNdYKaon1[i];

      ErdNchdY1[i] = TMath::Power(dNchdY[i],-2./3.)*ErdNchdY1[i]/3.;
      dNchdY1[i] = TMath::Power(dNchdY1[i],1./3.);
      //ErdNchdY1[i] = 0.1;
      
      //cout<<ErdNchdY1[i]<<endl;
      
      //printf("%d phi/K  %1.3f ± %1.4f (stat)  ± %1.4f (syst) \t %f\n",i,PhiByK[i],ErPhiByKStat[i],ErPhiByKSys[i],dNchdY1[i]);
      printf("%d K*/K  %1.3f ± %1.4f (stat)  ± %1.4f (syst)\n",i,KStarByK[i],ErKStarByKStat[i],ErKStarByKSys[i]);
    }
  TGraphErrors *phibykaStat1 = PlotGraph(kCenks,4,20,dNchdY1,ErdNchdY1,PhiByK,ErPhiByKStat);
  TGraphErrors *phibykaSys1 = PlotGraph(kCenks,4,20,dNchdY1,ErdNchdY1,PhiByK,ErPhiByKSys);
  phibykaStat1->SetName("gPhiByK_stat");

  TGraphErrors *KstarbykaStat1 = PlotGraph(kCenks,2,21,dNchdY1,ErdNchdY1,KStarByK,ErKStarByKStat);
  TGraphErrors *KstarbykaSys1 = PlotGraph(kCenks,2,21,dNchdY1,ErdNchdY1,KStarByK,ErKStarByKSys);
  KstarbykaStat1->SetName("gKStarByK_stat");
  
  //PrintValueGraph(phibykaStat1,phibykaSys1,phibykaSys1);

  //Color_t color[13] = { kMagenta+2, kMagenta, kRed, //K*/K value + sys
  //                   kMagenta-9,kAzure-9, kRed-9, //K*/K uncorr
  //                   kGreen+2, kMagenta, kRed,   //phi/K value+sys 
  //                   kSpring+2, kAzure-9, kRed-9,
  //                   kGreen+2}; //pp min bias INEL is black

  Color_t color[13] = { kBlue+2, kBlue, kRed, //K*/K value + sys
                       kMagenta-9,kAzure-9, kRed-9, //K*/K uncorr
                       kGreen+2, kBlue, kRed,   //phi/K value+sys 
                       kSpring+2, kAzure-9, kRed-9,
                       kGreen+2}; //pp min bias INEL is black
  
  TGraphErrors* ppt=new TGraphErrors(1);
  TBox* ppy=new TBox();
  get_phi_pp7(ppt,ppy);
  ppt->SetLineColor(1); ppt->SetMarkerColor(1); ppt->SetMarkerStyle(20); ppt->SetMarkerSize(1.3);
  ppy->SetLineColor(1); ppy->SetFillStyle(0);
  
  Int_t j;
  TGraphErrors* galkpPb=new TGraphErrors(5);
  TBox* balkpPb[5];
  TBox* balkpPb_uncorr[5];
  get_Kstar2Ka_pPb502(galkpPb,balkpPb,balkpPb_uncorr,0);
  galkpPb->SetLineColor(color[1]); 
  galkpPb->SetMarkerColor(color[1]); 
  galkpPb->SetMarkerStyle(25); 
  galkpPb->SetMarkerSize(1.5);
  for(j=0;j<5;j++){
    balkpPb[j]->SetFillStyle(0);
    balkpPb[j]->SetFillColor(color[1]);
    balkpPb[j]->SetLineColor(color[1]);

    balkpPb_uncorr[j]->SetFillStyle(1001);
    balkpPb_uncorr[j]->SetFillColor(color[4]);
    balkpPb_uncorr[j]->SetLineColor(color[4]);
  }
  
  TGraphErrors* gphipPb = 0x0;
  gphipPb = new TGraphErrors(7);
  TBox* bphipPb[7];
  TBox* bphipPb_uncorr[7];
  get_Phi2Ka_pPb502(gphipPb,bphipPb,bphipPb_uncorr,0);
  gphipPb->SetLineColor(color[7]);
  gphipPb->SetMarkerColor(color[7]); 
  gphipPb->SetMarkerStyle(24);
  gphipPb->SetMarkerSize(1.5);
  for(j=0;j<7;j++){
    bphipPb[j]->SetFillStyle(0);
    bphipPb[j]->SetFillColor(color[7]);
    bphipPb[j]->SetLineColor(color[7]);

    bphipPb_uncorr[j]->SetFillStyle(1001);
    bphipPb_uncorr[j]->SetFillColor(color[10]);
    bphipPb_uncorr[j]->SetLineColor(color[10]);
  }

  
  TGraphErrors* lnt=new TGraphErrors(6);
  TGraphErrors* lnu=new TGraphErrors(6);
  TBox* lny[6];
  get_Kstar_PbPb276_new(lnt,lnu,lny);
  //lnt->SetLineColor(kGreen+4); lnt->SetMarkerColor(kGreen+4); lnt->SetMarkerStyle(25); lnt->SetMarkerSize(1.3);
  lnt->SetLineColor(kRed); lnt->SetMarkerColor(kRed); lnt->SetMarkerStyle(25); lnt->SetMarkerSize(1.3);
  lnu->SetLineColor(33); lnu->SetFillColor(33);
  //for(Int_t j=0;j<6;j++){lny[j]->SetLineColor(kGreen+4); lny[j]->SetFillStyle(0);}
  for(Int_t j=0;j<6;j++){lny[j]->SetLineColor(kRed); lny[j]->SetFillStyle(0);}
  
  TGraphErrors* lkt=new TGraphErrors(4);//2010 data
  TGraphErrors* lku=new TGraphErrors(4);
  TBox* lky[4];
  get_Kstar_PbPb276(lkt,lku,lky);
  //lkt->SetLineColor(kGreen+4); lkt->SetMarkerColor(kGreen+4); lkt->SetMarkerStyle(25); lkt->SetMarkerSize(1.3);
  lkt->SetLineColor(kRed); lkt->SetMarkerColor(kRed); lkt->SetMarkerStyle(25); lkt->SetMarkerSize(1.3);
  lku->SetLineColor(33); lku->SetFillColor(33);
  //for(Int_t j=0;j<4;j++){lky[j]->SetLineColor(kGreen+4); lky[j]->SetFillStyle(0);}
  for(Int_t j=0;j<4;j++){lky[j]->SetLineColor(kRed); lky[j]->SetFillStyle(0);}

  TGraphErrors* lpt=new TGraphErrors(10);
  TGraphErrors* lpu=new TGraphErrors(10);
  TBox* lpy[10];

  get_phi_PbPb276(lpt,lpu,lpy);
  //lpt->SetLineColor(kGreen+4); lpt->SetMarkerColor(kGreen+4); lpt->SetMarkerStyle(24); lpt->SetMarkerSize(1.3);
  lpt->SetLineColor(kRed); lpt->SetMarkerColor(kRed); lpt->SetMarkerStyle(24); lpt->SetMarkerSize(1.3);
  lpu->SetLineColor(33); lpu->SetFillColor(33);
  //for(Int_t j=0;j<10;j++){lpy[j]->SetLineColor(kGreen+4); lpy[j]->SetFillStyle(0);}
  for(Int_t j=0;j<10;j++){lpy[j]->SetLineColor(kRed); lpy[j]->SetFillStyle(0);}
  
  TCanvas *c1 = DrawFrame("c1");
  phibykaSys1->SetMinimum(0.0);
  phibykaSys1->SetMaximum(0.5);
  phibykaSys1->GetXaxis()->SetLimits(0.,13);
  phibykaSys1->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT^{1/3}_{|#it{#eta}| < 0.5}");
  //phibykaSys1->GetYaxis()->SetTitle("#frac{2#phi}{(K^{+}+K^{-})}");
  phibykaSys1->GetYaxis()->SetTitle("Particle ratios");
  //phibykaSys1->SetFillColor(kGreen);
  phibykaSys1->SetFillStyle(0);
  phibykaSys1->SetMarkerColor(1);
  phibykaSys1->SetLineColor(1);
  phibykaSys1->Draw("5AP");
  phibykaStat1->SetMarkerColor(1);
  phibykaStat1->SetLineColor(1);
  phibykaStat1->Draw("P");
  
  KstarbykaSys1->SetFillStyle(0);
  KstarbykaSys1->SetMarkerColor(1);
  KstarbykaSys1->SetLineColor(1);
  KstarbykaSys1->Draw("5P");
  KstarbykaStat1->SetMarkerColor(1);
  KstarbykaStat1->SetLineColor(1);
  KstarbykaStat1->Draw("P");
  
  //lpu->Draw("e2same");
  lpt->Draw("pzsame");
  
  //lnu->Draw("e2same");
  lnt->Draw("pzsame");

  //lku->Draw("e2same");
  lkt->Draw("pzsame");

  //ppt->Draw("pzsame");
  //ppy->Draw("pzsame");

  //for(j=0;j<5;j++) balkpPb_uncorr[j]->Draw();//ALICE pPb
  for(j=0;j<5;j++) balkpPb[j]->Draw();

  //for(j=0;j<7;j++) bphipPb_uncorr[j]->Draw();//ALICE pPb
  for(j=0;j<7;j++) bphipPb[j]->Draw();
  
  for(j=0; j<4;  j++) lky[j]->Draw();
  for(j=0; j<6;  j++) lny[j]->Draw();
  for(j=0; j<10; j++) lpy[j]->Draw();
  
  galkpPb->Draw("pzsame");
  gphipPb->Draw("pzsame");
  KstarbykaStat1->Draw("PZ");
  phibykaStat1->Draw("PZ");

  TLegend *l1 = DrawLegend(0.38,0.89,0.4,0.9);
  l1->SetTextSize(0.04);
  l1->SetTextFont(42);
  l1->AddEntry((TObject*)0,"#bf{ALICE Preliminary}","");
  //l1->AddEntry((TObject*)0,"Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","");
  
  TLegend *l01 = DrawLegend(0.20,0.2,0.3,0.25);
  l01->SetTextSize(0.04);
  l01->SetTextFont(42);
  l01->AddEntry((TObject*)0,"#scale[1.0]{Uncertainties: stat. (bars), sys. (boxes)}", "");
  
  TLegend *l2 = DrawLegend(0.20,0.75,0.65,0.92);
  l2->SetTextSize(0.04);
  l2->SetTextFont(42);
  //l2->SetHeader("K*^{0}/K^{-}");
  l2->SetHeader(" K*^{0}/K");
  //l2->AddEntry((TObject*)0,"Pb-Pb","");
  l2->AddEntry(KstarbykaStat1,"","p");
  l2->AddEntry(lkt,"","p");
  l2->AddEntry(galkpPb,"","p");
  l2->Draw();

  TLegend *l3 = DrawLegend(0.30,0.75,0.65,0.92);
  l3->SetTextSize(0.038);
  l3->SetTextFont(42);
  //l3->SetHeader("#phi/K^{-}");
  l3->SetHeader(" #phi/K");
  //l3->AddEntry((TObject*)0,"Pb-Pb","");
  //l3->AddEntry(phibykaStat1,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","p");
  //l3->AddEntry(lpt,"Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV (PRC 95, 064606)","p");
  //l3->AddEntry(gphipPb,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (EPJC 76, 245)","p");
  l3->AddEntry(phibykaStat1,"Pb-Pb 5.02 TeV","p");
  l3->AddEntry(lpt,"Pb-Pb 2.76 TeV #scale[0.8]{(PRC 95 (2017) 064606)}","p");
  l3->AddEntry(gphipPb,"p-Pb 5.02 TeV #scale[0.8]{(EPJC 76 (2016) 245)}","p");
  l3->Draw();

  l1->Draw();
  l01->Draw();

  c1->SaveAs("ParticleRatioKStarPhiByK.gif");
  c1->SaveAs("ParticleRatioKStarPhiByK.eps");
  c1->SaveAs("ParticleRatioKStarPhiByK_pic.root");
  //c1->cd();
}
//------------------------------------------------
TGraphErrors * PlotGraph(Int_t NdataPoint, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, Double_t *X=0, Double_t *ErX=0,Double_t *Y=0, Double_t *ErY=0)
{
  TGraphErrors *gr = new TGraphErrors(NdataPoint, X ,Y , ErX, ErY);
  gr->SetTitle("");
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerColor(MarkerColor);
  gr->SetMarkerSize(1.5);
  //gr->GetXaxis()->SetTitle("p_{T} ( GeV/c)");
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(1.1);
  gr->GetXaxis()->SetTitleFont(42);
  //gr->GetXaxis()->CenterTitle(true);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetNdivisions(506);
  gr->GetXaxis()->SetTickLength(0.04);
  //gr->GetYaxis()->SetTitle("#Phi Mass (GeV/c^{2})");
  gr->GetYaxis()->SetTitleFont(42);
  //gr->GetYaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetTitleOffset(1.15);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetNdivisions(505);
  gr->GetYaxis()->SetTickLength(0.04);
  return gr;
}

//-----------------------------------------------------------
TCanvas *DrawFrame(TString o="c")
{
  TCanvas *c = new TCanvas(o.Data(),"c01",10,10,600,600);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.16);
  c->SetTicks(1,1);
  return c;
}
//--------------------------------------------------------------
Double_t RatioError(Double_t X,Double_t ErX,Double_t Y,Double_t ErY)
{
  return (X/Y)*TMath::Sqrt((ErX/X)*(ErX/X) + (ErY/Y)*(ErY/Y));
}
//==========================================
TLine *DrawLine(Double_t xmin = 0,Double_t ymin = 1, Double_t xmax = 1, Double_t ymax = 1,Int_t lStyle = 1, Int_t lColor = 1, Int_t lWidth = 1)
{
  TLine *line = new TLine(xmin,ymin,xmax,ymax);
  line->SetLineStyle(lStyle);
  line->SetLineColor(lColor);
  line->SetLineWidth(lWidth);
  //line->Draw();
  return line;
}
//==========================================
TLatex *DrawText(Double_t x = 0, Double_t y = 0,Int_t tColor = 2,TString name)
{
  TLatex* tex = new TLatex(x,y,name.Data());
  tex->SetTextSize(0.04);
  tex->SetTextColor(tColor);
  tex->SetTextFont(42);
  //tex->Draw();
  return tex;
}
//=================================
TLegend *DrawLegend(Double_t x1,Double_t y1,Double_t x2,Double_t y2)
{
  TLegend *legend = new TLegend(x1,y1,x2,y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  //legend->AddEntry(gr1,"(0 - 100) %","p");
  //legend->AddEntry(func1,"p_{0}[ 1 + 2 v_{2}^{obs} cos(2(#Phi - #Psi))]","l");
  return legend;
}
void get_phi_pp7(TGraphErrors* g,TBox* b){
  //dN_ch/deta: ALICE, Eur. Phys. J. C (2010) 68: 345–354
  //phi: ALICE, Eur. Phys. J. C (2012) 72:2183
  //K-: ALICE, in preparation
  double x=6.01,xT=0.01,xY=0.2;//sys.: +0.20/-0.12
  double xU=sqrt(xT*xT+xY*xY),X=pow(x,1),XU=2.*X*(xU/x)/3.;//factor of 2 for visibility
  double k=0.286,kT=0.0002,kY=0.016,r,rT,rL,rU;
  //double a=0.032,aT=0.0004,aL=0.0035,aU=0.004;
  double a=0.032,aT=0.0004,aL=sqrt(pow(0.0035,2)-pow(0.03*a,2)),aU=sqrt(pow(0.004,2)-pow(0.062*a,2));
  r=a/k;
  rT=r*sqrt(pow(aT/a,2)+pow(kT/k,2));
  rL=r*sqrt(pow(aL/a,2)+pow(kY/k,2));
  rU=r*sqrt(pow(aU/a,2)+pow(kY/k,2));
  g->SetPoint(0,X,r); g->SetPointError(0,XU,rT);
  b->SetX1(X-XU);
  b->SetY1(r-rL);
  b->SetX2(X+XU);
  b->SetY2(r+rU);
  Printf("pp 7 TeV phi/K Ratio = %6.4f +- %6.4f (stat) +- %6.4f (sys)", r, rT, rU);
  return;
}
void PrintValueGraph(TGraphErrors *grstat,TGraphErrors *grucorl,TGraphErrors *grsys)
{

  Double_t *y    = grstat->GetY();
  Double_t *Erystat  = grstat->GetEY();
  Double_t *Eryucol  = grucorl->GetEY();
  Double_t *Erysys   = grsys->GetEY();

  for (Int_t i = 0; i < grstat->GetN(); i++) {
    cout<<y[i]<<"\t"<<Erystat[i]<<"\t"<<Eryucol[i]<<"\t"<<Erysys[i]<<endl;
  }

  cout<<"Double_t y[7] = {";
  for(Int_t cen = 0; cen < grstat->GetN(); cen++) cout<<y[cen]<<", ";
  cout<<"};"<<endl;
  
  cout<<"Double_t ty[7] = {";
  for(Int_t cen = 0; cen < grstat->GetN(); cen++) cout<<Erystat[cen]<<", ";
  cout<<"};"<<endl;
  
  cout<<"Double_t yy[7] = {";
  for(Int_t cen = 0; cen < grstat->GetN(); cen++) cout<<Erysys[cen]<<", ";
  cout<<"};"<<endl;

  cout<<"Double_t uncor[7] = {";
  for(Int_t cen = 0; cen < grstat->GetN(); cen++) cout<<Eryucol[cen]<<", ";
  cout<<"};"<<endl;

}
void get_Kstar_PbPb276_new(TGraphErrors* t,TGraphErrors* u,TBox** b){
  double x,ux,y,ty,yy,uy;
  for(int cb=0;cb<6;cb++){                                                     // total syst            shaded (uncorr)         
    if(!cb){x=11.6956;        ux=1.461400016e-01; y=0.179479; ty=0.00855792338; yy=0.026541488340649; uy=0.0234859591397514;}
    else if(cb==1){x=10.8945; ux=1.375472698e-01; y=0.185667; ty=0.00727746905; yy=2.65747480695809994e-02; uy=2.4050786105929053996e-02;}
    else if(cb==2){x=9.88309; ux=1.262105672e-01; y=0.200729; ty=0.00932834347; yy=2.59677233646168036e-02; uy=2.3406326165695488000e-02;}
    else if(cb==3){x=8.65608; ux=8.831414558e-02; y=0.225400; ty=0.01102203695; yy=2.52619614709935006e-02; uy=2.32011900004163475e-02;}
    else if(cb==4){x=7.52285; ux=7.345605205e-02; y=0.245137; ty=0.00943639782; yy=2.51564762175866705e-02; uy=2.12017765202781316e-02;}
    else if(cb==5){x=6.38949; ux=7.116052607e-02; y=0.257756; ty=0.01105998460; yy=2.591203207009163994e-02; uy=2.2197753203779452037e-02;}

    t->SetPoint(cb,x,y); t->SetPointError(cb,ux,ty);
    u->SetPoint(cb,x,y); u->SetPointError(cb,ux,uy);
    b[cb]=new TBox(x-ux,y-yy,x+ux,y+yy);
  }
  return;
}
		    
void get_Kstar_PbPb276(TGraphErrors* t,TGraphErrors* u,TBox** b){
  double x,ux,y,ty,yy,uy;
  for(int cb=0;cb<4;cb++){
    if(!cb){x=10.6465; ux=0.13454; y=0.197871044776119; ty=0.007736119402985; yy=0.028320819037667; uy=0.026452072456978;}
    else if(cb==1){x=8.13; ux=0.0958; y=0.235555263157895; ty=0.020578947368421; yy=0.025157439187083; uy=0.022087105766683;}
    else if(cb==2){x=5.896; ux=0.0719; y=0.278084805653710; ty=0.022939929328622; yy=0.026850056476870; uy=0.022776313789396;}
    else if(cb==3){x=3.814; ux=0.06873; y=0.306378378378378; ty=0.023162162162162; yy=0.025113969030051; uy=0.019629448606742;}
    
    t->SetPoint(cb,x,y); t->SetPointError(cb,ux,ty);
    u->SetPoint(cb,x,y); u->SetPointError(cb,ux,uy);
    b[cb]=new TBox(x-ux,y-yy,x+ux,y+yy);
  }
  return;
}
void get_phi_PbPb276(TGraphErrors* t,TGraphErrors* u,TBox** b){
  t->SetPoint(0,1.169850713e+01,1.266071640e-01); t->SetPointError(0,1.461400016e-01,4.277173479e-03);
  u->SetPoint(0,1.169850713e+01,1.266071640e-01); u->SetPointError(0,1.461400016e-01,1.296099922e-02);
  b[0]=new TBox(1.155236713e+01,1.121204460e-01,1.184464713e+01,1.410938819e-01);

  t->SetPoint(1,1.089711227e+01,1.298899748e-01); t->SetPointError(1,1.375472698e-01,4.399779502e-03);
  u->SetPoint(1,1.089711227e+01,1.298899748e-01); u->SetPointError(1,1.375472698e-01,1.251583025e-02);
  b[1]=new TBox(1.075956500e+01,1.157222882e-01,1.103465954e+01,1.440576614e-01);

  t->SetPoint(2,9.885357396e+00,1.336099003e-01); t->SetPointError(2,1.262105672e-01,3.216668671e-03);
  u->SetPoint(2,9.885357396e+00,1.336099003e-01); u->SetPointError(2,1.262105672e-01,1.154274330e-02);
  b[2]=new TBox(9.759146829e+00,1.201982233e-01,1.001156796e+01,1.470215773e-01);

  t->SetPoint(3,8.657946522e+00,1.521116318e-01); t->SetPointError(3,1.022767180e-01,3.293738365e-03);
  u->SetPoint(3,8.657946522e+00,1.521116318e-01); u->SetPointError(3,1.022767180e-01,1.321885878e-02);
  b[3]=new TBox(8.555669804e+00,1.367758174e-01,8.760223240e+00,1.674474462e-01);
  
  t->SetPoint(4,7.524365204e+00,1.442048048e-01); t->SetPointError(4,8.831414558e-02,3.217651560e-03);
  u->SetPoint(4,7.524365204e+00,1.442048048e-01); u->SetPointError(4,8.831414558e-02,1.143516283e-02);
  b[4]=new TBox(7.436051058e+00,1.305999920e-01,7.612679349e+00,1.578096175e-01);

  t->SetPoint(5,6.390676528e+00,1.477874557e-01); t->SetPointError(5,7.345605205e-02,2.751514583e-03);
  u->SetPoint(5,6.390676528e+00,1.477874557e-01); u->SetPointError(5,7.345605205e-02,1.164438759e-02);
  b[5]=new TBox(6.317220476e+00,1.339075151e-01,6.464132580e+00,1.616673963e-01);

  t->SetPoint(6,5.301459192e+00,1.451979580e-01); t->SetPointError(6,7.116052607e-02,3.158255186e-03);
  u->SetPoint(6,5.301459192e+00,1.451979580e-01); u->SetPointError(6,7.116052607e-02,1.141614203e-02);
  b[6]=new TBox(5.230298666e+00,1.315815273e-01,5.372619718e+00,1.588143886e-01);

  t->SetPoint(7,4.235823584e+00,1.403154639e-01); t->SetPointError(7,7.431269446e-02,3.969243220e-03);
  u->SetPoint(7,4.235823584e+00,1.403154639e-01); u->SetPointError(7,7.431269446e-02,1.085985964e-02);
  b[7]=new TBox(4.161510890e+00,1.273011066e-01,4.310136279e+00,1.533298211e-01);

  t->SetPoint(8,3.271066310e+00,1.333000981e-01); t->SetPointError(8,6.230602496e-02,4.670836233e-03);
  u->SetPoint(8,3.271066310e+00,1.333000981e-01); u->SetPointError(8,6.230602496e-02,1.350670661e-02);
  b[8]=new TBox(3.208760285e+00,1.181722017e-01,3.333372335e+00,1.484279945e-01);
  
  t->SetPoint(9,2.375207716e+00,1.128118881e-01); t->SetPointError(9,9.453563457e-02,4.727025409e-03);
  u->SetPoint(9,2.375207716e+00,1.128118881e-01); u->SetPointError(9,9.453563457e-02,1.284109031e-02);
  b[9]=new TBox(2.280672081e+00,9.873558747e-02,2.469743350e+00,1.268881888e-01);

  return;
}
void get_Kstar2Ka_pPb502(TGraphErrors* g,TBox** b,TBox** buncorr, Bool_t plotPrelim = 0){
  //dN/deta
  // double x[5]={ 35.6, 23.2, 16.1, 9.8, 4.4};
  // double ux[5]={ 0.8,  0.5,  0.4, 0.2, 0.1};
  //(dNch/deta)^1/3
  double dndeta[5] = {35.6, 23.2, 16.1, 9.82, 4.41};
  double udndeta[5] = {0.8, 0.5, 0.4, 0.2, 0.1};
  double corrVisXsection[5] = {1.0, 1.0, 1.0, 0.998, 0.967};

  double x[5];
  double ux[5];
  for (int j=0;j<5;j++) {
    dndeta[j] *= corrVisXsection[j];
    udndeta[j]*= corrVisXsection[j];
    x[j]  = TMath::Power(dndeta[j], 1./3.);
    //ux[j] = TMath::Power(x[j], -2./3.) * udndeta[j] /3.;
    ux[j] = TMath::Power(dndeta[j], -2./3.) * udndeta[j] /3.;
  }
  

  //FINAL
  double y[5]  = {0.270104, 0.288596, 0.297597, 0.304737, 0.324889};
  double ty[5] = {0.003762, 0.003818, 0.003986, 0.004150, 0.005191};
  double yy[5] = {0.026760, 0.027251, 0.027792, 0.027855, 0.027892};
  double uncorr[5] = {0.025504, 0.025839, 0.026319, 0.026312, 0.026134};
  // 01/07/2015
  // double  y[5] = {0.2701, 0.2886, 0.2976, 0.3069, 0.3342};
  // double ty[5] = {0.0036, 0.0040, 0.0039, 0.0042, 0.0057};
  // double yy[5] = {0.0245, 0.0247, 0.0252, 0.0253, 0.0253};

  //intermediate 01/09
  // double  y[5] = {0.2825, 0.2920, 0.3017, 0.3027, 0.3413};
  // double ty[5] = {0.0048, 0.0049, 0.0050, 0.0054, 0.0068};
  // double yy[5] = {0.0257, 0.0250, 0.0256, 0.0249, 0.0258};
 
  //PRELIMINARY QM 2014
  double py[5] = {0.2841, 0.2918, 0.2998, 0.3026, 0.3411};
  double pty[5] = {0.0047, 0.0052, 0.0050, 0.0055, 0.0074};
  double pyy[5] = {0.0258, 0.0250, 0.0254, 0.0249, 0.0258};
  if (plotPrelim){
    for (Int_t j=0;j<5; j++) {
      y[j] = py[j];
      ty[j] = pty[j];
      yy[j] = pyy[j];
    }
  }
  //  for (int i=0;i<5; i++) { yy[i] = TMath::Sqrt(uncorr[i]*uncorr[i]+corr[i]*corr[i]);} //tot syst
  for(int j=0;j<5;j++){
    g->SetPoint(j,x[j],y[j]);
    g->SetPointError(j,ux[j],ty[j]);
    b[j]=new TBox(x[j]-ux[j],y[j]-yy[j],x[j]+ux[j],y[j]+yy[j]);
    buncorr[j]=new TBox(x[j]-ux[j],y[j]-uncorr[j],x[j]+ux[j],y[j]+uncorr[j]);
  }
  return;
}

//---------------------------------------------------------------
void get_Phi2Ka_pPb502(TGraphErrors* g,TBox** b,TBox** buncorr, Bool_t plotPrelim = 0){
 
  //(dNch/deta)^1/3
  double dndeta[7] = {45.15, 36.22, 30.46, 23.24, 16.08, 9.82, 4.41};
  double udndeta[7]= {1., 0.8, 0.67, 0.51, 0.35, 0.21, 0.1};
  double corrVisXsection[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.998, 0.967};
  
  double x[7];
  double ux[7];
  for (int j=0;j<7;j++) {
    dndeta[j] *= corrVisXsection[j];
    udndeta[j]*= corrVisXsection[j];
    x[j]  = TMath::Power(dndeta[j], 1./3.);
    //ux[j] = TMath::Power(x[j], -2./3.) * udndeta[j] /3.;
    ux[j] = TMath::Power(dndeta[j], -2./3.) * udndeta[j] /3.;
  }
  
  //final, 30/10/2014 by Ajay
  Double_t y[7] = {0.129025, 0.124117, 0.12537, 0.125046, 0.121316, 0.11525, 0.11931};//value
  Double_t ty[7] = {0.00128847, 0.00129688, 0.001028, 0.000796133, 0.00085392, 0.00115719, 0.00181553};//stat. err
  Double_t yy[7] = {0.0125894, 0.0112445, 0.0109847, 0.0107461, 0.0102444, 0.00893814, 0.0110205};//syst. err
  Double_t uncorr[7] = {0.00760133, 0.00570038, 0.00533615, 0.00532258, 0.00520986, 0.004639, 0.0080445};//uncorr
  
  //preliminary, QM 2014 09/05/2014 by Ajay
  //phi/0.5(K++K-)

// 0-5          0.129025        0.00128847      0.00670235      0.0125894
// 5-10 0.124117                0.00129688      0.00570038      0.0112445
// 10-20        0.12537         0.001028                0.00533615      0.0109847
// 20-40        0.125046        0.000796133     0.00532258      0.0107461
// 40-60        0.121316        0.00085392      0.00520986      0.0102444
// 60-80        0.11525         0.00115719      0.004639                0.00893814
   Double_t py[6] = {0.129025, 0.124117,        0.12537 , 0.125046, 0.121316, 0.11525 };
   Double_t pty[6] = {0.00128847, 0.00129688,0.001028,0.000796133, 0.00085392, 0.00115719};
   Double_t pyy[6] = {0.0125894, 0.0112445,0.0109847 ,0.0107461, 0.0102444,0.00893814 };
   Double_t puncorr[6] = {0.00670235,0.00570038, 0.00533615, 0.00532258, 0.00520986 ,0.004639}; 
   if (plotPrelim){
     for (Int_t j=0;j<6; j++) {
       y[j] = py[j];
       ty[j] = pty[j];
       yy[j] = pyy[j];
       uncorr[j]=puncorr[j]; 
     }
     y[6] = -10.;
     ty[6] = 0.;
     yy[6] = 0.;
   }
 
  //  for (int i=0;i<5; i++) { yy[i] = TMath::Sqrt(uncorr[i]*uncorr[i]+corr[i]*corr[i]);} //tot syst
  for(int j=0;j<7;j++){
    g->SetPoint(j,x[j],y[j]);
    g->SetPointError(j,ux[j],ty[j]);
    b[j]=new TBox(x[j]-ux[j],y[j]-yy[j],x[j]+ux[j],y[j]+yy[j]);
    buncorr[j]=new TBox(x[j]-ux[j],y[j]-uncorr[j],x[j]+ux[j],y[j]+uncorr[j]);
  }
  return;
}

void get_Kstar2Ka_PbPb276(TGraphErrors* t,TGraphErrors* u,TBox** b){
  double x,ux,y,ty,yy,uy;
  for(int cb=0;cb<4;cb++){
    if(!cb){x=10.6465; ux=0.13454; y=0.197871044776119; ty=0.007736119402985; yy=0.028320819037667; uy=0.026452072456978;}
    else if(cb==1){x=8.13; ux=0.0958; y=0.235555263157895; ty=0.020578947368421; yy=0.025157439187083; uy=0.022087105766683;}
    else if(cb==2){x=5.896; ux=0.0719; y=0.278084805653710; ty=0.022939929328622; yy=0.026850056476870; uy=0.022776313789396;}
    else if(cb==3){x=3.814; ux=0.06873; y=0.306378378378378; ty=0.023162162162162; yy=0.025113969030051; uy=0.019629448606742;}

    /*
    if(!cb){x=10.6465; ux=0.13454; y=0.19; ty=0.01; yy=0.03; uy=0.02;}
    else if(cb==1){x=8.13; ux=0.0958; y=0.23; ty=0.01; yy=0.03; uy=0.025;}
    else if(cb==2){x=5.896; ux=0.0719; y=0.27; ty=0.02; yy=0.04; uy=0.02;}
    else if(cb==3){x=3.814; ux=0.06873; y=0.30; ty=0.02; yy=0.04; uy=0.02;}
    */

    t->SetPoint(cb,x,y); t->SetPointError(cb,ux,ty);
    u->SetPoint(cb,x,y); u->SetPointError(cb,ux,uy);
    b[cb]=new TBox(x-ux,y-yy,x+ux,y+yy);
  }
  return;
}
