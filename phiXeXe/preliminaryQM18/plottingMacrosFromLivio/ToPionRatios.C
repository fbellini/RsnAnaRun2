const int n_particles = /*4*/8;
const int n_systems = 3;
const int n_htypes = 3; //3 types of histo: stat, tot syst e syst_unc
TGraphErrors *gr[n_particles][n_systems][n_htypes];
TGraphErrors *gr_STAR[n_particles][n_htypes];
TGraphErrors *gr_PHEN[n_particles][n_htypes];

enum{lpi, lk, lp, lk0s, llam, lphi, lxi, lom};
//enum{lk0s, llam, lxi, lom};
enum{lpp, lppb, lpbpb};
enum{lstat, lsyst, lsyst_unc};
bool presence[n_particles][n_systems][n_htypes];

TLegend *leg_th = 0x0;
TLegend *leg_PbPb_th = 0x0;
float leg_textsize_data = 0.04;
float leg_textsize_MC = 0.04;
float leg_textsize_Ther = 0.04;
int textfont = 42;

bool linlin=kTRUE;
TPad *padLeg_plotmode6 = 0x0;

void DefinePresence(){
  
  for(int i_part=lpi; i_part<n_particles; i_part++)
    for(int i_system=0; i_system<n_systems; i_system++)
      for(int i_htype=0; i_htype<3; i_htype++) presence[i_part][i_system][i_htype] = kFALSE;
  
  //ad-hoc selection
    
    presence [lk0s] [lpp] [lstat] 	= kTRUE;
    presence [lk0s] [lpp] [lsyst] 	= kTRUE;
    presence [lk0s] [lpp] [lsyst_unc] 	= kTRUE;
    presence [llam] [lpp] [lstat] 	= kTRUE;
    presence [llam] [lpp] [lsyst] 	= kTRUE;
    presence [llam] [lpp] [lsyst_unc] 	= kTRUE;
    presence [lxi]  [lpp] [lstat] 	= kTRUE;
    presence [lxi]  [lpp] [lsyst] 	= kTRUE;
    presence [lxi]  [lpp] [lsyst_unc] 	= kTRUE;
    presence [lom]  [lpp] [lstat] 	= kTRUE;
    presence [lom]  [lpp] [lsyst] 	= kTRUE;
    presence [lom]  [lpp] [lsyst_unc] 	= kTRUE;
    
    presence [lk0s] [lppb] [lstat] 	= kTRUE;
    presence [lk0s] [lppb] [lsyst] 	= kTRUE;
    presence [lk0s] [lppb] [lsyst_unc] 	= kTRUE;
    presence [llam] [lppb] [lstat] 	= kTRUE;
    presence [llam] [lppb] [lsyst] 	= kTRUE;
    presence [llam] [lppb] [lsyst_unc] 	= kTRUE;
    presence [lxi]  [lppb] [lstat] 	= kTRUE;
    presence [lxi]  [lppb] [lsyst] 	= kTRUE;
    presence [lxi]  [lppb] [lsyst_unc] 	= kTRUE;
    presence [lom]  [lppb] [lstat] 	= kTRUE;
    presence [lom]  [lppb] [lsyst] 	= kTRUE;
    presence [lom]  [lppb] [lsyst_unc] 	= kTRUE;
    
    presence [lk0s] [lpbpb] [lsyst] 	= kTRUE;
    presence [llam] [lpbpb] [lsyst] 	= kTRUE;
    presence [lxi]  [lpbpb] [lsyst] 	= kTRUE;
    presence [lom]  [lpbpb] [lsyst] 	= kTRUE;
    
}

void ToPionRatios(int plotmode=1 , bool linearlinear=kTRUE){
  
  //plotmode meanings:
  // 1 = all systems, as in arXiv.v0
  // 11= all systems, similar to arXiv.v0, but with grayed pPb and PbPb
  // 2 = pp and pPb + arrow showing fit to PbPb in 100-2000
  // 3 = ?
  // 4 = pp and pPb alone
  // 5 = Nch/R^2 where R is taken from Fig.1 right in http://cds.cern.ch/record/1971436/files/ALICE-PUBLIC-2014-003.pdf
  // 6 = divided x-axis -->  pp-pPb | PbPb (Npart)   linlin is valid for this mode. Both plots are in linear x-scale
  
  linlin = linearlinear;
  
  //style
    gROOT->LoadMacro("LoadPlotsStyle.C");
    LoadPlotsStyle(2);
  
  //building graphs
    DefinePresence();
    TFile *fin = new TFile("IdentifiedParticlesRatioPlots_pp_pA_AA_iterative.root","read");
    TFile *fin_2 = new TFile("IdentifiedParticlesRatioPlots_AA5.root","read");
    fin->cd();
    for(int i_part=lpi+1; i_part<n_particles; i_part++){
      for(int i_system=0; i_system<n_systems; i_system++){
        for(int i_htype=0; i_htype<3; i_htype++){
          if(presence[i_part][i_system][i_htype]) {
	    if(i_system!=lpbpb) gr[i_part][i_system][i_htype] = (TGraphErrors*) fin->Get(Form("gr[%d][%d][%d]",i_part,i_system,i_htype));
	    else gr[i_part][i_system][i_htype] = (TGraphErrors*) fin_2->Get(Form("gr[%d][%d][%d]",i_part,i_system,i_htype));
	  }
        }
      }
    }
    fin->Close();
//     if(plotmode==6 || plotmode ==11)  {
//       LoadStarGraphs(plotmode);
//       LoadPhenixGraphs(plotmode);
//     }
  
  //formatting graphs
    if(plotmode!=6 && plotmode !=11){
      float gralpha[n_systems] = {0.5, 1., 1.};
      //			                     lpi, lk, lp,   lk0s,    llam, lphi,      lxi,    lom
      int grcol[n_systems][n_particles] =    { { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //pp
                                             { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //pPb
					     { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1}   //PbPb
					   };
      int grmarker[n_systems][n_particles] = { { 0,  0,  0, 20,     20,   0,    20,      20},  //pp
                                             { 0,  0,  0, 27,     27,   0,    27,      27},  //pPb
					     { 0,  0,  0, 25,     25,   0,    25,      25}   //PbPb
					   };
      float grmarsiz[n_systems][n_particles] = { { 0., 0., 0., 2., 2., 0., 2., 2.},  //pp
					       { 0., 0., 0., 2.5, 2.5, 0., 2.5, 2.5},  //pPb
					       { 0., 0., 0., 2., 2., 0., 2., 2.}   //PbPb
					   };
    }
    else{
      float gralpha[n_systems] = {0.7, 0.7, 0.7};
      int grcol[n_systems][n_particles] =    { { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //pp
                                               { 0,  0,  0, kGray/*+1*/,   kGray/*+1*/,    0, kGray/*+1*/, kGray/*+1*/},  //pPb
  					     { 0,  0,  0, kGray/*+1*/,   kGray/*+1*/,    0, kGray/*+1*/, kGray/*+1*/}   //PbPb
  					   };
      int grmarker[n_systems][n_particles] = { { 0,  0,  0, 20,     20,   0,    20,      20},  //pp
                                               { 0,  0,  0, 21,     21,   0,    21,      21},  //pPb
  					     { 0,  0,  0, 33,     33,   0,    33,      33}   //PbPb
  					   };
      float grmarsiz[n_systems][n_particles] = { { 0., 0., 0., 2., 2., 0., 2., 2.},  //pp
  					       { 0., 0., 0., 2., 2., 0., 2., 2.},  //pPb
  					       { 0., 0., 0., 2.5, 2.5, 0., 2.5, 2.5}   //PbPb
  					   };
    }
    
    char drawop[n_htypes][10] = {"P","E5","E5"};
    
    for(int i_part=lk0s; i_part<n_particles; i_part++){
      for(int i_system=0; i_system<n_systems; i_system++){
	for(int i_htype=0; i_htype<3; i_htype++){
	  if(!presence[i_part][i_system][i_htype]) continue;
	  gr[i_part][i_system][i_htype]->SetLineColor(grcol[i_system][i_part]);
	  gr[i_part][i_system][i_htype]->SetMarkerColor(grcol[i_system][i_part]);
	  gr[i_part][i_system][i_htype]->SetMarkerColorAlpha(grcol[i_system][i_part], gralpha[i_system]);
	  gr[i_part][i_system][i_htype]->SetMarkerStyle(grmarker[i_system][i_part]);
	  gr[i_part][i_system][i_htype]->SetMarkerSize(grmarsiz[i_system][i_part]);
	  gr[i_part][i_system][i_htype]->SetFillStyle(0);
	  gr[i_part][i_system][i_htype]->SetFillColor(grcol[i_system][i_part]);
	  if(i_htype==lsyst_unc) {
	    gr[i_part][i_system][i_htype]->SetFillStyle(1001);
	    gr[i_part][i_system][i_htype]->SetFillColorAlpha(grcol[i_system][i_part],0.75);
	  }
	}
      }
    }
    if(plotmode==6 /*|| plotmode==11*/){
      gr_STAR[llam][lsyst]->SetLineColor(grcol[lpbpb][llam]+1);
      gr_STAR[llam][lsyst]->SetMarkerColor(grcol[lpbpb][llam]+1);
      gr_STAR[llam][lsyst]->SetMarkerColorAlpha(grcol[lpbpb][llam]+1, gralpha[lpbpb]);
      gr_STAR[llam][lsyst]->SetMarkerStyle(/*grmarker[lpbpb][llam]*/30);
      gr_STAR[llam][lsyst]->SetMarkerSize(2);
      gr_STAR[llam][lsyst]->SetFillStyle(0);
      gr_STAR[llam][lsyst]->SetFillColor(grcol[lpbpb][llam]+1);
      gr_STAR[lxi][lsyst]->SetLineColor(grcol[lpbpb][lxi]+1);
      gr_STAR[lxi][lsyst]->SetMarkerColor(grcol[lpbpb][lxi]+1);
      gr_STAR[lxi][lsyst]->SetMarkerColorAlpha(grcol[lpbpb][lxi]+1, gralpha[lpbpb]);
      gr_STAR[lxi][lsyst]->SetMarkerStyle(/*grmarker[lpbpb][lxi]*/30);
      gr_STAR[lxi][lsyst]->SetMarkerSize(2);
      gr_STAR[lxi][lsyst]->SetFillStyle(0);
      gr_STAR[lxi][lsyst]->SetFillColor(grcol[lpbpb][lxi]+1);
      gr_STAR[lom][lsyst]->SetLineColor(grcol[lpbpb][lom]+1);
      gr_STAR[lom][lsyst]->SetMarkerColor(grcol[lpbpb][lom]+1);
      gr_STAR[lom][lsyst]->SetMarkerColorAlpha(grcol[lpbpb][lom]+1, gralpha[lpbpb]);
      gr_STAR[lom][lsyst]->SetMarkerStyle(/*grmarker[lpbpb][lom]*/30);
      gr_STAR[lom][lsyst]->SetMarkerSize(2);
      gr_STAR[lom][lsyst]->SetFillStyle(0);
      gr_STAR[lom][lsyst]->SetFillColor(grcol[lpbpb][lom]+1);
      gr_PHEN[llam][lsyst]->SetLineColor(grcol[lpbpb][llam]+1);
      gr_PHEN[llam][lsyst]->SetMarkerColor(grcol[lpbpb][llam]+1);
      gr_PHEN[llam][lsyst]->SetMarkerColorAlpha(grcol[lpbpb][llam]+1, gralpha[lpbpb]);
      gr_PHEN[llam][lsyst]->SetMarkerStyle(/*grmarker[lpbpb][llam]*/24);
      gr_PHEN[llam][lsyst]->SetMarkerSize(1.5);
      gr_PHEN[llam][lsyst]->SetFillStyle(0);
      gr_PHEN[llam][lsyst]->SetFillColor(grcol[lpbpb][llam]+1);
//       gr_PHEN[lxi][lsyst]->SetLineColor(grcol[lpbpb][lxi]+1);
//       gr_PHEN[lxi][lsyst]->SetMarkerColor(grcol[lpbpb][lxi]+1);
//       gr_PHEN[lxi][lsyst]->SetMarkerColorAlpha(grcol[lpbpb][lxi]+1, gralpha[lpbpb]);
//       gr_PHEN[lxi][lsyst]->SetMarkerStyle(/*grmarker[lpbpb][lxi]*/25);
//       gr_PHEN[lxi][lsyst]->SetMarkerSize(1);
//       gr_PHEN[lxi][lsyst]->SetFillStyle(0);
//       gr_PHEN[lxi][lsyst]->SetFillColor(grcol[lpbpb][lxi]+1);
//       gr_PHEN[lom][lsyst]->SetLineColor(grcol[lpbpb][lom]+1);
//       gr_PHEN[lom][lsyst]->SetMarkerColor(grcol[lpbpb][lom]+1);
//       gr_PHEN[lom][lsyst]->SetMarkerColorAlpha(grcol[lpbpb][lom]+1, gralpha[lpbpb]);
//       gr_PHEN[lom][lsyst]->SetMarkerStyle(/*grmarker[lpbpb][lom]*/25);
//       gr_PHEN[lom][lsyst]->SetMarkerSize(1);
//       gr_PHEN[lom][lsyst]->SetFillStyle(0);
//       gr_PHEN[lom][lsyst]->SetFillColor(grcol[lpbpb][lom]+1);
    }
    
  
  // introduce arrows at the PbPb limit (useful for plotmode==2)  
    if(plotmode==2){
      double PbPblim_lam = 0;
      double PbPblim_k0s = 0;
      double PbPblim_xi = 0;
      double PbPblim_om = 0;
      TF1 *f_line_fit = new TF1("f_line_fit","[0]+[1]*x",100,2000);
      f_line_fit->SetParameters(100,0);
      f_line_fit->FixParameter(1,0);
      gr[lk0s][lpbpb][lsyst]->Fit(f_line_fit,"QR0","",100,2000);
      PbPblim_k0s = f_line_fit->GetParameter(0);
      gr[llam][lpbpb][lsyst]->Fit(f_line_fit,"QR0","",100,2000);
      PbPblim_lam = f_line_fit->GetParameter(0);
      gr[lxi][lpbpb][lsyst]->Fit(f_line_fit,"QR0","",100,2000);
      PbPblim_xi = f_line_fit->GetParameter(0);
      gr[lom][lpbpb][lsyst]->Fit(f_line_fit,"QR0","",100,2000);
      PbPblim_om = f_line_fit->GetParameter(0);
    }

  // for drawop==6 we change x values for PbPb  
    if(plotmode==6){
      double Npart_k0S_lam[10] = {382.7,329.4, 260.1, 157.2,68.56,22.5,7.485};
      double eNpart_k0S_lam[10] = {3,4.3,3.9,3.1,2.0,0.8,0.225};
      double Npart_xi_om[10]   = {356.1,260.1,157.2,68.6,22.5};
      double eNpart_xi_om[10]  = {3.6,3.9,3.1,2.0,0.8};
      // 2La/pi PbPb lambda is a special case!
      TGraphAsymmErrors *gr_laPbPb = 0x0;
      TFile *fin_L = new TFile("./cLambda2Pi.root","read");
      TCanvas *cLambda2Pi = fin_L->Get("cLambda2Pi");
      TList *list = cLambda2Pi->GetListOfPrimitives();
      gr_laPbPb = (TGraphAsymmErrors*) list->At(5);
      fin_L->Close();
      for(int i_point=0; i_point<gr_laPbPb->GetN(); i_point++){
        double x_gr, y_gr;
        double e_y_gr;
        gr_laPbPb->GetPoint(i_point,x_gr,y_gr);
        e_y_gr = gr_laPbPb->GetErrorYlow(i_point);
        gr[llam][lpbpb][lsyst]->SetPoint(i_point,Npart_k0S_lam[i_point],y_gr);
        ((TGraphAsymmErrors*)gr[llam][lpbpb][lsyst])->SetPointError(i_point,eNpart_k0S_lam[i_point],eNpart_k0S_lam[i_point],e_y_gr,e_y_gr);
      }
      for(int i_part=lk0s; i_part<n_particles; i_part++){
        for(int i_system=0; i_system<n_systems; i_system++){
  	for(int i_htype=0; i_htype<3; i_htype++){
  	  if(!presence[i_part][i_system][i_htype]) continue;
  	  if(i_part==llam) continue;
  	  if(i_system!=lpbpb) continue;
            for(int i_point=0; i_point<gr[i_part][i_system][i_htype]->GetN(); i_point++){
              double x_gr, y_gr;
              double e_y_gr;
              gr[i_part][i_system][i_htype]->GetPoint(i_point,x_gr,y_gr);
              e_y_gr = gr[i_part][i_system][i_htype]->GetErrorYlow(i_point);
              if(i_part==lk0s){
  	      gr[i_part][i_system][i_htype]->SetPoint(i_point,Npart_k0S_lam[i_point],y_gr);
                gr[i_part][i_system][i_htype]->SetPointError(i_point,eNpart_k0S_lam[i_point],e_y_gr);
  	    }
  	    else{
  	      if(4-i_point<0) continue;
  	      gr[i_part][i_system][i_htype]->SetPoint(i_point,Npart_xi_om[4-i_point],y_gr);
                gr[i_part][i_system][i_htype]->SetPointError(i_point,eNpart_xi_om[4-i_point],e_y_gr);
  	    }
            }
  	}
        }
      }
    }
  
  //scale
    if(plotmode==1 || plotmode==6 || plotmode==11) {float scalefactor[n_particles] = {0,0,0,1,2,0,6,16};}
    else if(plotmode==2 || plotmode==4 || plotmode==5) {float scalefactor[n_particles] = {0,0,0,1,2,0,6,16};}
    else if(plotmode==3) {float scalefactor[n_particles] = {0,0,0, 1,2.5,0,12,45};}
    for(int i_part=lk0s; i_part<n_particles; i_part++){
     if(presence[i_part][lpp][lsyst]) ScaleGraph(i_part,scalefactor[i_part]);
    }
    if(plotmode==6 /*|| plotmode==11*/) {
      ScaleGraph_STAR(llam,scalefactor[llam]);
      ScaleGraph_STAR(lxi,scalefactor[lxi]);
      ScaleGraph_STAR(lom,scalefactor[lom]);
      ScaleGraph_PHEN(llam,scalefactor[llam]);
      //ScaleGraph_PHEN(lxi,scalefactor[lxi]);
      //ScaleGraph_PHEN(lom,scalefactor[lom]);
    }
printf("here\n");

  //frame
    if(plotmode==1 || plotmode==11) {
      //double xrangeframe[2] = {2.,2000};
      double xrangeframe[2] = {2.,3000};
      double yrangeframe[2] = {0.9e-3,1.9e-1};
    }
    else if(plotmode==2) {
      double xrangeframe[2] = {1.5,150};
      double yrangeframe[2] = {9e-4,0.2};
    }
    else if(plotmode==3) {
      double xrangeframe[2] = {1,9e04};
      double yrangeframe[2] = {0,0.21};
    }
    else if(plotmode==4) {
      double xrangeframe[2] = {1.5,70};
      double yrangeframe[2] = {9e-4,0.2};
    }
    else if(plotmode==5) {
      double xrangeframe[2] = {2.,30};
      double yrangeframe[2] = {0.9e-3,1.9e-1};
    }
    else if(plotmode==6) {
      if(!linlin) {double xrangeframe[2] = {1.5,70};}
      else {double xrangeframe[2] = {0,49};}
      double yrangeframe[2] = {0.9e-3,1.9e-1};
    }
    TH1D *frame = new TH1D("frame","",1,xrangeframe[0],xrangeframe[1]);
    frame->GetYaxis()->SetRangeUser(yrangeframe[0],yrangeframe[1]);
    frame->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|< 0.5}");
    if(plotmode==5) frame->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|< 0.5} / R^{2} (1/fm^{2})");
    frame->GetYaxis()->SetTitle("Ratio of yields to (#pi^{#plus}+#pi^{#minus})");
    frame->GetXaxis()->SetLabelOffset(-0.025);
    if(plotmode==6){
      if(linlin){
        double sizeleft = 0.76;
//         frame->GetXaxis()->SetTickLength(0.02);
//         frame->GetYaxis()->SetTickLength(0.025);
        frame->GetXaxis()->CenterTitle();
        frame->GetXaxis()->SetNdivisions(206);
        frame->GetXaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetLabelSize(0.06);
        frame->GetXaxis()->SetTitleOffset(0.6);
        frame->GetXaxis()->SetLabelOffset(-0.02);
        frame->GetYaxis()->SetLabelSize(0.062);
        frame->GetYaxis()->SetLabelOffset(0.01);
        frame->GetYaxis()->SetTitleSize(0.07);
        frame->GetYaxis()->SetTitleOffset(1.43);
        TH1D *frame_PbPb = new TH1D("frame_PbPb","",1,-70,470);   //good if lin x
        frame_PbPb->GetYaxis()->SetRangeUser(yrangeframe[0],yrangeframe[1]);
        frame_PbPb->GetXaxis()->SetTickLength(0.02*0.8);
        frame_PbPb->GetYaxis()->SetTickLength(0.025/(1-sizeleft)*sizeleft);
        frame_PbPb->GetXaxis()->SetNdivisions(103/*,kFALSE*/);
        frame_PbPb->GetXaxis()->CenterTitle();
        frame_PbPb->GetXaxis()->SetTitle("#LT#it{N}_{part}#GT");
        frame_PbPb->GetXaxis()->SetTitleSize(0.2);
        frame_PbPb->GetXaxis()->SetLabelSize(0.19);
        frame_PbPb->GetXaxis()->SetTitleOffset(0.18);
        frame_PbPb->GetXaxis()->SetLabelOffset(-0.125);
      }
      else{
        double sizeleft = 0.76;
//         frame->GetXaxis()->SetTickLength(0.02);
//         frame->GetYaxis()->SetTickLength(0.025);
        frame->GetXaxis()->CenterTitle();
        frame->GetXaxis()->SetNdivisions(206);
        frame->GetXaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetLabelSize(0.06);
        frame->GetXaxis()->SetTitleOffset(0.6);
        if(linlin) frame->GetXaxis()->SetLabelOffset(-0.02);
	else frame->GetXaxis()->SetLabelOffset(-0.04);
        frame->GetYaxis()->SetLabelSize(0.062);
        frame->GetYaxis()->SetLabelOffset(0.01);
        frame->GetYaxis()->SetTitleSize(0.07);
        frame->GetYaxis()->SetTitleOffset(1.43);
        TH1D *frame_PbPb = new TH1D("frame_PbPb","",1,-70,470);   //good if lin x
        frame_PbPb->GetYaxis()->SetRangeUser(yrangeframe[0],yrangeframe[1]);
        frame_PbPb->GetXaxis()->SetTickLength(0.02*0.8);
        frame_PbPb->GetYaxis()->SetTickLength(0.025/(1-sizeleft)*sizeleft);
        frame_PbPb->GetXaxis()->SetNdivisions(103/*,kFALSE*/);
        frame_PbPb->GetXaxis()->CenterTitle();
        frame_PbPb->GetXaxis()->SetTitle("#LT#it{N}_{part}#GT");
        frame_PbPb->GetXaxis()->SetTitleSize(0.2);
        frame_PbPb->GetXaxis()->SetLabelSize(0.19);
        frame_PbPb->GetXaxis()->SetTitleOffset(0.18);
        frame_PbPb->GetXaxis()->SetLabelOffset(-0.125);
      }
    }

  // introduce arrows at the PbPb limit (useful for reply to Referee#1)  
    if(plotmode==2){
      TArrow *arr_PbPblim_k0s = new TArrow(70,PbPblim_k0s,150,PbPblim_k0s,0.03,"|>");
      arr_PbPblim_k0s->SetLineWidth(3);
      TArrow *arr_PbPblim_lam = new TArrow(70,PbPblim_lam,150,PbPblim_lam,0.03,"|>");
      arr_PbPblim_lam->SetLineWidth(3);
      TArrow *arr_PbPblim_xi = new TArrow(70,PbPblim_xi,150,PbPblim_xi,0.03,"|>");
      arr_PbPblim_xi->SetLineWidth(3);
      TArrow *arr_PbPblim_om = new TArrow(70,PbPblim_om,150,PbPblim_om,0.03,"|>");
      arr_PbPblim_om->SetLineWidth(3);
    }
  
  //make conversion Nch --> Nch/R^3 (plotmode==5)
    if(plotmode==5){
      for(int i_part=lk0s; i_part<n_particles; i_part++){
        for(int i_system=0; i_system<n_systems; i_system++){
	  for(int i_htype=0; i_htype<3; i_htype++){
	    if(!presence[i_part][i_system][i_htype]) continue;
	    for(int i_point = 0; i_point<gr[i_part][i_system][i_htype]->GetN(); i_point++){
	      double x_gr, y_gr;
	      double e_x_gr, e_y_gr;
	      gr[i_part][i_system][i_htype]->GetPoint(i_point,x_gr,y_gr);
	      double radius = EvaluateRadius(i_system,x_gr);
	      if(i_part==llam && i_system==lpbpb){
	        e_x_gr = gr[i_part][i_system][i_htype]->GetErrorXlow(i_point);
	        e_y_gr = gr[i_part][i_system][i_htype]->GetErrorYlow(i_point);
	      }
	      else{
	        e_x_gr = gr[i_part][i_system][i_htype]->GetErrorX(i_point);
	        e_y_gr = gr[i_part][i_system][i_htype]->GetErrorY(i_point);
	      }
	      printf("i_part %d  i_system %d  i_htype %d  radius = %f\n",i_part,i_system,i_htype,radius);
	      gr[i_part][i_system][i_htype]->SetPoint(i_point,x_gr/TMath::Power(radius,2),y_gr);
	      /*if(i_part==llam && i_system==lpbpb) gr[i_part][i_system][i_htype]->SetPointError(i_point,e_x_gr/TMath::Power(radius,3),e_x_gr/TMath::Power(radius,3),e_y_gr,e_y_gr);
	      else*/ 
	      if(i_part!=llam/* && i_system!=lpbpb*/) gr[i_part][i_system][i_htype]->SetPointError(i_point,e_x_gr/TMath::Power(radius,2),e_y_gr);
	      else if(i_system!=lpbpb) gr[i_part][i_system][i_htype]->SetPointError(i_point,e_x_gr/TMath::Power(radius,2),e_y_gr);
	      else ((TGraphAsymmErrors*)gr[i_part][i_system][i_htype])->SetPointError(i_point,e_x_gr/TMath::Power(radius,3),e_x_gr/TMath::Power(radius,3),e_y_gr,e_y_gr);
	    }
	  }
        }
      }
    }
    
    
  
  //canvasing
    TCanvas *can = new TCanvas("can","can",1.5*400,1.5*650);
    can->cd();
    if(plotmode!=6) can->SetLogx();
    else{
      can->cd();
      TPad *pad_l = new TPad("pad_l","",0,0,sizeleft,1,0,0,0);
      pad_l->SetMargin(0.18/sizeleft,0,0.08/0.5475*0.55,0.015/0.5475*0.55);
      //pad_l->SetMargin(0.17,0,0.12/0.5475*0.55,0.015/0.5475*0.55);
      if(!linlin) pad_l->SetLogx();
      pad_l->SetLogy();
      pad_l->Draw();
      can->cd();
      TPad *pad_r = new TPad("pad_r","",sizeleft,0,1,1,0,0,0);
      pad_r->SetMargin(0.,0.02/(1-sizeleft),0.08/0.5475*0.55,0.015/0.5475*0.55);
      //pad_r->SetMargin(0.,0.02/(1-sizeleft),0.12/0.5475*0.55,0.015/0.5475*0.55);
      pad_r->SetLogy();
      pad_r->Draw();
//       can->cd();
//       padLeg_plotmode6 = new TPad("padLeg_plotmode6","",0.26,0.13,0.75,0.35,0,0,0);
//       padLeg_plotmode6->SetFillColorAlpha(kWhite,0.5);
//       padLeg_plotmode6->Draw();
//       padLeg_plotmode6->cd();
//       TBox *b_legs = new TBox(0.01,0.01,0.99,0.99);
//       b_legs->SetFillColorAlpha(kWhite,0.7);
//       b_legs->SetLineColor(kBlack);
//       b_legs->SetLineWidth(1);
//       b_legs->SetLineStyle(1);
//       b_legs->Draw();
      pad_l->cd();
    }
    if(plotmode!=3) can->SetLogy();
    frame->Draw();
    int nsyststodraw = n_systems; 
    if(plotmode==2 || plotmode==4 || plotmode==6) nsyststodraw--;
    if(plotmode==6 || plotmode==11){
      if(plotmode==6) {pad_r->cd();
      frame_PbPb->Draw();}
      gr[lk0s][lpbpb][lsyst]->Draw("PE5");
      gr[llam][lpbpb][lsyst]->Draw("PE5");
      gr[lxi][lpbpb][lsyst]->Draw("PE5");
      gr[lom][lpbpb][lsyst]->Draw("PE5");
      //gr_STAR[lk0s][lsyst]->Draw("PE5");
//       gr_STAR[llam][lsyst]->Draw("P");
//       gr_STAR[lxi][lsyst]->Draw("P");
//       gr_STAR[lom][lsyst]->Draw("P");
//       gr_PHEN[llam][lsyst]->Draw("P");
      //gr_PHEN[lxi][lsyst]->Draw("P");
      //gr_PHEN[lom][lsyst]->Draw("P");
    }
    for(int i_part=lk0s; i_part<n_particles; i_part++){
      //for(int i_system=0; i_system<nsyststodraw; i_system++){
      for(int i_system=nsyststodraw-1; i_system>-1; i_system--){
	for(int i_htype=0; i_htype<3; i_htype++){
	  if(!presence[i_part][i_system][i_htype]) continue;
	  //if(i_htype==lstat) continue;
	  if (i_system != lpbpb)
	    gr[i_part][i_system][i_htype]->Draw(drawop[i_htype]);
	  else
	    gr[i_part][i_system][i_htype]->Draw("PE5");
	}
      }
    }
    
    if(plotmode!=5){
      BuildLegend_MC(plotmode);
      BuildLegend_thPbPb();
      if(plotmode==6) pad_l->cd();
  //     PlotMC("K0S","Pythia6-Perugia2011-NoCR-Retuned","Pythia6-P2011-NoCR",grcol[lpp][lk0s],scalefactor[lk0s],9);
  //     PlotMC("K0S","Pythia6-Perugia2011-WithCR","Pythia6-P2011-WithCR",grcol[lpp][lk0s],scalefactor[lk0s],9);
      PlotMC("K0S","Pythia8-Monash-WithCR","PYTHIA8",grcol[lpp][lk0s]/*+1*/,scalefactor[lk0s],1);
      PlotMC("K0S","DIPSY-ColorRopes","DIPSY ",grcol[lpp][lk0s]/*+1*/,scalefactor[lk0s],2);
      PlotMC("K0S","EPOSLHC","EPOS LHC",grcol[lpp][lk0s]/*+1*/,scalefactor[lk0s],3);
      //    PlotMC("K0S","Phojet","PHOJET",grcol[lpp][lk0s]+1,scalefactor[lk0s],4);
      //    PlotMC("K0S","Herwig","HERWIG",grcol[lpp][lk0s]+1,scalefactor[lk0s],9);
          
  //     PlotMC("Lambda","Pythia6-Perugia2011-NoCR-Retuned","Pythia6-P2011-NoCR",grcol[lpp][llam],scalefactor[llam],9);
  //     PlotMC("Lambda","Pythia6-Perugia2011-WithCR","Pythia6-P2011-WithCR",grcol[lpp][llam],scalefactor[llam],9);
      PlotMC("Lambda","Pythia8-Monash-WithCR","PYTHIA8",grcol[lpp][llam]+1,scalefactor[llam],1);
      PlotMC("Lambda","DIPSY-ColorRopes","DIPSY",grcol[lpp][llam]+1,scalefactor[llam],2);
      PlotMC("Lambda","EPOSLHC","EPOS LHC",grcol[lpp][llam]+1,scalefactor[llam],3);
      //    PlotMC("Lambda","Phojet","PHOJET",grcol[lpp][llam]+1,scalefactor[llam],4);
      //    PlotMC("Lambda","Herwig","HERWIG",grcol[lpp][llam]+1,scalefactor[llam],9);
      
  //     PlotMC("Xi","Pythia6-Perugia2011-NoCR-Retuned","Pythia6-P2011-NoCR",grcol[lpp][lxi],scalefactor[lxi],1);
  //     PlotMC("Xi","Pythia6-Perugia2011-WithCR","Pythia6-P2011-WithCR",grcol[lpp][lxi],scalefactor[lxi],2);
      PlotMC("Xi","Pythia8-Monash-WithCR","PYTHIA8",grcol[lpp][lxi]+1,scalefactor[lxi],1);
      PlotMC("Xi","DIPSY-ColorRopes","DIPSY",grcol[lpp][lxi]+1,scalefactor[lxi],2);
      PlotMC("Xi","EPOSLHC","EPOS LHC",grcol[lpp][lxi]+1,scalefactor[lxi],3);
      //    PlotMC("Xi","Phojet","PHOJET",grcol[lpp][lxi]+1,scalefactor[lxi],4);
      //    PlotMC("Xi","Herwig","HERWIG",grcol[lpp][lxi]+1,scalefactor[lxi],9);
  
  //     PlotMC("Omega","Pythia6-Perugia2011-NoCR-Retuned","Pythia6-Perugia2011-NoCR",grcol[lpp][lom],scalefactor[lom],1);
  //     PlotMC("Omega","Pythia6-Perugia2011-WithCR","Pythia6-Perugia2011-WithCR",grcol[lpp][lom],scalefactor[lom],2);
      PlotMC("Omega","Pythia8-Monash-WithCR","PYTHIA8",grcol[lpp][lom],scalefactor[lom],1);
      PlotMC("Omega","DIPSY-ColorRopes","DIPSY",grcol[lpp][lom],scalefactor[lom],2);
      PlotMC("Omega","EPOSLHC","EPOS LHC",grcol[lpp][lom],scalefactor[lom],3);
      //    PlotMC("Omega","Phojet","PHOJET",grcol[lpp][lom],scalefactor[lom],4);
      //    PlotMC("Omega","Herwig","HERWIG",grcol[lpp][lom],scalefactor[lom],9);
    }
    
    DrawParticles(scalefactor,plotmode);
    
    DrawLegend_data(plotmode);
    if(plotmode!=5) leg_th->Draw();
  
    if(plotmode==2){
      arr_PbPblim_k0s->SetLineColor(grcol[lpp][lk0s]);
      arr_PbPblim_lam->SetLineColor(grcol[lpp][llam]);
      arr_PbPblim_xi->SetLineColor(grcol[lpp][lxi]);
      arr_PbPblim_om->SetLineColor(grcol[lpp][lom]);
      
      arr_PbPblim_k0s->SetFillColor(grcol[lpp][lk0s]);
      arr_PbPblim_lam->SetFillColor(grcol[lpp][llam]);
      arr_PbPblim_xi->SetFillColor(grcol[lpp][lxi]);
      arr_PbPblim_om->SetFillColor(grcol[lpp][lom]);
      
      arr_PbPblim_k0s->Draw();
      arr_PbPblim_lam->Draw();
      arr_PbPblim_xi->Draw();
      arr_PbPblim_om->Draw();
      
      TArrow *legarr = new TArrow(4,2.2e-3,7,2.2e-3,0.03,"|>");
      legarr->SetLineWidth(3);
      legarr->Draw();
      
      TLatex lat_legarr;
      lat_legarr.SetTextFont(textfont);
      lat_legarr.SetTextAlign(12);
      lat_legarr.SetTextSize(0.04);
      lat_legarr.DrawLatex(7.5,2.2e-3,"Pb-Pb fit 100 < d#it{N}_{ch}/d#it{#eta} < 2000");
    }
    
//    if(plotmode==1 || plotmode==11) can->SaveAs("img_ToPionRatios_iterative.pdf");
//     if(plotmode==1 || plotmode==11) can->SaveAs("ToPionRatios.eps");
    //else if(plotmode==11)can->SaveAs("img_ToPionRatios_asArXiv_newcolcode.pdf");
//     else if(plotmode==2) can->SaveAs("img_ToPionRatios_arrowPbPb.pdf");
//     else if(plotmode==4) can->SaveAs("img_ToPionRatios_pppPbonly.pdf");
//     else if(plotmode==5) can->SaveAs("img_ToPionRatios_scaletoR2.pdf");
//     else if(plotmode==6 && linlin) can->SaveAs("img_ToPionRatios_separPbPb_linlin.pdf");
//     else if(plotmode==6 && !linlin) can->SaveAs("img_ToPionRatios_separPbPb_loglin.pdf");
    
  
}

void PlotMC(TString CascType = "Xi", TString lModel = "PythiaP0", TString LegOptional = "opt", int color, float scale=1, int linestyle=1){
  
  TFile *fin_model = new TFile(Form("./McPredictions_%s.root",lModel.Data()),"read");
  if (!fin_model || !fin_model->IsOpen()) return;
  TGraphErrors *g_MC = 0x0;
  if(CascType.Contains("Xi")) g_MC = (TGraphErrors*) fin_model->Get("g_XiToPion");
  else if(CascType.Contains("Omega")) g_MC = (TGraphErrors*) fin_model->Get("g_OmToPion");
  else if(CascType.Contains("Lambda")) g_MC = (TGraphErrors*) fin_model->Get("g_LaToPion");
  else if(CascType.Contains("K0S")) g_MC = (TGraphErrors*) fin_model->Get("g_K0SToPion");
  g_MC->SetName("g_MC");

  double x_point_MC = 0;
  double x_error_MC = 0;
  double y_point_MC = 0;
  double y_error_MC = 0;
  for(int ipoint=0; ipoint<g_MC->GetN(); ipoint++){
    g_MC->GetPoint(ipoint,x_point_MC,y_point_MC);
    x_error_MC = g_MC->GetErrorX(ipoint);
    y_error_MC = g_MC->GetErrorY(ipoint);
    g_MC->SetPoint(ipoint,x_point_MC,y_point_MC*scale);
    g_MC->SetPointError(ipoint,x_error_MC,y_error_MC*scale);
  }

  g_MC->SetLineColor(color);
  g_MC->SetLineWidth(3);
//   g_MC->SetLineStyle(9);
  g_MC->SetLineStyle(linestyle);
  
  g_MC->Draw("Lsame");
  
  if(CascType.Contains("K0S")) leg_th->AddEntry(g_MC,Form("#scale[1]{%s}",LegOptional.Data()),"l");
//   if(CascType.Contains("Omega") && LegOptional.Contains("Dipsy")) {
//     //leg_th->AddEntry(g_MC,Form("#scale[0.95]{%s}","EPOS LHC"),"l");
//     leg_th->AddEntry(g_MC,Form("#scale[0.95]{%s}","HERWIG++"),"l");
//   }
  
}

void DrawLegend_data(int plotmode = 1){
  
   if(plotmode==1 || plotmode==11) {
//     TLegend *leg = new TLegend(0.22,0.78,0.57,0.96,NULL,"brNDC");
    //TLegend *leg = new TLegend(0.41,0.27,0.76,0.37,NULL,"brNDC");
    TLegend *leg = new TLegend(0.19,0.13,0.79,0.25,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(leg_textsize_data*0.8);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
//     leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[0.9]{#splitline{ALICE pp #sqrt{#it{s}} = 7 TeV}{V0M multiplicity classes}}","p");
//     leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[0.9]{#splitline{ALICE p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}{V0A multiplicity classes}}","p");
//     leg->AddEntry(gr[lom][lpbpb][lsyst],"#scale[0.9]{#splitline{ALICE Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV}{centrality classes}}","p");
    leg->AddEntry((TObject *)NULL, "#scale[1]{ALICE}", "");
    leg->AddEntry(gr[lk0s][lpp][lstat]  ,"#scale[1]{pp, #sqrt{#it{s}} = 7 TeV  Nat. Phys. 13 (2017) 535-539}","p");
    leg->AddEntry(gr[lk0s][lppb][lstat] ,"#scale[1]{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV  PLB 728 (2014) 25-38}","p");
    leg->AddEntry(gr[lk0s][lpbpb][lsyst],"#scale[1]{Preliminary Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV}","p");
    leg->Draw();
  }
  else if(plotmode==2){
    TLegend *leg = new TLegend(0.48,0.25,0.83,0.35,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(leg_textsize_data);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject *)NULL, "#scale[1]{ALICE}", "");
    leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[1]{pp, #sqrt{#it{s}} = 7 TeV}","p");
    leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[1]{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV}","p");
//     leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[0.9]{#splitline{ALICE pp #sqrt{#it{s}} = 7 TeV}{V0M multiplicity classes}}","pef");
//     leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[0.9]{#splitline{ALICE p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}{V0A multiplicity classes}}","pef");
    leg->Draw();
    
  }
  else if(plotmode==3){
    TLegend *leg = new TLegend(0.22,0.77,0.59,0.96,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(leg_textsize_data);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[1]{#splitline{ALICE pp #sqrt{#it{s}} = 7 TeV}{V0M multiplicity classes}}","pef");
    leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[1]{#splitline{ALICE p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}{V0A multiplicity classes}}","pef");
    leg->AddEntry(gr[lom][lpbpb][lsyst],"#scale[1]{#splitline{ALICE Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV}{centrality classes}}","pef");
    leg->Draw();
  }
  else if(plotmode==4){
    TLegend *leg = new TLegend(0.5,0.14,0.84,0.24,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(leg_textsize_data);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject *)NULL, "#scale[1]{ALICE}", "");
    leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[1]{pp, #sqrt{#it{s}} = 7 TeV}","p");
    leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[1]{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV}","p");
//     leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[0.9]{#splitline{ALICE pp #sqrt{#it{s}} = 7 TeV}{V0M multiplicity classes}}","pef");
//     leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[0.9]{#splitline{ALICE p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}{V0A multiplicity classes}}","pef");
    leg->Draw();
    
  }
   else if(plotmode==5) {
//     TLegend *leg = new TLegend(0.22,0.78,0.57,0.96,NULL,"brNDC");
    //TLegend *leg = new TLegend(0.41,0.27,0.76,0.37,NULL,"brNDC");
    TLegend *leg = new TLegend(0.46,0.14,0.86,0.27,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(leg_textsize_data);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
//     leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[0.9]{#splitline{ALICE pp #sqrt{#it{s}} = 7 TeV}{V0M multiplicity classes}}","p");
//     leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[0.9]{#splitline{ALICE p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}{V0A multiplicity classes}}","p");
//     leg->AddEntry(gr[lom][lpbpb][lsyst],"#scale[0.9]{#splitline{ALICE Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV}{centrality classes}}","p");
    leg->AddEntry((TObject *)NULL, "#scale[1]{ALICE}", "");
    leg->AddEntry(gr[lom][lpp][lstat]  ,"#scale[1]{pp, #sqrt{#it{s}} = 7 TeV}","p");
    leg->AddEntry(gr[lom][lppb][lstat] ,"#scale[1]{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV}","p");
    leg->AddEntry(gr[lom][lpbpb][lsyst],"#scale[1]{Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV}","p");
    leg->Draw();
  }
   else if(plotmode==6) {
    //TLegend *leg = new TLegend(0.41,0.14,0.81,0.26,NULL,"brNDC");
    if(linlin) TLegend *leg = new TLegend(0.41,0.23,0.81,0.35,NULL,"brNDC");
    else TLegend *leg = new TLegend(0.42,0.13,0.82,0.25,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(leg_textsize_data*1.2);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject *)NULL, "#scale[1]{ALICE}", "");
    leg->AddEntry(gr[lk0s][lpp][lstat]  ,"#scale[1]{pp, #sqrt{#it{s}} = 7 TeV}","p");
    leg->AddEntry(gr[lk0s][lppb][lstat] ,"#scale[1]{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV}","p");
    leg->AddEntry(gr[lk0s][lpbpb][lsyst],"#scale[1]{Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV}","p");
    leg->Draw();
  }
  
}

void BuildLegend_MC(int plotmode = 1){
   //if(plotmode==1 || plotmode==11) leg_th = new TLegend(0.45,0.16-0.02,0.85,0.25-0.02,NULL,"brNDC");
   if(plotmode==1 || plotmode==11) leg_th = new TLegend(0.60,0.25,0.99,0.34,NULL,"brNDC");
   else if(plotmode==2) leg_th = new TLegend(0.68,0.10,0.98,0.20,NULL,"brNDC");
   else if(plotmode==3) leg_th = new TLegend(0.5,0.20,0.90,0.30,NULL,"brNDC");
   else if(plotmode==4) leg_th = new TLegend(0.7,0.25,0.99,0.35,NULL,"brNDC");
   else if(plotmode==5) leg_th = new TLegend(0.45,0.16-0.02,0.85,0.25-0.02,NULL,"brNDC");
   //else if(plotmode==6) {leg_th = new TLegend(0.6,0.28,0.99,0.37,NULL,"brNDC");}
   else if(plotmode==6 & linlin) {leg_th = new TLegend(0.42,0.13,0.81,0.22,NULL,"brNDC");}
   else if(plotmode==6 & !linlin) {leg_th = new TLegend(0.61,0.24,0.99,0.33,NULL,"brNDC");}
   leg_th->SetBorderSize(0);
   if(plotmode!=6) leg_th->SetTextSize(leg_textsize_MC*0.95*0.8);
   else leg_th->SetTextSize(leg_textsize_MC*1.15*0.8);
   leg_th->SetLineColor(0);
   leg_th->SetLineStyle(1);
   leg_th->SetLineWidth(1);
   leg_th->SetFillColor(0);
   leg_th->SetFillStyle(0);
}

void BuildLegend_thPbPb(){
    leg_PbPb_th = new TLegend(0.64,0.26,0.93,0.63);
    leg_PbPb_th->SetHeader("#bf{Thermal models (Pb-Pb)}"); 
    leg_PbPb_th->SetTextSize(leg_textsize_Ther);
    //leg_PbPb_th->SetNColumns(2);
    leg_PbPb_th->SetBorderSize(0);
    leg_PbPb_th->SetFillColor(0);
    leg_PbPb_th->SetFillStyle(0);
    leg_PbPb_th->SetLineColor(0);
    leg_PbPb_th->SetLineStyle(1);
    leg_PbPb_th->SetLineWidth(1);
}

void ScaleGraph(int part, float nscale){
  
  if((nscale-1)<1e-5) return;
  
  double x_point = 0;
  double x_error = 0;
  double y_point = 0;
  double y_error = 0;
  
  for(int i_system=0; i_system<n_systems; i_system++){
    for  (int i_htype=0; i_htype<3; i_htype++){
      if(!presence[part][i_system][i_htype]) continue;
      for(int ipoint=0; ipoint<gr[part][i_system][i_htype]->GetN(); ipoint++){
        gr[part][i_system][i_htype]->GetPoint(ipoint,x_point,y_point);
        x_error = gr[part][i_system][i_htype]->GetErrorX(ipoint);
        y_error = gr[part][i_system][i_htype]->GetErrorY(ipoint);
	gr[part][i_system][i_htype]->SetPoint(ipoint,x_point,y_point*nscale);
        if(!(part==llam && i_system==lpbpb)) gr[part][i_system][i_htype]->SetPointError(ipoint,x_error,y_error*nscale);
      }
    }    
  }
  
}

void ScaleGraph_STAR(int part, float nscale){
  
  if((nscale-1)<1e-5) return;
  
  double x_point = 0;
  double x_error = 0;
  double y_point = 0;
  double y_error = 0;
  
  for(int ipoint=0; ipoint<gr_STAR[part][lstat]->GetN(); ipoint++){
    gr_STAR[part][lstat]->GetPoint(ipoint,x_point,y_point);
    x_error = gr_STAR[part][lstat]->GetErrorX(ipoint);
    y_error = gr_STAR[part][lstat]->GetErrorY(ipoint);
    gr_STAR[part][lstat]->SetPoint(ipoint,x_point,y_point*nscale);
    gr_STAR[part][lstat]->SetPointError(ipoint,x_error,y_error*nscale);
  }
  for(int ipoint=0; ipoint<gr_STAR[part][lsyst]->GetN(); ipoint++){
    gr_STAR[part][lsyst]->GetPoint(ipoint,x_point,y_point);
    x_error = gr_STAR[part][lsyst]->GetErrorX(ipoint);
    y_error = gr_STAR[part][lsyst]->GetErrorY(ipoint);
    gr_STAR[part][lsyst]->SetPoint(ipoint,x_point,y_point*nscale);
    gr_STAR[part][lsyst]->SetPointError(ipoint,x_error,y_error*nscale);
  }
  
  
}

void ScaleGraph_PHEN(int part, float nscale){
  
  if((nscale-1)<1e-5) return;
  
  double x_point = 0;
  double x_error = 0;
  double y_point = 0;
  double y_error = 0;
  
  for(int ipoint=0; ipoint<gr_PHEN[part][lstat]->GetN(); ipoint++){
    gr_PHEN[part][lstat]->GetPoint(ipoint,x_point,y_point);
    x_error = gr_PHEN[part][lstat]->GetErrorX(ipoint);
    y_error = gr_PHEN[part][lstat]->GetErrorY(ipoint);
    gr_PHEN[part][lstat]->SetPoint(ipoint,x_point,y_point*nscale);
    gr_PHEN[part][lstat]->SetPointError(ipoint,x_error,y_error*nscale);
  }
  for(int ipoint=0; ipoint<gr_PHEN[part][lsyst]->GetN(); ipoint++){
    gr_PHEN[part][lsyst]->GetPoint(ipoint,x_point,y_point);
    x_error = gr_PHEN[part][lsyst]->GetErrorX(ipoint);
    y_error = gr_PHEN[part][lsyst]->GetErrorY(ipoint);
    gr_PHEN[part][lsyst]->SetPoint(ipoint,x_point,y_point*nscale);
    gr_PHEN[part][lsyst]->SetPointError(ipoint,x_error,y_error*nscale);
  }
  
  
}

Double_t ErrRat ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ){
  //Error in a Ratio
        if(B!=0){
                Double_t errorfromtop = Aerr*Aerr / (B*B) ;
                Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
                return TMath::Sqrt( errorfromtop + errorfrombottom );
        }
        return 1;
}

Double_t ErrRatSyst ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ){
  //Error in a Ratio
        if(B!=0){
                Double_t errorfromA = Aerr/A*Aerr/A ;
                Double_t errorfromB = Berr/B*Berr/B;
                return TMath::Sqrt( errorfromA + errorfromB )*A/B;
        }
        return 1;
}

void DrawParticles(float *scalef, int plotmode=1){
  
  double y_lat[n_particles];
  double dum = 0;
  
  TLatex lat;
  //lat.SetNDC();
  lat.SetTextFont(textfont);
  //   lat.SetTextAlign(12);
  lat.SetTextAlign(22);
  lat.SetTextSize(0.04);
  if(plotmode==1 || plotmode==11){
    for(int i=lk0s+1; i<n_particles; i++){
      if(presence[i][lpbpb][lsyst]) gr[i][lpbpb][lsyst]->GetPoint(gr[i][lpbpb][lsyst]->GetN()-1,dum,y_lat[i]);
    }
    //     lat.DrawLatex(3000,y_lat[lom],Form("#Omega^{-}+#bar{#Omega}^{+} (#times%d)",(int)scalef[lom]));
    //     lat.DrawLatex(3000,y_lat[lxi],Form("#Xi^{-}+#bar{#Xi}^{+} (#times%d)",(int)scalef[lxi]));
    //     lat.DrawLatex(3000,y_lat[llam],Form("#Lambda+#bar{#Lambda} (#times%d)",(int)scalef[llam]));
    //     lat.DrawLatex(3000,1.6e-1,"2K^{0}_{S}");
    lat.DrawLatex(500,0.008,Form("#Omega^{#minus}+#bar{#Omega}^{#plus} (#times%d)",(int)scalef[lom]));
    lat.DrawLatex(500,0.9*y_lat[lxi],Form("#Xi^{#minus}+#bar{#Xi}^{#plus} (#times%d)",(int)scalef[lxi]));
    lat.DrawLatex(500,0.065,Form("#Lambda+#bar{#Lambda} (#times%d)",(int)scalef[llam]));
    lat.DrawLatex(500,0.12,"2K^{0}_{S}");
  }
  else if(plotmode==2){
    for(int i=lk0s; i<n_particles; i++){
      if(presence[i][lpbpb][lsyst]) gr[i][lppb][lsyst]->GetPoint(/*gr[i][lppb][lsyst]->GetN()-1*/0,dum,y_lat[i]);
    }
    lat.DrawLatex(50,y_lat[lom]*0.65,Form("#Omega^{#minus}+#bar{#Omega}^{#plus} (#times%d)",(int)scalef[lom]));
    lat.DrawLatex(50,y_lat[lxi]*0.7,Form("#Xi^{#minus}+#bar{#Xi}^{#plus} (#times%d)",(int)scalef[lxi]));
    lat.DrawLatex(50,y_lat[llam]*0.75,Form("#Lambda+#bar{#Lambda} (#times%d)",(int)scalef[llam]));
    lat.DrawLatex(50,y_lat[lk0s]*1,"2K^{0}_{S}");
  }
  if(plotmode==3){
    for(int i=lk0s+1; i<n_particles; i++){
      if(presence[i][lpbpb][lsyst]) gr[i][lppb][lsyst]->GetPoint(/*gr[i][lppb][lsyst]->GetN()-1*/0,dum,y_lat[i]);
      if(i==lom) gr[i][lpbpb][lsyst]->GetPoint(gr[i][lpbpb][lsyst]->GetN()-1,dum,y_lat[i]);
    }
    lat.DrawLatex(3000,y_lat[lom],Form("#frac{#Omega^{#minus}+#bar{#Omega}^{#plus}}{#pi^{#minus}+#pi^{#plus}} (#times%d)",(int)scalef[lom]));
    lat.DrawLatex(3000,y_lat[lxi],Form("#frac{#Xi^{#minus}+#bar{#Xi}^{#plus}}{#pi^{#minus}+#pi^{#plus}} (#times%d)",(int)scalef[lxi]));
    lat.DrawLatex(3000,y_lat[llam],Form("#frac{#Lambda+#bar{#Lambda}}{#pi^{#minus}+#pi^{#plus}} (#times%d)",(int)scalef[llam]));
    lat.DrawLatex(3000,1.4e-1,"#frac{2K^{0}_{S}}{#pi^{#minus}+#pi^{#plus}}");
  }
  else if(plotmode==4){
    for(int i=lk0s; i<n_particles; i++){
      if(presence[i][lpbpb][lsyst]) gr[i][lppb][lsyst]->GetPoint(/*gr[i][lppb][lsyst]->GetN()-1*/0,dum,y_lat[i]);
    }
    lat.DrawLatex(40,y_lat[lom]*0.65,Form("#Omega^{#minus}+#bar{#Omega}^{#plus} (#times%d)",(int)scalef[lom]));
    lat.DrawLatex(40,y_lat[lxi]*0.7,Form("#Xi^{#minus}+#bar{#Xi}^{#plus} (#times%d)",(int)scalef[lxi]));
    lat.DrawLatex(40,y_lat[llam]*0.75,Form("#Lambda+#bar{#Lambda} (#times%d)",(int)scalef[llam]));
    lat.DrawLatex(40,y_lat[lk0s]*0.8,"2K^{0}_{S}");
  }
  else if(plotmode==5){
    for(int i=lk0s+1; i<n_particles; i++){
      if(presence[i][lpbpb][lsyst]) gr[i][lpbpb][lsyst]->GetPoint(gr[i][lpbpb][lsyst]->GetN()-1,dum,y_lat[i]);
    }
    //     lat.DrawLatex(3000,y_lat[lom],Form("#Omega^{-}+#bar{#Omega}^{+} (#times%d)",(int)scalef[lom]));
    //     lat.DrawLatex(3000,y_lat[lxi],Form("#Xi^{-}+#bar{#Xi}^{+} (#times%d)",(int)scalef[lxi]));
    //     lat.DrawLatex(3000,y_lat[llam],Form("#Lambda+#bar{#Lambda} (#times%d)",(int)scalef[llam]));
    //     lat.DrawLatex(3000,1.6e-1,"2K^{0}_{S}");
    lat.DrawLatex(500,0.0095,Form("#Omega^{#minus}+#bar{#Omega}^{#plus} (#times%d)",(int)scalef[lom]));
    lat.DrawLatex(500,0.8*y_lat[lxi],Form("#Xi^{#minus}+#bar{#Xi}^{#plus} (#times%d)",(int)scalef[lxi]));
    lat.DrawLatex(500,0.055,Form("#Lambda+#bar{#Lambda} (#times%d)",(int)scalef[llam]));
    lat.DrawLatex(500,0.11,"2K^{0}_{S}");
  }
  else if(plotmode==6){
    for(int i=lk0s+1; i<n_particles; i++){
      if(presence[i][lpbpb][lsyst]) gr[i][lpbpb][lsyst]->GetPoint(gr[i][lpbpb][lsyst]->GetN()-1,dum,y_lat[i]);
    }
    lat.SetTextSize(leg_textsize_data*1.2);
    lat.DrawLatex(38,0.007,Form("#Omega^{#minus}+#bar{#Omega}^{#plus} (#times%d)",(int)scalef[lom]));
    lat.DrawLatex(38,0.8*y_lat[lxi],Form("#Xi^{#minus}+#bar{#Xi}^{#plus} (#times%d)",(int)scalef[lxi]));
    lat.DrawLatex(38,0.06,Form("#Lambda+#bar{#Lambda} (#times%d)",(int)scalef[llam]));
    lat.DrawLatex(38,0.112,"2K^{0}_{S}");
  }
  
}

double EvaluateRadius(int system, double Nch){
  
  TF1 *f_conversion = new TF1("f_conversion","[0]+[1]*x",0,200);
  if(system==lpp) f_conversion->SetParameters(0.332,0.405);
  else if (system==lppb) f_conversion->SetParameters(-0.054,0.585);
  else f_conversion->SetParameters(0.049,0.772);
  
  return f_conversion->Eval(TMath::Power(Nch,1./3));
  
}

void LoadStarGraphs(int plotmode = 6){
  
//
// http://arxiv.org/pdf/nucl-ex/0606014.pdf
//
// 0-5%, 10-20%, 20-40%, 40-60%, 60-80%
double npart[5]      = { 352., 235., 141.,  62.,  21. };
double nparte[5]     = {   0.,   0.,   0.,   0.,   0. };
double npartE[5]     = {   3.,   9.,   8.,   9.,   6. };
//
double lambda[5]     = { 16.7, 10.0, 5.53, 2.07, 0.58 };
double lambdae[5]    = {  0.2,  0.1, 0.05, 0.03, 0.01 };
double lambdaE[5]    = {  1.1,  0.7, 0.39, 0.14, 0.04 };
//
double lambdabar[5]  = { 12.7,  7.7, 4.30, 1.64, 0.48 };
double lambdabare[5] = {  0.2,  0.1, 0.04, 0.03, 0.01 };
double lambdabarE[5] = {  0.9,  0.5, 0.30, 0.11, 0.03 };
//
double xi[5]         = { 2.17, 1.41, 0.72, 0.26, 0.063 };
double xie[5]        = { 0.06, 0.04, 0.02, 0.01, 0.004 };
double xiE[5]        = { 0.19, 0.08, 0.02, 0.02, 0.003 };
//
double xibar[5]      = { 1.83, 1.14, 0.62, 0.23, 0.061 };
double xibare[5]     = { 0.05, 0.04, 0.02, 0.01, 0.004 };
double xibarE[5]     = { 0.20, 0.08, 0.03, 0.02, 0.002 };
//
double omegasum[5]   = { 0.53, 0.  , 0.17, 0.063, 0.   };
double omegasume[5]  = { 0.04, 0.  , 0.02, 0.008, 0.   };
double omegasumE[5]  = { 0.04, 0.  , 0.01, 0.004, 0.   };
//
// http://journals.aps.org/prc/pdf/10.1103/PhysRevC.79.034909
//
double dndeta[5]     = { 691., 421., 241.,  102,  33.5 };
double dndetae[5]    = {   0.,   0.,   0.,    0,   0   };
double dndetaE[5]    = {  49.,  30.,  17.,  7.5,   2.5 };
//
double piminus[5]    = { 327., 196., 112.8, 47.6, 16. };
double piminusE[5]   = {  25.,  15.,   8.4, 3.65, 1.2 };
//
double piplus[5]     = { 322., 194., 112.1, 47.45, 15.95 };
double piplusE[5]    = {  25.,  15.,   8.4,   3.6,  1.2 };
//

double lamrat[5];
double xirat[5];
double omrat[5];
double estat_lamrat[5];
double estat_xirat[5];
double estat_omrat[5];
double esyst_lamrat[5];
double esyst_xirat[5];
double esyst_omrat[5];

const int numpoints = 5;

for(int i=0; i<numpoints; i++){
  lamrat[i] = (lambda[i]+lambdabar[i]) / (piminus[i]+piplus[i]);
  xirat[i] = (xi[i]+xibar[i]) / (piminus[i]+piplus[i]);
  if(i!=1 && i!=4) omrat[i] = (omegasum[i]) / (piminus[i]+piplus[i]);
  estat_lamrat[i] = ErrRat(lambda[i]+lambdabar[i],TMath::Sqrt(lambdae[i]*lambdae[i]+lambdabare[i]*lambdabare[i]),piminus[i]+piplus[i],0);
  estat_xirat[i]  = ErrRat(xi[i]+xibar[i],TMath::Sqrt(xie[i]*xie[i]+xibare[i]*xibare[i]),piminus[i]+piplus[i],0);
  if(i!=1 && i!=4) estat_omrat[i]  = ErrRat(omegasum[i],omegasume[i],piminus[i]+piplus[i],0);
  esyst_lamrat[i] = ErrRat(lambda[i]+lambdabar[i],(lambdaE[i]+lambdabarE[i])/2,piminus[i]+piplus[i],(piminusE[i]+piplusE[i])/2);
  esyst_xirat[i]  = ErrRat(xi[i]+xibar[i],(xie[i]+xibare[i])/2,piminus[i]+piplus[i],(piminusE[i]+piplusE[i])/2);
  if(i!=1 && i!=4) esyst_omrat[i]  = ErrRat(omegasum[i],omegasume[i],piminus[i]+piplus[i],(piminusE[i]+piplusE[i])/2);
  esyst_lamrat[i] = TMath::Sqrt(esyst_lamrat[i]*esyst_lamrat[i]+estat_lamrat[i]*estat_lamrat[i]);
  esyst_xirat[i] = TMath::Sqrt(esyst_xirat[i]*esyst_xirat[i]+estat_xirat[i]*estat_xirat[i]);
  esyst_omrat[i] = TMath::Sqrt(esyst_omrat[i]*esyst_omrat[i]+estat_omrat[i]*estat_omrat[i]);
}

  if(plotmode==6){
    gr_STAR[llam][lstat] = new TGraphErrors(numpoints,npart,lamrat,nparte,estat_lamrat);
    gr_STAR[lxi][lstat] = new TGraphErrors(numpoints,npart,xirat,nparte,estat_xirat);
    gr_STAR[lom][lstat] = new TGraphErrors(numpoints,npart,omrat,nparte,estat_omrat);
    
    gr_STAR[llam][lsyst] = new TGraphErrors(numpoints,npart,lamrat,npartE,esyst_lamrat);
    gr_STAR[lxi][lsyst] = new TGraphErrors(numpoints,npart,xirat,npartE,esyst_xirat);
    gr_STAR[lom][lsyst] = new TGraphErrors(numpoints,npart,omrat,npartE,esyst_omrat);
  }
  else if(plotmode==11){
    gr_STAR[llam][lstat] = new TGraphErrors(numpoints,dndeta,lamrat,dndetae,estat_lamrat);
    gr_STAR[lxi][lstat] = new TGraphErrors(numpoints,dndeta,xirat,dndetae,estat_xirat);
    gr_STAR[lom][lstat] = new TGraphErrors(numpoints,dndeta,omrat,dndetae,estat_omrat);
    
    gr_STAR[llam][lsyst] = new TGraphErrors(numpoints,dndeta,lamrat,dndetaE,esyst_lamrat);
    gr_STAR[lxi][lsyst] = new TGraphErrors(numpoints,dndeta,xirat,dndetaE,esyst_xirat);
    gr_STAR[lom][lsyst] = new TGraphErrors(numpoints,dndeta,omrat,dndetaE,esyst_omrat);
    
  }
  
//   TCanvas *can = new TCanvas("can","can");
//   can->cd();
//   TH1D *frame_PbPb = new TH1D("frame_PbPb","",1,-70,470);   //good if lin x
//   frame_PbPb->GetYaxis()->SetRangeUser(0.9e-3,1.9e-1);
// //   frame_PbPb->GetXaxis()->SetTickLength(0.02*0.8);
// //   frame_PbPb->GetYaxis()->SetTickLength(0.025/(1-sizeleft)*sizeleft);
//   frame_PbPb->GetXaxis()->SetNdivisions(103/*,kFALSE*/);
//   frame_PbPb->GetXaxis()->CenterTitle();
//   frame_PbPb->GetXaxis()->SetTitle("#LT#it{N}_{part}#GT");
//   frame_PbPb->GetXaxis()->SetTitleSize(0.2);
//   frame_PbPb->GetXaxis()->SetLabelSize(0.19);
//   frame_PbPb->GetXaxis()->SetTitleOffset(0.18);
//   frame_PbPb->GetXaxis()->SetLabelOffset(-0.125);
//   frame_PbPb->Draw();
//   gr_STAR[llam][lstat]->Draw("P");
//   gr_STAR[lxi][lstat]->Draw("P");
//   //gr_STAR[lom][lstat]->Draw("P");
//   gr_STAR[llam][lsyst]->Draw("E5");
//   gr_STAR[lxi][lsyst]->Draw("E5");
//   //gr_STAR[lom][lsyst]->Draw("E5");
  
  
}

void LoadPhenixGraphs(int plotmode = 6){
  
//
// http://arxiv.org/pdf/nucl-ex/0606014.pdf
//
// 0-5%, 10-20%, 20-40%, 40-60%, 60-80%
double npart[5]      = { 348., 271., 180.,  79.,  14. };
double nparte[5]     = {   0.,   0.,   0.,   0.,   0. };
double npartE[5]     = {  10.,  8.4,  6.6,  4.6,  3.3 };
//
double lambda[5]     = { 17.9, 13.0, 10.1, 5.9, 1.61 };
double lambdae[5]    = {  0.4,  0.3,  0.2, 0.2, 0.05 };
double lambdaE[5]    = {  1.8,  1.3, 1.01,0.59, 0.16 };
//
double lambdabar[5]  = { 12.0,  9.6, 7.40, 4.60, 1.26 };
double lambdabare[5] = {  0.3,  0.3, 0.20, 0.10, 0.04 };
double lambdabarE[5] = {  1.2, 0.96, 0.74, 0.46, 0.13 };
//
// double xi[5]         = { 2.17, 1.41, 0.72, 0.26, 0.063 };
// double xie[5]        = { 0.06, 0.04, 0.02, 0.01, 0.004 };
// double xiE[5]        = { 0.19, 0.08, 0.02, 0.02, 0.003 };
// //
// double xibar[5]      = { 1.83, 1.14, 0.62, 0.23, 0.061 };
// double xibare[5]     = { 0.05, 0.04, 0.02, 0.01, 0.004 };
// double xibarE[5]     = { 0.20, 0.08, 0.03, 0.02, 0.002 };
// //
// double omegasum[5]   = { 0.53, 0.  , 0.17, 0.063, 0.   };
// double omegasume[5]  = { 0.04, 0.  , 0.02, 0.008, 0.   };
// double omegasumE[5]  = { 0.04, 0.  , 0.01, 0.004, 0.   };
//
// http://journals.aps.org/prc/pdf/10.1103/PhysRevC.79.034909
//
double dndeta[5]     = { 691., 421., 241.,  102,  33.5 };
double dndetae[5]    = {   0.,   0.,   0.,    0,   0   };
double dndetaE[5]    = {  49.,  30.,  17.,  7.5,   2.5 };
//
double piminus[5]    = { 270., 200., 129.0, 53.3, 8.6 };
double piminusE[5]   = {  0.13*270,  0.13*200., 0.13*129.0, 0.13*53.3, 0.13*8.6 };
//
double piplus[5]     = { 276., 216., 141, 57, 9.6 };
double piplusE[5]    = { 0.13*276., 0.13*216., 0.13*141, 0.13*57, 0.13*9.6 };
//

double lamrat[5];
double xirat[5];
double omrat[5];
double estat_lamrat[5];
double estat_xirat[5];
double estat_omrat[5];
double esyst_lamrat[5];
double esyst_xirat[5];
double esyst_omrat[5];

const int numpoints = 5;

for(int i=0; i<numpoints; i++){
  lamrat[i] = (lambda[i]+lambdabar[i]) / (piminus[i]+piplus[i]);
//   xirat[i] = (xi[i]+xibar[i]) / (piminus[i]+piplus[i]);
//   if(i!=1 && i!=4) omrat[i] = (omegasum[i]) / (piminus[i]+piplus[i]);
  estat_lamrat[i] = ErrRat(lambda[i]+lambdabar[i],TMath::Sqrt(lambdae[i]*lambdae[i]+lambdabare[i]*lambdabare[i]),piminus[i]+piplus[i],0);
//   estat_xirat[i]  = ErrRat(xi[i]+xibar[i],TMath::Sqrt(xie[i]*xie[i]+xibare[i]*xibare[i]),piminus[i]+piplus[i],0);
//   if(i!=1 && i!=4) estat_omrat[i]  = ErrRat(omegasum[i],omegasume[i],piminus[i]+piplus[i],0);
  esyst_lamrat[i] = ErrRat(lambda[i]+lambdabar[i],(lambdaE[i]+lambdabarE[i])/2,piminus[i]+piplus[i],(piminusE[i]+piplusE[i])/2);
//   esyst_xirat[i]  = ErrRat(xi[i]+xibar[i],(xie[i]+xibare[i])/2,piminus[i]+piplus[i],(piminusE[i]+piplusE[i])/2);
//   if(i!=1 && i!=4) esyst_omrat[i]  = ErrRat(omegasum[i],omegasume[i],piminus[i]+piplus[i],(piminusE[i]+piplusE[i])/2);
  esyst_lamrat[i] = TMath::Sqrt(esyst_lamrat[i]*esyst_lamrat[i]+estat_lamrat[i]*estat_lamrat[i]);
//   esyst_xirat[i] = TMath::Sqrt(esyst_xirat[i]*esyst_xirat[i]+estat_xirat[i]*estat_xirat[i]);
//   esyst_omrat[i] = TMath::Sqrt(esyst_omrat[i]*esyst_omrat[i]+estat_omrat[i]*estat_omrat[i]);
}

  if(plotmode==6){
    gr_PHEN[llam][lstat] = new TGraphErrors(1,npart,lamrat,nparte,estat_lamrat);
//     gr_PHEN[lxi][lstat] = new TGraphErrors(numpoints,npart,xirat,nparte,estat_xirat);
//     gr_PHEN[lom][lstat] = new TGraphErrors(numpoints,npart,omrat,nparte,estat_omrat);
    
    gr_PHEN[llam][lsyst] = new TGraphErrors(1,npart,lamrat,npartE,esyst_lamrat);
//     gr_PHEN[lxi][lsyst] = new TGraphErrors(numpoints,npart,xirat,npartE,esyst_xirat);
//     gr_PHEN[lom][lsyst] = new TGraphErrors(numpoints,npart,omrat,npartE,esyst_omrat);
  }
  else if(plotmode==11){
    gr_PHEN[llam][lstat] = new TGraphErrors(1,dndeta,lamrat,dndetae,estat_lamrat);
//     gr_PHEN[lxi][lstat] = new TGraphErrors(numpoints,dndeta,xirat,dndetae,estat_xirat);
//     gr_PHEN[lom][lstat] = new TGraphErrors(numpoints,dndeta,omrat,dndetae,estat_omrat);
    
    gr_PHEN[llam][lsyst] = new TGraphErrors(1,dndeta,lamrat,dndetaE,esyst_lamrat);
//     gr_PHEN[lxi][lsyst] = new TGraphErrors(numpoints,dndeta,xirat,dndetaE,esyst_xirat);
//     gr_PHEN[lom][lsyst] = new TGraphErrors(numpoints,dndeta,omrat,dndetaE,esyst_omrat);
    
  }
  
//   TCanvas *can = new TCanvas("can","can");
//   can->cd();
//   TH1D *frame_PbPb = new TH1D("frame_PbPb","",1,-70,470);   //good if lin x
//   frame_PbPb->GetYaxis()->SetRangeUser(0.9e-3,1.9e-1);
// //   frame_PbPb->GetXaxis()->SetTickLength(0.02*0.8);
// //   frame_PbPb->GetYaxis()->SetTickLength(0.025/(1-sizeleft)*sizeleft);
//   frame_PbPb->GetXaxis()->SetNdivisions(103/*,kFALSE*/);
//   frame_PbPb->GetXaxis()->CenterTitle();
//   frame_PbPb->GetXaxis()->SetTitle("#LT#it{N}_{part}#GT");
//   frame_PbPb->GetXaxis()->SetTitleSize(0.2);
//   frame_PbPb->GetXaxis()->SetLabelSize(0.19);
//   frame_PbPb->GetXaxis()->SetTitleOffset(0.18);
//   frame_PbPb->GetXaxis()->SetLabelOffset(-0.125);
//   frame_PbPb->Draw();
//   gr_STAR[llam][lstat]->Draw("P");
//   gr_STAR[lxi][lstat]->Draw("P");
//   //gr_STAR[lom][lstat]->Draw("P");
//   gr_STAR[llam][lsyst]->Draw("E5");
//   gr_STAR[lxi][lsyst]->Draw("E5");
//   //gr_STAR[lom][lsyst]->Draw("E5");
  
  
}