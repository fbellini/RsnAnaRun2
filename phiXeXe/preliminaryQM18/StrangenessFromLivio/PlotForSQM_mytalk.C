const int n_exp = 2;
const int n_particles = 8;
const int n_systems = 15;
const int n_htypes = 4; //3 types of histo: stat, tot syst e syst_unc
TGraphErrors *gr[n_exp][n_particles][n_systems][n_htypes];

enum{lALICE, lSTAR};
enum{lpi, lk, lp, lk0s, llam, lphi, lxi, lom};
enum{lpp7, lppb5, lpbpb276, lpp13, lpbpb5, lppb8, lauau02, lpp02, lauau624, lauau77, lauau115, lauau196, lauau27, lauau39, lcucu02};
enum{lstat, lsyst, lsyst_unc, ltot};
bool presence[n_exp][n_particles][n_systems][n_htypes];
char expname[n_exp][10] = {"ALICE","STAR"};
char partname[n_particles][10] = {"Pion","Kaon","Proton","K0S","Lambda","Phi","Xi","Omega"};
char leg_partname[n_particles][10] = {"#pi","K","p","K0S","#Lambda","#phi","#Xi","#Omega"/*,"#phi"*/};
char systname[n_systems][10] = {"pp_7","pPb_5","PbPb_2.76","pp_13","PbPb_5","pPb_8","AuAu_0.2","pp_0.2","AuAu_62.4","AuAu_7.7","AuAu_11.5","AuAu_19.6","AuAu_27.0","AuAu_39.0","CuCu_0.2"};
char systen_exp[n_systems][50] = {"pp 7 TeV","p-Pb 5.02 TeV","Pb-Pb 2.76 TeV","pp 13 TeV","PbPb 5 TeV","pPb 8 TeV","200 GeV","pp 200 GeV","62.4 GeV","7.7 GeV","11.5 GeV","19.6 GeV","27 GeV","39 GeV"};



void DefinePresence(){
  
  for(int i_exp=lALICE; i_exp<n_exp; i_exp++)
    for(int i_part=lpi; i_part<n_particles; i_part++)
      for(int i_system=0; i_system<n_systems; i_system++)
        for(int i_htype=0; i_htype<n_htypes; i_htype++) presence[i_exp][i_part][i_system][i_htype] = kFALSE;
  
  //ad-hoc selection
    
//     presence [lALICE] [lk0s] [lpp7] [ltot] 	= kTRUE;
    presence [lALICE] [llam] [lpp7] [ltot] 	= kTRUE;
    presence [lALICE] [lxi]  [lpp7] [ltot] 	= kTRUE;
    presence [lALICE] [lom]  [lpp7] [ltot] 	= kTRUE;
    
//     presence [lALICE] [lk0s] [lppb5] [ltot] 	= kTRUE;
    presence [lALICE] [llam] [lppb5] [ltot] 	= kTRUE;
    presence [lALICE] [lxi]  [lppb5] [ltot] 	= kTRUE;
    presence [lALICE] [lom]  [lppb5] [ltot] 	= kTRUE;
    
//     presence [lALICE] [lk0s] [lpbpb276] [ltot] 	= kTRUE;
//     presence [lALICE] [llam] [lpbpb276] [ltot] 	= kTRUE;
//     presence [lALICE] [lxi]  [lpbpb276] [ltot] 	= kTRUE;
//     presence [lALICE] [lom]  [lpbpb276] [ltot] 	= kTRUE;
    
//     presence [lSTAR] [lk0s] [lauau02] [ltot] 	= kTRUE;
//     presence [lSTAR] [lp] [lauau02] [ltot] 	= kTRUE;
//     presence [lSTAR] [llam] [lauau02] [ltot] 	= kTRUE;
//     presence [lSTAR] [lxi]  [lauau02] [ltot] 	= kTRUE;
//     presence [lSTAR] [lom]  [lauau02] [ltot] 	= kTRUE;
    
//     presence [lSTAR] [llam] [lcucu02] [ltot] 	= kTRUE;
//     presence [lSTAR] [lxi]  [lcucu02] [ltot] 	= kTRUE;
//     presence [lSTAR] [lom]  [lcucu02] [ltot] 	= kTRUE;
    
//     presence [lALICE] [lk0s] [lpbpb5] [ltot] 	= kTRUE;
    presence [lALICE] [llam] [lpbpb5] [ltot] 	= kTRUE;
    presence [lALICE] [lxi]  [lpbpb5] [ltot] 	= kTRUE;
    presence [lALICE] [lom]  [lpbpb5] [ltot] 	= kTRUE;

//     presence [lSTAR] [lp] [lpp02] [ltot] 	= kTRUE;
//     presence [lSTAR] [llam] [lpp02] [ltot] 	= kTRUE;
//     presence [lSTAR] [lxi]  [lpp02] [ltot] 	= kTRUE;
//     presence [lSTAR] [lom]  [lpp02] [ltot] 	= kTRUE;
    
//     presence [lALICE] [lp] [lpp7] [ltot] 	= kTRUE;
//     presence [lALICE] [lp] [lpbpb5] [ltot] 	= kTRUE;
//     presence [lALICE] [lp] [lppb5] [ltot] 	= kTRUE;

  
}


void PlotForSQM_mytalk(bool bands=kFALSE, int sumplusminus=0, bool addpp = kFALSE){
  
  //style
    gROOT->LoadMacro("LoadPlotsStyle.C");
    LoadPlotsStyle(2);
  
  //build graphs from txt file
    DefinePresence();
    gROOT->LoadMacro("GenerateGraphFromTxt.C");
    float scalefactor[n_particles] = {1,1,2,1,1.,0.5,2,4};
    //float scalefactor[n_particles] = {0,0,0,1,1.5,0,6,16};
    //float scalefactor[n_particles] = {0,0,0,1,1,1,1,1};
    for(int i_exp=lALICE; i_exp<n_exp; i_exp++){
      for(int i_part=lpi; i_part<n_particles; i_part++){
        for(int i_system=0; i_system<n_systems; i_system++){
          for(int i_htype=0; i_htype<n_htypes; i_htype++){
            if(presence[i_exp][i_part][i_system][i_htype]) {
	      gr[i_exp][i_part][i_system][i_htype] = (TGraphErrors*) GenerateGraphFromTxt(expname[i_exp],systname[i_system],partname[i_part],kFALSE,sumplusminus,addpp);
	      ScaleGraph(i_exp,i_part,i_system,i_htype,scalefactor[i_part]);
            }
          }
        }
      }
    }
  
  
  
      float gralpha[n_systems] = {0.5, 1., 1.};
	  //generate color gradients
	  const int NRGBs = 5;
	  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	  Double_t red[NRGBs]   = { 0.00, 0.00, 0.9*0.87, 1.00, 0.51 };
	  Double_t green[NRGBs] = { 0.00, 0.81, 0.9*1.00, 0.20, 0.00 };
	  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	  int steps = 7; //number of multiplicity bins for that system and those fitted particles
	  int FI = TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,steps+1);
	
//       //			               lpi, lk, lp,  		lk0s,		llam, 		lphi,	lxi,   		lom
//       int grcol[n_systems][n_particles] =    { { 0,  0,  kGray+3,	kGray+3,	kGray+3,	0, 	kGray+3,	kGray+3},  //pp7
//                                                { 0,  0,  kGray+2, 	kGray+2,  	kGray+2,	0, 	kGray+2, 	kGray+2},  //pPb5
// 					       { 0,  0,  kGray+1, 	kGray+1,   	kGray+1, 	0, 	kGray+1, 	kGray+1},  //PbPb276
// 					       { 0,  0,  kBlack, 	kBlack,		kBlue,		0, 	kGreen+2, 	kRed+1},  //pp13
// 					       { 0,  0,  kGray+1, 	kGray+1,	kGray+1,	0, 	kGray+1, 	kGray+1},  //PbPb5
// 					       { 0,  0,  0, 	kBlack,		kBlue,		0, 	kGreen+2, 	kRed+1},  //pPb8
// 					       { 0,  0,  FI+6, 	FI+6,		FI+6,		0, 	FI+6, 		FI+6},  //AuAu02
// 					       { 0,  0,  kRed, 	kRed,		kRed,		0, 	kRed, 		kRed},  //pp02
// 					       { 0,  0,  0, 	FI+5,		FI+5,		0, 	FI+5, 		FI+5},  //AuAu624
// 					       { 0,  0,  0, 	FI,		FI,		0, 	FI, 		FI},  //AuAu77
// 					       { 0,  0,  0, 	FI+1,		FI+1,		0, 	FI+1, 		FI+1},  //AuAu115
// 					       { 0,  0,  0, 	FI+2,		FI+2,		0, 	FI+2, 		FI+2},  //AuAu196
// 					       { 0,  0,  0, 	FI+3,		FI+3,		0, 	FI+3, 		FI+3},  //AuAu27
// 					       { 0,  0,  0, 	FI+4,		FI+4,		0, 	FI+4, 		FI+4}   //AuAu39
// 					   };
      //			               lpi, lk, lp,   lk0s,    llam, lphi,      lxi,    lom
      int grcol[n_systems][n_particles] =    { { 0,  0,  0, kBlack,   kBlack,    0, 	kBlack, 	kBlack},  //pp7
                                               { 0,  0,  0, kBlack,   kBlack,    0, 	kBlack, 	kBlack},  //pPb5
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //PbPb276
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //pp13
					       { 0,  0,  0, kBlack,   kBlack,    0, 	kBlack, 	kBlack},  //PbPb5
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //pPb8
					       { 0,  0,  kRed, 	kRed,		kRed,		0, 	kRed, 		kRed},  //AuAu02
					       { 0,  0,  kRed, 	kRed,		kRed,		0, 	kRed, 		kRed},  //pp02
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //AuAu624
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //AuAu77
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //AuAu115
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //AuAu196
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},  //AuAu27
					       { 0,  0,  0, kBlack,   kBlue,    0, kGreen+2, kRed+1},   //AuAu39
					       { 0,  0,  0, kGreen+2,   kGreen+2,    0, kGreen+2, kGreen+2}   //CuCu200
					   };
      int grmarker[n_systems][n_particles] = { { 0,  0,  20,     20,      20,    0,       20,     20},  //pp7
                                               { 0,  0,  27,     27,      27,    0,       27,     27},  //pPb5
					       { 0,  0,  25,     25,      25,    0,       25,     25},  //PbPb276
					       { 0,  0,  25,     25,      25,    0,       25,     25},  //pp13
					       { 0,  0,  25,     25,      25,    0,       25,     25},  //PbPb5
					       { 0,  0,  25,     25,      25,    0,       25,     25},  //pPb8
					       { 0,  0,  30,     30,      30,    0,       30,     30},  //AuAu02
					       { 0,  0,  29,     29,      29,    0,       29,     29},  //pp02
					       { 0,  0,  28,     28,      28,    0,       28,     28},  //AuAu624
					       { 0,  0,  21,     21,      21,    0,       21,     21},  //AuAu77
					       { 0,  0,  0,     22,      22,    0,       22,     22},  //AuAu115
					       { 0,  0,  0,     23,      23,    0,       23,     23},  //AuAu196
					       { 0,  0,  0,     26,      26,    0,       26,     26},  //AuAu27
					       { 0,  0,  0,     32,      32,    0,       32,     32},   //AuAu39
					       { 0,  0,  30,     30,      30,    0,       30,     30}   //CuCu200
					   };
      float grmarsiz[n_systems][n_particles]={ { 0,  0,  2.,     2.,      2.,    0,        2,      2},  //pp7
					       { 0,  0,  2.5,    2.5,     2.5,    0,      2.5,    2.5},  //pPb5
					       { 0,  0,  1.5,    1.5,     1.5,    0,      1.5,    1.5},  //PbPb276
					       { 0,  0,  2.,     2.,      2.,    0,        2,      2},  //pp13
					       { 0,  0,  2.,     2.,      2.,    0,        2,      2},  //PbPb5
					       { 0,  0,  2.,     2.,      2.,    0,        2,      2},  //pPb8
					       { 0,  0,  2.,     2.,      2.,    0,        2,      2},  //AuAu02
					       { 0,  0,  2.,     2.,      2.,    0,        2,      2},  //pp02
					       { 0,  0,  0,     2.,      2.,    0,        2,      2},  //AuAu624
					       { 0,  0,  0,     2.,      2.,    0,        2,      2},  //AuAu77
					       { 0,  0,  0,     2.,      2.,    0,        2,      2},  //AuAu115
					       { 0,  0,  0,     2.,      2.,    0,        2,      2},  //AuAu196
					       { 0,  0,  0,     2.,      2.,    0,        2,      2},  //AuAu27
					       { 0,  0,  0,     2.,      2.,    0,        2,      2},   //AuAu39
					       { 0,  0,  0,     2.,      2.,    0,        2,      2}   //CuCu200
					   };
   
    for(int i_exp=lALICE; i_exp<n_exp; i_exp++){
      for(int i_part=lp; i_part<n_particles; i_part++){
        for(int i_system=0; i_system<n_systems; i_system++){
  	  for(int i_htype=0; i_htype<n_htypes; i_htype++){
  	    if(!presence[i_exp][i_part][i_system][i_htype]) continue;
  	    gr[i_exp][i_part][i_system][i_htype]->SetLineColor(grcol[i_system][i_part]);
  	    gr[i_exp][i_part][i_system][i_htype]->SetMarkerColor(grcol[i_system][i_part]);
  	    //gr[i_exp][i_part][i_system][i_htype]->SetMarkerColorAlpha(grcol[i_system][i_part], gralpha[i_system]);
  	    gr[i_exp][i_part][i_system][i_htype]->SetMarkerStyle(grmarker[i_system][i_part]);
  	    gr[i_exp][i_part][i_system][i_htype]->SetMarkerSize(grmarsiz[i_system][i_part]);
  	    if(!bands) gr[i_exp][i_part][i_system][i_htype]->SetFillStyle(0);
  	    gr[i_exp][i_part][i_system][i_htype]->SetFillStyle(1001);
  	    //gr[i_exp][i_part][i_system][i_htype]->SetFillStyle(0);
  	    gr[i_exp][i_part][i_system][i_htype]->SetFillColorAlpha(grcol[i_system][i_part],0.1);
  	  }
        }
      }
    }



    int nsyststodraw = n_systems; 
    char drawop[n_htypes][10] = {"P","E5","E5","P"};
    
    double xrangeframe[2] = {1.5,7000};
    double yrangeframe[2] = {1.8e-4,0.07};
    
    TH1D *frame = new TH1D("frame","",1,xrangeframe[0],xrangeframe[1]);
    frame->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|< 0.5}");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetTitleOffset(1.);
    frame->GetXaxis()->SetNdivisions(206);
    frame->GetYaxis()->SetRangeUser(yrangeframe[0],yrangeframe[1]);
    frame->GetYaxis()->SetTitle("ratio of yields to (#pi^{#plus}+#pi^{#minus})");
    //frame->GetYaxis()->CenterTitle();
    //frame->GetXaxis()->CenterTitle();
    
    TCanvas *can = new TCanvas("can","can",700,950);
    can->cd(0)->SetLogx();
    can->cd(0)->SetLogy();
    can->SetMargin(0.16,0.01,0.12,0.01);
    frame->DrawClone();
    for(int i_exp=lALICE; i_exp<n_exp; i_exp++){
      for(int i_part=lp; i_part<n_particles; i_part++){
        for(int i_system=0; i_system<nsyststodraw; i_system++){
	  for(int i_htype=0; i_htype<n_htypes; i_htype++){
	    if(!presence[i_exp][i_part][i_system][i_htype]) continue;
	    if(i_system==lauau02 || i_system==lcucu02) gr[i_exp][i_part][i_system][i_htype]->DrawClone("3");
	    gr[i_exp][i_part][i_system][i_htype]->Draw(drawop[i_htype]);
	    //leg->AddEntry(gr[i_exp][llam][i_system][i_htype],Form("%s , %s , %s",expname[i_exp],leg_partname[llam],systname[i_system]),"P");
	  }
        }
      }
    } 
    TLegend *leg_ALICE = new TLegend(0.30,0.15,0.69,0.3);
    leg_ALICE-> SetHeader("ALICE");
//     leg_ALICE-> SetNColumns(4);
    leg_ALICE->SetBorderSize(0);
    leg_ALICE->SetTextSize(0.03);
    leg_ALICE->SetFillColor(0);
    leg_ALICE->SetFillStyle(0);
    leg_ALICE->AddEntry(gr[lALICE][lom][lpbpb5][ltot],Form("Preliminary Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"p");
    //leg_ALICE->AddEntry(gr[lALICE][lom][lppb5][ltot],Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV , PLB 728 (2014) 25-38"),"p");
    leg_ALICE->AddEntry(gr[lALICE][lom][lppb5][ltot],Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"p");
    //leg_ALICE->AddEntry(gr[lALICE][lom][lpp7][ltot],Form("pp #sqrt{#it{s}} = 7 TeV , Nat. Phys. 13 (2017) 535-539"),"p");
    leg_ALICE->AddEntry(gr[lALICE][lom][lpp7][ltot],Form("pp #sqrt{#it{s}} = 7 TeV"),"p");
    leg_ALICE->Draw();
    TLegend *leg_STAR = new TLegend(0.61,0.30,0.82,0.43);
    leg_STAR-> SetHeader("STAR");
//     leg_STAR-> SetNColumns(2);
    leg_STAR->SetBorderSize(0);
    leg_STAR->SetTextSize(0.03);
    leg_STAR->SetFillColor(0);
    leg_STAR->SetFillStyle(0);
    leg_STAR->AddEntry(gr[lSTAR][lom][lauau02][ltot],Form("Au-Au #sqrt{#it{s}_{NN}} = 200 GeV"),"p");
    leg_STAR->AddEntry(gr[lSTAR][lom][lcucu02][ltot],Form("Cu-Cu #sqrt{#it{s}_{NN}} = 200 GeV"),"p");
    leg_STAR->AddEntry(gr[lSTAR][lom][lpp02][ltot],Form("pp #sqrt{#it{s}} = 200 GeV"),"p");
    //leg_STAR->Draw();
    
    TLatex lat_part;
    lat_part.SetTextAlign(11);
    lat_part.SetTextFont(42);
    lat_part.DrawLatex(2500,2.8e-3,"#Omega#times4");
    lat_part.DrawLatex(2500,1.2e-2,"#Xi#times2");
    lat_part.DrawLatex(2500,3.5e-2,"#Lambda");
    //lat_part.DrawLatex(3000,8.0e-2,"p");
    
    
  
  TString *canfilename = new TString();
  canfilename->Append("MyTalkPlot_noSTAR");
  canfilename->Append(".pdf");
  
  can->SaveAs(canfilename->Data());
}

void ScaleGraph(int exp, int part, int system, int type, float nscale){
  
  if((nscale-1)<1e-5) return;
  
  double x_point = 0;
  double x_error = 0;
  double y_point = 0;
  double y_error = 0;
  
//   for(int i_system=0; i_system<n_systems; i_system++){
//     for  (int i_htype=0; i_htype<n_htypes; i_htype++){
//       if(!presence[exp][part][i_system][i_htype]) continue;
      for(int ipoint=0; ipoint<gr[exp][part][system][type]->GetN(); ipoint++){
        gr[exp][part][system][type]->GetPoint(ipoint,x_point,y_point);
        x_error = gr[exp][part][system][type]->GetErrorX(ipoint);
        y_error = gr[exp][part][system][type]->GetErrorY(ipoint);
	gr[exp][part][system][type]->SetPoint(ipoint,x_point,y_point*nscale);
        if(!(part==llam && system==lpbpb276)) gr[exp][part][system][type]->SetPointError(ipoint,x_error,y_error*nscale);
      }
//     }    
//   }
  
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