TGraphErrors *GenerateGraphFromTxt(TString experiment, TString system_energy, TString particle, bool onlyunc = kFALSE, int sumplusminus = 0, bool addSTARpp=kFALSE){
  
  int npoints = 0;
  
  FILE * infile;
  infile = fopen("Data_sumPARToverPI.txt","r");
//   if(sumplusminus==0) infile = fopen("Data_sumPARToverPI_v2.txt","r");
//   else if(sumplusminus==1) infile = fopen("Data_antiPARToverPI_v2.txt","r");
//   else if(sumplusminus==2)infile = fopen("Data_PARToverPI_v2.txt","r");
  
  if(addSTARpp && system_energy.Contains("AuAu_0.2")) system_energy.Append("_pp_0.2");
  
  char exp[20];
  char syst[20];
  char part[20];
  float x, ex, y, estaty, esysty, esystuncy;
  char comm[10000];
  
  char dum[10000];
  fgets(dum, 10000, infile);
  fgets(dum, 10000, infile);
  
  //arrays for graph storing
    double x_var[100], y_var[100], ex_var[100], ey_var[100];
  
  while(!feof(infile)){
    fscanf(infile,"%s %s %s %f %f %f %f %f %f %s",exp,syst,part,&x,&ex,&y,&estaty,&esysty,&esystuncy,comm);
    if(experiment.Contains(exp) && system_energy.Contains(syst) && particle.Contains(part)){
      //printf("here\n");
      x_var[npoints]  =  x;
      ex_var[npoints] = ex;
      y_var[npoints]  = y;
      if(!onlyunc) ey_var[npoints] = TMath::Sqrt(estaty*estaty+esysty*esysty);
      else ey_var[npoints] = TMath::Sqrt(estaty*estaty+esystuncy*esystuncy);
      if(experiment.Contains("STAR") && particle.Contains("Proton")){ //roughly correct for feed-down
	y_var[npoints] *= 0.6;
	ey_var[npoints] *= 0.6;
      }
      if(experiment.Contains("STAR") && system_energy.Contains("AuAu_0.2")){ //roughly correct for feed-down
        if(particle.Contains("Lambda")) ey_var[npoints] = y_var[npoints]*7./100;
	else if(particle.Contains("Xi")) ey_var[npoints] = y_var[npoints]*8./100;
      }
      npoints++;
    }
  }

  fclose(infile);
  
  const int n = npoints;
  
  TGraphErrors *gr = new TGraphErrors(n,x_var,y_var,ex_var,ey_var);
  
  return gr;
  
}



// void ReadTxtData(){
//   
//   FILE * infile;
//   infile = fopen("Data_PARToverPI.txt","r");
//   
//   char exp[20];
//   char syst[20];
//   char part[20];
//   float x, ex, y, estaty, esysty, esystuncy;
//   char comm[100];
//   
//   char dum[10000];
//   fgets(dum, 10000, infile);
//   fgets(dum, 10000, infile);
//   
//   while(!feof(infile)){
//     fscanf(infile,"%s %s %s %f %f %f %f %f %f %s",exp,syst,part,&x,&ex,&y,&estaty,&esysty,&esystuncy,comm);
//     printf("%s %s %s %f %f %f %f %f %f %s\n",exp,syst,part,x,ex,y,estaty,esysty,esystuncy,comm);
//   }
//   
//   
//   
// }