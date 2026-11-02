/*=========================================================================
       This code uses gridded visibility data in 3D for and computes the maps.
       Authors : Sukhdeep Singh Gill (IIT Kharagpur)
===========================================================================*/
#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<unistd.h>
#include<omp.h>
#include<fitsio.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>

/*---- Global variables ----*/


int Nr, NRi, nc10;     
double Umax, Umin;                  
double *ULow, *UHigh, *Uavg;
char data_path[512];//="../3_vcg/output_data/GV_data_1";




long NN,N1,N2,nc,Ng1,Ng2; // Dimension
int Nthreads, F_Task;
double Du, pi=M_PI,nu_chan0,del_chan,chan0; //Du is grid spacing in uv plane
double **ro, **ro_u; // for power spectrum
double **b_ro, **b_ro_u; // for bispectrum
fftw_plan p_ro, q_ro;


int Nm,Nu;        
double fx,FWHM,dUx;

double *nua;                   
double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL,**NTiG;
char GV[256];//="../3_vcg/output_data/GV_data_1";
/*---------------------------*/

int main(int argc, char *argv[])
{
  int Nbin, Rlzn; 
  long index;
  double sigma1, sigma2, ff1, LL, Area,ff2,Av=1.,beta=-1.; 
  unsigned long int TotN=0;
  char  OPT;
  char file[512],file1[512],file2[512];
  FILE *fp,*fp1,*fp2;
  double t,T=omp_get_wtime();
  
  int i, j, l, NB, *num;
  double *binu,*binps,*logU,*logPU;
  
  int ii, jj;
  double Bispec, Biana,Denom,Bierr, Num, tt, mu, PU1, PU2, PU3;
  double U1, U2, U3, dU;
  double Umax1, Umin1;
  double **I1, **I2, **I3, **D1, **D2, **D3, **O1, **O2, **O3; ;

  
  // void read_uvimage(char filename1[], char filename2[]);
  void read_uvdata(char filename1[], char filename2[]);
  void read_inputs(char *in_file);
  void power_spec_GV(int NB, double umax, double umin, double *usum, double *pksum, int *no);
  void power_spec_GV_Ntig(int nn,int NB, double umax, double umin, double *usum,double *pksum,int *no);
  void cl_dnu(int nn1, int nn2, int NB, double umax, double umin, double *usum,double *pksum,int *no);
  double **allocate_2d_double(int N1,int N2);
  double P_T(double UU);
  
  int initialize();
  void ReadGVnAC_gridding();
  
  /*------ reading parameters from runtime arguments --------*/

  
  sscanf(argv[1],"%d",&Rlzn);
  sscanf(argv[2],"%s",GV); 

  NRi = Rlzn;

  
  read_inputs("input_params");
  nc = 768; //remove this line in general

  
  //Du = 1.0254384902642049;//((180.*60.)/(3.14159265359*2048.*1.719));//dUx; 
  N1 = NN;
  N2 = NN;
  

  initialize(); 
  printf("\nArrays initialized\n");
  ReadGVnAC_gridding();
  
  
  
    /*---------Data normalization-------------------*/
  
  
    for(ii=0; ii<N1;ii++)
     for(jj=0; jj<N2/2+1;jj++)
      {
	index = (long)ii*(N2/2+1)+jj; 
	for(i=0; i<nc; i++)
	  {
	    if (NTiG[index][i]>0)
	    {
	    GVRRre[index][i] /= NTiG[index][i];
	    GVRRim[index][i] /= NTiG[index][i];
	    }
	  }
      }
  
    for(ii=0; ii<N1;ii++)
     for(jj=0; jj<N2/2+1;jj++)
      {
	index = (long)ii*(N2/2+1)+jj; 
	for(i=0; i<nc; i++)
	  {
	  if (NTiG[index][i]>0)
	   {
	    NTiG[index][i]=1;
	   }
	  }
      }
  
  /*----------------------------------------------------------*/
  
  
  
  printf("\nRead gridded visibilities\n");
  printf("\nTotal no. of channels (nc) = %ld",nc);
  printf("\nFWHM = %.2f arcmin.",FWHM);
  printf("\nGrid spacing, dU = %.2f wavelengths\n",dUx);
  printf("\nFreq. nu_0 = %e",nua[0]);
  printf("\nNm = %d",Nm);
  printf("\nNu = %d",Nu);
  printf("\nN1 = %ld",N1);
  


  
 /*---------------------- Binning done --------------------------*/
 
 
  fprintf(stdout,"\nComputing the MAPS. The result will be in results/cldnu.\n");  
  
  /*------------ Initializing fft plan and openmp ------------*/
  fftw_init_threads();
  fftw_plan_with_nthreads(Nthreads);
  omp_set_num_threads(Nthreads);
  

  /*================ binned power spectrum computation ===============*/ 
  NB = 15;
  Umax1 = 500.;
  Umin1 =1.;
  int plane_nn,nn1,nn2;
  
  num=(int*)calloc(NB,sizeof(int));
  binu=(double*)calloc(NB,sizeof(double));
  binps=(double*)calloc(NB,sizeof(double));
  
  double **cldnu,*numdnu;
  cldnu = (double**)calloc(NB,sizeof(double*));  
  numdnu = (double*)calloc(nc,sizeof(double));  
  for(i=0;i<NB;i++)
    {cldnu[i]=(double*)calloc(nc,sizeof(double)); }
    
    for (j=0;j<nc;j++)
    {numdnu[j] = 0.;
    for(i=0;i<NB;i++){
    cldnu[i][j]=0.;}}
  

  
 for(nn1=0;nn1<nc;++nn1){  // for cldnu loops
 printf("Plane_1=%d\n",nn1);
  for(nn2=nn1;nn2<nc;++nn2){
  

  
 // t=omp_get_wtime();
  
   /*power_spec_GV_Ntig(plane_nn,NB, Umax1, Umin1, binu, binps, num); // computes cl(dnu=0) for single plane.
   snprintf(file2, 512*sizeof(char),"results/psout_p%d_%s",plane_nn,GV); 
   fp2=fopen(file2,"w");

  for(i=0;i<NB;++i)
    if(num[i]>0)
      {
  	fprintf(fp2,"%e\t%e\t%d\n",binu[i],binps[i],num[i]);
  	//logU[i]=log10(binu[i]); // log is done for the interpolation used in B_Ana
  	//logPU[i]=log10(binps[i]);
      }
  fclose(fp2);
  fprintf(stdout,"Power spectrum time = %e sec\n",omp_get_wtime()-t);*/
  
  
  
  cl_dnu( nn1, nn2, NB, Umax1, Umin1, binu, binps, num);
  numdnu[nn2-nn1] += 1.;
  for(i=0;i<NB;++i)
      {cldnu[i][nn2-nn1] += binps[i];}

  
  }
  }
  
   snprintf(file2, 512*sizeof(char),"results/cldnu"); // open file to write bispectrum
   fp2=fopen(file2,"w");

for(i=0;i<NB;++i){
 for(j=0;j<nc;++j){
  
  	fprintf(fp2,"%d\t%e\t%d\t%e\t%e\t%d\n",i,binu[i],j,cldnu[i][j],numdnu[j],num[i]);
      
      }}
  fclose(fp2);
  /*===================================================================*/
  
 
 
 
fprintf(stdout,"Total time taken = %d hr %d min %d sec\n",(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);
 

}
