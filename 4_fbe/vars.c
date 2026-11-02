#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<unistd.h>
#include<omp.h>
#include<fitsio.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include <string.h>
#include "vars.h"

/*---- Global variables ----*/

long NN,N1,N2,nc,Ng1,Ng2; // Dimension
int Nthreads, F_Task, Nr, NRi, nc10, uf_nc, *uf_chans;
double Du, pi=M_PI,nu_chan0,del_chan,chan0; //Du is grid spacing in uv plane
double **ro, **ro_u; // for power spectrum
double **b_ro, **b_ro_u; // for bispectrum
fftw_plan p_ro, q_ro;


int Nm,Nu, flag_key;        
double fx,FWHM,dUx, Umax, Umin;

double *nua;                   
double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL,**NTiG;
double *ULow, *UHigh, *Uavg;
char GV[256],data_path[512], res_path[512];//="../3_vcg/output_data/GV_data_1";
