#ifndef VARS_H
#define VARS_H


/*---- Global variables ----*/

extern long NN,N1,N2,nc,Ng1,Ng2; // Dimension
extern int Nthreads, F_Task, Nr, NRi, nc10, uf_nc, *uf_chans;
extern double Du, pi,nu_chan0,del_chan,chan0; //Du is grid spacing in uv plane
extern double **ro, **ro_u; // for power spectrum
extern double **b_ro, **b_ro_u; // for bispectrum
extern fftw_plan p_ro, q_ro;


extern int Nm,Nu, flag_key;        
extern double fx,FWHM,dUx, Umax, Umin;

extern double *nua;                   
extern double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL,**NTiG;
extern double *ULow, *UHigh, *Uavg;
extern char GV[256],data_path[512], res_path[512];//="../3_vcg/output_data/GV_data_1";


#endif // VARS_H
