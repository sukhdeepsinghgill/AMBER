#ifndef FUNCS_H
#define FUNCS_H

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

  // void read_uvimage(char filename1[], char filename2[]);
  void read_uvdata(char filename1[], char filename2[]);
  void read_inputs(char *in_file);
  void read_inputs_ess(char *in_file);
  void fft_init(int Nbin, double *uL,double *uH);
  void fft_init_orig(int nn, int Nbin, double *uL,double *uH);
  void calCR();
  void read_rings(int *rcountPtr, double **rdataPtr);
  void Dfg(int plane,int ring, double **DD, fftw_plan pd)  ;
  void wrapid();
  void get_D(int plane, int ring, int rlz, double **D);
  void wrapidg();
  void get_Dg(int plane, int ring, int rlz, double **D);
  void shapes(int nbin, double *uavg);
  double **allocate_2d_double(long N1,long N2);
  double **allocate_2d_double_aligned(long N1, long N2);
  float **allocate_2d_float_aligned(long N1, long N2);
  double ***allocate_3d_double(long N1,long N2,long N3);
  double ***allocate_3d_array(long nc, long Ng1, long Ng2);
  double ****allocate_4d_array(long nc, long Nr, long Ng1, long Ng2);
  void free_4d_array(double ****array);
  void free_3d_array(double ***array, long nc, long Ng1);
  void free_2d_double(double **array, int rows);
  void free_2d_double_v2(double **array);
  void free_2d_float_aligned(float **array);
  void fill_s4d(double ****sd, int rlz, int read_or_cal);
  void fill_s3d(double ***sd, int ring, int rlz, int read_or_cal);
  void fill_s2d(double **sd, int ring, int rlz, int read_or_cal);
  void fill_s2d_float(float **sd, int ring, int rlz, int read_or_cal);
  void fill_s2d_float_adaptive(float **sd, int ring, int rlz, int read_or_cal, int iplane);
  void fill_s2d_float_unflg(float **sd, int ring, int rlz, int read_or_cal);
  double P_T(double UU);
  int precompute_vcp(int nc, int nc10, int **vcp);
  int precompute_vcp_dnu(int nc, int nc10, int **vcp);
  int precompute_vcp_dnu_expt0(int nc, int nc10, int **vcp);
  int precompute_vcp_adaptive(int nc10, int **vcp);
  int precompute_unflg_vcp(int nc10, int **vcp);
  void fill_s2I(float **sI, int ring);
  void Ifg(int plane,int ring, double **II, fftw_plan pd);
  int get_uf_chans(const char *filename);
  
  int initialize();
  void ReadGVnAC_gridding();
  void normalize_data();
  void annular_rings();
  
  
 #endif // FUNCS_H
