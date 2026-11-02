/*=========================================================================
This code uses gridded visibility data in 3D and computed MABS (Final code: 20/02/2025)//
       Authors : Sukhdeep Singh Gill (IIT Kharagpur)
                 Prof. Somnath Bharadwaj (IIT Kharagpur)
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
#include <string.h>
#include "mabs_funcs.h"
#include "vars.h"
#include <stdalign.h>


int main(int argc, char *argv[])
{
  int Nbin, Rlzn,i, j, l, NB, *num,rcount, r1,r2,r3,rp;
  int nclen, ncsize, *vcp, total_vcp;
  long Nmax, sd_size;
  double LL, Area,t,T=omp_get_wtime(), Biana,Denom,Bierr, Num, tt, mu, PU1, PU2, PU3, Ntri; 
  double U1, U2, U3, dU,Umax1, Umin1, gfact, gfact_sq,*rdata;
  char file[512],file1[512],file2[512];
  FILE *fp,*fp1,*fp2;
 
  
    if(argc!=2)
   {
      printf("Usage: %s <Rlzn>\n", argv[0]);
      return 1;
    }
  sscanf(argv[1],"%d",&Rlzn);
  NRi = Rlzn;
  
  read_inputs_ess("input_params_ess");
  
  
  Umin          =  1.;
  Nbin   	=  22;
  Nr		=  Nbin;
  nc10          =  383;//95 (basically nc//2-1); //130;  //// 768*0.04/(4/17) // nc10 is included, ie tot_nc10=(2nc10+1). nc10_max=Nchan-1
  nclen         =  2*nc10+1;
  ncsize        =  nclen*nclen;
  uf_nc         =  552;//552; 144
  
  int max_threads = omp_get_max_threads();
  printf("Maximum available threads: %d\n", max_threads);
  Nthreads	=  max_threads/2 - 2;
  printf("Utilizing threads: %d\n", Nthreads);
  
  
  /*----------------------------------------------------------*/
  
  initialize(); 
  ReadGVnAC_gridding();
  if (flag_key == 1){printf("\nData is flagged\n");} else if (flag_key == 0){printf("\nData is Not flagged\n");} else {printf("\nGive correct flag key\n");}
  
  printf("NN=%ld, nc=%ld, Du=%lf, Flag=%d, nc10=%d, Unflagged Channles=%d\n", NN,nc,Du,flag_key,nc10,uf_nc);
  
  printf("data_path: %s\n", data_path);
  
  snprintf(res_path, sizeof(res_path), "results/v1");
  printf("result will be saved in %s\n", res_path);

  
  /*---------Normalize Data-------------------*/

  normalize_data();

  printf("\nData is ready\n");
  /*----------------------------------------------------------*/
  
  
  printf("\nTotal no. of channels (nc) = %ld",nc);
  printf("\nOriginal Dimensions N1 = %ld\n\n",N1);
  //printf("\nComputing on Reduced Dimensions N1 = %ld\n\n",Ng1);
  

  /*------- runtime arguments reading done. ---------*/
  
  /*------------ Checking the limit of U1 ------------*/
  annular_rings();

  printf("\t%s\t\t%s\t\t\t%s\n","ULow","UHigh","Uavg") ;
  for(i=0;i<Nr;++i)
  {
  printf("\t%f\t%f\t%f\n",ULow[i],UHigh[i],Uavg[i]) ;}
  
  
  
  /*------------ Use this if you want to define arrays outside rp loop ------------*/
  /*-- (1) Remove free(sd) in rp loop----*/
    
    Umax = UHigh[Nr-1]+3;
    Nmax = (long) (Umax/Du);
    Ng1  = (long) (3*Nmax);
    Ng1  = (Ng1 % 2 == 0) ? Ng1+2 : Ng1 + 3;
    Ng2  = Ng1;
    printf("\n\nNew grids Ng1 = %ld\n\n",Ng1);
    sd_size = Ng1*Ng2;
    /* -- This has been used to get final initial results without flaging. Ultimate --*/
    /*float **sd1 = allocate_2d_float_aligned(nc10+1, sd_size);
    float **sd2 = allocate_2d_float_aligned(nc10+1, sd_size);
    float **sd3 = allocate_2d_float_aligned(nc10+1, sd_size);*/
    
    float **sd1 = allocate_2d_float_aligned(nc, sd_size);
    float **sd2 = allocate_2d_float_aligned(nc, sd_size);
    float **sd3 = allocate_2d_float_aligned(nc, sd_size);

  
  /*------------ Initializing fft plan and openmp ------------*/

  
  fftw_init_threads();
  fftw_plan_with_nthreads(Nthreads);
  omp_set_num_threads(Nthreads);
  

  
  
  fprintf(stdout,"\nComputing triangles for all possible combination of rings.\n"); 
  calCR();
  
  read_rings(&rcount, &rdata);
  
  /*fprintf(stdout,"\nComputing Ds.\n"); 
  //wrapid();
  wrapidg();
  fprintf(stdout,"\nComputing Ds Done!!.\n"); */
  
  

/*----  (i) precompute_vcp_adaptive: Used for initial ultimate results, wo flagging and wo collapse (LL=4, nc10=130,nc=768)
   (ii) precompute_unflg_vcp: Used for flagged+collapsed chans (LL_para=0.68, nc10=95,nc=192, n_coll = 4)
   (iii)precompute_vcp_dnu: for unflagged collapsed channels (LL_para=0.68, nc10=95,nc=192, n_coll = 0) ---- */
    
    /*----  Unflagged Channel Index Read ---- */
    if (flag_key == 1){
    get_uf_chans("uf_chans.txt");
    printf("Unflagged Channels List:\n"); for ( i = 0; i < uf_nc; i++) {printf("%d ", uf_chans[i]);}   printf("\n\n");
    vcp = NULL;
    total_vcp = precompute_unflg_vcp(nc10, &vcp);
    printf("Total Viable Combinations of Planes: %d, reduced by:%f\n", total_vcp,(1.*total_vcp)/(1.*nc*nc*nc));
    }
   
    else if (flag_key == 0){
    vcp = NULL;
    //int total_vcp = precompute_vcp_adaptive(nc10, &vcp); // //precompute_vcp(nc, nc10, &vcp);
    total_vcp = precompute_vcp_dnu(nc,nc10, &vcp);
    ////int total_vcp = precompute_unflg_vcp(nc10, &vcp);
    printf("Total Viable Combinations of Planes: %d, reduced by:%f\n", total_vcp,(1.*total_vcp)/(1.*nc*nc*nc));
   }
   
   else {printf("\n\nERROR: Flag key is not given\n\n");}

    
  
  /*=====================  Bispectrum computation ====================*/
 
   fprintf(stdout,"\nComputing BS!!.\n"); 
   double start_time=omp_get_wtime();
   double inv_Ng1_Ng2 ;
   
   int iplane   = 0;                                    //0;150;300;450;600
   int out_Rlzn = Rlzn  + 0;//+ 240+ 0; 40; 80;120;160//0; 60; 120;180;240
   printf("\n\nInitial Plane=%d, Saved in R%d\n\n",iplane,out_Rlzn);
   

   float *bldnu = (float *)aligned_alloc(64, ncsize * sizeof(float));
   int *numdnu = (int *)aligned_alloc(64, ncsize * sizeof(int));

    
    fprintf(stdout,"\nStarting rp loop!!.\n");
    
    snprintf(file2, 512*sizeof(char),"results/time_R%d",out_Rlzn); // open file to write time
    
    /*fp2=fopen(file2,"w");
    fprintf(fp2,"Time for Realization = %d\n",out_Rlzn);
    fclose(fp2);*/
	   
   //for(rp=0;rp<rcount;++rp)
   for(rp=365;rp<rcount;++rp)   // rp=365(r1=15) ;561 is last r1=19 (with Umax=200)
   {
   r1=(int)rdata [rp * 4 + 0] ; 
   r2=(int)rdata [rp * 4 + 1] ; 
   r3=(int)rdata [rp * 4 + 2] ;  
    
    Umax = UHigh[r1]+3;
    Nmax = (long) (Umax/Du);
    Ng1  = (long) (3*Nmax);
    Ng1  = (Ng1 % 2 == 0) ? Ng1+2 : Ng1 + 3;
    Ng2  = Ng1;
    gfact= (1.*Ng1)/(1.*N1);
    gfact_sq = gfact*gfact;
    Num=gfact_sq*rdata [rp * 4 + 3] ;
    inv_Ng1_Ng2 = 1.0 / (Ng1*Ng2*1.);
    sd_size = Ng1*Ng2;
    
    printf("\n\nNew grids for r1=%d: Ng1 = %ld (rp=%d)\n\n",r1,Ng1,rp);
   
    /* --- This has not been used yet ---*/
    /*float **sd1 = allocate_2d_float_aligned(nc10+1, sd_size);   //remember: free(sd) in rp loop
    float **sd2 = allocate_2d_float_aligned(nc10+1, sd_size);
    float **sd3 = allocate_2d_float_aligned(nc10+1, sd_size);*/
    
    /* -- This has been used to get final initial results without flaging. Ultimate --*/
    /*fill_s2d_float_adaptive(sd1, r1, 1, 1, iplane); //Used for initial ultimate results, wo flagging and wo collapse (LL=4, nc10=130,nc=768)
    fill_s2d_float_adaptive(sd2, r2, 1, 1, iplane);
    fill_s2d_float_adaptive(sd3, r3, 1, 1, iplane);  */
    
    fill_s2d_float(sd1, r1, 1, 1); // for unflagged collapsed channels (LL_para=0.68, nc10=95,nc=192, n_coll = 0) 
    fill_s2d_float(sd2, r2, 1, 1);
    fill_s2d_float(sd3, r3, 1, 1);
    
    /*
    fill_s2d_float_unflg(sd1, r1, 1, 1); // Used for flagged+collapsed chans (LL_para=0.68, nc10=95,nc=192, n_coll = 4); fill_s2d_float_float can also used
    fill_s2d_float_unflg(sd2, r2, 1, 1);
    fill_s2d_float_unflg(sd3, r3, 1, 1);  */
    
 
   memset(bldnu, 0, ncsize * sizeof(float));
   memset(numdnu, 0, ncsize * sizeof(int));
	    
    //fprintf(stdout,"\nStarting p1 loop!!.\n");
    start_time=omp_get_wtime();
    
   
   omp_set_num_threads(Nthreads);
   #pragma omp parallel firstprivate(sd_size,ncsize,vcp)
{
    int ipc, base, p1, p2, p3, idx; //, dnu2, dnu3
    long index;
    float Bispec;
    
   float *private_bldnu = (float *)aligned_alloc(64, ncsize * sizeof(float));
   int *private_numdnu = (int *)aligned_alloc(64, ncsize * sizeof(int));
   
   memset(private_bldnu, 0, ncsize * sizeof(float));
   memset(private_numdnu, 0, ncsize * sizeof(int));
   
   #pragma omp for schedule(dynamic) nowait
   for(ipc=0;ipc<total_vcp;++ipc){
     base = 6 * ipc;
     p1 = vcp[base];
     p2 = vcp[base+1];
     p3 = vcp[base+2];
     //dnu2 = vcp[base+3];
     //dnu3 = vcp[base+4];
     idx = vcp[base+5];
     Bispec=0.;
     
     float *restrict ptr1 = sd1[p1];
     float *restrict ptr2 = sd2[p2];
     float *restrict ptr3 = sd3[p3];
        
     for(index=0;index<sd_size;++index){  
         Bispec+=ptr1[index]*ptr2[index]*ptr3[index];
         }
  
     private_bldnu[idx] += Bispec;
     private_numdnu[idx] += 1;
     
        
  }
    
    //fprintf(stdout,"Before Critical Area; time taken = %d hr %d min %d sec\n",(int)((omp_get_wtime()-start_time)/3600), (int)((omp_get_wtime()-start_time)/60)%60, (int)(omp_get_wtime()-start_time)%60);
    
       #pragma omp critical
    {
        for ( i = 0; i < ncsize; ++i) {
            bldnu[i] += private_bldnu[i];
            numdnu[i] += private_numdnu[i];
        }
    }
    free(private_bldnu);
    free(private_numdnu);
    
  }
  
  fprintf(stdout,"Nested loops time taken = %d hr %d min %d sec\n",(int)((omp_get_wtime()-start_time)/3600), (int)((omp_get_wtime()-start_time)/60)%60, (int)(omp_get_wtime()-start_time)%60);
          
          Ntri=Num*inv_Ng1_Ng2;
          
          snprintf(file1, 512 * sizeof(char), "%s/bsout_final_ra%d_rb%d_rc%d_R%d", res_path, r1, r2, r3, out_Rlzn);
          //snprintf(file1, 512*sizeof(char),"results/final_pblsh/flagged/coll/bsout_final_ra%d_rb%d_rc%d_R%d",r1,r2,r3,out_Rlzn); // open file to write bispectrum
   	  fp1=fopen(file1,"wb");
          //fprintf(fp1,"%.0lf\n",Ntri);  
          
          fwrite(&Ntri, sizeof(double), 1, fp1);
          fwrite(&Num, sizeof(double), 1, fp1);
	  if (fwrite(bldnu, sizeof(float), ncsize, fp1) != ncsize) {
	    perror("fwrite bldnu");
	    exit(EXIT_FAILURE);}

	  if (fwrite(numdnu, sizeof(int), ncsize, fp1) != ncsize) {
	    perror("fwrite numdnu");
	    exit(EXIT_FAILURE);}
	    
	  fclose(fp1);
          
          /*for (i = 0; i < nclen; i++) {
           for (j = 0; j < nclen; j++) {
	    fprintf(fp1,"%d\t%d\t%e\t%d\n",i-nc10,j-nc10,bldnu[i*nclen+j]/Num,numdnu[i*nclen+j]); 
	    }}*/
	    
	    
	    
	    
	  
	   /* fprintf(stdout,"rp=%d in Total time taken = %d hr %d min %d sec\n",rp,(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);
	    
	    fp2=fopen(file2,"a");
	   fprintf(fp2,"rp=%d in Total time taken = %d hr %d min %d sec\n",rp,(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);
	   fclose(fp2);*/
  
  
 /* free_2d_float_aligned(sd1);
  free_2d_float_aligned(sd2);
  free_2d_float_aligned(sd3);*/
  
  }
 free_2d_float_aligned(sd1);
 free_2d_float_aligned(sd2);
 free_2d_float_aligned(sd3);
 free(vcp);
 free(bldnu);
 free(numdnu);
 fprintf(stdout,"Total time taken = %d hr %d min %d sec\n",(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);

}

