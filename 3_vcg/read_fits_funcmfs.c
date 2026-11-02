//UVFITS reading functions
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <fitsio.h>
# include "read_fits_func.h"

void printerror(int status)
{
  if (status)
    {
      fits_report_error(stderr, status);
      exit( status );
    }
  return;
}

static fitsfile *fptr;
long el1,nel;
int status,anynul;
  
// extern in main
long gcount,pcount,nstokes,nchan,ncmplx;
float *data;
float *randpar; 
float  chan0,nu_chan0,del_chan;

int read_fits_header(char* in_file)
{
  char key_simple[FLEN_KEYWORD]="SIMPLE";
  char key_naxis[FLEN_KEYWORD]="NAXIS";
  char key_naxis2[FLEN_KEYWORD]="NAXIS2";
  char key_naxis3[FLEN_KEYWORD]="NAXIS3";
  char key_naxis4[FLEN_KEYWORD]="NAXIS4";
  char key_gcount[FLEN_KEYWORD]="GCOUNT";
  char key_pcount[FLEN_KEYWORD]="PCOUNT";
  char keyvalue[FLEN_VALUE];
  char comment[FLEN_COMMENT];
  
  int ii,jj,kk,st,ttest;

  status=0;
  anynul=0;
  printf("\n-----------Fits header-------------------\n");
  fits_open_file(&fptr, in_file, READONLY, &status);
  if(status==1){fprintf(stderr,"unable to open %s", in_file);	return 1;}// input file pointer is fptr

  fits_read_keyword(fptr,key_simple,keyvalue,comment,&status);
  if(*keyvalue!='T'){ fprintf(stderr,"Input File is NOT a SIMPLE FITS file\n");return 1;}

  if(fits_read_key_lng(fptr,key_naxis2,&ncmplx,comment,&status))
    printerror( status );
  printf("NAXIS2=%ld\n",ncmplx);

  if(fits_read_key_lng(fptr,key_naxis3,&nstokes,comment,&status))
    printerror( status );
  printf("NAXIS3=%ld\n",nstokes);
  
  if(fits_read_key_lng(fptr,key_naxis4,&nchan,comment,&status))
    printerror( status );
  printf("NAXIS4=%ld\n",nchan);
  
  if(fits_read_key_lng(fptr,key_gcount,&gcount,comment,&status))
    printerror( status );
  printf("GCOUNT=%ld\n",gcount);
  
  if(fits_read_key_lng(fptr,key_pcount,&pcount,comment,&status))
    printerror( status );
  printf("PCOUNT=%ld\n",pcount);
   
  if(fits_read_key_flt(fptr, "CRVAL4", &nu_chan0, comment,&status))
    printerror( status );
  printf("CRVAL4=%e\n",nu_chan0);
  
  if(fits_read_key_flt(fptr, "CDELT4", &del_chan, comment,&status))
    printerror( status );
  printf("CDELT4=%e\n",del_chan);
  
  if(fits_read_key_flt(fptr, "CRPIX4", &chan0, comment,&status))
    printerror( status );
  printf("CRPIX4=%f\n",chan0);//ref. pixel; CRPIXn can be a floating point no. x for which the physical value is CRVALn & increment CDELTn . We should interpret CRPIXn as the location of a FORTRAN index. For any given channel i the freq. is nu[i]=CRVAL4 + CDELT4*(i+1-CRPIX4).
  
   
  printf("chan0=%f\tnu_chan0=%f\n",chan0,nu_chan0);
  
  randpar=(float*) calloc((pcount), sizeof(float));

  nel=(long)ncmplx*nstokes*nchan; // (re,im,wt)*(RR*LL)*(n_ave)
  printf("nel=%ld\n",nel);

  if((data = (float*)calloc(nel, sizeof(float)))==NULL)
      { fprintf(stderr,"malloc failure\n"); return 1;}
  
  return 0;
  
}

int read_ranpar(long grp)
{
  status=0;
  el1=1;  
  if(fits_read_grppar_flt(fptr,grp,el1,pcount,randpar,&status))
  printerror( status );
  return 0;
 }

int read_data(long group)
{
  float nulval;

  status=0;
  anynul=0;
  nulval=0.;
  el1=1;

  if(fits_read_img_flt(fptr,group,el1,nel,nulval,data,&anynul,&status))
  printerror( status );
  return 0;
}

int close_fits() 
{
  if(fits_close_file(fptr,&status))
    printerror( status );
  return(0);
}
