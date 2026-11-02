# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include <fitsio.h>
# include <time.h>
# include "read_fits_func.h"
# include "tgefuncs.h"

extern long gcount,nstokes,nchan;
extern float *randpar,*data;
extern float del_chan,nu_chan0,chan0; 
extern char INFITS[128],input[128],GV[128],file_ntig[128];
extern long NN, n_collapse;
extern double dU,theta_0,FWHM,Umax;
extern int flag_key;

int nc, nc_coll,Nm,Ng,Nu,chan1=1;        
double f,Uval,theta_w,theta_eff;//dU

long group;
double *nua;                   
double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL,**NTiG;
double **GVRRre_coll, **GVRRim_coll ,**NTiG_coll;


void read_inputs(char *in_file)
{
  FILE *fp;
  fp=fopen(in_file,"r"); 
  if( fscanf(fp,"%d%lf%lf%lf",&nc,&Umax,&FWHM,&f) == 4){
    printf("\nInput parameters: nc=%d, Umax=%.1lf, FWHM=%.1lf arcmin, f=%.1lf",nc,Umax,FWHM,f);
  }
  fclose(fp);
}

double winf(double udif)
{
  double y;
  y=M_PI*pow(theta_w,2.)*exp(-1.*M_PI*M_PI*pow(theta_w,2.)*udif); // udif = (U_g -U_i)^2
  return(y);
}

int initialize()
{
  int ii,jj,nn,sf,Ngrid; 

  read_fits_header(INFITS);  

  theta_0=0.6*FWHM*M_PI/(180.*60.);         
  theta_w=f*theta_0; 
  theta_eff=f*theta_0/sqrt(1+f*f);

  sf = 2; //sampling factor
  dU=sqrt(log(2))/(M_PI*theta_eff*sf); // sampling the effective PB twice.
  Nm=3*sf; //sqrt(ln(500)/ln(2) = 3)

  Ng = (int) ceil(Umax/dU) + Nm; // padding for convolution. 
  Nu = 2*Ng +1;                        
  Ngrid = Nu*Nu;

  nua = (double*)calloc(nc,sizeof(double));
  for(nn=0;nn<nc;nn++)
    {
      nua[nn]= (nu_chan0 + del_chan*(nn - (chan0-1)));
    }
                    
  
  GVRRre = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRim = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLre = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLim = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      GVRRre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVRRim[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLim[ii]=(double*)calloc(nc,sizeof(double)); 
    }

  SCRR = (double**)calloc(Ngrid,sizeof(double*));  
  SCLL = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      SCRR[ii]=(double*)calloc(nc,sizeof(double)); 
      SCLL[ii]=(double*)calloc(nc,sizeof(double)); 
    }

  printf("\nTotal no. of channels (nc) = %d",nc);
  printf("\nFWHM = %.2f arcmin.",FWHM);
  printf("\nGrid spacing, dU = %.2f wavelengths",dU);
  printf("\n Ng = %d",Ng);
}



void group_loop()
{
  float u1,v1,signv;
  int nn,ii,jj,kk,ll,nx,ny,index,ch1,ch2,dnu;
  double diff,wt,gRRr,gRRi,gLLr,gLLi,sRRr,sLLr;
  float *RRre,*RRim,*LLre,*LLim; 
  double *V2RR, *V2LL;
  RRre = (float*)calloc(nc,sizeof(float));
  RRim = (float*)calloc(nc,sizeof(float));
  LLre = (float*)calloc(nc,sizeof(float));
  LLim = (float*)calloc(nc,sizeof(float));
  V2RR = (double*)calloc(nc,sizeof(double));
  V2LL = (double*)calloc(nc,sizeof(double));

  typedef struct visibility_type {float re,im,wt;} visibility; 
  visibility *v,*RR,*LL;
  
  int count=0;
  time_t time_s = 0,time_e,time_d = 0,time_d1=0;
  for(group=1;group<=gcount;group++)
    {    
      read_ranpar(group);
      signv = (randpar[1]<0.) ? -1. :1. ;       // flip uv
      u1 = nu_chan0*signv*randpar[0];
      v1 = nu_chan0*signv*randpar[1];

      Uval=sqrt(u1*u1+v1*v1);

      if(Uval<=Umax) 
	{	  
	  read_data(group);
	  v=(visibility *)data;
	  v+=(chan1-1)*nstokes;
	  RR=v;
	  LL=v+1;
	    
	  for(nn=0;nn<nc;nn++)
	    {	         
	      if(RR->wt>0.)
		{
		  RRre[nn]=RR->re; 
		  RRim[nn]=signv*RR->im;
		}
	      else 
		{
		  RRre[nn]=0.; RRim[nn]=0.;
		}

	      if(LL->wt>0.)
		{
		  LLre[nn]=LL->re; 
		  LLim[nn]=signv*LL->im;
		}
	      else 
		{
		  LLre[nn]=0.; LLim[nn]=0.;
		}
	      v+=nstokes; 
	      RR=v;
	      LL=v+1; 
	    }



	  time(&time_s);
	  for(nn=0;nn<nc;nn++)
	    { 
	      V2RR[nn]=0.; V2LL[nn]=0.;
	    }
	  for(ch1=0;ch1<nc;ch1++)	  
	    for(ch2=ch1;ch2<nc;ch2++)
	      {
		dnu = ch2-ch1;
		V2RR[dnu] += (RRre[ch1]*RRre[ch2] + RRim[ch1]*RRim[ch2]);		
		V2LL[dnu] += (LLre[ch1]*LLre[ch2] + LLim[ch1]*LLim[ch2]);
	      }
	  time(&time_e);
	  time_d += difftime(time_e,time_s);



	  time(&time_s);	  
	  ii = (int) round(u1/dU); 
	  jj = (int) round(v1/dU); 

	  for (kk=-Nm;kk<=Nm;kk++)           
	    for(ll=-Nm;ll<=Nm;ll++)
	      {
		nx = ii+kk;
		ny = jj+ll;
		diff = pow(nx*dU-u1,2.) +  pow(ny*dU-v1,2.);  
		wt = winf(diff);
		    
		for(nn=0;nn<nc;nn++) // gridding & self-correlation.
		  {
		    gRRr = wt*RRre[nn]; 
		    gRRi = wt*RRim[nn];
		    gLLr = wt*LLre[nn]; 
		    gLLi = wt*LLim[nn];
		    sRRr = wt*wt*V2RR[nn];
		    sLLr = wt*wt*V2LL[nn];

		    		  
		    index = (Ng+nx)*Nu + (Ng+ny);		  
		    GVRRre[index][nn] += gRRr; 
		    GVRRim[index][nn] += gRRi;
		    GVLLre[index][nn] += gLLr; 
		    GVLLim[index][nn] += gLLi;
		    SCRR[index][nn] += sRRr;		    	  
		    SCLL[index][nn] += sLLr;		    	  

		    index = (Ng-nx)*Nu + (Ng-ny); // inversion points.
		    GVRRre[index][nn] += gRRr; 
		    GVRRim[index][nn] -= gRRi;
		    GVLLre[index][nn] += gLLr; 
		    GVLLim[index][nn] -= gLLi;
		    SCRR[index][nn] += sRRr;		    	  
		    SCLL[index][nn] += sLLr;		    	  
		  }
	      } 
	  time(&time_e);
	  time_d1 += difftime(time_e,time_s);
	  
	  count += 1;
	}
    }
  printf("\nTime taken in self-correlation: %d hr %d min %d sec\n", (int)(time_d/3600), (int)((time_d/60)%60),(int)(time_d%60));
  printf("\nTime taken in gridding: %d hr %d min %d sec\n", (int)(time_d1/3600), (int)((time_d1/60)%60),(int)(time_d1%60));
  printf("\nThere are %d baselines within the uv range.",count);       
  close_fits();
}

void WriteGVnAC()
{
  int ii,jj,nn, index;
  FILE *fp;
  int MM = 2*Nm;
  printf("\nSize of the output GV: %d*%d*%d.",Nu-2*MM,Nu-2*MM,nc);  

  fp=fopen(GV,"wb"); 
  fwrite(&FWHM,sizeof(FWHM),1,fp);
  fwrite(&dU,sizeof(dU),1,fp);
  fwrite(&f,sizeof(f),1,fp);
  for(ii=0;ii<nc;ii++)
    fwrite(&nua[ii],sizeof(nua[ii]),1,fp);

  for(ii=MM; ii<Nu-MM;ii++)
    for(jj=MM; jj<Nu-MM;jj++)
      {
	index=ii*Nu+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVRRre[index][nn],sizeof(GVRRre[index][jj]),1,fp);
	    fwrite(&GVRRim[index][nn],sizeof(GVRRim[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",GVRRre[index][nn],GVRRim[index][nn]);
	  }
      }
    
  for(ii=MM; ii<Nu-MM;ii++)
    for(jj=MM; jj<Nu-MM;jj++)
      {
	index=ii*Nu+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVLLre[index][nn],sizeof(GVLLre[index][jj]),1,fp);
	    fwrite(&GVLLim[index][nn],sizeof(GVLLim[index][jj]),1,fp);
	  }
      }

  //Self-correlation. (real) 
  for(ii=MM; ii<Nu-MM;ii++)
    for(jj=MM; jj<Nu-MM;jj++)
      {
	index=ii*Nu+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&SCRR[index][nn],sizeof(SCRR[index][jj]),1,fp);
	    fwrite(&SCLL[index][nn],sizeof(SCLL[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",SCRR[index][nn],SCRR[index][nn]);
	  }
      }
  fclose(fp); 
}









int initialize_gridding()
{
  int ii,jj,nn,sf;
  long Ngrid; 

  read_fits_header(INFITS);  

 // theta_0=0.6*FWHM*M_PI/(180.*60.);         

  //sf = 2; //sampling factor
  //dU=sqrt(log(2))/(M_PI*theta_eff*sf); // sampling the effective PB twice.
  //Nm=3*sf; //sqrt(ln(500)/ln(2) = 3)

  //Ng = NN;//(int) ceil(Umax/dU) ; // padding for convolution. 
  //Nu = NN;//2*Ng +1;     
  
  theta_0=0.6*FWHM*M_PI/(180.*60.); 
  dU= sqrt(log(2))/(M_PI*theta_0);//((180.*60.)/(3.14159265359*2048.*1.719));
  Ng = (int) ceil(Umax/dU); 
  NN = (long) 2*Ng ; //2*Ng +1; 
  nc_coll = nc/n_collapse;   
  printf("\n Grids after collapse = %ld",nc_coll);
                   
  Ngrid = (long) NN*(NN/2+1);//Nu*Nu;
  printf("\n Ngrid = %ld",Ngrid);

  nua = (double*)calloc(nc,sizeof(double));
  for(nn=0;nn<nc;nn++)
    {
      nua[nn]= (nu_chan0 + del_chan*(nn - (chan0-1)));
    }
                    
  
  GVRRre = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRim = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLre = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLim = (double**)calloc(Ngrid,sizeof(double*));  
  NTiG   = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRre_coll = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRim_coll = (double**)calloc(Ngrid,sizeof(double*)); 
  NTiG_coll   = (double**)calloc(Ngrid,sizeof(double*));
  for(ii=0;ii<Ngrid;ii++)
    {
      GVRRre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVRRim[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLim[ii]=(double*)calloc(nc,sizeof(double)); 
      NTiG[ii]=(double*)calloc(nc,sizeof(double));
      GVRRre_coll[ii]=(double*)calloc(nc_coll,sizeof(double)); 
      GVRRim_coll[ii]=(double*)calloc(nc_coll,sizeof(double));
      NTiG_coll[ii]=(double*)calloc(nc_coll,sizeof(double));
    }

  /*SCRR = (double**)calloc(Ngrid,sizeof(double*));  
  SCLL = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      SCRR[ii]=(double*)calloc(nc,sizeof(double)); 
      SCLL[ii]=(double*)calloc(nc,sizeof(double)); 
    }*/

  printf("\nTotal no. of channels (nc) = %d",nc);
  printf("\nFWHM = %.2f arcmin.",FWHM);
  printf("\nGrid spacing, dU = %.2f wavelengths",dU);
  printf("\n NN = %ld",NN);
}





















void group_loop_gridding()
{
  float u1,v1,signv;
  int nn,ii,jj,kk,ll,nx,ny,ch1,ch2,dnu;
  long index;
  double diff,wt,gRRr,gRRi,gLLr,gLLi,sRRr,sLLr;
  float *RRre,*RRim,*LLre,*LLim, *Ng_RR, *Ng_LL; 
  double *V2RR, *V2LL;
  RRre = (float*)calloc(nc,sizeof(float));
  RRim = (float*)calloc(nc,sizeof(float));
  LLre = (float*)calloc(nc,sizeof(float));
  LLim = (float*)calloc(nc,sizeof(float));
  //V2RR = (double*)calloc(nc,sizeof(double));
  //V2LL = (double*)calloc(nc,sizeof(double));
  Ng_RR = (float*)calloc(nc,sizeof(float));
  Ng_LL = (float*)calloc(nc,sizeof(float));

  typedef struct visibility_type {float re,im,wt;} visibility; 
  visibility *v,*RR,*LL;
  
  int count=0;
  time_t time_s = 0,time_e,time_d = 0,time_d1=0;
  for(group=1;group<=gcount;group++)
    {    
      read_ranpar(group);
      signv = (randpar[1]<0.) ? -1. :1. ;       // flip uv
      u1 = nu_chan0*signv*randpar[0];
      v1 = nu_chan0*signv*randpar[1];

      Uval=sqrt(u1*u1+v1*v1);

      if(Uval<=Umax) 
	{	  
	  read_data(group);
	  v=(visibility *)data;
	  v+=(chan1-1)*nstokes;
	  RR=v;
	  LL=v+1;
	  
	 
	    
	  for(nn=0;nn<nc;nn++)
	    {	         
	      if(RR->wt>0.)
		{
		  RRre[nn]=RR->re; 
		  RRim[nn]=signv*RR->im;
		  Ng_RR[nn] = 1.;
		}
	      else 
		{
		  RRre[nn]=0.; RRim[nn]=0.; Ng_RR[nn] = 0.;
		}

	      if(LL->wt>0.)
		{
		  LLre[nn]=LL->re; 
		  LLim[nn]=signv*LL->im;
		  Ng_LL[nn] = 1.;
		}
	      else 
		{
		  LLre[nn]=0.; LLim[nn]=0.; Ng_LL[nn] = 0.;
		}
	      v+=nstokes; 
	      RR=v;
	      LL=v+1; 
	    }



	  time(&time_s);
	  /*for(nn=0;nn<nc;nn++)
	    { 
	      V2RR[nn]=0.; V2LL[nn]=0.;
	    }
	  for(ch1=0;ch1<nc;ch1++)	  
	    for(ch2=ch1;ch2<nc;ch2++)
	      {
		dnu = ch2-ch1;
		V2RR[dnu] += (RRre[ch1]*RRre[ch2] + RRim[ch1]*RRim[ch2]);		
		V2LL[dnu] += (LLre[ch1]*LLre[ch2] + LLim[ch1]*LLim[ch2]);
	      }*/
	  time(&time_e);
	  time_d += difftime(time_e,time_s);



	  time(&time_s);	  
	  ii = (int) round(u1/dU); 
	  ii = (ii<0) ? NN+ii : ii ;
	  jj = (int) round(v1/dU); 
	  
  /*if (group<50)
      {nn=0;
      //printf("(%f,%f,%f) \n",RR->re,RR->im,RR->wt);
      //printf("(%f,%f) \n",RRre[nn],RRim[nn]);
      printf("(%d,%d) \n",ii,jj);
      //printf("(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f) \n",data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11]);
      }*/
	  

	    
	for(nn=0;nn<nc;nn++) // gridding & self-correlation.
	  {		  
	    //index = (Ng+ii)*Nu + (Ng+jj);
	    index = (long) ii*(NN/2+1)+jj; 		  
	    GVRRre[index][nn] += RRre[nn]; 
	    GVRRim[index][nn] += RRim[nn];
	    GVLLre[index][nn] += LLre[nn]; 
	    GVLLim[index][nn] += LLim[nn];
	    //SCRR[index][nn] += V2RR[nn];		    	  
	    //SCLL[index][nn] += V2LL[nn];
	    
	    //NTiG[index][nn] += 1.; // use it if flagged data also contribute to Ng;
	    NTiG[index][nn] += Ng_RR[nn];		    	  	    	  
	  }
	      
	  time(&time_e);
	  time_d1 += difftime(time_e,time_s);
	  
	  count += 1;
	}
    }
  printf("\nTime taken in self-correlation: %d hr %d min %d sec\n", (int)(time_d/3600), (int)((time_d/60)%60),(int)(time_d%60));
  printf("\nTime taken in gridding: %d hr %d min %d sec\n", (int)(time_d1/3600), (int)((time_d1/60)%60),(int)(time_d1%60));
  printf("\nThere are %d baselines within the uv range.",count);       
  close_fits();
}





void group_loop_gridding_without_flag()
{
  float u1,v1,signv;
  int nn,ii,jj,kk,ll,nx,ny,ch1,ch2,dnu;
  long index;
  double diff,wt,gRRr,gRRi,gLLr,gLLi,sRRr,sLLr;
  float *RRre,*RRim,*LLre,*LLim; 
  double *V2RR, *V2LL;
  RRre = (float*)calloc(nc,sizeof(float));
  RRim = (float*)calloc(nc,sizeof(float));
  LLre = (float*)calloc(nc,sizeof(float));
  LLim = (float*)calloc(nc,sizeof(float));
  V2RR = (double*)calloc(nc,sizeof(double));
  V2LL = (double*)calloc(nc,sizeof(double));

  typedef struct visibility_type {float re,im,wt;} visibility; 
  visibility *v,*RR,*LL;
  
  int count=0;
  time_t time_s = 0,time_e,time_d = 0,time_d1=0;
  for(group=1;group<=gcount;group++)
    {    
      read_ranpar(group);
      signv = (randpar[1]<0.) ? -1. :1. ;       // flip uv
      u1 = nu_chan0*signv*randpar[0];
      v1 = nu_chan0*signv*randpar[1];

      Uval=sqrt(u1*u1+v1*v1);

      if(Uval<=Umax) 
	{	  
	  read_data(group);
	  v=(visibility *)data;
	  v+=(chan1-1)*nstokes;
	  RR=v;
	  LL=v+1;
	  
	 
	    
	  for(nn=0;nn<nc;nn++)
	    {	         
	      RRre[nn]=RR->re; 
	      RRim[nn]=signv*RR->im;

	      LLre[nn]=LL->re; 
	      LLim[nn]=signv*LL->im;
		
	      v+=nstokes; 
	      RR=v;
	      LL=v+1; 
	    }



	  /*time(&time_s);
	  for(nn=0;nn<nc;nn++)
	    { 
	      V2RR[nn]=0.; V2LL[nn]=0.;
	    }
	  for(ch1=0;ch1<nc;ch1++)	  
	    for(ch2=ch1;ch2<nc;ch2++)
	      {
		dnu = ch2-ch1;
		V2RR[dnu] += (RRre[ch1]*RRre[ch2] + RRim[ch1]*RRim[ch2]);		
		V2LL[dnu] += (LLre[ch1]*LLre[ch2] + LLim[ch1]*LLim[ch2]);
	      }*/
	  time(&time_e);
	  time_d += difftime(time_e,time_s);



	  time(&time_s);	  
	  ii = (int) round(u1/dU); 
	  ii = (ii<0) ? NN+ii : ii ;
	  jj = (int) round(v1/dU); 
	  
  /*if (group<50)
      {nn=0;
      //printf("(%f,%f,%f) \n",RR->re,RR->im,RR->wt);
      //printf("(%f,%f) \n",RRre[nn],RRim[nn]);
      printf("(%d,%d) \n",ii,jj);
      //printf("(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f) \n",data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11]);
      }*/
	  

	    
	for(nn=0;nn<nc;nn++) // gridding & self-correlation.
	  {		  
	    //index = (Ng+ii)*Nu + (Ng+jj);
	    index = (long) ii*(NN/2+1)+jj; 		  
	    GVRRre[index][nn] += RRre[nn]; 
	    GVRRim[index][nn] += RRim[nn];
	    //GVLLre[index][nn] += LLre[nn]; 
	    //GVLLim[index][nn] += LLim[nn];
	    //SCRR[index][nn] += V2RR[nn];		    	  
	    //SCLL[index][nn] += V2LL[nn];
	    
	    NTiG[index][nn] += 1.;		    	  	    	  
	  }
	      
	  time(&time_e);
	  time_d1 += difftime(time_e,time_s);
	  
	  count += 1;
	}
    }
  printf("\nTime taken in self-correlation: %d hr %d min %d sec\n", (int)(time_d/3600), (int)((time_d/60)%60),(int)(time_d%60));
  printf("\nTime taken in gridding: %d hr %d min %d sec\n", (int)(time_d1/3600), (int)((time_d1/60)%60),(int)(time_d1%60));
  printf("\nThere are %d baselines within the uv range.",count);       
  close_fits();
}






void WriteGVnAC_gridding()
{
  int ii,jj,nn;
  long index;
  FILE *fp;
  int MM = 2*Nm;
  long nc_long = (long)nc;
  printf("\nWriting Output\n");
  printf("\nSize of the output GV: %d*%d*%d.",NN,NN/2+1,nc);  
  printf("\nDu = %lf, flag = %d.\n",dU,flag_key); 

 fp=fopen(GV,"wb"); 
 fwrite(&NN,sizeof(NN),1,fp);
 fwrite(&nc_long, sizeof(nc_long), 1, fp);
 fwrite(&dU,sizeof(dU),1,fp);
 fwrite(&flag_key,sizeof(flag_key),1,fp);

  for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = (long) ii*(NN/2+1)+jj; 
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVRRre[index][nn],sizeof(GVRRre[index][nn]),1,fp);
	    fwrite(&GVRRim[index][nn],sizeof(GVRRim[index][nn]),1,fp);
	    //printf("\nRe= %e, Im= %e",GVRRre[index][nn],GVRRim[index][nn]);
	  }
      }
      
       printf("GV[0][0]=%f\n",GVRRre[0][0]);
       fclose(fp); 
       
      
       
       
    
 /* for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = ii*(NN/2+1)+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVLLre[index][nn],sizeof(GVLLre[index][jj]),1,fp);
	    fwrite(&GVLLim[index][nn],sizeof(GVLLim[index][jj]),1,fp);
	  }
      }

  //Self-correlation. (real) 
  for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = ii*(NN/2+1)+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&SCRR[index][nn],sizeof(SCRR[index][jj]),1,fp);
	    fwrite(&SCLL[index][nn],sizeof(SCLL[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",SCRR[index][nn],SCRR[index][nn]);
	  }
      }*/
      
      
 // fclose(fp); 
}


void WriteGVnAC_coll()
{
  int ii,jj,nn;
  long index;
  FILE *fp;
  int MM = 2*Nm;
  long nc_long = (long)nc_coll;
  printf("\nWriting Output\n");  
  printf("\nSize of the output GV: %d*%d*%d.\n",NN,NN/2+1,nc_coll);  
  printf("\nDu = %lf, flag = %d.\n",dU,flag_key);  
  
  
 fp=fopen(GV,"wb"); 
 fwrite(&NN,sizeof(NN),1,fp);
 fwrite(&nc_long, sizeof(nc_long), 1, fp);
 fwrite(&dU,sizeof(dU),1,fp);
 fwrite(&flag_key,sizeof(flag_key),1,fp);

  for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = (long) ii*(NN/2+1)+jj; 
	for(nn=0; nn<nc_coll; nn++)
	  {
	    fwrite(&GVRRre_coll[index][nn],sizeof(GVRRre_coll[index][nn]),1,fp);
	    fwrite(&GVRRim_coll[index][nn],sizeof(GVRRim_coll[index][nn]),1,fp);
	  }
      }
      
       fclose(fp); 
       
}



void Write_NTiG()
{
  int ii,jj,nn;
  long index;
  FILE *fp;
       
    //fp=fopen("output_data/sampling","wb");    
    fp=fopen(file_ntig,"wb"); 
   for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = (long) ii*(NN/2+1)+jj; 
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&NTiG[index][nn],sizeof(NTiG[index][nn]),1,fp);
	  }
      }
      
       fclose(fp); 
       
}



void Write_NTiG_coll()
{
  int ii,jj,nn;
  long index;
  FILE *fp;
       
    //fp=fopen("output_data/sampling","wb");    
    fp=fopen(file_ntig,"wb"); 
   for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = (long) ii*(NN/2+1)+jj; 
	for(nn=0; nn<nc_coll; nn++)
	  {
	    fwrite(&NTiG_coll[index][nn],sizeof(NTiG_coll[index][nn]),1,fp);
	  }
      }
      
       fclose(fp); 
       
}






void WriteGVnAC_gridding_realpart()
{
  int ii,jj,nn;
  long index;
  FILE *fp;
  int MM = 2*Nm;
  printf("\nSize of the output GV: %d*%d*%d.",NN,NN/2+1,nc);  

  fp=fopen(GV,"wb"); 
 // fwrite(&FWHM,sizeof(FWHM),1,fp);
 // fwrite(&dU,sizeof(dU),1,fp);
 // fwrite(&f,sizeof(f),1,fp);
 // for(ii=0;ii<nc;ii++)
   // fwrite(&nua[ii],sizeof(nua[ii]),1,fp);

  for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = (long) ii*(NN/2+1)+jj; 
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVRRre[index][nn],sizeof(GVRRre[index][nn]),1,fp);
	  }
      }
      
       printf("GV[0][0]=%f\n",GVRRre[0][0]);
       fclose(fp); 
       
}








void collapse_channels(long n_collapse)
{	
	double sumre, sumim, sumNg;
	long Ngrid, idx, col, ki, nn_col;
	Ngrid = (long)NN*(NN/2+1);
        nn_col = nc/n_collapse;
	for (idx = 0; idx < Ngrid; idx++) {
	    for (col = 0; col < nn_col; col++) {
		sumre = 0.0;
		sumim = 0.0;
		sumNg = 0.0;
		for (ki = 0; ki < n_collapse; ki++) {
		    sumre += GVRRre[idx][col * n_collapse + ki];
		    sumim += GVRRim[idx][col * n_collapse + ki];
		    sumNg += NTiG[idx][col * n_collapse + ki];
		}
		 GVRRre_coll[idx][col] = sumre;
		 GVRRim_coll[idx][col] = sumim;
		 NTiG_coll[idx][col] = sumNg;
	    }}


}











