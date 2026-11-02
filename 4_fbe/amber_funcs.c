/*=========================================================================
       This supports in computing the power spectra and bispectra from
       a convolved gridded visiblity data of a 2D non-Gaussian field generated
       from simulations. Bispectrum computation uses FFT method.
       Authors : SS Gill 
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
//#include <mpi.h>
#include "vars.h"
#include <stdalign.h>


/*-------------- Array Allocation functions -----------------*/
double **allocate_2d_double(long N1,long N2)
{
  int ii;
  long asize;
  double **phia, *phi;

  phia=(double **)calloc(N1,sizeof(double *));
  
  asize = (long)N1*N2;

  if(!(phi = (double *) calloc(asize,sizeof(double))))
    {
      fprintf(stderr,"Error in allocate_2d_double.\n");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
	phia[ii]=phi+ N2*ii;

  return(phia);
}



double **allocate_2d_double_aligned(long N1, long N2) {
    int ii;
    long asize;
    double **phia, *phi;
    int CACHE_LINE_SIZE = 64;

    // Allocate the row pointers
    if (posix_memalign((void **)&phia, CACHE_LINE_SIZE, N1 * sizeof(double *)) != 0) {
        fprintf(stderr, "Error: posix_memalign failed for row pointers.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate the 2D array with cache-line alignment
    asize = N1 * N2 * sizeof(double);
    if (posix_memalign((void **)&phi, CACHE_LINE_SIZE, asize) != 0) {
        fprintf(stderr, "Error: posix_memalign failed for array data.\n");
        exit(EXIT_FAILURE);
    }

    // Assign row pointers
    for (ii = 0; ii < N1; ++ii)
        phia[ii] = phi + N2 * ii;

    return phia;
}





float **allocate_2d_float_aligned(long N1, long N2) {
    int ii;
    long asize;
    float **phia, *phi;
    int CACHE_LINE_SIZE = 64;

    // Allocate the row pointers
    if (posix_memalign((void **)&phia, CACHE_LINE_SIZE, N1 * sizeof(float *)) != 0) {
        fprintf(stderr, "Error: posix_memalign failed for row pointers.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate the 2D array with cache-line alignment
    asize = N1 * N2 * sizeof(float);
    if (posix_memalign((void **)&phi, CACHE_LINE_SIZE, asize) != 0) {
        fprintf(stderr, "Error: posix_memalign failed for array data.\n");
        exit(EXIT_FAILURE);
    }

    // Assign row pointers
    for (ii = 0; ii < N1; ++ii)
        phia[ii] = phi + N2 * ii;

    return phia;
}


void free_2d_float_aligned(float **array) {
    if (array != NULL) {
        free(array[0]); // Free the contiguous data block
        free(array);    // Free the array of row pointers
    }
}
/*-------------------------------------------------------*/

double ***allocate_3d_double(long N1,long N2,long N3)   // dynamicL memory allocation to 3D array using pointers
{
  long ii,jj;
  long asize,index;
  double ***phia, *phi;

  phia=(double ***)calloc(N1,sizeof(double **));


  for(ii=0;ii<N1;++ii)
      phia[ii]=(double **)calloc(N2,sizeof(double *));

  asize = N1*N2;
  asize = asize*N3;

  if(!(phi = (double *) calloc(asize,sizeof(double))))
    {
      printf("error in allocate_3d_double.\n");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
    for(jj=0;jj<N2;++jj)
      {
	index = N2*N3;
	index = index*ii + N3*jj;
	phia[ii][jj]=phi+ index;   // it is address for contiguous (elements of the array are adjacent to one another in memory with no gaps between them)  memory allottment// adress increments with loop// phi is any random pointer to location, adding index to phi increment in location
      }						// it is allocting memory location to 0th element of [i][j] array (={1,2,3,4,...}, here 0th element is 1; see notebook
  return(phia);
}

/*-------------------------------------------------------*/

double P_T(double UU)// return input brightness temperature PS in K^2
{
  double ell,pu;  
  double Av=1.0, betav=-1.0; // amplitude and exponent of PS

  ell=2.0*M_PI*UU;// calculate ell  
  //pu=Av*pow(ell/1000.0,betav);  // use this if you are reading 2D_ngrf_field file
  pu=Av*pow(ell/1000.0,betav)*exp(-1.*pow(UU/500.,2));
  //pu=1.;   // use this if you are reading 2D_Ungrf_field file
  return(pu);
}


void read_inputs(char *in_file)
{
  FILE *fp;
  fp=fopen(in_file,"r"); 
  
   if (fp == NULL) {
     perror("Error opening file: input_params");
     exit(EXIT_FAILURE); }
  
  if( fscanf(fp,"%ld%lf%lf%lf%lf%ld%511s",&NN,&Du,&nu_chan0,&del_chan,&chan0,&nc,data_path) == 7){
    printf("\nInput parameters:\n\nNN=%ld\ndu=%0.2f\nnu0=%.2lf\ndeltanu=%.2lf\nchan0=%.2lf\nNchan=%ld\nData Path=%511s\n\n",NN,Du,nu_chan0,del_chan,chan0,nc,data_path);
  }
  fclose(fp);
}


void read_inputs_ess(char *in_file)
{
  FILE *fp;
  fp=fopen(in_file,"r"); 
  
   if (fp == NULL) {
     perror("Error opening file: input_params_ess");
     exit(EXIT_FAILURE); }
  
  if( fscanf(fp,"%lf%lf%lf%511s",&nu_chan0,&del_chan,&chan0,data_path) == 4){
    printf("\nInput parameters:\n\nnu0=%.2lf\ndeltanu=%.2lf\nchan0=%.2lf\nData Path=%511s\n\n",nu_chan0,del_chan,chan0,data_path);
  }
  fclose(fp);
}

void fft_init(int Nbin, double *uL,double *uH)
{
  double uu;
  double **I, **D;
  int i,j,k,m,index,ia,ja,index_1;
  unsigned long int count;
  FILE *pp,*pf;
  char fname[128];
  int nn=0;
  //int index1;
  
  I=allocate_2d_double(N1,N2+2); // for denominator
  D=allocate_2d_double(N1,N2+2); // for numerator
  
  fftw_complex *i1,*d1, *comp_ro;
  
  i1=(fftw_complex*)&(I[0][0]);
  d1=(fftw_complex*)&(D[0][0]);
 
  /*---------------------------------------------------------*/
  
  fftw_plan p1, pd1,po1;
  p1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(I[0][0]), &(I[0][0]), FFTW_ESTIMATE);
  pd1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(D[0][0]), &(D[0][0]), FFTW_ESTIMATE);
 
 //comp_ro = (fftw_complex *)&(b_ro[0][0]); // storing data visibilities 
 
  
  pf=fopen("modes","w"); // contains modes within U bins
  for(m=0;m<Nbin;++m)
  {
   count=0;
   for(i=0;i<N1;++i)
      for(j=0;j<N2/3;++j)//for(j=(N2-1)/2;j<N2;++j) 
      {
	index=i*N1+j;
	
	
  	//ia=((N1-1)/2-i);
  	//ja=((N2-1)/2-j);
  	
        ia=(i>N1/2) ? (N1-i) : i ;
        index_1=i*(N2/2+1) + j;
  	
  	//index_1=ia*(N2/2-1) + ja;
  
  	uu = sqrt(1.*(ia*ia+j*j))*Du; 
  	//uu = sqrt(1.*(ia*ia+ja*ja))*Du; 

  	if(uu>=uL[m] && uu<uH[m] && uu>0.) // checking which grid is within m-th ring
  	{
  	 d1[index_1][0]=GVRRre[index][nn]; d1[index_1][1]=GVRRim[index][nn]; // filling NG data into array
	 i1[index_1][0]=1.; i1[index_1][1]=0.;
	
  	 count++; // count modes within the bin
  	}
      }
    fftw_execute(p1); // perform FT to obtain I
    fftw_execute(pd1); // perform FT to obtain D
 
    /*--------- Writing D and I for every rings in files ----------*/
    snprintf(fname,128*sizeof(char),"masked_%d",m);
    pp=fopen(fname,"w");
    
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
           fwrite(&D[i][j],sizeof(double),1,pp);
           
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
          fwrite(&I[i][j],sizeof(double),1,pp);

    fclose(pp);
    /*-----------------------------------------------------------*/
    
    /*------- Initializing for the next bin --------*/
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
          {
            I[i][j]=0.;
            D[i][j]=0.;
	  }
    /*----------------------------------------------*/
    fprintf(pf,"%d\t%ld\n",m,count);
   }
 fclose(pf);
 fftw_destroy_plan(p1);
 fftw_destroy_plan(pd1);
}








void fft_init_orig(int nn,int Nbin, double *uL,double *uH)
{
  double uu;
  double **I, **D;
  int i,j,k,m,ia,ja;
  long index,index_1;
  int ic, jc, iagv;
  unsigned long int count;
  FILE *pp,*pf;
  char fname[128];

  I=allocate_2d_double(N1,N2+2); 
  D=allocate_2d_double(N1,N2+2);
  
  fftw_complex *i1,*d1, *comp_ro;
  i1=(fftw_complex*)&(I[0][0]);
  d1=(fftw_complex*)&(D[0][0]);
  
  fftw_plan p1, pd1;
  p1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(I[0][0]), &(I[0][0]), FFTW_ESTIMATE);
  pd1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(D[0][0]), &(D[0][0]), FFTW_ESTIMATE);

  //comp_ro = (fftw_complex *)&(ro[0][0]);
 
  pf=fopen("modes","w"); // contains modes within k bins
  for(m=0;m<Nbin;++m)
  {
   count=0;
   for(i=0;i<N1;++i)
    for(j=0;j<N2/3;++j)    //(j=0;j<N2/2+1;++j)
      {
  	ia=(i>N1/2) ? (N1-i) : i ;
        index=(long)i*(N2/2+1) + j;
  	uu=sqrt(1.*(ia*ia+j*j));
  	
  	/*ic = i+N2/2;
        jc = j+N2/2;
        iagv=(ic>=N1) ? (ic-N1) : ic ;
        index_1 = iagv*N1+jc;*/
  	
  	//igv = i+(N2/2+1);
  	//jgv = j+(N2/2+1);
  	//iagv=(igv>N1/2) ? (N1-igv) : igv ;
  	//index_1=iagv*N1+jgv;
  	
  	
  	if(uu>=uL[m] && uu<uH[m] && uu>0.)
  	{
  	 //i1[index][0]=1.; i1[index][1]=0.;
  	 i1[index][0]=NTiG[index][nn]; i1[index][1]=0.;
  	 //d1[index][0]=GVRRre[index_1][nn]; d1[index][1]=GVRRim[index_1][nn];
  	 d1[index][0]=GVRRre[index][nn]; d1[index][1]=GVRRim[index][nn];// use if center at 0,0
  	 count++; // count modes within the bin
  	}
      }
    fftw_execute(p1); // perform FT
    fftw_execute(pd1); // perform FT
    
    snprintf(fname,128*sizeof(char),"masked_%d",m);
    pp=fopen(fname,"w");
    
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
           fwrite(&D[i][j],sizeof(double),1,pp);
           
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
          fwrite(&I[i][j],sizeof(double),1,pp);
    
    fclose(pp);
    /*------- Initializing for the next bin --------*/
    for(i=0;i<N1;++i)
      for(j=0;j<N2+2;++j)
          {
            I[i][j]=0.;
            D[i][j]=0.;
          }
    /*----------------------------------------------*/
    fprintf(pf,"%d\t%ld\n",m,count);
   }
 fclose(pf);
 fftw_destroy_plan(p1);
 fftw_destroy_plan(pd1);
}








/*-------------------------------------------------------*/







/*-------------------------------------------------------*/

void calCR() {
    int r, i, j, ia, r1, r2, r3, ii, jj;
    long index;
    double tt, mu, uu, Num;
    char filename[128], filename2[128];
    double data[4];
    int rcount = 0;

    // Allocate memory
    double **AA = allocate_2d_double(N1, N2 + 2);
    fftw_complex *to_AA = (fftw_complex *)&(AA[0][0]);
    double ***II = allocate_3d_double(Nr, NN, NN);

    fftw_plan p = fftw_plan_dft_c2r_2d(N1, N2, to_AA, &(AA[0][0]), FFTW_ESTIMATE);

    // Prepare output file for writing ring data
    //snprintf(filename, sizeof(filename), "results/data_files/closed_rings");
    snprintf(filename, sizeof(filename), "%s/closed_rings",res_path);
    FILE *file_rings = fopen(filename, "wb");
    if (!file_rings) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Loop over all rings
    for (r = 0; r < Nr; ++r) {
        // Reset AA array using memset for efficiency
        memset(&(AA[0][0]), 0, sizeof(double) * N1 * (N2 + 2));

        // Populate AA with valid values
        for (i = 0; i < N1; ++i) {
            ia = (i > N1 / 2) ? (N1 - i) : i;
            for (j = 0; j < N2 / 3; ++j) {
                index = (long)i * (N2 / 2 + 1) + j;
                uu = sqrt(1.0 * (ia * ia + j * j));
                if (uu >= ULow[r] && uu < UHigh[r] && uu > 0.0) {
                    to_AA[index][0] = NTiG[index][5];
                    to_AA[index][1] = 0.0;
                }
            }
        }

        // Perform FFT
        fftw_execute(p);

        // Store result in II[r]
        for (i = 0; i < N1; ++i) {
            for (j = 0; j < N2; ++j) {
                II[r][i][j] = AA[i][j];
            }
        }
    }

    // Loop over combinations of rings
    for (r1 = 0; r1 < Nr; ++r1) {
        for (r2 = 0; r2 <= r1; ++r2) {
            tt = Uavg[r2] / Uavg[r1];
            if (tt >= 0.5 && tt <= 1.0) {
                for (r3 = 0; r3 <= r2; ++r3) {
                    mu=0.5*(tt+1./tt-(Uavg[r3]/Uavg[r1])*(Uavg[r3]/Uavg[r2]));
                    if (mu >= 0.5 / tt && mu <= 1.0) {
                        // Compute Num
                        Num = 0.0;
                        for (ii = 0; ii < N1; ++ii) {
                            for (jj = 0; jj < N2; ++jj) {
                                Num += II[r1][ii][jj] * II[r2][ii][jj] * II[r3][ii][jj];
                            }
                        }

                        // Write to file
                        data[0] = r1;
                        data[1] = r2;
                        data[2] = r3;
                        data[3] = Num;
                        fwrite(data, sizeof(double), 4, file_rings);
                        rcount++;
                        //printf("%lf %lf %lf %lf\n", data[0], data[1], data[2], Num/(N1*N2));
                    }
                }
            }
        }
    }

    fclose(file_rings);

    // Write ring count to file
    //snprintf(filename2, sizeof(filename2), "results/data_files/rcount");
    snprintf(filename2, sizeof(filename2), "%s/rcount", res_path);
    FILE *file_count = fopen(filename2, "w");
    if (!file_count) {
        fprintf(stderr, "Error opening file: %s\n", filename2);
        exit(EXIT_FAILURE);
    }
    fprintf(file_count, "%d", rcount);
    fclose(file_count);

    printf("\nTotal combinations of rings on a plane forming closed triangle = %d\n", rcount);

    // Cleanup
    fftw_destroy_plan(p);
    free(II);
    free(AA);
}







void read_rings(int *rcountPtr, double **rdataPtr) {
    char filename_count[128], filename_rings[128];

    // Read total number of configurations of closed rings
    //snprintf(filename_count, sizeof(filename_count), "results/data_files/rcount");
    snprintf(filename_count, sizeof(filename_count), "%s/rcount", res_path);
    FILE *file_count = fopen(filename_count, "r");
    if (!file_count) {
        fprintf(stderr, "Error: Unable to open file %s\n", filename_count);
        exit(EXIT_FAILURE);
    }
    if (fscanf(file_count, "%d", rcountPtr) != 1) {
        fprintf(stderr, "Error: Failed to read ring count from %s\n", filename_count);
        fclose(file_count);
        exit(EXIT_FAILURE);
    }
    fclose(file_count);

    // Read ring combination data
    //snprintf(filename_rings, sizeof(filename_rings), "results/data_files/closed_rings");
    snprintf(filename_rings, sizeof(filename_rings), "%s/closed_rings",res_path);
    FILE *file_rings = fopen(filename_rings, "rb");
    if (!file_rings) {
        fprintf(stderr, "Error: Unable to open file %s\n", filename_rings);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for ring data
    *rdataPtr = (double *)malloc(*rcountPtr * 4 * sizeof(double));
    if (!*rdataPtr) {
        fprintf(stderr, "Error: Memory allocation failed for ring data\n");
        fclose(file_rings);
        exit(EXIT_FAILURE);
    }

    // Read the ring data into the allocated memory
    if (fread(*rdataPtr, sizeof(double), *rcountPtr * 4, file_rings) != *rcountPtr * 4) {
        fprintf(stderr, "Error: Failed to read ring data from %s\n", filename_rings);
        fclose(file_rings);
        free(*rdataPtr);
        exit(EXIT_FAILURE);
    }
    fclose(file_rings);

    // Confirmation (optional, for debugging)
    printf("Successfully read %d ring configurations.\n", *rcountPtr);
}



/*-------------------  Df ----------------------*/

void Df(int plane,int ring, double **DD, fftw_plan pd)   
  	{
  	int i, j, ia;
  	double uu;
  	long index;
  	fftw_complex *to_DD, *comp_ro;
	to_DD=(fftw_complex*)&(DD[0][0]);
	memset(&(DD[0][0]), 0, sizeof(double) * N1 * (N2 + 2));
	for (i=0;i<N1;++i)
  		for (j=0;j<N2/3;++j)
  		{	
			ia=(i>N1/2) ? (N1-i) : i ;
        		index=(long)i*(N2/2+1) + j;
  			uu=sqrt(1.*(ia*ia+j*j));	
			if (uu >= ULow[ring] && uu < UHigh[ring] && uu > 0.0) {
			to_DD[index][0]=GVRRre[index][plane]; //*RR[ring][index]
			to_DD[index][1]=GVRRim[index][plane];}
		}

	fftw_execute(pd);	     
	}
	
void Dfg(int plane,int ring, double **DD, fftw_plan pd)   
  	{
  	int iN, jN, ia, j, iNg, jNg;
  	double uu;
  	long indexN,indexNg;
  	fftw_complex *to_DD, *comp_ro;
	to_DD=(fftw_complex*)&(DD[0][0]);
	memset(&(DD[0][0]), 0, sizeof(double) * Ng1 * (Ng2 + 2));
	for (ia=-Ng1/2+1;ia<Ng1/2+1;++ia)
  		for (j=0;j<Ng2/3;++j)   //for (j=0;j<Ng2/2+1;++j)-> this also works but loop is unnencessarily large
  		{	
			iN   = (ia<0) ? (N1+ia) : ia ;
			iNg  = (ia<0) ? (Ng1+ia) : ia ;
        		indexN= (long)iN*(N2/2+1) + j;
        		indexNg= (long)iNg*(Ng2/2+1) + j;
  			uu=sqrt(1.*(ia*ia+j*j));	
			if (uu >= ULow[ring] && uu < UHigh[ring] && uu > 0.0) {
			to_DD[indexNg][0]=GVRRre[indexN][plane]; 
			to_DD[indexNg][1]=GVRRim[indexN][plane];}
		}

	fftw_execute(pd);	     
	}
	
	
	
void Ifg(int plane,int ring, double **II, fftw_plan pd)   
  	{
  	int iN, jN, ia, j, iNg, jNg;
  	double uu;
  	long indexN,indexNg;
  	fftw_complex *to_II, *comp_ro;
	to_II=(fftw_complex*)&(II[0][0]);
	memset(&(II[0][0]), 0, sizeof(double) * Ng1 * (Ng2 + 2));
	for (ia=-Ng1/2+1;ia<Ng1/2+1;++ia)
  		for (j=0;j<Ng2/3;++j)   //for (j=0;j<Ng2/2+1;++j)-> this also works but loop is unnencessarily large
  		{	
			iN   = (ia<0) ? (N1+ia) : ia ;
			iNg  = (ia<0) ? (Ng1+ia) : ia ;
        		indexN= (long)iN*(N2/2+1) + j;
        		indexNg= (long)iNg*(Ng2/2+1) + j;
  			uu=sqrt(1.*(ia*ia+j*j));	
			if (uu >= ULow[ring] && uu < UHigh[ring] && uu > 0.0) {
			to_II[indexNg][0]=NTiG[indexN][plane]; //*RR[ring][index]
			to_II[indexNg][1]=0.;}
		}

	fftw_execute(pd);	     
	}
	
/*------------------- Wrapid ----------------------*/	
void wrapid()  // writes rings and planes ift of data
{
// writing IFT of each rings on each plane field in binary files i.e D

int p1,r1,ii,jj;
long index;
double **D;
double *buffer = (double *)malloc(N1 * N2 * sizeof(double));

D=allocate_2d_double(N1,N2+2); 
fftw_complex *to_D=(fftw_complex*)&(D[0][0]);
fftw_plan pd = fftw_plan_dft_c2r_2d(N1, N2, to_D, &(D[0][0]), FFTW_ESTIMATE);
char dfilename[128];

for(p1=0;p1<nc;++p1){
 for(r1=0;r1<Nr;++r1){
   Df(p1,r1,D,pd); 
   snprintf(dfilename,sizeof(dfilename), "results/data_files/dataPlaneRing_%d_%d_R%d",p1,r1,NRi);
   FILE *dfile= fopen(dfilename, "wb");
   
   for (ii = 0; ii < N1; ++ii) { // this gives same as element-wise buffer filling. But it is fast.
    memcpy(buffer + ii * N2, &(D[ii][0]), N2 * sizeof(double));}
	
   fwrite(buffer, sizeof(double), N1*N2, dfile);
   fclose(dfile);}}
   
fftw_destroy_plan(pd);
free(D);
free(buffer);

/*------------For future reference --------------*/
/* for(ii=0;ii<NN;++ii){ // this is the correct way to write Ds in file. But it is bit slow.
     for(jj=0;jj<NN;++jj){
	index=(long)ii*(NN)+ jj;
	buffer[index] = D[ii][jj];//to_D[index][0];
        //buffer[index + 1] = to_D[index][1];
	}}

//fwrite(&(D[0][0]), sizeof(double), N1 * N2, dfile);// this is wrong

for (ii = 0; ii < N1; ++ii) {  // it is fine and gives correct results. But fwrite operations are many
    fwrite(&(D[ii][0]), sizeof(double), N2, dfile);}*/
}

/*------------------- Wrapid grid reduced ----------------------*/	
void wrapidg()  // writes rings and planes ift of data
{
// writing IFT of each rings on each plane field in binary files i.e D

int p1,r1,ii,jj;
long index;
double **D;
double *buffer = (double *)malloc(Ng1 * Ng2 * sizeof(double));

D=allocate_2d_double(Ng1,Ng2+2); 
fftw_complex *to_D=(fftw_complex*)&(D[0][0]);
fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_D, &(D[0][0]), FFTW_ESTIMATE);
char dfilename[128];

for(p1=0;p1<nc;++p1){
 for(r1=0;r1<Nr;++r1){
   Dfg(p1,r1,D,pd); 
   snprintf(dfilename,sizeof(dfilename), "results/data_files/dataPlaneRing_%d_%d_R%d",p1,r1,NRi);
   FILE *dfile= fopen(dfilename, "wb");
   
   for (ii = 0; ii < Ng1; ++ii) { // this gives same as element-wise buffer filling. But it is fast.
    memcpy(buffer + ii * Ng2, &(D[ii][0]), Ng2 * sizeof(double));}
	
   fwrite(buffer, sizeof(double), Ng1*Ng2, dfile);
   fclose(dfile);}}
   
fftw_destroy_plan(pd);
free(D);
free(buffer);
}

/*------------read wrapid --------------*/
void get_D(int plane, int ring, int rlz, double **D) {
    char filename[128];
    snprintf(filename,sizeof(filename), "results/data_files/dataPlaneRing_%d_%d_R%d",plane,ring,rlz);
    FILE *dfile = fopen(filename, "rb");  // Open the file for reading in binary mode
    if (!dfile) {
        fprintf(stderr, "Error: Could not open file %s for reading.\n", filename);
        exit(1);
    }

    // Read N1 * N2 doubles from the file into D
    size_t result = fread(&(D[0][0]), sizeof(double), N1 * N2, dfile);
    if (result != N1 * N2) {
        fprintf(stderr, "Error: Could not read enough data from the file.\n");
        fclose(dfile);
        exit(1);
    }

    fclose(dfile);  // Close the file after reading
}

/*------------read wrapid grid reduced --------------*/
void get_Dg(int plane, int ring, int rlz, double **D) {
    char filename[128];
    snprintf(filename,sizeof(filename), "results/data_files/dataPlaneRing_%d_%d_R%d",plane,ring,rlz);
    FILE *dfile = fopen(filename, "rb");  // Open the file for reading in binary mode
    if (!dfile) {
        fprintf(stderr, "Error: Could not open file %s for reading.\n", filename);
        exit(1);
    }

    // Read N1 * N2 doubles from the file into D
    size_t result = fread(&(D[0][0]), sizeof(double), Ng1 * Ng2, dfile);
    if (result != Ng1 * Ng2) {
        fprintf(stderr, "Error: Could not read enough data from the file.\n");
        fclose(dfile);
        exit(1);
    }

    fclose(dfile);  // Close the file after reading
}

// Function to allocate a 3D array dynamically
double ***allocate_3d_array(long nc, long Ng1, long Ng2) {
    double ***array = (double ***)malloc(nc * sizeof(double **));
    if (!array) {
        fprintf(stderr, "Memory allocation failed for array\n");
        exit(EXIT_FAILURE);
    }

    for (long i = 0; i < nc; i++) {
        array[i] = (double **)malloc(Ng1 * sizeof(double *));
        if (!array[i]) {
            fprintf(stderr, "Memory allocation failed for array[%ld]\n", i);
            exit(EXIT_FAILURE);
        }

        for (long j = 0; j < Ng1; j++) {
            array[i][j] = (double *)malloc(Ng2 * sizeof(double));
            if (!array[i][j]) {
                fprintf(stderr, "Memory allocation failed for array[%ld][%ld]\n", i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
    return array;
}


void free_3d_array(double ***array, long nc, long Ng1) {
    for (long i = 0; i < nc; i++) {
        for (long j = 0; j < Ng1; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}



// Function to allocate 4D array dynamically
double ****allocate_4d_array(long nc, long Nr, long Ng1, long Ng2) {
    double ****array = (double ****)malloc(nc * sizeof(double ***));
    for (int i = 0; i < nc; i++) {
        array[i] = (double ***)malloc(Nr * sizeof(double **));
        for (int j = 0; j < Nr; j++) {
            array[i][j] = (double **)malloc(Ng1 * sizeof(double *));
            for (int k = 0; k < Ng1; k++) {
                array[i][j][k] = (double *)malloc(Ng2 * sizeof(double));
            }
        }
    }
    return array;
}


void free_2d_double(double **array, int rows) {
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

void free_2d_double_v2(double **array) {
    if (array != NULL) {
        free(array[0]);  // Free the contiguous block of data
        free(array);     // Free the array of row pointers
    }
}



// Function to free 4D array
void free_4d_array(double ****array) {
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < Nr; j++) {
            for (int k = 0; k < Ng1; k++) {
                free(array[i][j][k]);
            }
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

void fill_s4d(double ****sd, int rlz, int read_or_cal) {
    int p1, r1, ii, jj;
    
    if (read_or_cal == 0) {
        double **D = allocate_2d_double(Ng1, Ng2);
        
        for (p1 = 0; p1 < nc; ++p1) {
            for (r1 = 0; r1 < Nr; ++r1) {
                get_Dg(p1, r1, rlz, D);
                
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        sd[p1][r1][ii][jj] = D[ii][jj];
                    }
                }
            }
        }
        free(D);
    } 
    else if (read_or_cal == 1) {
        double **D = allocate_2d_double(Ng1, Ng2 + 2);
        fftw_complex *to_D = (fftw_complex*)&(D[0][0]);
        fftw_plan_with_nthreads(Nthreads);
        fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_D, &(D[0][0]), FFTW_ESTIMATE);
        
        for (p1 = 0; p1 < nc; ++p1) {
            for (r1 = 0; r1 < Nr; ++r1) {
                Dfg(p1, r1, D, pd);
                
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        sd[p1][r1][ii][jj] = D[ii][jj];
                    }
                }
            }
        }
        
        fftw_destroy_plan(pd);
        free(D);
    } 
    else {
        fprintf(stderr, "Error: Give either 0 to read or 1 to calculate\n");
    }
}





void fill_s3d(double ***sd, int ring, int rlz, int read_or_cal) {
    int p1, ii, jj;
    
    if (read_or_cal == 0) {
        double **D = allocate_2d_double(Ng1, Ng2);
        
        for (p1 = 0; p1 < nc; ++p1) {
                get_Dg(p1, ring, rlz, D);
                
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        sd[p1][ii][jj] = D[ii][jj];
                    }
                }
        }
        free(D);
    } 
    else if (read_or_cal == 1) {
        double **D = allocate_2d_double(Ng1, Ng2 + 2);
        fftw_complex *to_D = (fftw_complex*)&(D[0][0]);
        fftw_plan_with_nthreads(Nthreads);
        fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_D, &(D[0][0]), FFTW_ESTIMATE);
        
        for (p1 = 0; p1 < nc; ++p1) {
                Dfg(p1, ring, D, pd);
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        sd[p1][ii][jj] = D[ii][jj];
                    }
                }
            
        }
        
        fftw_destroy_plan(pd);
        free(D);
    } 
    else {
        fprintf(stderr, "Error: Give either 0 to read or 1 to calculate\n");
    }
}

void fill_s2d(double **sd, int ring, int rlz, int read_or_cal) {
    int p1, ii, jj;
    long index;
    
    if (read_or_cal == 0) {
        double **D = allocate_2d_double(Ng1, Ng2);
        
        for (p1 = 0; p1 < nc; ++p1) {
                get_Dg(p1, ring, rlz, D);
                
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = D[ii][jj];
                    }
                }
        }
        free(D);
    } 
    else if (read_or_cal == 1) {
        double **D = allocate_2d_double(Ng1, Ng2 + 2);
        fftw_complex *to_D = (fftw_complex*)&(D[0][0]);
        fftw_plan_with_nthreads(Nthreads);
        fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_D, &(D[0][0]), FFTW_ESTIMATE);
        
        for (p1 = 0; p1 < nc; ++p1) {
                Dfg(p1, ring, D, pd);
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = D[ii][jj];
                    }
                }
            
        }
        
        fftw_destroy_plan(pd);
        free(D);
    } 
    else {
        fprintf(stderr, "Error: Give either 0 to read or 1 to calculate\n");
    }
}







void fill_s2d_float(float **sd, int ring, int rlz, int read_or_cal) {
    int p1, ii, jj;
    long index;
    
    if (read_or_cal == 0) {
        double **D = allocate_2d_double(Ng1, Ng2);
        
        for (p1 = 0; p1 < nc; ++p1) {
                get_Dg(p1, ring, rlz, D);
                
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = (float) D[ii][jj];
                    }
                }
        }
        free(D);
    } 
    else if (read_or_cal == 1) {
        double **D = allocate_2d_double(Ng1, Ng2 + 2);
        fftw_complex *to_D = (fftw_complex*)&(D[0][0]);
        fftw_plan_with_nthreads(Nthreads);
        fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_D, &(D[0][0]), FFTW_ESTIMATE);
        
        for (p1 = 0; p1 < nc; ++p1) {
                Dfg(p1, ring, D, pd);
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = (float) D[ii][jj];
                    }
                }
            
        }
        
        fftw_destroy_plan(pd);
        //free(D);
        free_2d_double_v2(D);
    } 
    else {
        fprintf(stderr, "Error: Give either 0 to read or 1 to calculate\n");
    }
}


void fill_s2d_float_adaptive(float **sd, int ring, int rlz, int read_or_cal, int iplane) {
    int p1, ii, jj;
    long index;
    //printf("######## Ng1=%d\n#######",Ng1);
    
    if (read_or_cal == 0) {
        double **D = allocate_2d_double(Ng1, Ng2);
        
        for (p1 = 0; p1 <= nc10; ++p1) {
                get_Dg(p1+iplane, ring, rlz, D);
                
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = (float) D[ii][jj];
                    }
                }
        }
        free(D);
    } 
    else if (read_or_cal == 1) {
        double **D = allocate_2d_double(Ng1, Ng2 + 2);
        fftw_complex *to_D = (fftw_complex*)&(D[0][0]);
        fftw_plan_with_nthreads(Nthreads);
        fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_D, &(D[0][0]), FFTW_ESTIMATE);
        
        for (p1 = 0; p1 <= nc10; ++p1) {
                Dfg(p1+iplane, ring, D, pd);
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = (float) D[ii][jj];
                    }
                }
            
        }
        
        fftw_destroy_plan(pd);
        //free(D);
        free_2d_double_v2(D);
    } 
    else {
        fprintf(stderr, "Error: Give either 0 to read or 1 to calculate\n");
    }
}



void fill_s2d_float_unflg(float **sd, int ring, int rlz, int read_or_cal) {
    int p1, ii, jj, p1_i;
    long index;
    //printf("######## Ng1=%d\n#######",Ng1);
    
    if (read_or_cal == 0) {
        double **D = allocate_2d_double(Ng1, Ng2);
        
        for (p1_i = 0; p1_i < uf_nc; ++p1_i) {
                p1 = uf_chans[p1_i];
                get_Dg(p1, ring, rlz, D);
                
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = (float) D[ii][jj];
                    }
                }
        }
        free(D);
    } 
    else if (read_or_cal == 1) {
        double **D = allocate_2d_double(Ng1, Ng2 + 2);
        fftw_complex *to_D = (fftw_complex*)&(D[0][0]);
        fftw_plan_with_nthreads(Nthreads);
        fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_D, &(D[0][0]), FFTW_ESTIMATE);
        
        for (p1_i = 0; p1_i < uf_nc; ++p1_i) {
        	p1 = uf_chans[p1_i];
                Dfg(p1, ring, D, pd);
                for (ii = 0; ii < Ng1; ++ii) {
                    for (jj = 0; jj < Ng2; ++jj) {
                        index = ii * Ng2 + jj;  // Flattening 2D to 1D
                        sd[p1][index] = (float) D[ii][jj];
                    }
                }
            
        }
        
        fftw_destroy_plan(pd);
        free_2d_double_v2(D);
    } 
    else {
        fprintf(stderr, "Error: Give either 0 to read or 1 to calculate\n");
    }
}


void fill_s2I(float **sI, int ring) {
    int p1, ii, jj;
    long index;

	double **I = allocate_2d_double(Ng1, Ng2 + 2);
	fftw_complex *to_I = (fftw_complex*)&(I[0][0]);
	fftw_plan_with_nthreads(Nthreads);
	fftw_plan pd = fftw_plan_dft_c2r_2d(Ng1, Ng2, to_I, &(I[0][0]), FFTW_ESTIMATE);

	for (p1 = 0; p1 < nc; ++p1) {
		Ifg(p1, ring, I, pd);
		for (ii = 0; ii < Ng1; ++ii) {
		    for (jj = 0; jj < Ng2; ++jj) {
		        index = ii * Ng2 + jj;  // Flattening 2D to 1D
		        sI[p1][index] = (float) I[ii][jj];
		    }
		}
	    
	}

	fftw_destroy_plan(pd);
	free_2d_double_v2(I);
}

//=======================================================================//
// Function to calculate the power spectrum from Fourier modes
void power_spec(int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
  fftw_complex *in, *in_u;
  
  in=(fftw_complex*)&b_ro[0][0];
  in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  double binsiz,uu, numer, denom;
  int i,j,ia,tmp;
  long index;

  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
  
  binsiz=log10(umax/umin)/(1.*NB);

  for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=(long)i*(N2/2+1)+j;
  	ia=(i>N1/2) ? (N1-i) : i ;
  
  	uu = sqrt(1.*(ia*ia+j*j))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	    denom = in_u[index][0]*in_u[index][0]+in_u[index][1]*in_u[index][1];
  	    if(denom >= 1.)
  	    {
  	      numer = in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	      tmp=(int) floor(log10(uu/umin)/binsiz);
  	      no[tmp]++;
	      pksum[tmp]+=(numer/denom);
  	      usum[tmp]+=uu;
  	    }
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/




void power_spec_GVRL(int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
 // fftw_complex *in, *in_u;
  
 // in=(fftw_complex*)&b_ro[0][0];
  //in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  
  double binsiz,uu, numer, denom, re_part, im_part;
  int i,j,ia,ja,tmp;
  long index;
  int nn=0;
  int MM = 2*Nm;
  
  int Nlim;
  
  Nlim=Nu-2*MM;
   
  printf("\n Nlim=%d", Nlim); 
  printf("\n hiii");
  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
  
  binsiz=log10(umax/umin)/(1.*NB);
  printf("\n %f",binsiz);
  for(i=0;i<Nlim;++i)
    for(j=0;j<Nlim;++j)
      {
  	index=(long)i*Nlim+j;//i*((Nu-2*MM)/2+1)+j;
  	//ia=(i>(Nlim-1)/2) ? ((Nlim-1)/2-i) : i ;
  	//ja=(j>(Nlim-1)/2) ? ((Nlim-1)/2-j) : j ;
  	
  	ia=((Nlim-1)/2-i);
  	ja=((Nlim-1)/2-j);
  
  	uu = sqrt(1.*(ia*ia+ja*ja))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	     re_part = GVRRre[index][nn] + GVLLre[index][nn];
  	     im_part = GVRRim[index][nn] + GVLLim[index][nn];
  	     numer =  re_part*re_part+im_part*im_part;//in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	     tmp=(int) floor(log10(uu/umin)/binsiz);
  	     no[tmp]++;
	     pksum[tmp]+=(numer);
	     //if (GVRRre[index][nn]!=0 ){
	     //printf("\n %d",tmp);}
	     //if (tmp<5){
	     //printf("\n %d",tmp);}
  	     usum[tmp]+=uu;
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/













void power_spec_GV(int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
 // fftw_complex *in, *in_u;
  
 // in=(fftw_complex*)&b_ro[0][0];
  //in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  
  double binsiz,uu, numer, denom;
  int i,j,ia,ja,tmp;
  long index;
  int nn=0;
  int MM = 2*Nm;
  
  int Nlim;
  
  Nlim=Nu-2*MM;
   
  printf("\n Nlim=%d", Nlim); 
  printf("\n hiii");
  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
  
  binsiz=log10(umax/umin)/(1.*NB);
  printf("\n %f",binsiz);
  for(i=0;i<Nlim;++i)
    for(j=0;j<Nlim;++j)
      {
  	index=(long)i*Nlim+j;//i*((Nu-2*MM)/2+1)+j;
  	//ia=(i>(Nlim-1)/2) ? ((Nlim-1)/2-i) : i ;
  	//ja=(j>(Nlim-1)/2) ? ((Nlim-1)/2-j) : j ;
  	
  	ia=((Nlim-1)/2-i);
  	ja=((Nlim-1)/2-j);
  
  	uu = sqrt(1.*(ia*ia+ja*ja))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	     numer =  GVRRre[index][nn]*GVRRre[index][nn]+GVRRim[index][nn]*GVRRim[index][nn];//in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	     tmp=(int) floor(log10(uu/umin)/binsiz);
  	     no[tmp]++;
	     pksum[tmp]+=(numer);
	     //if (GVRRre[index][nn]!=0 ){
	     //printf("\n %d",tmp);}
	     //if (tmp<5){
	     //printf("\n %d",tmp);}
  	     usum[tmp]+=uu;
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/

// following works very fine with gridded visibility data//
void power_spec_GV_Ntig(int nn, int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
 // fftw_complex *in, *in_u;
  
 // in=(fftw_complex*)&b_ro[0][0];
  //in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  
  double binsiz,uu, numer, denom;
  int i,j,index,ia,ja,tmp;
  int MM = 2*Nm;
  
  int Nlim;
  
  Nlim=Nu-2*MM;
   
  printf("\n Nlim=%d", Nlim); 
  printf("\n hiii");
  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
    
    
      /*for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=i*(N2/2+1)+j;
  	ia=(i>N1/2) ? (N1-i) : i ;*/
    
    
  
  binsiz=log10(umax/umin)/(1.*NB);
  printf("\n %f",binsiz);
  for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=i*(N2/2+1)+j;
  	
	ia=(i>N1/2) ? (N1-i) : i ;
  
  	uu = sqrt(1.*(ia*ia+j*j))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	     numer =  GVRRre[index][nn]*GVRRre[index][nn]+GVRRim[index][nn]*GVRRim[index][nn];//in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	     numer *= NTiG[index][nn];
  	     tmp=(int) floor(log10(uu/umin)/binsiz);
  	     no[tmp]+= NTiG[index][nn];
	     pksum[tmp]+=(numer);
	     //if (GVRRre[index][nn]!=0 ){
	     //printf("\n %d",tmp);}
	     //if (tmp<5){
	     //printf("\n %d",tmp);}
  	     usum[tmp]+=(uu*NTiG[index][nn]);
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/


void cl_dnu(int nn1, int nn2, int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
 // fftw_complex *in, *in_u;
  
 // in=(fftw_complex*)&b_ro[0][0];
  //in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  
  double binsiz,uu, numer, denom;
  int i,j,index,ia,ja,tmp;
  int MM = 2*Nm;
  
  int Nlim;
  
  Nlim=Nu-2*MM;
   
  //printf("\n Nlim=%d", Nlim); 

  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
    
    
      /*for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=i*(N2/2+1)+j;
  	ia=(i>N1/2) ? (N1-i) : i ;*/
    
    
  
  binsiz=log10(umax/umin)/(1.*NB);
  //printf("\n %f",binsiz);
  for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=i*(N2/2+1)+j;
  	
	ia=(i>N1/2) ? (N1-i) : i ;
  
  	uu = sqrt(1.*(ia*ia+j*j))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	     numer =  GVRRre[index][nn1]*GVRRre[index][nn2]+GVRRim[index][nn1]*GVRRim[index][nn2];//in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	     numer *= NTiG[index][nn1];
  	     tmp=(int) floor(log10(uu/umin)/binsiz);
  	     no[tmp]+= NTiG[index][nn1];
	     pksum[tmp]+=(numer);
	     //if (GVRRre[index][nn]!=0 ){
	     //printf("\n %d",tmp);}
	     //if (tmp<5){
	     //printf("\n %d",tmp);}
  	     usum[tmp]+=(uu*NTiG[index][nn1]);
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/







/*-------- Computing shapes for various sizes of closed triangles ---------*/
void shapes(int nbin, double *uavg)
{
  int i,j,l; 
  double MU,TT, uu1;
  char file[128];
  FILE *fp1;
 /*----- Evaluating all possible unique shapes of closed triangles -----*/  
  for(i=0;i<nbin;++i)
  {
   uu1=uavg[i];
   
   snprintf(file, 128*sizeof(char),"shape_u%d_%.3lf",i,uu1);
   fp1=fopen(file,"w");
   
   for(j=0;j<=i;++j)
    {
     TT=uavg[j]/uavg[i];
     if(TT>=0.5 && TT<=1.)
     {
      for(l=0;l<=j;++l)
      {
       MU=0.5*(TT+1./TT-(uavg[l]/uavg[i])*(uavg[l]/uavg[j]));
       if(MU>=0.5/TT && MU<=1.) // this incorporates 2*mu*t>=1
          fprintf(fp1,"%e\t%e\n",TT,MU); // writing output
      }
     }
    }
   fclose(fp1);
  } 
}

/*--------------------------------------------------------------------------*/

/*****************************************************/
/* Print out cfitsio error messages and exit program */
/*****************************************************/
void printerror( int status)
{
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}
/*--------------------------------------------------------------------------*/


/************************************************************************/
/*           Read a UV data in Multifrequency FITS format               */
/*         This module reads both the NG data and UAPS data             */
/************************************************************************/
void read_uvdata(char filename1[], char filename2[])
{
    fitsfile *fptr, *fptr1;       /* pointer to the FITS file, defined in fitsio.h */
    int status, nfound, anynull;
    long naxes[3], fpixel[3], lpixel[3], inc[3]={1,1,1}, nbuffer, npixels, ii, jj, index, index1;
    long Nu, Nv, Nz, Nvis;
    double *buffer1, *buffer2;
    double cdelt[3], crpix[3], crval[3];
    void *nullval;
    int ia, ja;
    
    fftw_complex *in, *in_u, *in1, *in_u1;
    
    if(access(filename1, F_OK)!=0)
    {
      printf("Input File %s does not exists\n",filename1);
      exit(0);
    }
  
    if(access(filename2, F_OK)!=0)
    {
      printf("Input File %s does not exists\n",filename2);
      exit(0);
    }
    
    /*------------- Reading NG data file --------------*/
    status = 0;

    if ( fits_open_file(&fptr, filename1, READONLY, &status) )
         printerror( status );

    /* read the NAXIS1, NAXIS2 and NAXIS3 keywords to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 3, naxes, &nfound, &status) )
         printerror( status );
    
    /* read the CDELT1, CDELT2 and CDELT3 keywords to get data size */
    if ( fits_read_keys_dbl(fptr, "CDELT", 1, 3, cdelt, &nfound, &status) )
         printerror( status );
         
    if ( fits_read_keys_dbl(fptr, "CRVAL", 1, 3, crval, &nfound, &status) )
         printerror( status );
    
    if ( fits_read_keys_dbl(fptr, "CRPIX", 1, 3, crpix, &nfound, &status) )
         printerror( status );     

    Nu = naxes[0]; Nv = naxes[1]; Nz = naxes[2];          /* number of grids along the three sides */
    Nvis  = Nu*Nv;         /* number of visiblities in one channel */
    //nullval  = 0;                /* don't check for null values in the image */
    
    printf("dU=%lf naxes vis={%ld,%ld,%ld}\n",cdelt[0],naxes[0],naxes[1],naxes[2]);
    
    //memory allocation for rho[][] array
    N1 = Nu;
    N2 = Nu;
   
    Du = fabs(cdelt[0]); // grid size in Fourier space in the unit of lambda
    
    ro = allocate_2d_double(N1,N2+1); //   ro = allocate_2d_double(N1,N2+1);,N2 is odd
    ro_u = allocate_2d_double(N1,N2+1);
    
    b_ro = allocate_2d_double(N1,N2+1);
    b_ro_u = allocate_2d_double(N1,N2+1);
    
    in = (fftw_complex*)&ro[0][0];
    in_u = (fftw_complex*)&ro_u[0][0];
    
    in1 = (fftw_complex*)&b_ro[0][0];
    in_u1 = (fftw_complex*)&b_ro_u[0][0];
    	
    buffer1 = (double*)calloc(Nvis,sizeof(double));
    buffer2 = (double*)calloc(Nvis,sizeof(double));
    
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=1;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=1;
    
    if( fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nullval,buffer1,&anynull,&status) )
         printerror( status );
         
         
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=2;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=2;
    
    if( fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nullval,buffer2,&anynull,&status) )
         printerror( status );
         
         
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2/2+1;jj++)
        {
         index = (long)jj*N1+ii;
         in[index][0] = buffer1[index];
         in[index][1] = buffer2[index];
        }
    
   /*------------------ Shifting the U-plane origin for FFTW ------------------*/
    for(ii=0;ii<N1;++ii)
    {
      ia = ii-N1/2;
      ja = (ia<0) ? (ia+N1) : ia;
      for(jj=0;jj<N2/2+1;++jj)
          {
           index = (long)jj*N1 + ii;
           index1 = (long)jj + ja*(N2/2+1);
           in1[index1][0]=in[index][0];
           in1[index1][1]=in[index][1];
          }
     }
  /*---------------------------------------------------------*/
        
    free(buffer1);
    free(buffer2);
    if ( fits_close_file(fptr, &status) )
         printerror( status );
   /*------------- Reading UAPS data file --------------*/
   
   status = 0;

    if ( fits_open_file(&fptr1, filename2, READONLY, &status) )
         printerror( status );

    
    buffer1 = (double*)calloc(Nvis,sizeof(double));
    buffer2 = (double*)calloc(Nvis,sizeof(double));
    
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=1;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=1;
    
    if( fits_read_subset(fptr1,TDOUBLE,fpixel,lpixel,inc,nullval,buffer1,&anynull,&status) )
         printerror( status );
         
         
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=2;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=2;
    
    if( fits_read_subset(fptr1,TDOUBLE,fpixel,lpixel,inc,nullval,buffer2,&anynull,&status) )
         printerror( status );
    
    
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2/2+1;jj++)
        {
         index = (long)jj*N1+ii;
         in_u[index][0] = buffer1[index];
         in_u[index][1] = buffer2[index];
        }
        
    /*------------------ Shifting the U-plane origin for FFTW------------------*/
    for(ii=0;ii<N1;++ii)
    {
      ia = ii-N1/2;
      ja = (ia<0) ? (ia+N1) : ia;
      for(jj=0;jj<N2/2+1;++jj)
          {
           index = (long)jj*N1 + ii;
	   index1 = (long)jj + ja*(N2/2+1);
           in_u1[index1][0]=in_u[index][0];
           in_u1[index1][1]=in_u[index][1];
          }
     }
  /*---------------------------------------------------------*/
    free(buffer1);
    free(buffer2);
    
    if ( fits_close_file(fptr1, &status) )
         printerror( status );
}
/*--------------------------------------------------------------------------*/


int initialize()
{
  int ii,jj,nn,sf;
  long Ngrid; 

  
  ULow=(double*)calloc(Nr,sizeof(double));
  UHigh=(double*)calloc(Nr,sizeof(double));
  Uavg=(double*)calloc(Nr,sizeof(double));
  uf_chans =(int*)calloc(uf_nc,sizeof(int));
  
  char file_GV[512];
  FILE *fp;
  snprintf(file_GV, 512*sizeof(char),"%s/GV_data_%d",data_path,NRi);
  fp=fopen(file_GV,"rb"); 
  if (fp == NULL) {perror("Error opening header file: GV_data");exit(1);}
  fread(&NN,sizeof(NN),1,fp);
  fread(&nc, sizeof(nc), 1, fp);
  fread(&Du,sizeof(Du),1,fp);
  fread(&flag_key,sizeof(flag_key),1,fp);
  fclose(fp);
  printf("\nReading Input parameters of GV: %ld*%ld*%ld\t Du = %lf, flag = %d.",NN,NN/2+1,nc,Du,flag_key);   
  
  N1 = NN;
  N2 = NN;                
  Ngrid = (long)N1*(N2/2+1);//Nu*Nu;
  printf("\n Ngrid = %d",Ngrid);

  nua = (double*)calloc(nc,sizeof(double));
  for(nn=0;nn<nc;nn++)
    {
      nua[nn]= (nu_chan0 + del_chan*(nn - (chan0-1)));
      //nua[nn]=nu0+(nn+1.-chan0)*deltanu;
    }
                    
  
  GVRRre = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRim = (double**)calloc(Ngrid,sizeof(double*));  
  //GVLLre = (double**)calloc(Ngrid,sizeof(double*));  
  //GVLLim = (double**)calloc(Ngrid,sizeof(double*));  
  NTiG   = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      GVRRre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVRRim[ii]=(double*)calloc(nc,sizeof(double)); 
      //GVLLre[ii]=(double*)calloc(nc,sizeof(double)); 
      //GVLLim[ii]=(double*)calloc(nc,sizeof(double)); 
      NTiG[ii]=(double*)calloc(nc,sizeof(double));
     
    }

  /*SCRR = (double**)calloc(Ngrid,sizeof(double*));  
  SCLL = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      SCRR[ii]=(double*)calloc(nc,sizeof(double)); 
      SCLL[ii]=(double*)calloc(nc,sizeof(double)); 
    }*/

}







/*
void ReadGVnAC()
{ 
  int ii,jj,nn, index;
  char file_GV[256];
  FILE *fp;
  
  int MM = 2*Nm, Nlim;
  
  Nlim=Nu-2*MM;
 // printf("\nSize of the output GV: %d*%d*%d.",Nu-2*MM,Nu-2*MM,nc);  



   snprintf(file_GV, 256*sizeof(char),"../3_vcg/output_data/%s",GV); // open file to write bispectrum
   fp=fopen(file_GV,"rb");
  //fp=fopen(GV,"rb"); 
  fread(&FWHM,sizeof(FWHM),1,fp);
  fread(&dUx,sizeof(dUx),1,fp);
  fread(&fx,sizeof(fx),1,fp);
  for(ii=0;ii<nc;ii++)
    fread(&nua[ii],sizeof(nua[ii]),1,fp);

  for(ii=0; ii<Nlim;ii++)
    for(jj=0; jj<Nlim;jj++)
      {
	index=ii*Nlim+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&GVRRre[index][nn],sizeof(GVRRre[index][jj]),1,fp);
	    fread(&GVRRim[index][nn],sizeof(GVRRim[index][jj]),1,fp);
	   // printf("\ni= %d, j= %d",ii,jj);
	   // printf("\nRe= %e, Im= %e",GVRRre[index][nn],GVRRim[index][nn]);
	  }
      }
      
      
   
  for(ii=0; ii<Nlim;ii++)
    for(jj=0; jj<Nlim;jj++)
      {
	index=ii*Nlim+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&GVLLre[index][nn],sizeof(GVLLre[index][jj]),1,fp);
	    fread(&GVLLim[index][nn],sizeof(GVLLim[index][jj]),1,fp);
	  }
      }

  //Self-correlation. (real) 
  for(ii=0; ii<Nlim;ii++)
    for(jj=0; jj<Nlim;jj++)
      {
	index=ii*Nlim+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&SCRR[index][nn],sizeof(SCRR[index][jj]),1,fp);
	    fread(&SCLL[index][nn],sizeof(SCLL[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",SCRR[index][nn],SCRR[index][nn]);
	  }
      }
  fclose(fp); 
}*/











void ReadGVnAC_gridding()
{
  int ii,jj,nn;
  long index;
  char file_GV[512],file_samp[512];
  FILE *fp;
  
   //snprintf(file_GV, 512*sizeof(char),"/scratch/sukhdeep/project_MABS/3VC_B/3_vcg/output_data/GV_data_1"); // open file to write bispectrum
   //snprintf(file_GV, 512*sizeof(char),"/scratch/sukhdeep/project_MABS/3VC_B/3_vcg/output_data/%s",GV);
   snprintf(file_GV, 512*sizeof(char),"%s/GV_data_%d",data_path,NRi);
   fp=fopen(file_GV,"rb"); 
   if (fp == NULL) {
    perror("Error opening file: GV_data");
    exit(1);
}
 
 fread(&NN,sizeof(NN),1,fp);
 fread(&nc, sizeof(nc), 1, fp);
 fread(&Du,sizeof(Du),1,fp);
 fread(&flag_key,sizeof(flag_key),1,fp);

  for(ii=0; ii<N1;ii++)
    for(jj=0; jj<N2/2+1;jj++)
      {
	index = (long)ii*(N2/2+1)+jj; 
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&GVRRre[index][nn],sizeof(GVRRre[index][nn]),1,fp);
	    fread(&GVRRim[index][nn],sizeof(GVRRim[index][nn]),1,fp);
	    //printf("\nRe= %e, Im= %e",GVRRre[index][nn],GVRRim[index][nn]);
	  }
      }
      
       //printf("GV[0][0]=%f\n",GVRRre[0][0]);
       fclose(fp); 
   printf("\n Successfully Read the Input GV Data: %ld*%ld*%ld\t Du = %lf, flag = %d.",NN,NN/2+1,nc,Du,flag_key);   
       
   snprintf(file_samp, 512*sizeof(char),"%s/sampling",data_path);
   fp=fopen(file_samp,"rb"); 
    //fp=fopen("../3_vcg/output_data/sampling","rb");    
   for(ii=0; ii<N1;ii++)
    for(jj=0; jj<N2/2+1;jj++)
      {
	index = (long)ii*(N2/2+1)+jj; 
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&NTiG[index][nn],sizeof(NTiG[index][nn]),1,fp);
	  }
      }
      
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




void normalize_data()
{
int ii,jj,i;
long index;

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

}


/*
void normalize_coll_data()
{
int ii,jj,i;
long index;

    for(ii=0; ii<N1;ii++)
     for(jj=0; jj<N2/2+1;jj++)
      {
	index = (long)ii*(N2/2+1)+jj; 
	for(i=0; i<nc_coll; i++)
	  {
	    if (NTiG_coll[index][i]>0)
	    {
	    GVRRre_coll[index][i] /= NTiG_coll[index][i];
	    GVRRim_coll[index][i] /= NTiG_coll[index][i];
	    }
	  }
      }
  
    for(ii=0; ii<N1;ii++)
     for(jj=0; jj<N2/2+1;jj++)
      {
	index = (long)ii*(N2/2+1)+jj; 
	for(i=0; i<nc_coll; i++)
	  {
	  if (NTiG_coll[index][i]>0)
	   {
	    NTiG_coll[index][i]=1;
	   }
	  }
      }

}*/







void annular_rings()
{

    int intervals[][3] =  {{0, 50, 5}, {50, 100, 10}, {100, 200, 20}, {200, 270, 35}};
    int num_intervals  =  4;
    
    //int intervals[][5] = {{0, 50, 5}, {50, 100, 10}, {100, 200, 20}, {200, 410, 35}, {410, 660, 50}};
    //int num_intervals = 5;

    // Index for ULow array
    long index = 0;
    int i,j;

    // Loop through intervals and fill ULow and UHigh
    for (i = 0; i < num_intervals; ++i) {
        int start = intervals[i][0];
        int end = intervals[i][1];
        int step = intervals[i][2];

        for (j = start; j < end; j += step) {
            if (index < Nr) {
                ULow[index] = j;
                if (index > 0) {
                    UHigh[index - 1] = ULow[index];
                }
                ++index;
            }
        }
    }
 
  ULow[0] = 1.;
  UHigh[Nr-1]=270.;
  //UHigh[Nr-1]=1.*Umax;
  
  
  for(i=0;i<Nr;++i)
  { //k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;

/*if (i<10)
   {k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;}
else
   {k_perp_avg[i]=(2./3.)*(pow(k_perp_High[i],3.)-pow(k_perp_Low[i],3.))/(pow(k_perp_High[i],2.)-pow(k_perp_Low[i],2.)); }// bin average value */
   
   Uavg[i]=(2./3.)*(pow(UHigh[i],3.)-pow(ULow[i],3.))/(pow(UHigh[i],2.)-pow(ULow[i],2.));
   
   //printf("\t%f\t%f\t%f\n",ULow[i],UHigh[i],Uavg[i]) ;
  }
  
   /*  if(Umax>Du*N1/3)
  {
   printf("Umax is out of the bound dU*(N1/3).\n");
   return 1; 
   } */

  /*----------- Dividing U1 plane for bispectra -----------*/

  
 /*------------- Dividing k1 space into linear bins -------------*/
 /* dU=(Umax-Umin)/(1.*Nbin); // goes max upto dU*(N/3) to avoid modes outside boxes
   fprintf(stdout,"Binsize= %e\n", dU);

  for(i=0;i<Nbin;++i)
  {
   ULow[i]=i*dU+Umin;
   UHigh[i]=(i+1)*dU+Umin;
   Uavg[i]=(2./3.)*(pow(UHigh[i],3.)-pow(ULow[i],3.))/(pow(UHigh[i],2.)-pow(ULow[i],2.)); // bin average value
   
   printf("Sl=%d \t ULow=%e\t UHigh=%e \t Uave=%e\n",i+1,ULow[i],UHigh[i], Uavg[i]);

  }
  */

}




// Function to count and store valid combinations
int precompute_vcp(int nc, int nc10, int **vcp) {
    int count = 0;

    // First pass: Count the number of valid combinations
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int p2 = p1 - nc10; p2 <= p1 + nc10; ++p2) {
            if (p2 < 0 || p2 >= nc) continue; // Skip if p2 is out of bounds
            int dnu2 = abs(p2 - p1);
            for (int p3 = p1 - nc10; p3 <= p1 + nc10; ++p3) {
                if (p3 < 0 || p3 >= nc) continue; // Skip if p3 is out of bounds
                int dnu3 = abs(p3 - p1);

                // Condition (1): Skip if both p2-p1 and p3-p1 are negative
                //if ((p2 - p1) < 0 && (p3 - p1) < 0) continue;

                // If all conditions are met, increment the count
                count++;
            }
        }
    }

    // Allocate memory for the combinations array
    *vcp = (int *)malloc(6 * count * sizeof(int));
    if (*vcp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Second pass: Store the valid combinations
    int index = 0;
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int p2 = p1 - nc10; p2 <= p1 + nc10; ++p2) {
            if (p2 < 0 || p2 >= nc) continue; // Skip if p2 is out of bounds
            int dnu2 = p2 - p1 + nc10;
            for (int p3 = p1 - nc10; p3 <= p1 + nc10; ++p3) {
                if (p3 < 0 || p3 >= nc) continue; // Skip if p3 is out of bounds
                int dnu3 = p3 - p1 + nc10;

                // Condition (1): Skip if both p2-p1 and p3-p1 are negative
                //if ((p2 - p1) < 0 && (p3 - p1) < 0) continue;

                // Store the combination
                (*vcp)[index++] = p1;
                (*vcp)[index++] = p2;
                (*vcp)[index++] = p3;
                (*vcp)[index++] = dnu2;
                (*vcp)[index++] = dnu3;
                (*vcp)[index++] = dnu2 * (2*nc10+1) + dnu3;
            }
        }
    }

    return count;
}



int precompute_vcp_adaptive(int nc10, int **vcp) {
    int count = 0;

    // First pass: Count the number of valid combinations
    for (int p1 = 0; p1 <= nc10; ++p1) {
        for (int p2 = 0; p2 <= nc10; ++p2) {
            int dnu2 = abs(p2 - p1);
            for (int p3 = 0; p3 <= nc10; ++p3) {
                int dnu3 = abs(p3 - p1);

                // Condition (1): Skip if both p2-p1 and p3-p1 are negative
                //if ((p2 - p1) < 0 && (p3 - p1) < 0) continue;

                // If all conditions are met, increment the count
                count++;
            }
        }
    }

    // Allocate memory for the combinations array
    *vcp = (int *)malloc(6 * count * sizeof(int));
    if (*vcp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Second pass: Store the valid combinations
    int index = 0;
    for (int p1 = 0; p1 <= nc10; ++p1) {
        for (int p2 = 0; p2 <= nc10; ++p2) {
            int dnu2 = p2 - p1 + nc10;
            for (int p3 = 0; p3 <= nc10; ++p3) {
                int dnu3 = p3 - p1 + nc10;

                // Store the combination
                (*vcp)[index++] = p1;
                (*vcp)[index++] = p2;
                (*vcp)[index++] = p3;
                (*vcp)[index++] = dnu2;
                (*vcp)[index++] = dnu3;
                (*vcp)[index++] = dnu2 * (2*nc10+1) + dnu3;
            }
        }
    }

    return count;
}



int precompute_unflg_vcp(int nc10, int **vcp) {
    int count = 0;

    // First pass: Count the number of valid combinations
    for (int p1_i = 0; p1_i < uf_nc; ++p1_i) {
    	int p1 = uf_chans[p1_i];
        for (int p2_i = 0; p2_i < uf_nc; ++p2_i) {
            int p2 = uf_chans[p2_i];
            int dnu2 = abs(p2 - p1);
            if ( dnu2 > nc10 ) continue;
            for (int p3_i = 0; p3_i < uf_nc; ++p3_i) {
                int p3 = uf_chans[p3_i];
                int dnu3 = abs(p3 - p1);
                if ( dnu3 > nc10 ) continue;

                count++;
            }
        }
    }

    // Allocate memory for the combinations array
    *vcp = (int *)malloc(6 * count * sizeof(int));
    if (*vcp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Second pass: Store the valid combinations
    int index = 0;
    for (int p1_i = 0; p1_i < uf_nc; ++p1_i) {
        int p1 = uf_chans[p1_i];
        for (int p2_i = 0; p2_i < uf_nc; ++p2_i) {
            int p2 = uf_chans[p2_i];
            if ( abs(p2-p1) > nc10 ) continue;
            int dnu2 = p2 - p1 + nc10;
            for (int p3_i = 0; p3_i < uf_nc; ++p3_i) {
                int p3 = uf_chans[p3_i];
                if ( abs(p3-p1) > nc10 ) continue;
                int dnu3 = p3 - p1 + nc10;

                // Store the combination
                (*vcp)[index++] = p1;
                (*vcp)[index++] = p2;
                (*vcp)[index++] = p3;
                (*vcp)[index++] = dnu2;
                (*vcp)[index++] = dnu3;
                (*vcp)[index++] = dnu2 * (2*nc10+1) + dnu3;
            }
        }
    }

    return count;
}




int precompute_vcp_dnu(int nc, int nc10, int **vcp) {
    int count = 0;

    // First pass: Count the number of valid combinations
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int p2 = 0; p2 < nc; ++p2) {
            int dnu2 = abs(p2 - p1);
            if ( dnu2 > nc10 ) continue;
            for (int p3 = 0; p3 < nc; ++p3) {
                int dnu3 = abs(p3 - p1);
                if ( dnu3 > nc10 ) continue;

                count++;
            }
        }
    }

    // Allocate memory for the combinations array
    *vcp = (int *)malloc(6 * count * sizeof(int));
    if (*vcp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Second pass: Store the valid combinations
    int index = 0;
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int p2 = 0; p2 < nc; ++p2) {
            if ( abs(p2-p1) > nc10 ) continue;
            int dnu2 = p2 - p1 + nc10;
            for (int p3 = 0; p3 < nc; ++p3) {
                if ( abs(p3-p1) > nc10 ) continue;
                int dnu3 = p3 - p1 + nc10;

                // Store the combination
                (*vcp)[index++] = p1;
                (*vcp)[index++] = p2;
                (*vcp)[index++] = p3;
                (*vcp)[index++] = dnu2;
                (*vcp)[index++] = dnu3;
                (*vcp)[index++] = dnu2 * (2*nc10+1) + dnu3;
            }
        }
    }

    return count;
}

 

// Function to precompute valid combinations
int precompute_vcp_dnu_old(int nc, int nc10, int **vcp) {
    int count = 0;

    // First pass: Count the number of valid combinations
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int dnu2 = 0; dnu2 <= nc10; ++dnu2) {
            int p1_plus_dnu2 = p1 + dnu2;
            int p1_minus_dnu2 = p1 - dnu2;
            if (p1_minus_dnu2 < 0 || p1_plus_dnu2 >= nc) continue;
            
            for (int dnu3 = 0; dnu3 <= nc10; ++dnu3) {
                int p1_plus_dnu3 = p1 + dnu3;
                int p1_minus_dnu3 = p1 - dnu3;
                if (p1_minus_dnu3 < 0 || p1_plus_dnu3 >= nc) continue;
                
                count++;
            }
        }
    }

    // Allocate memory for the combinations array
    *vcp = (int *)malloc(7 * count * sizeof(int));
    if (*vcp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Second pass: Store the valid combinations
    int index = 0;
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int dnu2 = 0; dnu2 <= nc10; ++dnu2) {
            int p1_plus_dnu2 = p1 + dnu2;
            int p1_minus_dnu2 = p1 - dnu2;
            if (p1_minus_dnu2 < 0 || p1_plus_dnu2 >= nc) continue;
            
            for (int dnu3 = 0; dnu3 <= nc10; ++dnu3) {
                int p1_plus_dnu3 = p1 + dnu3;
                int p1_minus_dnu3 = p1 - dnu3;
                if (p1_minus_dnu3 < 0 || p1_plus_dnu3 >= nc) continue;
                
                (*vcp)[index++] = p1;
                (*vcp)[index++] = dnu2;
                (*vcp)[index++] = dnu3;
                (*vcp)[index++] = p1_plus_dnu2;
                (*vcp)[index++] = p1_minus_dnu2;
                (*vcp)[index++] = p1_plus_dnu3;
                (*vcp)[index++] = p1_minus_dnu3;
            }
        }
    }

    return count;
}



// Function to precompute valid combinations
int precompute_vcp_dnu_expt0(int nc, int nc10, int **vcp) {
    int count = 0;

    // First pass: Count the number of valid combinations
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int dnu2 = 1; dnu2 <= nc10; ++dnu2) {
            int p1_plus_dnu2 = p1 + dnu2;
            int p1_minus_dnu2 = p1 - dnu2;
            if (p1_minus_dnu2 < 0 || p1_plus_dnu2 >= nc) continue;
            
            for (int dnu3 = 1; dnu3 <= nc10; ++dnu3) {
                int p1_plus_dnu3 = p1 + dnu3;
                int p1_minus_dnu3 = p1 - dnu3;
                if (p1_minus_dnu3 < 0 || p1_plus_dnu3 >= nc) continue;
                
                count++;
            }
        }
    }

    // Allocate memory for the combinations array
    *vcp = (int *)malloc(7 * count * sizeof(int));
    if (*vcp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Second pass: Store the valid combinations
    int index = 0;
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int dnu2 = 1; dnu2 <= nc10; ++dnu2) {
            int p1_plus_dnu2 = p1 + dnu2;
            int p1_minus_dnu2 = p1 - dnu2;
            if (p1_minus_dnu2 < 0 || p1_plus_dnu2 >= nc) continue;
            
            for (int dnu3 = 1; dnu3 <= nc10; ++dnu3) {
                int p1_plus_dnu3 = p1 + dnu3;
                int p1_minus_dnu3 = p1 - dnu3;
                if (p1_minus_dnu3 < 0 || p1_plus_dnu3 >= nc) continue;
                
                (*vcp)[index++] = p1;
                (*vcp)[index++] = dnu2;
                (*vcp)[index++] = dnu3;
                (*vcp)[index++] = p1_plus_dnu2;
                (*vcp)[index++] = p1_minus_dnu2;
                (*vcp)[index++] = p1_plus_dnu3;
                (*vcp)[index++] = p1_minus_dnu3;
            }
        }
    }

    return count;
}




int precompute_vcp_old(int nc, int nc10, int **vcp) {
    int count = 0;

    // First pass: Count the number of valid combinations
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int p2 = p1 - nc10; p2 <= p1 + nc10; ++p2) {
            if (p2 < 0 || p2 >= nc) continue; // Skip if p2 is out of bounds
            int dnu2 = abs(p2 - p1);
            for (int p3 = p1 - nc10; p3 <= p1 + nc10; ++p3) {
                if (p3 < 0 || p3 >= nc) continue; // Skip if p3 is out of bounds
                int dnu3 = abs(p3 - p1);

                // Condition (1): Skip if both p2-p1 and p3-p1 are negative
                if ((p2 - p1) < 0 && (p3 - p1) < 0) continue;

                // If all conditions are met, increment the count
                count++;
            }
        }
    }

    // Allocate memory for the combinations array
    *vcp = (int *)malloc(6 * count * sizeof(int));
    if (*vcp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Second pass: Store the valid combinations
    int index = 0;
    for (int p1 = 0; p1 < nc; ++p1) {
        for (int p2 = p1 - nc10; p2 <= p1 + nc10; ++p2) {
            if (p2 < 0 || p2 >= nc) continue; // Skip if p2 is out of bounds
            int dnu2 = abs(p2 - p1);
            for (int p3 = p1 - nc10; p3 <= p1 + nc10; ++p3) {
                if (p3 < 0 || p3 >= nc) continue; // Skip if p3 is out of bounds
                int dnu3 = p3 - p1 + nc10;

                // Condition (1): Skip if both p2-p1 and p3-p1 are negative
                if ((p2 - p1) < 0 && (p3 - p1) < 0) continue;

                // Store the combination
                (*vcp)[index++] = p1;
                (*vcp)[index++] = p2;
                (*vcp)[index++] = p3;
                (*vcp)[index++] = dnu2;
                (*vcp)[index++] = dnu3;
                (*vcp)[index++] = dnu2 * (2*nc10+1) + dnu3;
            }
        }
    }

    return count;
}


int get_uf_chans(const char *filename) {
    char file_ufchans[512];
    snprintf(file_ufchans, sizeof(file_ufchans),"%s/%s",data_path, filename);

    FILE *fp = fopen(file_ufchans, "r");
    if (!fp) {
        perror("Error opening file");
        return -1;
    }

    int i = 0;
    while (i < uf_nc && fscanf(fp, "%d", &uf_chans[i]) == 1) {
        i++;
    }

    if (i < uf_nc) {
        fprintf(stderr, "Error: expected %zu integers, but only read %zu\n", uf_nc, i);
        fclose(fp);
        return -1;
    }

    fclose(fp);
    return 0;
}


