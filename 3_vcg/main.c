# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include "tgefuncs.h"

char INFITS[128], input[128], GV[128],file_ntig[128];
long NN,n_collapse; 
double dU,theta_0,FWHM,Umax;
int flag_key;

int main(int argc, char *argv[])
{ 
  int key_ntig;
  int coll_key;
  
  if(argc!=8) 
    {
      printf("Usage: %s <input FITS file> <input parameters file> <output GV file> <flagging?> <Collapse?> <want sampling_fle?> <output sampling_file>\n", argv[0]);
      return 1;
    }
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",input); 
  sscanf(argv[3],"%s",GV); 
  sscanf(argv[4],"%d",&flag_key); 
  sscanf(argv[5],"%d",&coll_key);
  sscanf(argv[6],"%d",&key_ntig); 
  sscanf(argv[7],"%s",file_ntig);  
  
  n_collapse = 4;
  
  printf("\n##################\n");
  printf("\nkey_ntig=%d\n",key_ntig);
  printf("\nFlag Key=%d\n",flag_key);
  printf("\nCollapse Key=%d\n",coll_key);
  printf("\n##################\n");
  
  read_inputs(input); 
 
  
  initialize_gridding();
  printf("\nArrays Initialized\n");

 /* -- Flagging Info ---   */
 if (flag_key==1)      {group_loop_gridding();}
 else if (flag_key==0) {group_loop_gridding_without_flag();}
 else {printf("\nWrong Flag Key\n");}
  
 /* -- Collapsing Info ---   */
  if (coll_key==1)
    {collapse_channels(n_collapse);
    WriteGVnAC_coll();}
  else if (coll_key==0)
    {printf("\nNo Channel Collapsing\n");
     WriteGVnAC_gridding();}
  else {printf("\nWrong Collapse Key\n");}
  
  //WriteGVnAC_coll();
  //WriteGVnAC_gridding();
  //WriteGVnAC_gridding_realpart();
  
  if (key_ntig==1)
   {
   if (coll_key==0)
   {Write_NTiG();}
   else if (coll_key==1)
   {Write_NTiG_coll();}
   else {printf("\nWrong Collapse Key\n");}
   }
  

  
  printf("\nEnd.\n");
  return 0;   
} 
