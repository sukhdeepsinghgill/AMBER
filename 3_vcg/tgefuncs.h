void read_inputs(char *in_file); 
double winf(double udif); 
int initialize(); 
void group_loop(); // for a random group, reads u,v and data from a fits file.
void free_func(); // frees all functions in the end.
void WriteGVnAC();

int initialize_gridding(); 
void group_loop_gridding();
void group_loop_gridding_without_flag();
void WriteGVnAC_gridding();
void Write_NTiG();
void WriteGVnAC_gridding_realpart();

void collapse_channels(long n_collapse);
void Write_NTiG_coll();
void WriteGVnAC_coll();
