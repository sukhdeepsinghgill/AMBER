void printerror( int status); // prints error messages, if any.
int read_fits_header(char* in_file); // reads the header
int read_ranpar(long grp); // reads random parameters u and v
int read_data(long group); // reads visibility data
int close_fits(); // closes an opened fits file
