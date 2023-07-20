/*
 Header file to supply the definitions of reimplemented lapack functionality

 These two functions are used in PSFEx solver, but have been reimplemented to
 avoid a dependancy on anything fortran. the source can be found in src/clapack.c
*/

void dposv_internal(char *, long *, long *, double *, long *, double *, long *, long *);

void dtrtri_internal(char *, char *, long *, double *, long *, long *);
