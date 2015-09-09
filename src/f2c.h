/*
 These compiler directives are used to control various types that may be
 encountered in PSFEx code, they are reproduced here, a subset of a complete
 f2c.h, to ensure the stub of lapack which we re-implemented will function
 correctly with the existing PSFEx code
*/
#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
#ifdef INTEGER_STAR_8	/* Adjust for integer*8. */
typedef long long longint;		/* system-dependent */
typedef unsigned long long ulongint;	/* system-dependent */
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif
#endif
