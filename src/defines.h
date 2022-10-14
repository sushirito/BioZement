#ifndef DEFINES_H
#define DEFINES_H

#define GOOD 		1
#define BAD			0

#define FLUID	    0
#define WALL 	    1
#define INLET       2
#define OUTLET      3
#define SOLID	    4
#define SOLID_INERT 5
#define PERIODIC    6
#define GHOST       12
#define DISCARD     13
#define NUM_MODE    14

#define TRUE 		1
#define FALSE 		0
#define u_0 		0. /* low speed 0.1 */
#define faktor		1.0
#define _EPS_ 		1E-14
#define real double
#define MAXSTR 1000



#define GKS			0.04

#define ISQRT2		0.707106781186550
#define SQRT2		1.414213562373095
#define ISQRT3		0.577350269189625
#define PI		3.14159265358979

#define M2_TO_DARCY     1.013250E12
//#ifndef abs
//#define abs(a)     ( ((a) >= 0) ? (a) : -(a) )
//#endif

//#ifndef max
//#define max(a,b)    (((a) > (b)) ? (a) : (b))
//#endif

//#ifndef min
//#define min(a,b)    (((a) < (b)) ? (a) : (b))
//#endif

#ifndef dirac_delta
#define dirac_delta(a,b)  (((a) == (b)) ? 1 : 0)
#endif

#ifndef pow10
#define pow10(a)  exp(a*2.3025850929940459)
//#define pow10(a)  exp(a*2.30258509299404590109361379)
//#define pow10(a)  pow(10,a)
#endif

// control the datatype used in the output VTK files
// uncomment two next lines for DOUBLE output (8 bytes)
//#define VTK_DTYPE  "double"
//#define WRITE_TO_VTK(a, fp)  fwrite_reverse_endian_double(a, fp)
// uncomment two next lines for FLOAT output (4 bytes)
#define OUTPUT_DTYPE  float
#define VTK_DTYPE  "float"
#define WRITE_TO_VTK(a, fp)  fwrite_reverse_endian_float(a, fp)


/************************ CHEM *************************
#define SURFACE 1
#define BULK 0
#define COMPLEXNAME 20
#define TINY 1.0e-20 
#define EPSILON 1.0e-20 
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;} */


#endif // DEFINES_H

