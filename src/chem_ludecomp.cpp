#include "chem_global.h"

void ludcmp(double **a, int n, int *indx, double *d)
     /*Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise*/
     /*permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;*/
     /*indx[1..n] is an output vector that records the row permutation eected by the partial*/
     /*pivoting; d is output as ±1 depending on whether the number of row interchanges was even*/
     /*or odd, respectively. This routine is used in combination with lubksb to solve linear equations*/
     /*or invert a matrix.*/
{
int i,imax,j,k;
double big,dum,sum,temp;
 double *vv; /*vv stores the implicit scaling of each row.*/
 vv= (double *) calloc(n+1, sizeof(double)); 
 *d=1.0; /*No row interchanges yet.*/
 for (i=1;i<=n;i++) { /*Loop over rows to get the implicit scaling information.*/
   big=0.0;
   for (j=1;j<=n;j++)
     if ((temp=fabs(a[i][j])) > big) big=temp;
   if (big == 0.0) {
     printf("Singular matrix in routine ludcmp\n");
     for (j=1;j<=n;j++) {
       for (k=1;k<=n;k++) {
	 printf("%g, ", a[j][k]);
       }
       printf("\n");
     }
     printf("\n");
     exit(1);
   }
   /*No nonzero largest element.*/
   vv[i]=1.0/big; /*Save the scaling.*/
 }
 for (j=1;j<=n;j++) { 
   for (i=1;i<j;i++) { 
     sum=a[i][j];
     for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
     a[i][j]=sum;
   }
   big=0.0;/* Initialize for the search for largest pivot element.*/
   for (i=j;i<=n;i++) { /*This is i = j of equation (2.3.12) and i = j+1. . .N*/
     /*of equation (2.3.13). */
     sum=a[i][j];
     for (k=1;k<j;k++)
       sum -= a[i][k]*a[k][j];
     a[i][j]=sum;
     if ( (dum=vv[i]*fabs(sum)) >= big) {
       /*Is the gure of merit for the pivot better than the best so far?*/
       big=dum;
       imax=i;
     }
   }
   if (j != imax) { /*Do we need to interchange rows?*/
     for (k=1;k<=n;k++) { /*Yes, do so...*/
       dum=a[imax][k];
       a[imax][k]=a[j][k];
       a[j][k]=dum;
     }
     *d = -(*d);/* ...and change the parity of d.*/
     vv[imax]=vv[j]; /*Also interchange the scale factor.*/
   }
   indx[j]=imax;
   if (a[j][j] == 0.0) a[j][j]=CHEM_EPSILON_;
   /*If the pivot element is zero the matrix is singular (at least to the precision of the */
   /*algorithm). For some applications on singular matrices, it is desirable to substitute */
   /*TINY for zero.*/
   if (j != n) { /*Now, nally, divide by the pivot element.*/
     dum=1.0/(a[j][j]);
     for (i=j+1;i<=n;i++) a[i][j] *= dum;
   }
 } /*Go back for the next column in the reduction.*/
 free(vv);
}

void lubksb(double **a, int n, int *indx, double b[])
     /*Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix*/
     /*A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input*/
     /*as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector*/
     /*B, and returns with the solution vector X. a, n, and indx are not modied by this routine*/
     /*and can be left in place for successive calls with dierent right-hand sides b. This routine takes*/
     /*into account the possibility that b will begin with many zero elements, so it is ecient for use*/
     /*in matrix inversion.*/
{
  int i,ii=0,ip,j;
  double sum;
  for (i=1;i<=n;i++) { /*When ii is set to a positive value, it will become the*/
    /*index of the rst nonvanishing element of b. Wenow*/
    /*do the forward substitution, equation (2.3.6). The*/
    /*only new wrinkle is to unscramble the permutation*/
    /*as we go.*/
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i; /*A nonzero element was encountered, so from now on we*/
    /*will have to do the sums in the loop above. */
    b[i]=sum;
}
  for (i=n;i>=1;i--) { /*Now we do the backsubstitution, equation (2.3.7).*/
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];/* Store a component of the solution vector X.*/
  }/* All done!*/
}

