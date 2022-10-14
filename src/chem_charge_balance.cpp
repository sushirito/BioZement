#include "chem_global.h"



double f_pH(double x, int no_call, struct BasVec *Vchem, void (*mbal)(int, struct BasVec *))
{
  int i;
  real f_tmp;
  struct ChemTable *SM_basis, *SM_all;

  SM_basis = Vchem->ICS->SM_basis;
  SM_all   = Vchem->ICS->SM_all;
  /* H+ */
  Vchem->log_a[Vchem->pos_pH] = x;
  Vchem->log_m[Vchem->pos_pH] = x-Vchem->log_g[Vchem->pos_pH];


  /*--------------------- apply mass balance -------------------------*/
  if (Vchem ->size_mass>0)
    {
      //if (_DEBUG_FLAG_) printf("\mbal(no_call=%d, Vchem);\n",no_call);
      mbal(no_call, Vchem);
      //if (_DEBUG_FLAG_) printf("\nafter\n");
    }
  else /* calculate the concentration of the rock buffered species and complexes */
    {
      if(Vchem->size_rock>0) 
	calc_rock_spec_conc(Vchem);
      calc_complex(SM_all, Vchem);
    }


  /* sum up basis species */

  f_tmp = 0.;
  if(Vchem->DL)
    {
      for(i=0;i<Vchem->size;++i)
	{
	  f_tmp +=  pow10(Vchem->log_m[i])*SM_basis->charge[i];
	}
    }
  else
    {
      for(i=0;i<Vchem->size;++i)
	{
	  f_tmp +=  pow10(Vchem->log_m[i])*(SM_basis->scharge[i]+SM_basis->charge[i]);
	}
    }
	
	
  /* SM_all->charge[i]  = 0 for all surface complexes */
  /* SM_all->scharge[i] = 0 for all aqueous complexes */
  if(Vchem->DL)
    for(i=0;i<SM_all->size[0];++i) f_tmp +=  pow10(SM_all->log_m[i])*SM_all->charge[i];
  else
    for(i=0;i<SM_all->size[0];++i) f_tmp +=  pow10(SM_all->log_m[i])*(SM_all->charge[i]+SM_all->scharge[i]);
	
  /* add the charge of the cation exchanger */
  if(Vchem->pos_X >= 0) f_tmp -= Vchem->ctot[Vchem->pos_X];
  return (double) f_tmp;
}


/*
 ************************************************************************
 *	    		    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,tol)
 *	double ax; 			Root will be seeked for within
 *	double bx;  			a range [ax,bx]
 *	double (*f)(double x);		Name of the function whose zero
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

//double zeroin( x_init, indx, ax, bx,f,tol,Vchem)		/* An estimate to the root	*/
//double ax;				/* Left border | of the range	*/
//double bx;  				/* Right border| the root is seeked*/
//double x_init;       /* previous solution */
//int indx ;           /* use initial guess ? */
//double (*f)(double x, int no_call, struct BasVec *Vchem);			/* Function under investigation	*/
//double tol;				/* Acceptable tolerance		*/
//struct BasVec *Vchem;
//{
//  double a,b,c;				/* Abscissae, descr. see above	*/
//  double fa;				/* f(a)				*/
//  double fb;				/* f(b)				*/
//  double fc;				/* f(c)				*/
//  double crit;
//  int    icrit;
//  int iter;
//
//  crit  = 2.;
//  icrit = 0;
//  iter =0;
//  if(indx)
//  {
//  a=x_init-2.;
//  b=x_init+1.;
//  fa = (*f)(a,icrit,Vchem );
//  fb = (*f)(b,icrit,Vchem );
//  c=a; fc=fa;
//  } else
//  {
//	  fa=fb=1.;
//  }
//  if(fa*fb>0) /* if init guess is bad */
//  {
//  a = ax;  b = bx;  fa = (*f)(a,icrit, Vchem );  fb = (*f)(b,icrit,Vchem);
//  c = a;   fc = fa;
//  }
//  for(;;)		/* Main iteration loop	*/
//  {
//    double prev_step = b-a;		/* Distance from the last but one*/
//					/* to the last approximation	*/
//    double tol_act;			/* Actual tolerance		*/
//    double p;      			/* Interpolation step is calcu- */
//    double q;      			/* lated in the form p/q; divi- */
//  					/* sion operations is delayed   */
// 					/* until the last moment	*/
//    double new_step;      		/* Step at this iteration       */
//
//    if( fabs(fc) < fabs(fb) )
//    {                         		/* Swap data for b to be the 	*/
//	a = b;  b = c;  c = a;          /* best approximation		*/
//	fa=fb;  fb=fc;  fc=fa;
//    }
//    tol_act = 2*CHEM_EPSILON_*fabs(b) + tol/2;
//    new_step = (c-b)/2;
//
//    if( fabs(new_step) <= tol_act || fb == (double)0 )
//	{
//		if(Vchem->ICS->PRINT_DEBUG_CHEM) printf("ZEROIN: No iter: %d \n", iter);
//		return b;				/* Acceptable approx. is found	*/
//	}
//
//    			/* Decide if the interpolation can be tried	*/
//    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
//	&& fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
//    {					/* Interpolatiom may be tried	*/
//	register double t1,cb,t2;
//	cb = c-b;
//	if( a==c )			/* If we have only two distinct	*/
//	{				/* points linear interpolation 	*/
//	  t1 = fb/fa;			/* can only be applied		*/
//	  p = cb*t1;
//	  q = 1.0 - t1;
// 	}
//	else				/* Quadric inverse interpolation*/
//	{
//	  q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
//	  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
//	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
//	}
//	if( p>(double)0 )		/* p was calculated with the op-*/
//	  q = -q;			/* posite sign; make p positive	*/
//	else				/* and assign possible minus to	*/
//	  p = -p;			/* q				*/
//
//	if( p < (0.75*cb*q-fabs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
//	    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
//	  new_step = p/q;			/* it is accepted	*/
//					/* If p/q is too large then the	*/
//					/* bissection procedure can 	*/
//					/* reduce [b,c] range to more	*/
//					/* extent			*/
//    }
//
//    if( fabs(new_step) < tol_act ) {/* Adjust the step to be not less*/
//      if( new_step > (double)0 ) {	/* than tolerance		*/
//        new_step = tol_act;
//      } else {
//        new_step = -tol_act;
//      }
//    }
//    iter++;
//    a = b;  fa = fb;			/* Save the previous approx.	*/
//    b += new_step;
//	if(fabs(a-b)<crit) icrit = 1;
//	fb = (*f)(b, icrit, Vchem);	/* Do step to a new approxim.	*/
//    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
//    {                 			/* Adjust c for it to have a sign*/
//      c = a;  fc = fa;                  /* opposite to that of b	*/
//    }
//
//  }
//
//}

/*
 ************************************************************************
 *	    		    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,tol)
 *	double ax; 			Root will be seeked for within
 *	double bx;  			a range [ax,bx]
 *	double (*f)(double x);		Name of the function whose zero
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

//double zeroin2( ax, bx,f,tol,Vchem)		/* An estimate to the root	*/
//double ax;				/* Left border | of the range	*/
//double bx;  				/* Right border| the root is seeked*/
//double (*f)(double x, struct BasVec *Vchem);			/* Function under investigation	*/
//double tol;				/* Acceptable tolerance		*/
//struct BasVec *Vchem;
//{
//  double a,b,c;				/* Abscissae, descr. see above	*/
//  double fa;				/* f(a)				*/
//  double fb;				/* f(b)				*/
//  double fc;				/* f(c)				*/
//  int iter;
//  iter =0;
//  a = ax;  b = bx;  fa = (*f)(a,Vchem );  fb = (*f)(b,Vchem);
//  c = a;   fc = fa;
//  iter = 0;
//  if(fa*fb>0) printf("ZERION2 : fa and fb have same sign!\n");
//  while(fabs(a-b) > tol)
//  {
//	  c =(a+b)*.5;fc=(*f)(c,Vchem);
//	  if(fc*fa > 0)
//	  {
//		  a = c;
//		  fa=fc;
//	  }
//	  else
//	  {
//		  b=c;
//	  }
//
//	  iter ++;
//  }
//
//
//}
/*
void print_fpH(real x1, real x2, struct BasVec *Vchem)
{
	int i, pos_ph;
	int N=100;
	real dx, x_init, value;
	FILE *fp;

	fp=my_fopen("pH_debug.out", "w");
	fclose(fp);
	dx = (x2-x1)/((real) N);

	x_init = x1;
	for(i=0;i<N;++i)
	{
		fp=my_fopen("pH_debug.out", "a");
		value = f_pH(x_init,0,Vchem);
		fprintf(fp,"%lf\t%lf\n", -x_init, value);
		if(Vchem->ICS->PRINT_DEBUG_CHEM)
			printf("%lf\t%lf\n", -x_init, value);
		x_init +=dx;

		fclose(fp);
	}

return;
}
*/
