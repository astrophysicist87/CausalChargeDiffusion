/* complex variable auxilliary routines
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

sign  returns sign of type double argument, i.e. -1,+1, or 0
argmt	argument of complex number (i.e., angle in polar representation)
		returns angle between -pi to pi
Argmt as above, but returns 0 to 2pi
clog	complex logarithm
csqrt	complex square root
polarxy convert from polar to rectangular
cexp	complex exponential
ctrig	complex cosine and sine
		complex trig functions:
csin	sine
ccos	cosine
ctan  tangent
ccot  cotangent
csinh	sinh
ccosh	cosh
casin arcsine
cacos arccosine
catan arctangent
printc	print a complex number

*/

#include "stdio.h"
#include "stdlib.h"
#include "cmlib.h"
#include "complex.h"
#include "protom.h"
/*	#define DEFATAN2*/

/*static unsigned int btm,top;*/

namespace SFL
{
	double sign(double x)
	{
		if(x==0.)
			return 0.;
		return ( x >= 0. ) ? 1. : -1.;
	}

	double argmt(double y,double x)
	{
		/* returns answer between -pi to pi, in conformity with atan2
		 and the complex log branch cut along negative real axis.
		This way, the logs of complex conjugates are conjugates */
		/* caveat- Aztec C returns 0,not + or -halfpi if x=0., y nonzero*/
		double ans,ratio,/*twopi=6.283185307,*/halfpi=.5*pi;
		double undef=0.0;/* change if desired*/
		/* use atan2 if it exists:*/
		#ifdef DEFATAN2
		printf(" using system atan2\n");
		return atan2(y,x);
		/*otherwise, hand code it*/
		#else
		if(y==0.)
		{
			if(x>0.)return 0.;
			else if(x<0.)return pi;
			else return undef;
		}
			
		if (x==0.)
		{
			if(y>0.) return (halfpi);
			if(y<0.) return(-halfpi);
			/*return (halfpi*sign(y));*/
			return(undef);
		};
		ratio=(y/x);ans=atan(ratio);
		/*atan returns answer between -halfpi and halfpi*/
		if(x>0.) return ans;/* Quadrants I and IV*/
		/* x<0.*/
		if( y<0.) return ans-pi;/*Quadrant III y/x and ans >0*/
		return  ans+pi;/*Quadrant II x<0, y>=0 ans<0*/
		#endif
	}

	double Argmt(double y,double x)
	{
		/* returns answer between 0 and twopi, as needed in
		complex principal argument*/
		/* caveat- Aztec C returns 0,not + or -halfpi if x=0., y nonzero*/
		double ans,ratio,twopi=6.283185307,halfpi=1.570796327;
		double undef=0.0;/* change if desired*/
		if(y==0.)
		{
			if(x>0.)return 0.;
			else if(x<0.)return pi;
			else return undef;
		}
		if (x==0.)
		{
			if(y>0.) return(halfpi);
			if(y<0.) return(pi+halfpi);
			return(undef);
		};
		/* if -pi to pi: if(x==0.)return (halfpi*sign(y));*/
		ratio=(y/x);ans=atan(ratio);
		/*atan returns answer between -halfpi and halfpi*/
		/* now move to correct quadrant  between -pi and pi*/
		if (ratio>0.){/* ratio, ans>0. */
				  if (x>0.) return(ans);/* quadrant I*/
				/* else x<0.,y<0. quadrant III*/
				return (pi+ans);
				      };
		/* else ratio,ans<0.*/
		if(x>0.)return(twopi+ans);/*quadrant IV*/
		/*else x<0.,y>0., quadrant II*/
		return(pi+ans);

		/* if answer bwtn -pi and pi desired: if ans<=pi accept ans unchanged
			else   ans-twopi
		   this will affect quadrant III,IV only. change:
		   	III: pi+ans to ans-pi
		   	 IV: twopi+ans to ans
		*/
	}

	/*------------- complex logarithm function ---------------*/

	int clog(struct complex *x,struct complex *ans)
	{
		double r,angle;
		r= sqrt( CNORM(*x) );
		angle=argmt(x->y,x->x);
		ans->x=log(r);
		ans->y=angle;
		return 0;
	}

	/* ----------------- complex square root --------------------*/

	int csqrt(struct complex *z,struct complex *ans)
	{
		double x,y,r,angle;
		r=  sqrt(sqrt( CNORM(*z) ) );
		angle=.5*argmt(z->y,z->x);
		polarxy(r,angle,&x,&y);
		ans->x=x;
		ans->y=y;
		return 0;
	}

	/* ------- convert from polar to rectangular coordinates---- */

	int polarxy(double r,double angle,double * x,double * y)
	{
		*x=r*cos(angle);
		*y=r*sin(angle);
		return 0;
	}

	/* complex exponential function */

	int cexp(struct complex *x,struct complex *ans)
	{
		double y;
		y = exp ( x->x);
		ans->x= y*cos (x->y);
		ans->y= y*sin (x->y);
		return 0;
	}

	/* cosine, sine of complex arguments
	 .5(cexp(i*z)+(cexp(-i*z)),etc.*/

	int ctrig(struct complex *z,struct complex *ccos,struct complex *csin)
	{
		double si,co,real,imag,e,ei,sinh,cosh;
		real=z->x;imag=z->y;
		e = exp (imag);ei=1./e;
		sinh= .5*(e-ei);
		cosh=.5*(e+ei);
		co=cos(real);
		si=sin(real);
		ccos->x=co*cosh;
		ccos->y=-si*sinh;
		csin->x=si*cosh;
		csin->y=co*sinh;
		return 0;
	}


	int csin(struct complex *z,struct complex *ans)
	{
		double x,y;
		x=z->x;y=z->y;
		ans->x= sin(x)*cosh(y);ans->y= cos(x)*sinh(y);
		return 0;
	}

	int ccos(struct complex *z,struct complex *ans)
	{
		double x,y;
		x=z->x;y=z->y;
		ans->x= cos(x)*cosh(y);ans->y= -sin(x)*sinh(y);
		return 0;
	}


	int ctan(struct complex *x, struct complex *ans)
	{
		struct complex c,s;
		ctrig(x,&c,&s);
		if( cabs( c) != 0.)
			{CDIV( (*ans),s,c);return 0;}
		CMPLX( *ans, errorcode,errorcode);return 1;
	}

	int ccot(struct complex *x, struct complex *ans)
	{
		struct complex ccos,csin;
		ctrig(x,&ccos,&csin);
		CDIV((*ans),ccos,csin);
		return 0;
	}


	int csinh(struct complex *x, struct complex *ans)
	{
		struct complex i,y,z; CMPLX(i,0.,1.);
		CMULT(y,i,*x);
		ctrig(&y,&z,ans);
		CMULT(z,*ans,i);CTREAL(*ans,z,-1.);
		return 0;
	}

	int ccosh(struct complex *x, struct complex *ans)
	{
		struct complex i,y,z; CMPLX(i,0.,1.);
		CMULT(y,i,*x);ctrig(&y,ans,&z);
		return 0;
	}

	/* asin z= -i log(iz+sqrt(1-z^2))
	   acos z= -i log( z+sqrt(z^2-1)) can add any multiple of twopi*/

	int casin(struct complex *x, struct complex *ans)
	{
		struct complex ci,z2,arg,sum,z;
		CLET(z,*x);
		CMULT(z2,z,z);
		CMPLX(ci,0.,1.);CMPLX(arg,1.,0.);CSUB(arg,arg,z2);
		clog(&arg,&sum);CTREAL(sum,sum,.5);cexp(&sum,&arg);
		CMULT(sum,z,ci);CADD(sum,sum,arg);clog(&sum,&arg);
		CMULT(sum,arg,ci);CTREAL((*ans),sum,-1.);
		return 0;
	}

	int cacos(struct complex *x, struct complex *ans)
	{
		struct complex ci,z2,arg,sum,z;
		CLET(z,*x);
		CMULT(z2,z,z);
		CMPLX(ci,0.,1.);CMPLX(arg,1.,0.);CSUB(arg,z2,arg);
		clog(&arg,&sum);CTREAL(sum,sum,.5);cexp(&sum,&arg);
		CADD(sum,z,arg);clog(&sum,&arg);
		CMULT(sum,arg,ci);CTREAL((*ans),sum,-1.);
		return 0;
	}

	int catan(struct complex *x, struct complex *ans)
	{
		struct complex i,n,d,r;CMPLX(i,0.,1.);
		CADD(n,i,*x);CSUB(d,i,*x);CDIV(r,n,d); clog(&r,&d);
		CMULT(r,i,d); CTREAL(*ans,d, .5);return 0;
	}

	int printc(struct complex *z)
	{
		char *s;
		if(z->y <0.) s="";
		else s="+";
		printf(" %le %s %le i ",z->x,s,z->y);
		return 0;
	}}
