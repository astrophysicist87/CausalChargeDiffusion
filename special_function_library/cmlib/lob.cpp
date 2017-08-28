/*  Lobachevsky function

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

static double
	lcoef[25]= {-1./6.,  -1./60.,  -1./315.,  -17./22680.,  -31./155925.,
   -691./12162150.,  -10922./638512875.,  -929569./173675502000.,  
	-3202291./1856156927625.,  -221930581./389792954801250.,
   -9444233042./49308808782358125.,  -56963745931./870155449100437500.,  
   -29435334228302./1298054391195577640625.,  
   -2093660879252671./263505041412702261046875., 
   -344502690252804724./122529844256906551386796875., 
   -129848163681107301953./129391515535293318264457500000., 
   -868320396104950823611./2405873491984360136479756640625., 
   -209390615747646519456961./1602311745661583850895517922656250., 
   -28259319101491102261334882./593656501767616816756789390344140625., 
	-16103843159579478297227731./923715998955305103872044212679687500.,
   -705105746914492614138024917342./
	 109894723324712387033933067993555591796875.,
   -129024520859926228378837238913451./
    54397888045732631581796868656810017939453125., 
   -51768752002756374733644471311053204./
    58804116977436974739922415018011629392548828125., 
   -22680552792491997823522126468923904459./
    69153641565465882294148760061181676165637421875000., 
   -43588258945274750397872591971875276602./
    355527794338584677117095439830671923835434326171875.};


/* for x near pi/2, accuracy of sum is only about .001*/
/* should do expansion abt x=pi/2, where L(pi/2)= ln2*/
double Lob(double x)
	{int k,mult=0;double y,sum,l1,l2;
	if(x==0.)return 0.;
	if(x==pi*.5)return pi*.5*log(2.);
	if(x>pi)return pi*log(2.)+Lob(x-pi);
	if(x<0.){x=-x;mult=1;}
	if(x>pi*.25 && x<pi*.5)/* pi/4<x<pi/2*/
		{y=pi*.5-x;/*   0< y < pi*.25 */
		sum=pi*.5-2.*y;/*  0 <sum <pi*.5 */
		if(sum<=0. || y<= 0.){printf(" bad Lob %le %le\n",sum,y);
						return errorcode;}
		l1=Lob(y);
		l2=Lob(sum);
		sum= l1+.5*l2-(pi*.25-x)*log(2.);
		return mult?-sum:sum;
		}
	if(x>pi*.5 && x<=pi)
		{return pi*log(2.)-Lob(pi-x);}
	/* at this point, 0<x<pi/4*/
	y=x*x;sum=lcoef[24];
	for(k=23;k>=0;k--){sum= sum*y+lcoef[k];	}
	sum*=x*y;
	return mult? sum:-sum;
	}