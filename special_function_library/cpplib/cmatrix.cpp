/*  C++ vector/matrix support routines 
	(& main test driver)

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include <stdlib.h>
#include "cmatrix.hpp"

error(char* msg,int index,int s)
	{
	printf("%s%d %d\n",msg,index,s);
	exit(1);
	}

void Cvector::check( int index)
	{int loc;
	loc=index-(*this).base;
	if (loc<0)
		error(" Cvector error: index-base<0",index,(*this).base);
	if (loc>= (*this).size)
		error(" Cvector error: index too large",index,(*this).size);
	}

void Cmatrix::check( int i,int j)
	{int loc;
	loc=i-(*this).base;
	if (loc<0)
		error(" Cmatrix error: row index-base<0",i,(*this).base);
	if (loc>= (*this).rowkt)
		error(" Cmatrix error: row index too large",i,(*this).rowkt);

	loc=j-(*this).base;
	if (loc<0)
		error(" Cmatrix error: column index-base<0",j,(*this).base);
	if (loc>= (*this).colkt)
		error(" Cmatrix error: column index too large",j,(*this).colkt);

	}


void Cmatrix::init(int row, int col, int b)
	{
	if(row<1||col<1)error(" matrix needs at least one row/col ",row,col);
	rowkt=row;colkt=col;base=b;
	m= new complex *[rowkt];
	for(int j=0;j<rowkt;j++)m[j]= new complex[colkt];
	complex zero;complex one(1.,0.);
	//initialize to identity matrix. 
	for(int i=0;i<rowkt;i++){for(j=0;j<colkt;j++)m[i][j]=zero;
		m[i][i]= one;// omit this line for init. to zero matrix
		}
	}

Cmatrix::Cmatrix(Cmatrix& a)
	{ init(a.rowkt,a.colkt,a.base);
	 for(int i=0;i<a.rowkt;i++)for(int j=0;j<a.colkt;j++)m[i][j]=a.m[i][j];
	}

Cmatrix::~Cmatrix()
	{for(int i=0;i<rowkt;i++)delete m[i];
	delete m;
	}


complex& Cvector::element(int i)
	{
	check(i);return head[i-base];
	}

void Cvector::setelement( int i, complex& value)
	{
	check(i);
	head[i]=value;
	}
complex& Cmatrix::element(int i,int j)
	{
	check(i,j);return m[i-base][j-base];
	}

void Cmatrix::setelement( int i, int j, complex value)
	{
	check(i,j);m[i][j]=value;
	}


void Cmatrix::operator=(Cmatrix& a)
	{
	if(this==&a)return;
	for(int i=0;i<rowkt;i++)delete m[i];
	delete m;
	init(a.rowkt,a.colkt,a.base);
	for(i=0;i<a.rowkt;i++)for(int j=0;j<a.colkt;j++)m[i][j]=a.m[i][j];
	}

 Cvector::Cvector( Cvector& orig) //copy
	{
	size=orig.size;base=orig.base;head=new complex[orig.size];
	for(int i=0;i<size;i++)head[i]=orig.head[i];
	}

Cmatrix operator*(Cmatrix&a, Cmatrix& b)
	{
	int i,j,k;complex zero;
	if(a.colkt!=b.rowkt)
		error(" cannot multiply matrices colkt,rowkt ",a.colkt,b.rowkt);
	Cmatrix prod(a.rowkt,b.colkt,a.base);
	for(i=0;i<a.rowkt;i++)for(j=0;j<b.colkt;j++)
		{prod.m[i][j]=zero;
		for(k=0;k<a.colkt;k++)prod.m[i][j] += a.m[i][k]*b.m[k][j];
		}
	return prod;
	}	


Cvector operator*(Cmatrix& a, Cvector& b)
	{
	int k;complex zero;
	if(a.colkt!=b.size)error(" cannot mult. matrix,vector ",a.colkt,b.size);
	Cvector prod(b.size,b.base,zero);
	for(int i=0;i<a.rowkt;i++)
		{prod.head[i]=zero;
		for(k=0;k<a.colkt;k++)prod.head[i] += (a.m[i][k])*(b.head[k]);
		}
	return prod;
	}	

complex Cvector::operator*(Cvector& rvalue) //dot product
	{int i; complex *element,*relement; complex sum(0.,0.);// sum 
	if(rvalue.size!= size)
		error(" dot product unequal length vectors ",size,rvalue.size);
	for(i=0,element=head,relement=rvalue.head;i< size;
		i++,element++,relement++)
		sum += *element * *relement;
	return sum;
	}

void Cvector::operator=(  Cvector& rhs)
	{
	if( this== &rhs)return;
	if(size!=rhs.size)
		{
		delete head;
		head= new complex[rhs.size];
		}
	for(int j=0;j<size;j++){head[j]=rhs.head[j];}
	}

main() 
{
complex a,b(1.,0.),c(0.,1.); double dot; int i,j;
printf(" size of complex=%d\n",sizeof(complex));
a=b+c;
a.print("summand=","\n");
printf(" magnitude=%e\n", a.abs());
b=a/c;
b.print("quotient=","\n");
a=b.cexp();
a.print(" exponential=","\n");
a= complex(1.,0.);b=complex(3.,2.);
printf(" size of Cvector=%d\n",sizeof(Cvector));
//
	{Cvector A(4,0,a);Cvector B(4,0,b);Cvector D(4,0,b);
for(j=0;j<B.size;j++) B.element(j).print(" B=","\n");	
	a= A*B;
	B=A;
for(j=0;j<B.size;j++) B.element(j).print(" B=","\n");	
	a.print(" dot product=","\n");
		{Cmatrix m(4,4,0),n(4,4,0),o(4,4,0);
		c=complex(0.,1.);
		for (i=0;i<4;i++)
			for(j=0;j<4;j++){
				m.m[i][j].print(" matrix element=","\n");}
		for (i=0;i<4;i++)
			for(j=0;j<4;j++){
				n.setelement( i,j, complex((double)(i+j),0.));
				n.m[i][j].print(" matrix element=","\n");}
		D = (m * A);
for(j=0;j<B.size;j++) B.element(j).print(" B=","\n");	
		B = (m * A);
for(j=0;j<B.size;j++) B.element(j).print(" B=","\n");	
		o= m * n;

		for (i=0;i<4;i++)
			for(j=0;j<4;j++){
				o.m[i][j].print(" matrix element=","\n");}
		printf(" end of inner block\n");
		}
printf(" end of block\n");
	}
printf(" end of program\n");
}


