// complex.hpp
#include <math.h>
#include <stdio.h>

class complex {
   protected: 
      double x,y;
   public:
//      complex( double xx = 0., double yy = 0.) //create
//         { x=xx;y=yy;}
        complex( double xx , double yy = 0.) //create
         { x=xx;y=yy;}
        complex( ) //create
         { x=0.;y=0.;}
      inline void operator=(complex rvalue)
         {x=rvalue.x;y=rvalue.y;}
      inline void operator-=(complex rvalue)
         {x-=rvalue.x;y-=rvalue.y;}
      inline void operator+=(complex rvalue)
         {x+=rvalue.x;y+=rvalue.y;}
      inline void operator*=(complex rvalue)
         {
          *this=complex(rvalue.x*x-rvalue.y*y, 
            rvalue.x*y+rvalue.y*x);
          //return *this;
         }
      inline void operator*=(double rvalue)
         {
         *this=complex(rvalue*x, rvalue*y);//return *this;
         }
      inline complex operator+(complex rvalue)
         {return complex(x+rvalue.x,y+rvalue.y);}
      inline complex operator-(complex rvalue)
         {return complex(x-rvalue.x,y-rvalue.y);}
      inline complex operator-() //unary minus
         {return complex(-x,-y);}
      inline complex operator*(complex rvalue)
         {return complex(
          rvalue.x*x-rvalue.y*y,
          rvalue.x*y+rvalue.y*x);}
      inline friend complex operator/(double dividend,complex divisor)
         { double temp;
         temp=1./(divisor.x*divisor.x+divisor.y*divisor.y);
           return complex((dividend*divisor.x)*temp,
           (-dividend*divisor.y)*temp);
         }
      inline complex operator/(complex divisor)
         { 
         double temp;
         temp=1./(divisor.x*divisor.x+divisor.y*divisor.y);
               return complex((divisor.x*x+divisor.y*y)*temp,
            (divisor.x*y-divisor.y*x)*temp);
         }
      inline int operator==(complex rvalue)
         {return (x==rvalue.x && y==rvalue.y);}
      inline double real() {return x;}
      inline double imaginary() {return y;}
      inline complex conjugate()
         {return complex(x,-y);}
      inline friend complex operator*(complex num,double real)
         {return complex(num.x*real,num.y*real);}
      inline friend complex operator*(double real,complex num)
         {return complex(num.x*real,num.y*real);}
      inline friend complex operator+(complex num,double real)
         {return complex(num.x+real,num.y);}
      inline friend complex operator+(double real,complex num)
         {return complex(num.x+real,num.y);}
      inline complex operator+=(double real)
         {return complex(x+=real,y);}
      inline complex operator-=(double real)
         {
         return complex(x-=real,y);}
      inline complex operator++()
         {x+=1.;return *this;}
      inline friend complex operator/(complex num,double real)
         {return complex(num.x/real,num.y/real);}
      inline friend complex operator-(complex num,double real)
         {return complex(num.x-real,num.y);}
      inline friend complex operator-(double real,complex num)
         {return complex(real-num.x,-num.y);}
      double abs();
      complex cexp();
      complex clog();
		complex csin(); complex ccos();
      complex operator^(double expon);
      complex operator^(complex expon);
      friend complex operator^(double base, complex expon);
      void print( char *ahead="",char *behind="");
      complex hurwitz(complex );
   };
