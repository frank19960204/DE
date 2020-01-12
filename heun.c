/*============================*/
/* Heun method for ODE (IVP) */
/*    x' = f(t,x)
/*    x(T0)=X0
/*============================*/
#include <stdio.h>
#include <math.h>

/*--- Calculate from t=T0 to t=T1 ---*/
#define T0   0.0
#define T1   5.0

/*--- Initial data ---*/
#define X0   2.0

/*--- parameter(s) ---*/
#define A  1.0

/*--- division number ---*/
#define N    50

/*---  f(t,x)  ---*/
double f(double t, double x) { return A*x; }
double f_2(double t,double x) {return x - x*x;}

/*--- Exact solution ---*/
double ExactSol(double t) { return exp(A*t); }
double ExactSol_2(double t, double c) { 
  return exp(t)/(exp(t) + c);
}

/*--- Heun's scheme ---*/
double heun(double h, double t, double x)
{
  double k1,k2,y_next;

  k1 = f_2(t,x) ;
  y_next = x+h*f_2(t,x) ;
  k2 = f_2(t+h,y_next) ;

 return(x+h*(k1+k2)/2.0);  // x: X_n, t: t_n
}

/*--- main ---*/
/*
     h: mesh size
     x=X_n: approx. sol.
     t: time
     en=|e_n|: error
     maxen=max_n |e_n|
*/
int main(void)
{
 int n;
  double h=(T1-T0)/N;  
  double x=X0;  
  double t;  
  double maxen=0.0, en;
  double c = 0.0;

  c = (exp(T0) - exp(T0)*X0)/X0 ; 
  
  printf("heun method\n");
 
  for(n=0; n<N; n++) {
    t=T0+n*h;
    en=fabs(x-ExactSol_2(t,c));
    maxen=(maxen>en)? maxen: en ;

    /* display t_n, X_n, x(t_n) and |e_n| */
    printf("%6.4f %16.14f %16.14f %16.10e\n",
           t, x, ExactSol_2(t,c), en); 
    x=heun(h,t,x);  // X_{n+1}=euler(h, t_n, X_n)
    
  }
  
  /* max of |e_n| */
  printf("max_n |e_n|=%16.10e\n",maxen);

}

