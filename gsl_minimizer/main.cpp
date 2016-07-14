#include <iostream>
#include "gsl_minimizer.h"

enum {count = 30};

double xsquare(const gsl_vector *v, void *params) 
{
   double result = 0;
   vector x(count,0);
   for(int i=0;i<count;++i) x[i] = gsl_vector_get(v, i);
   for(int i=0;i<count;++i) result += pow((x[i]-i),2);
   return result;   
}

void dxsquare(const gsl_vector *v, void *params, gsl_vector *df) 
{
   vector x(count,0);
   for(int i=0;i<count;++i) x[i] = gsl_vector_get(v, i);
   double* p = static_cast<double*>(params);

   for(int i=0;i<count;++i) gsl_vector_set(df,i,2*(x[i]-i));
}

void square_fdf(const gsl_vector *x, void *params,double *f, gsl_vector *df) {
  *f = xsquare(x, params); 
  dxsquare(x, params, df);
}

int main(int argc,char* argv[]) {
   vector arr(count,0);
   vector par(1,0);
   gsl_minimizer solver(xsquare,dxsquare,square_fdf,par,count);

   solver.solve("steepest_descent",arr,1000);
   return 0;
}

