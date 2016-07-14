#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <gsl/gsl_multimin.h>

typedef std::vector<double> vector;
typedef std::string string;
typedef double (*my_f_ptr)  (const gsl_vector *v, void *params);
typedef void   (*my_df_ptr) (const gsl_vector *v, void *params, gsl_vector *df);
typedef void   (*my_fdf_ptr)(const gsl_vector *x, void *params,double *f, gsl_vector *df);

class gsl_minimizer {
public:
   gsl_minimizer(const my_f_ptr   my_f,
                 const my_df_ptr  my_df,
                 const my_fdf_ptr my_fdf,
                 const vector&    params,
                 const int        dim);
   ~gsl_minimizer();
   void solve (const string  solver,
                     vector& arr,
               const int     iteration,
               const double  tol=1e-5,
               const double  step=1e-2);
private:
   const my_f_ptr   _my_f;
   const my_df_ptr  _my_df;
   const my_fdf_ptr _my_fdf;
   double*    par;
   const int  _dim;
};

gsl_minimizer::gsl_minimizer (const my_f_ptr   my_f,
                              const my_df_ptr  my_df,
                              const my_fdf_ptr my_fdf,
                              const vector&    params,
                              const int        dim)
:_my_f(my_f),
 _my_df(my_df),
 _my_fdf(my_fdf),
 par(NULL),
 _dim(dim)
{
   par = new double [params.size()];
   for(int i=0;i<params.size();++i) par[i] = params[i];
}

gsl_minimizer::~gsl_minimizer ()
{
   delete [] par; par = NULL;
}


void gsl_minimizer::solve (const string solver,
                                 vector& arr,
                           const int iteration,
                           const double tol,
                           const double step) 
{
   int iter = 0;
   int status;

   const gsl_multimin_fdfminimizer_type *T; 
   gsl_multimin_fdfminimizer *s; 

   gsl_vector *x;
   gsl_multimin_function_fdf my_func;

   my_func.n      = _dim;
   my_func.f      = _my_f;
   my_func.df     = _my_df;
   my_func.fdf    = _my_fdf;
   my_func.params = par;

   x = gsl_vector_alloc(_dim);
   for(int i=0;i<_dim;++i) gsl_vector_set (x, i, arr[i]);

   if("bfgs2" == solver) 
      T = gsl_multimin_fdfminimizer_vector_bfgs2;
   else if("bfgs" == solver) 
      T = gsl_multimin_fdfminimizer_vector_bfgs;
   else if("steepest_descent" == solver)
      T = gsl_multimin_fdfminimizer_steepest_descent;
   else if("conjugate_pr" == solver)
      T = gsl_multimin_fdfminimizer_conjugate_pr;
   else if("conjugate_fr" == solver)
      T = gsl_multimin_fdfminimizer_conjugate_fr;
   else {
      std::cout<<"there are no corresponded to algorithms what you entered."<<std::endl;
      abort();
   }

   s = gsl_multimin_fdfminimizer_alloc (T, _dim);

   gsl_multimin_fdfminimizer_set (s, &my_func, x, step, tol);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        break;

      status = gsl_multimin_test_gradient (s->gradient, tol);

      if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");

        printf("%5d ",iter);
        for(int i=0;i<_dim;++i) printf("%.5f ",gsl_vector_get (s->x, i));
        printf("%10.5f\n",s->f);

    }
  while (status == GSL_CONTINUE && iter < iteration);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
}
