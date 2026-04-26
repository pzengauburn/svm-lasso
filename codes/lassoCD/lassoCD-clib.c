/********************************************************************
 * Coordinate Descent Algorithm for lasso  
 * 05-20-2011
 ********************************************************************/

#include "mydef.h"

/********************************************************************
 * sum(x * y)
 ********************************************************************/

double sum2(const int *n, const double *x, const double *y)
{
  double ans = 0.0;
  for(int i = 0; i < (*n); i++) ans += (x[i] * y[i]);
  return(ans);
}

/********************************************************************
 * calculate sign(a) * (abs(a) - lambda)_+ / s
 ********************************************************************/

double uni_lasso_CD(double s, double a, double lambda)
{
  double ans = fabs(a) - lambda;

  if(s < EPS) error("uni_lasso_CD(): s < EPS. stop.\n");

  if(ans > 0.0) ans *= (-sign(a) / s);
  else ans = 0.0;

  return(ans);
}

/********************************************************************
 * update the k-th component of x
 ********************************************************************/

double lasso_CDk(int k, const int *p, 
     const double *Smat, const double *avec, const double *lambda, 
     double *x)
{
  double a;

  x[k] = 0.0;
  a = sum2(p, Smat + k * (*p), x) - avec[k];
  x[k] = uni_lasso_CD(Smat[k + k * (*p)], a, *lambda);

  return(x[k]);
}


void lasso_CD_core(const int *p, 
     const double *Smat, const double *avec, const double *lambda,
     const int *maxiter, const double *tol,
     double *x, int *iter, double *conv)
{
  double oldx, newx;
  int k;

  (*iter) = 0; (*conv) = (*tol) + 1.0;
  while(((*iter) < (*maxiter)) & ((*conv) > (*tol)))
  {
    (*conv) = 0.0;
    for(k = 0; k < (*p); k++)
    { 
      oldx = x[k];
      newx = lasso_CDk(k, p, Smat, avec, lambda, x);
      (*conv) = fmax(fabs(oldx - newx), *conv);
    }

    (*iter) ++;
  }
}

void lasso_CD_core2(const int *p, 
     const double *Smat, const double *avec, const double *lambda,
     const int *maxiter, const double *tol,
     double *x, int *iter, double *conv)
{
  double oldx, newx;
  int k;

  (*iter) = 0; (*conv) = (*tol) + 1.0;
  while(((*iter) < (*maxiter)) & ((*conv) > (*tol)))
  {
    (*conv) = 0.0;
    for(k = 0; k < (*p); k++)
    { 
      oldx = x[k];
      newx = lasso_CDk(k, p, Smat, avec, lambda + k, x);
      (*conv) = fmax(fabs(oldx - newx), *conv);
    }

    (*iter) ++;
  }
}

/********************************************************************
 * END OF THE FILE
 ********************************************************************/
