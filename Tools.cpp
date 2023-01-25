#include "Tools.h"
#include <math.h>
#define pi acos(-1)
#define eps 2.2204e-16
using namespace std;

double Func_source(const double &x)
{
  return exp(x);
}
void Gauss(int N, double a, double b, double * qp, double * wq)
{
  N -= 1;
  const int N1 = N + 1;
  const int N2 = N + 2;
  double xu[N1];
  for (int ii = 0; ii < N1; ii++)
  {
    xu[ii] = -1.0 + double(ii) * 2.0/(double(N1)-1.0);
  }
  double y[N1];
  for (int ii = 0; ii < N1; ii++)
  {
    y[ii] = cos((2.0*double(ii)+1)*pi/(2.0*double(N)+2.0))+(0.27/double(N1))*sin(pi*xu[ii]*double(N)/double(N2));
  } 
  double L[N1][N2];
  double * Lp = new double[N1];
  double y0[N1];
  double error = 1.0;

  while (error > eps)
  {
    for (int ii = 0; ii < N1; ii++)
    {
      L[ii][0] = 1.0;
      L[ii][1] = y[ii]; 
    }
    for (int ii = 1; ii < N1; ii++)
    {
      for (int jj = 0; jj < N1; jj++)
      {
        L[jj][ii+1] = ( (2.0*double(ii)+1.0) * y[jj] * L[jj][ii] - double(ii)*L[jj][ii-1] ) / double(ii+1);
      }
    }
    for (int ii = 0; ii < N1; ii++)
    {
      Lp[ii] = double(N2) * (L[ii][N1-1] - y[ii]*L[ii][N2-1] ) / (1.0 - y[ii]*y[ii]);
    }
    for (int ii = 0; ii < N1; ii++)
    {
      y0[ii] = y[ii];
      y[ii]  =y0[ii] - L[ii][N2-1] / Lp[ii];
    }

    double error0 = 0.0;
    for (int ii = 0; ii < N1; ii++)
    {
      error = (error0 > abs( y[ii]-y0[ii] )) ? error : abs(y[ii]-y0[ii]);
    }
    error0 = error;
  }
  for (int ii = 0; ii < N1; ii++)
  {
    qp[ii] = (a*(1.0-y[ii])+b*(1.0+y[ii]))/2.0;
    wq[ii] = (b-a) / ((1.0-y[ii]*y[ii]) * Lp[ii] * Lp[ii]) * (double(N2)/double(N1))*(double(N2)/double(N1));
  }
  return;
}
// Cubic Hermite-shape function
double HermiteBasis(const double &x1, const double &x2, int i, int der, double &x)
{
  const double hh = abs(x2 - x1);
  const double h3 = hh * hh * hh;
  const double h2 = hh * hh;
  double hermite = 0.0;
  switch (i) {
    case 1:
      // N1(x) -- u1
      switch (der) {
        case 0:
          hermite = (1.0/h3) * (x-x2) * (x-x2) * ( 2.0*(x-x1) + hh );
          break;
        case 1:
          hermite = (2.0 * (x-x2) / h3) * ( (x-x2) + 2.0*(x-x1) + hh );
          break;
        case 2:
          hermite = ( 2.0 / h3) * ( 4.0 * (x-x2) + 2.0 * (x-x1) + hh );
          break;
        case 3:
          hermite = 12.0 / h3;
          break;
        case 4:
          hermite = 0.0; 
          break;
      }
      break;
    case 2:
      // N2(x) -- u1'
      switch (der) {
        case 0:
          hermite = (x-x1) * (x-x2) * (x-x2) / h2;
          break;
        case 1:
          hermite = ( (x-x2) / h2) * ( 2.0*(x-x1) + (x-x2) );
          break;
        case 2:
          hermite = (1.0/h2) * (4.0*(x-x2) + 2.0*(x-x1));
          break;
        case 3:
          hermite = 6.0 / h2;
          break;
        case 4:
          hermite = 0.0;
          break;
      }
      break;
    case 3:
      // N3(x) -- u2
      switch (der) {
        case 0:
          hermite = ( 1.0/h3 ) * (x-x1) * (x-x1) * ( 2.0*(x2-x) + hh );
          break;
        case 1:
          hermite = ( 2.0 * (x-x1) / h3 ) * ( (x1-x) + 2.0 * (x2-x) + hh );
          break;
        case 2:
          hermite = ( 2.0/h3 ) * ( 4.0 * (x1-x) + 2.0*(x2-x) +hh );
          break;
        case 3:
          hermite = -12.0 / h3;
          break;
        case 4:
          hermite = 0.0;
          break;
      }
      break;
    case 4:
      // N4(x) -- h2'
      switch (der) {
        case 0:
          hermite = (x-x2) * (x-x1) * (x-x1) / h2;
          break;
        case 1:
          hermite = ( (x-x1) / h2 ) * ( (x-x1) + 2.0 * (x-x2) );
          break;
        case 2:
          hermite = ( 2/h2) * ( 2.0*(x-x1) + (x-x2) );
          break;
        case 3:
          hermite = 6.0 / h2;
          break;
        case 4:
          hermite = 0.0;
          break;
      }
      break;
  }
  return hermite;
}
