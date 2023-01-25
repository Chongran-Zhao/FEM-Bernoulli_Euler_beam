#ifndef _TOOL_H_
#define _TOOL_H_

void Gauss(int N, double a, double b, double * qp, double * wq);

double HermiteBasis(const double &x1, const double &x2, int i, int der, double &x);

double Func_source(const double &x);
#endif
