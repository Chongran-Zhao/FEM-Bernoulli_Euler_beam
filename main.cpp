#include <iostream>
#include "Tools.h"
#include <Eigen/Eigen>
#include <iomanip>
using namespace Eigen;
using namespace std;

int main()
{
  // Define the physical problem
  // u_xxxx + ll = 0; x in (0,L) 
  // u is displacement
  // u(L) = 0;
  // u_x(L) = 0;
  // EI u_xx(0) = M;
  // EI u_xxx(0) = Q;
  const double EI = 1.0;
  const double M = 0.0/EI;
  const double Q = 0.0/EI;

  const double omega_l = 0.0;
  const double omega_r = 1.0;
  // exact solution can be solved

  const int nElem = 10;
  const int nLocBas = 4;
  const int nFunc = 4 + 2 * (nElem - 1);

  const double hh = (omega_r - omega_l) / nElem;

  double x_coor[nElem+1];
  for (int ii = 0; ii < nElem+1; ii++)
  {
    x_coor[ii] = hh * double(ii);
  }

  int IEN[nLocBas][nElem];
  for (int ee = 0; ee < nElem; ee++)
  {
    for (int aa = 0; aa < nLocBas; aa++)
    {
      IEN[aa][ee] = 2 * ee + aa; 
    } 
  }

  int ID[nFunc];
  for (int ii = 0; ii < nFunc; ii++)
  {
    ID[ii] = ii;
  }
  // assign the ID for the Dirichlet node to be -1
  ID[nFunc-1] = -1;
  ID[nFunc-2] = -1;

  SparseMatrix <double,RowMajor> K(nFunc,nFunc);
  VectorXd F(nFunc);
  VectorXd uu(nFunc);

  for (int ee = 0; ee < nElem; ee++)
  {
    double k_ele[nLocBas][nLocBas];
    double f_ele[nLocBas];
    double x_ele[nLocBas];

    for (int ii = 0; ii < nLocBas; ii++)
    {
      for (int jj = 0; jj < nLocBas; jj++)
      {
        k_ele[ii][jj] = 0.0;
      }
      f_ele[ii] = 0.0;
      x_ele[ii] = 0.0;
    }

    const double x1 = x_coor[ee];
    const double x2 = x_coor[ee+1];

    // quadrature rule
    const int nqp = 10;
    double * qp = new double[nqp];
    double * wq = new double[nqp];
    Gauss(nqp, x1, x2, qp, wq);

    for (int qua = 0; qua < nqp; qua++)
    {
      for (int aa = 0; aa < nLocBas; aa++)
      {
        const double Na    = HermiteBasis(x1, x2, aa+1, 0, qp[qua]);
        const double Na_xx = HermiteBasis(x1, x2, aa+1, 2, qp[qua]);

        f_ele[aa] += wq[qua] * Func_source(qp[qua]) * Na * EI;

        for (int bb = 0; bb < nLocBas; bb++)
        {
          const double Nb_xx = HermiteBasis(x1, x2, bb+1, 2, qp[qua]);
          k_ele[aa][bb] += wq[qua] * Na_xx * Nb_xx * EI;
        }
      }
    }
    // k_ele goes into K and f_ele goes into F, which is called global assembly
    for (int aa = 0; aa < nLocBas; aa++)
    {
      const int AA = ID[IEN[aa][ee]];
      if (AA >= 0)
      {
        F(AA) += f_ele[aa];
        for (int bb = 0; bb < nLocBas; bb++)
        {
          const int BB = IEN[bb][ee];
          K.coeffRef(AA,BB) += k_ele[aa][bb];
        }
      }
      else
      { K.coeffRef(IEN[aa][ee],IEN[aa][ee]) = 1.0;
        F(IEN[aa][ee]) = 0.0;
      }
    }
  }
  F(0) += Q;
  F(1) -= M;
  SparseLU< SparseMatrix<double> > solver;
  solver.compute(K);
  uu = solver.solve(F);

  cout << "The number of elements = " << nElem << endl;
  cout << "The number of nodes = " << nElem+1 << endl;
  cout << "The degree of Hermite-shape function = " << 2 << endl;
  cout << "Displacement     Slope" << endl;
  for (int ii = 0; ii < nElem+1; ii++)
  {
    cout << "u = "; 
    cout << fixed <<setprecision(6) << uu(2*ii);
    cout << "     u_x = " << uu(2*ii+1) << endl;
  }
}
