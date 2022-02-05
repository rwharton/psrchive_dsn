/***************************************************************************
 *
 *   Copyright (C) 2011 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Dirac.h"
#include "Jacobi.h"

using namespace std;

int main () 
{
  // test that Dirac matrices are orthogonal to each other

  for (unsigned i=0; i<4; i++)
   {
    for (unsigned j=0; j<4; j++)
    {
      Dirac::type M1 = Dirac::matrix (i, j);

      for (unsigned k=0; k<4; k++)
      {
        for (unsigned l=0; l<4; l++)
        {
          Dirac::type M2 = Dirac::matrix (k, l);

          // Dirac::type R = M1 * M2;

          complex<double> tr = trace(M1 * M2);

          complex<double> expect (0,0);
          if (i==k && j==l)
            expect = complex<double> (4,0);

          if (tr != expect)
          {   
            cerr << "fail " << i << " " << j << " " << k << " " << l << " " << tr << endl;
            return -1;
          }
        }
      }
    }
  }


  for (unsigned i=0; i<4; i++)
    for (unsigned j=0; j<4; j++)
    {
      Dirac::type M1 = Dirac::matrix (i, j);

      // test the mixed product property of the Kronecker product
      Dirac::type M2 = Dirac::matrix(i,0) * Dirac::matrix(0,j);
      if (M1 != M2)
      {
        cout << i << " " << j << endl << M1 << endl;
        cout << "!=\n" << M2 << endl << endl;
        return -1;
      }

      // test for symmetry
      bool symmetric = true;

      for (unsigned k=0; k<4 && symmetric; k++)
        for (unsigned l=k+1; l<4 && symmetric; l++)
          if ( M1[k][l] != M1[l][k] )
          {
            cerr << "not symmetric i=" << i << " j=" << j << endl;
            cerr << M1 << endl;
            symmetric = false;
          }
    }

  for (unsigned i=0; i<2; i++)
  {
    for (unsigned j=1; j<4; j++)
    {
      Dirac::type covar (0);
      covar[0][0] = covar[j][j] = 1.0;

      if (i == 0)
	covar[j][0] = covar[0][j] = 1.0;

    cout << "covar=\n" << covar << endl;

    Dirac::type N (0);
    for (unsigned i=0; i<4; i++)
      for (unsigned j=0; j<4; j++)
	N += Dirac::type(covar[i][j]) * Dirac::matrix (i, j);

    cout << "N=\n" << N << endl << endl;

    Dirac::type evec;
    Vector<4,double> eval;
    Jacobi (N, evec, eval);

    for (unsigned i=0; i<4; i++)
      cout << eval[i] << " " << evec[i] << endl;
    }
  }

  return 0;
}
