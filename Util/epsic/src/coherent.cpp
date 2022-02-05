/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "sample.h"
#include "Pauli.h"

using namespace std;

Spinor<double> spinor (const Stokes<double>& stokes)
{
  cerr << "Stokes=" << stokes << endl;
  Quaternion<double, Hermitian> H = natural(stokes);
  cerr << "Natural=" << H << endl;
  Quaternion<double, Unitary> U = eigen (H);
  cerr << "Unitary=" << U << endl;

  Jones<double> jones = herm( convert(U) );
  cerr << "Jones=" << jones << endl;
  
  Spinor<double> result (jones.j00, jones.j10);
  result *= sqrt(stokes[0]);
  return result;
}

coherent::coherent (double _coh)
{
  built = false;
  coupling = new mode;

  coherence = _coh;
  Stokes<double> S (2,0,2*coherence,0);
  coupling->set_Stokes (S);

  a_xform = b_xform = 0;
}

void coherent::set_normal (BoxMuller* n)
{
  coupling->set_normal(n);
  combination::set_normal (n);
}

Stokes<double> coherent::get_Stokes ()
{
  if (!built)
    build ();

#if 1
  double phi = drand48() * 2*M_PI;
  Stokes<double> S (1,0, coherence*cos(phi), coherence*sin(phi));
  coupling->set_Stokes (2.0*S);
#endif

  Stokes<double> result;
  for (unsigned i=0; i<sample_size; i++)
  {
    Spinor<double> amps = coupling->get_field();
    Spinor<double> a_e = amps.x * a;
    if (a_xform)
      a_e = a_xform->transform(a_e);

    Spinor<double> b_e = amps.y * b;
    if (b_xform) 
      b_e = b_xform->transform(b_e); 

    Vector<4, double> tmp;
    compute_stokes (tmp, a_e + b_e);
    result += tmp;
  }
  result /= sample_size;
  return result;
}

void coherent::build ()
{
  a = spinor (A->get_Stokes());
  b = spinor (B->get_Stokes());

  a_xform = dynamic_cast<field_transformer*> (A);
  b_xform = dynamic_cast<field_transformer*> (B);

  built = true;
}

Vector<4, double> coherent::get_mean ()
{
  return A->get_mean() + B->get_mean();
}

//! Implements Equation (42) of van Straten & Tiburzi (2017)
Matrix<4,4, double> coherent::get_covariance ()
{
  Matrix<4,4, double> result = sample::get_covariance (A, sample_size);
  result += sample::get_covariance (B, sample_size);
  
  Stokes<double> mean_A = A->get_mean();
  Stokes<double> mean_B = B->get_mean();
  
  /* Minkowski::outer implements A \otimes B - 0.5 \eta A \cdot B
     such that Minkowski::outer(A,B) +  Minkowski::outer(B,A)
     yields Equation (43) of van Straten & Tiburzi (2017) */
  Matrix<4,4, double> xcovar = Minkowski::outer(mean_A, mean_B);
  xcovar /= sample_size;
  result += xcovar + transpose(xcovar);
  
  return result;
}
