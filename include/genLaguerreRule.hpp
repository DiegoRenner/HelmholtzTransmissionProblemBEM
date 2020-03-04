//****************************************************************************
//
//  Purpose:
//
//   Computation of nodes and weights for generalized Laguerre rules
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
// Source: http://people.sc.fsu.edu/~jburkardt/cpp_src/gen_laguerre_rule/gen_laguerre_rule.html
//****************************************************************************

# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] );

void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double t[], double wts[] );

double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] );

void imtqlx ( int n, double d[], double e[], double z[] );

void parchk ( int kind, int m, double alpha, double beta );

double r8_epsilon ( );

double r8_huge ( );

double r8_sign ( double x );

void r8mat_write ( string output_filename, int m, int n, double table[] );

void rule_write ( int order, string filename, double x[], double w[],  double r[] );

void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b );

void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] );

void timestamp ( );
