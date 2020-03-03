#define _USE_MATH_DEFINES
#include "BoundaryMesh.hpp"
#include "buildK.hpp"
#include "buildM.hpp"
#include "buildV.hpp"
#include "buildW.hpp"
#include "evaluateK.hpp"
#include "evaluateV.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>
using namespace Eigen;

namespace IndirectSecondKind {

// 1.9.k -> solving dirichlet problem using indirect second kind BIE

template <typename FUNC>
Eigen::VectorXd solveDirichlet(const BoundaryMesh &mesh, const FUNC &g) {
  double eta = 1e-5;
  MatrixXd K;
  computeK00(K, mesh, eta);

  SparseMatrix<double> M00(mesh.numElements(), mesh.numElements());
  MatrixXd M0;
  computeM00(M00, mesh);
  M0 = MatrixXd(M00);

  SparseMatrix<double> M01(mesh.numElements(), mesh.numElements());
  MatrixXd M1;
  computeM01(M01, mesh);
  M1 = MatrixXd(M01);

  int N = mesh.numVertices();
  VectorXd gNcoeff(N);
  for (int i = 0; i < N; ++i)
    gNcoeff(i) = (g(mesh.getVertex(i)));

  VectorXd sol;
  VectorXd rhs = M1 * gNcoeff;
  MatrixXd lhs = -.5 * M0 + K;
  sol = lhs.lu().solve(rhs);
  return sol;
}

// 1.9.i -> reconstructing the solution
double reconstructSolution(const Vector2d &X, const VectorXd &f,
                           const BoundaryMesh &mesh) {
  double eta = 1e-5;
  Eigen::VectorXd DLf_x(1);
  evaluateK(DLf_x, mesh, f, X.transpose(), eta);

  return DLf_x(0);
}
} // namespace IndirectSecondKind

namespace IndirectFirstKind {

// 1.9.h -> solving dirichlet problem using indirect first kind BIE

template <typename FUNC>
Eigen::VectorXd solveDirichlet(const BoundaryMesh &mesh, const FUNC &g) {
  double eta = 1e-5;
  MatrixXd V;
  computeV(V, mesh, eta);

  SparseMatrix<double> M01(mesh.numElements(), mesh.numElements());
  MatrixXd M;
  computeM01(M01, mesh);
  M = MatrixXd(M01);

  int N = mesh.numVertices();
  VectorXd gNcoeff(N);
  for (int i = 0; i < N; ++i)
    gNcoeff(i) = g(mesh.getVertex(i));

  VectorXd sol;
  VectorXd rhs = M * gNcoeff;
  sol = V.lu().solve(rhs);
  return sol;
}

// 1.9.i -> reconstructing the solution
double reconstructSolution(const Vector2d &X, const VectorXd &phi,
                           const BoundaryMesh &mesh) {
  double eta = 1e-5;
  Eigen::VectorXd SLphi_x(1);
  evaluateV(SLphi_x, mesh, phi, X.transpose(), eta);

  return SLphi_x(0);
}
} // namespace IndirectFirstKind

namespace DirectSecondKind {

// 1.9.e -> solving dirichlet problem using direct second kind BIE

template <typename FUNC>
Eigen::VectorXd solveDirichlet(const BoundaryMesh &mesh, const FUNC &g) {
  double eta = 1e-5;
  MatrixXd W;
  computeW(W, mesh, eta);

  MatrixXd K;
  computeK(K, mesh, eta);

  SparseMatrix<double> M01(mesh.numElements(), mesh.numElements());
  MatrixXd M;
  computeM01(M01, mesh);
  M = MatrixXd(M01);

  int N = mesh.numVertices();
  VectorXd gNcoeff(N);
  for (int i = 0; i < N; ++i)
    gNcoeff(i) = g(mesh.getVertex(i));

  MatrixXd MT = M.transpose();
  MatrixXd KT = K.transpose();

  VectorXd sol;
  MT = 0.5 * MT;
  MT -= KT;
  VectorXd rhs = W * gNcoeff;
  sol = MT.lu().solve(rhs);
  return sol;
}
} // namespace DirectSecondKind

// definition of the function gamma
struct Kite { // kite function for a scalar input
  VectorXd operator()(double &t) const {
    assert(t <= 1 && t >= -1);
    VectorXd value(2);
    value << 0.25 * std::cos(M_PI * t) + 0.1625 * std::cos(2 * M_PI * t),
        0.375 * std::sin(M_PI * t);
    return value;
  }
  // kite function for a vector input
  MatrixXd operator()(VectorXd &input) const {
    int N = input.size();
    MatrixXd output(N, 2);
    for (int i = 0; i < N; ++i)
      output.row(i) = this->operator()(input(i));
    return output;
  }
};

double gfunc(VectorXd input) {
  assert(input.size() == 2);
  double val = std::sin(input(0) - input(1)) * std::sinh(input(0) + input(1));
  return val;
}

/*
//construct the BoundaryMesh object
struct BoundaryMesh{

        public:
        int nE; //no. of panels
        int nC; //no. of vertices
        MatrixXd Vertices;

        //VectorXd<double,8,1> y;
        //constructor
        BoundaryMesh(const int& N,MatrixXd&
inpvert):nE(N),nC(N),Vertices(inpvert)
        {
        //MatrixXd Vertices(N,2);
        //Vertices=MatrixXd::Zero(4*N,2);
        //x=VectorXd::Zero(4*N);
        //y==VectorXd::Zero(4*N);
        }

};
*/

// 1.9.a -> sets up a uniform mesh of the boundary of the square [0,1/2]^2
BoundaryMesh createMiniSquareMesh(const int &N) {
  assert(N > 0);
  VectorXd sq_side = VectorXd::LinSpaced(N + 1, 0, 1. / 2.);
  VectorXd zero_side = VectorXd::Zero(N + 1);
  VectorXd const_side = VectorXd::Constant(N + 1, 1. / 2.);

  // anticlockwise order for boundaries gamma(i). segment on x axis is gamma(1)

  VectorXd Verticesx(4 * N);
  VectorXd Verticesy(4 * N);
  VectorXd reverse = sq_side.reverse();
  Verticesx << sq_side, const_side.segment(1, N), reverse.segment(1, N),
      zero_side.segment(1, N - 1);
  Verticesy << zero_side, sq_side.segment(1, N), const_side.segment(1, N),
      reverse.segment(1, N - 1);
  MatrixXd temp(4 * N, 2);
  temp << Verticesx, Verticesy;
  MatrixXi temp1(4 * N, 2);
  temp1 << VectorXi::LinSpaced(4 * N, 0, 4 * N - 1),
      VectorXi::LinSpaced(4 * N, 1, 4 * N);
  temp1(4 * N - 1, 1) = 0;
  BoundaryMesh mesh(temp, temp1);
  return mesh;
}

// 1.9.b -> creating mesh with a given parametrization
template <typename PARAM>
BoundaryMesh createMeshwithGamma(const PARAM &gamma, const int &N) {
  assert(N > 0);
  VectorXd t = VectorXd::LinSpaced(N + 1, -1., 1.);
  VectorXd t1 = t.segment(0, N);
  MatrixXd temp(N, 2);
  for (unsigned i = 0; i < N; ++i)
    temp.row(i) = gamma(t1(i));

  MatrixXi temp1(N, 2);
  // std::cout << "Test" <<std::endl;
  temp1 << VectorXi::LinSpaced(N, 0, N - 1), VectorXi::LinSpaced(N, 1, N);
  temp1(N - 1, 1) = 0;
  // std::cout << temp1 << std::endl;
  BoundaryMesh mesh(temp, temp1);
  return mesh;
}

// 1.9.c -> solving dirichlet boundary value problem

template <typename FUNC>
Eigen::VectorXd solveDirichlet(const BoundaryMesh &mesh, const FUNC &g) {
  double eta = 1e-5;
  MatrixXd V;
  computeV(V, mesh, eta);

  MatrixXd K;
  computeK(K, mesh, eta);

  SparseMatrix<double> M01(mesh.numElements(), mesh.numElements());
  MatrixXd M;
  computeM01(M01, mesh);
  M = MatrixXd(M01);

  int N = mesh.numVertices();
  VectorXd gNcoeff(N);
  for (int i = 0; i < N; ++i)
    gNcoeff(i) = g(mesh.getVertex(i));

  VectorXd sol;
  M = 0.5 * M;
  M += K;
  VectorXd rhs = M * gNcoeff;
  sol = V.lu().solve(rhs);
  return sol;
}

double TNu_pt(VectorXd &point, VectorXd &a, VectorXd &b) {
  // assert(point.size()==2);
  VectorXd grad_u(2); // gradient of u
  double X1 = point(0), X2 = point(1);
  grad_u << std::sin(X1 - X2) * std::cosh(X1 + X2) +
                std::cos(X1 - X2) * std::sinh(X1 + X2),
      std::sin(X1 - X2) * std::cosh(X1 + X2) -
          std::cos(X1 - X2) * std::sinh(X1 + X2);
  VectorXd n(2); // normal vector
  VectorXd temp = VectorXd(b - a);
  n(0) = temp(1);
  n(1) = -temp(0);
  // n << 0.375* std::cos(M_PI*t) , 0.25*cos(M_PI*t)+2*.1625*cos(2*M_PI*t);
  n /= n.norm();
  // n=-n;
  double traceval = grad_u.dot(n);
  return traceval;
}

// 1.9.d -> error calculation
void error_calc(int &N, Kite &gamma, double (*fun)(VectorXd)) {
  BoundaryMesh mesh = createMeshwithGamma(gamma, N);
  VectorXd solcoeffs = solveDirichlet(mesh, fun); // the MU vector
  // build a vector for the coefficients of neumann trace of given solution.
  // Using the trace value at the midpoint of each panel int
  // N=mesh.numVertices();
  VectorXd TNu(N);
  // VectorXd t=VectorXd::LinSpaced(N+1,-1.,1.);
  for (int i = 0; i < N; ++i) {
    VectorXd a = mesh.getVertex(i % N), b = mesh.getVertex((i + 1) % N);
    VectorXd point = (a + b) / 2;
    TNu(i) = TNu_pt(point, a, b);
  }
  // for computing the Av norm, V is required and for L2 norm, M is required
  double eta = 1e-5;
  MatrixXd V;
  computeV(V, mesh, eta);

  SparseMatrix<double> M00(mesh.numElements(), mesh.numElements());
  MatrixXd M;
  computeM00(M00, mesh);
  M = MatrixXd(M00);
  // std::cout << "solcoeffs: \n" << solcoeffs <<std::endl;
  // std::cout << "Neumann trace: \n" << TNu <<std::endl;
  VectorXd err_coeff = VectorXd(solcoeffs - TNu);
  // std::cout << "Error coeffs: \n" << err_coeff <<std::endl;
  std::cout << "Error calculated by Av norm: "
            << err_coeff.transpose().dot(V * err_coeff)
            << " Error calculated by L2 norm: "
            << (err_coeff.transpose()).dot((M * err_coeff)) << " for N: " << N
            << std::endl;

  // std::cout << solcoeffs << std::endl;
}

// 1.9.g -> test mass matrix
void testMassMatrixSVD(const BoundaryMesh &mesh) {
  SparseMatrix<double> M01(mesh.numElements(), mesh.numElements());
  MatrixXd M;
  computeM01(M01, mesh);
  M = MatrixXd(M01);
  MatrixXd MT = M.transpose();
  Eigen::JacobiSVD<MatrixXd> svd(MT, ComputeThinU | ComputeThinV);
  VectorXd sv = svd.singularValues();
  int N = mesh.numVertices();
  if (sv(N - 1) < 1e-5)
    std::cout << "Singular Matrix!" << std::endl;
}

// 1.9.d -> error calculation
void ind_error_calc(int N, Kite &gamma, double (*fun)(VectorXd)) {
  Vector2d X;
  X << 0., 0.3;
  BoundaryMesh mesh = createMeshwithGamma(gamma, N);
  VectorXd solcoeffs = solveDirichlet(mesh, fun);
  // double uapprx = IndirectSecondKind::reconstructSolution(X,solcoeffs,mesh);
  // double uex = fun(X);
  // double error = uex-uapprx;
  // if (error<0)
  //	error = - error;
  std::cout << "solcoeffs: \n" << solcoeffs << std::endl;
}

int mainfunc() {
  Kite gamma;
  // int N=15;
  //	BoundaryMesh temp = createMiniSquareMesh(4);
  // oundaryMesh mesh=createMeshwithGamma(gamma,N);
  // VectorXd solcoeffs=solveDirichlet(mesh,gfunc);
  // for (int N=10;N<3201;N=N*2)
  //{
  // BoundaryMesh mesh = createMeshwithGamma(gamma,N);
  // testMassMatrixSVD(mesh);
  // ind_error_calc(10,gamma,gfunc);
  // error_calc(N,gamma,gfunc);
  //}
  int N = 10;
  BoundaryMesh mesh = createMeshwithGamma(gamma, N);
  Eigen::VectorXd sol = IndirectSecondKind::solveDirichlet(mesh, gfunc);
  std::cout << sol << std::endl;
  return 0;
}
