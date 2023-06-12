#include "single_layer.hpp"
#include "cbessel.hpp"

    namespace single_layer_helmholtz {

        typedef std::complex<double> complex_t;
        complex_t ii = complex_t (0.0,1.0);
        double epsilon = std::numeric_limits<double>::epsilon();
        Eigen::MatrixXcd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                           const AbstractParametrizedCurve &pi_p,
                                           const AbstractBEMSpace &space,
                                           const QuadRule &GaussQR,
                                           const QuadRule &CGaussQR,
                                           const complex_t k,
                                           const double c) {
            if (&pi == &pi_p) { // Same Panels case
                return ComputeIntegralCoinciding(pi, pi_p, space, CGaussQR, k, c);
            }
            else if ((pi(1) - pi_p(-1)).norm() / 100. < epsilon ||
                     (pi(-1) - pi_p(1)).norm() / 100. < epsilon) {// Adjacent Panels case
                return ComputeIntegralAdjacent(pi, pi_p, space, CGaussQR, k, c);
            }
            else {// Disjoint panels case
                return ComputeIntegralGeneral(pi, pi_p, space, GaussQR, k, c);
            }
        }

        Eigen::MatrixXcd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                                 const AbstractParametrizedCurve &pi_p,
                                                 const AbstractBEMSpace &space,
                                                 const QuadRule &GaussQR,
                                                 const complex_t k,
                                                 const double c) {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = space.getQ();
            // Interaction matrix with size Qtest x Qtrial
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Computing the (i,j)th matrix entry
            bool swap = (pi(1) - pi_p(-1)).norm() / 100. > epsilon;
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    // Lambda expression for functions F and G
                    auto F = [&](double t) {
                        if (swap) {
                            return space.evaluateShapeFunction_01_swapped(j, t) * pi_p.Derivative_01_swapped(t).norm();
                        } else {
                            return space.evaluateShapeFunction(j, t) * pi_p.Derivative_01(t).norm();
                        }
                    };
                    auto G = [&](double s) {
                        if (swap) {
                            return space.evaluateShapeFunction(i, s) * pi.Derivative_01(s).norm();
                        } else {
                            return space.evaluateShapeFunction_01_swapped(i, s) * pi.Derivative_01_swapped(s).norm();
                        }
                    };
                    // Lambda expression for the integrand
                    auto integrand = [&](double s, double t) {
                        complex_t result = complex_t(0.,0.);
                        if (swap) {
                            if ( abs(k*sqrt(c))*(pi[s]-pi_p.swapped_op(t)).norm() > epsilon) {
                                result = ii*complex_bessel::H1(0,k*sqrt(c)*(pi[s]-pi_p.swapped_op(t)).norm())/4.;
                            } else if ((pi[s]-pi_p.swapped_op(t)).norm() > epsilon){
                                result = -1/(2*M_PI)*log((pi[s]-pi_p.swapped_op(t)).norm());
                            }
                        } else {
                            if ( abs(k*sqrt(c))*(pi.swapped_op(s)-pi_p[t]).norm() > epsilon) {
                                result = ii*complex_bessel::H1(0,k*sqrt(c)*(pi.swapped_op(s)-pi_p[t]).norm())/4.;
                            } else if ((pi.swapped_op(s)-pi_p[t]).norm() > epsilon){
                                result = -1/(2*M_PI)*log((pi.swapped_op(s)-pi_p[t]).norm());
                            }
                        }
                        return result * F(t) * G(s);
                    };
                    complex_t integral = complex_t(0.,0.);
                    // Tensor product quadrature for double integral
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            double s = GaussQR.x(k)*GaussQR.x(l);
                            double t = GaussQR.x(k);
                            double w = GaussQR.x(k)*GaussQR.w(k)*GaussQR.w(l);
                            integral += w*integrand(s,t);
                            integral += w*integrand(t,s);
                        }
                    }
                    // Filling the matrix entry
                    interaction_matrix(i, j) = integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                                   const AbstractParametrizedCurve &pi_p,
                                                   const AbstractBEMSpace &space,
                                                   const QuadRule &GaussQR,
                                                   const complex_t k,
                                                   const double c) {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // Calculating the quadrature order for stable evaluation of integrands for
            // disjoint panels
            // No. of Reference Shape Functions in trial/test space
            int Q = space.getQ();
            // Interaction matrix with size Q x Q
            Eigen::MatrixXcd interaction_matrix(Q, Q);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Q; ++i) {
                for (int j = 0; j < Q; ++j) {
                    // Lambda expression for functions F and G
                    // Single Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunction(j, t) * pi_p.Derivative_01(t).norm();
                    };
                    auto G = [&](double s) { // Function associated with panel pi
                        return space.evaluateShapeFunction(i, s) * pi.Derivative_01(s).norm();
                    };
                    auto integrand = [&](double s, double t) {
                        complex_t result = complex_t(0.0,0.0);
                        //std::cout << (pi[s]-pi_p[t]).norm() << std::endl;
                        if ( abs(k*sqrt(c))*(pi[s]-pi_p[t]).norm() > epsilon) {
                            result = ii*complex_bessel::H1(0,k*sqrt(c)* (pi[s] - pi_p[t]).norm())/4.;
                        } else if ((pi[s]-pi_p[t]).norm() > epsilon){
                            result = -1/(2*M_PI)*log((pi[s]-pi_p[t]).norm());
                        }
                        return result*F(t)*G(s);
                    };
                    // Integrating by tensor product quadrature rule
                    complex_t integral = complex_t(0.,0.);
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
			    // Generate points and weights from existing QR
                            double s = GaussQR.x(l)*(1.-GaussQR.x(k));
                            double t = GaussQR.x(l);
                            double w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
			    // Evaluate integrand symetrically
                            integral += w*integrand(s,t);
                            integral += w*integrand(t,s);
                        }
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                                const AbstractParametrizedCurve &pi_p,
                                                const AbstractBEMSpace &space,
                                                const QuadRule &GaussQR,
                                                const complex_t k,
                                                const double c) {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // Calculating the quadrature order for stable evaluation of integrands for
            // disjoint panels
            // No. of Reference Shape Functions in trial/test space
            int Q = space.getQ();
            // Interaction matrix with size Q x Q
            Eigen::MatrixXcd interaction_matrix(Q, Q);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Q; ++i) {
                for (int j = 0; j < Q; ++j) {
                    // Lambda expression for functions F and G
                    // Single Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunction(j, t) * pi_p.Derivative_01(t).norm();
                    };
                    auto G = [&](double s) { // Function associated with panel pi
                        return space.evaluateShapeFunction(i, s) * pi.Derivative_01(s).norm();
                    };
                    auto integrand = [&](double s, double t) {
                        complex_t result = complex_t(0.0,0.0);
                        if ( abs(k*sqrt(c))*(pi[s]-pi_p[t]).norm() > epsilon) {
                            result = ii*complex_bessel::H1(0,k*sqrt(c)*(pi[s]-pi_p[t]).norm())/4.;
                        } else if ((pi[s]-pi_p[t]).norm() > epsilon){
                            result = -1/(2*M_PI)*log((pi[s]-pi_p[t]).norm());
                        }
                        return result*F(t)*G(s);
                    };
                    complex_t integral = complex_t(0.,0.);
                    // Tensor product quadrature rule
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            double s = GaussQR.x(k);
                            double t = GaussQR.x(l);
                            double w = GaussQR.w(k)*GaussQR.w(l);
                            integral += w*integrand(s,t);
                        }
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &space,
                                        const unsigned int &N,
                                        const complex_t k,
                                        const double c) {
            // Getting the number of panels in the mesh
            unsigned int numpanels = mesh.getNumPanels();
            // Getting dimensions of trial/test space
            unsigned int dims = space.getSpaceDim(numpanels);
            // Getting the panels from the mesh
            PanelVector panels = mesh.getPanels();
            // Getting the number of local shape functions in the trial/test space
            unsigned int Q = space.getQ();
            // Initializing the Galerkin matrix with zeros
            Eigen::MatrixXcd output = Eigen::MatrixXcd::Zero(dims, dims);
            // Panel oriented assembly
            QuadRule GaussQR = getGaussQR(N,0.,1.);
            QuadRule CGaussQR = getCGaussQR(N);
            for (unsigned int i = 0; i < numpanels; ++i) {
                for (unsigned int j = 0; j < numpanels; ++j) {
                    // Getting the interaction matrix for the pair of panels i and j
                    Eigen::MatrixXcd interaction_matrix =
                            InteractionMatrix(*panels[i], *panels[j], space, GaussQR, CGaussQR, k, c);
                    // Local to global mapping of the elements in interaction matrix
                    for (unsigned int I = 0; I < Q; ++I) {
                        for (unsigned int J = 0; J < Q; ++J) {
                            int II = space.LocGlobMap(I + 1, i + 1, numpanels) - 1;
                            int JJ = space.LocGlobMap(J + 1, j + 1, numpanels) - 1;
                            // Filling the Galerkin matrix entries
                            output(II, JJ) += interaction_matrix(I, J);
                        }
                    }
                }
            }
            return output;
        }

    } //namespace single_layer_helmholtz



