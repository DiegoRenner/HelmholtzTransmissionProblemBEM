/**
 * \file transmission_problem_verification.cpp
 * \brief This target builds a script that computes solutions to
 * the analytically solvable case of the Helmholtz transmission
 * problem where the scatterer is a circle using second-kind direct
 * BIEs and Galerkin BEM.
 * The results are written to file.
 * The script can be run as follows:
 * <tt>
 * /path/to/transmission_problem_verification \<radius of circle\>
 *      \<\#coeffs for series expansion of solution\> \<refraction inside\>
 *      \<refraction outside\> \<wavenumber\>
 *      \<order of quadrature rule\> \<outputfile\>
 * </tt>
 * This output file will contain two columns.
 * The first will contain the current panel size.
 * The second will contain the residual error in the euclidean 
 * norm of the computed FEM-space interpolation coefficients to the 
 * known FEM-space interpolation coefficients for the current number of panels.
 * Then the columns will contain the computed derivative, the 
 * extrapolated derivative, the computed second derivative and the extrapolated second 
 * derivative in this order.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "solvers.hpp"
#include "gen_sol.hpp"
#include "continuous_space.hpp"
#include "mass_matrix.hpp"
#include <chrono>

using namespace std::chrono;

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main(int argc, char** argv) {

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[3]);
    double c_o = atof(argv[4]);
    double k = atof(argv[5]);

    //set boundary
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);

    // define order of quadrature rule used to compute matrix entries
    unsigned order = atoi(argv[6]);


    // set number of panels for which to compute BIOs
    unsigned n_runs = 7;
    int numpanels[n_runs];
    numpanels[0] = 50;
    for (int i=1; i<n_runs; i++){
        numpanels[i] = 2*numpanels[i-1];
    }
    // define incoming wave and resulting wave
    int l = atoi(argv[2]);
    double a_n[2*l+1];
    for( int i = 0; i<2*l+1; i++) {
        a_n[i] = 1.;///((k*k*(c_o-c_i))*sqrt((2*l+1)*M_PI*eps*eps*(jn(i-l,k)*jn(i-l,k)-jn(i-l-1,k)*jn(i-l+1,k))));
    }
    auto u_i_dir = [&] (double x1, double x2) {
        return sol::u_i(x1, x2, l, a_n, k);
    };
    auto u_t_dir = [&] (double x1, double x2) {
        return sol::u_t(x1, x2, l, eps, a_n, k, c_i);
    };
    auto u_i_neu = [&] (double x1, double x2) {
        return sol::u_i_neu(x1, x2, l, a_n, k);
    };
    auto u_t_neu = [&] (double x1, double x2) {
        return sol::u_t_neu(x1, x2, l, eps, a_n, k, c_i);
    };

    // set FEM-sapces of lowest order for validation
    ContinuousSpace<1> cont_space;

    // generate outputfilename
    std::string base_order = "../data/file_order_";
    std::string suffix = ".dat";
    std::string divider = "_";
    std::string file_order = base_order.append(argv[6])
                             + divider.append(argv[2]) + suffix;
    // clear existing file
    std::ofstream file_out;
    file_out.open(file_order, std::ofstream::out | std::ofstream::trunc);
    file_out.close();

#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Computing convergence rate for Helmholtz TP BIEs." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
#endif

    // Loop over number of panels
    for (unsigned i = 0; i < n_runs; i++) {

        // generate mesh on boudary
        ParametrizedMesh mesh(curve.split(numpanels[i]));

        // compute interpolation coefficients
        // in FEM-sapces for resulting waves
        // DEBUG auto start = high_resolution_clock::now();
        Eigen::VectorXcd sol = tp::direct_second_kind::solve(
                mesh, u_i_dir, u_i_neu, order, k, c_o, c_i);
        // DEBUG auto end = high_resolution_clock::now();
        // DEBUG auto duration = duration_cast<milliseconds>(end - start);

        // compute mass matrix for projection on to orthonromal basis functions
        //Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
        //Eigen::MatrixXcd M(2*numpanels[i],2*numpanels[i]);
        //M = Eigen::MatrixXcd::Zero(2*numpanels[i],2*numpanels[i]);
        //M.block(0,0,numpanels[i],numpanels[i]) = M_cont;
        //M.block(numpanels[i],numpanels[i],numpanels[i],numpanels[i]) = M_cont;

        // compute interpolation coefficients in FEM-spaces of known solution
        Eigen::VectorXcd u_t_dir_N = cont_space.Interpolate_helmholtz(u_t_dir,mesh);
        Eigen::VectorXcd u_t_neu_N = cont_space.Interpolate_helmholtz(u_t_neu,mesh);
        Eigen::VectorXcd u_t_N(2*numpanels[i]);
        u_t_N << u_t_dir_N, u_t_neu_N;
        //std::cout << sol.transpose() << std::endl;
        //std::cout << u_t_N.transpose() << std::endl;
        PanelVector panels_coarse = mesh.getPanels();
        unsigned N = 20;
        QuadRule GaussQR = getGaussQR(N,0.,1.);
        //for (int j = 0; j < numpanels[i]; j++){
        //    std::cout << u_t_dir(panels_coarse[(j)%numpanels[i]]->operator[](0.0).x(), panels_coarse[(j) % numpanels[i]]->operator[](0.0).y());
        //}
        std::cout << std::endl;

        //generate finer mesh for error approximation
        //use panels_fine of computation mesh for evaluating shape functions
        //generate finer mesh for error approximation
        //ParametrizedMesh mesh_refined(curve.split(4*numpanels[n_runs-1]));
        //PanelVector panels_fine = mesh_refined.getPanels();
        //generate QR of higher order for error approximation
        //double res1 = 0.0;
        //for (unsigned j=0; j < 4*numpanels[n_runs-1]; j++){
        //    //compute index for shapefunction on coarser mesh
        //    //dividing j by the number of points the finer mesh has within one panel of the coarser mesh
        //    int j_coarse = j / 4 / pow(2, (n_runs - 1 - i));
        //    for (unsigned k=0; k < N; k++){
        //        //rescaling quadrature rule to finer mesh
        //        double x =
        //                (GaussQR.x(k)
        //                 +(j%(4*int(pow(2,(n_runs-1-i))))))
        //                *1.0/double(4*pow(2,n_runs-1-i));
        //        //contribution of first shape fct.
        //        complex_t temp = (sol[j_coarse] * cont_space.evaluateShapeFunction(1, x)
        //                          //*panels_fine[i]->Derivative_01(k*0.25).norm()
        //                          +
        //                          //contribution of second shape fct.
        //                          sol[(j_coarse + 1) % numpanels[i]] * cont_space.evaluateShapeFunction(0, x)
        //                          //*panels_fine[i]->Derivative_01(k*0.25).norm()
        //                          -
        //                          //exact solution
        //                          u_t_dir(panels_coarse[j_coarse]->operator[](x).x(), panels_coarse[j_coarse]->operator[](x).y()));
        //        //squared norm multiplied by scaling factors
        //        res1 += (temp*(temp.real()-ii*temp.imag())).real()*panels_fine[j]->length() * GaussQR.w(k);
        //    }
        //}
        //std::cout << "Test1: " << sqrt(res1) << std::endl;

        complex_t res2 = 0.0;
        for (int j=0; j < numpanels[i]; j++){
            //compute index for shapefunction on coarser mesh
            //dividing j by the number of points the finer mesh has within one panel of the coarser mesh
            for (int m=0; m < N; m++){
                //rescaling quadrature rule to finer mesh
                //contribution of first shape fct.
                //std::cout << (sol[j] * cont_space.evaluateShapeFunction(1, GaussQR.x(m))) << std::endl;
                //std::cout << (sol[(j+1)%numpanels[i]] * cont_space.evaluateShapeFunction(1, GaussQR.x(m))) << std::endl;
                //std::cout << u_t_dir(panels_coarse[j]->operator[](GaussQR.x(m)).x(), panels_coarse[j]->operator[](GaussQR.x(m)).y())<< std::endl;
                complex_t temp = (sol[j] * cont_space.evaluateShapeFunction(1, GaussQR.x(m))
                                   +
                                   //contribution of second shape fct.
                                   sol[(j + 1) % numpanels[i]] * cont_space.evaluateShapeFunction(0, GaussQR.x(m))
                                   -
                                   //exact solution
                                   u_t_dir(panels_coarse[(j)%numpanels[i]]->operator[](GaussQR.x(m)).x(), panels_coarse[(j) % numpanels[i]]->operator[](GaussQR.x(m)).y()));
                //squared norm multiplied by scaling factors
                res2 += (temp*(temp.real()-ii*temp.imag())).real() * GaussQR.w(m) * panels_coarse[j]->Derivative_01(GaussQR.x(m)).norm();
            }
        }
        std::cout << "Test2: " << sqrt(res2) << std::endl;
        complex_t res3 = 0.0;
        for (int j=0; j < numpanels[i]; j++){
            //compute index for shapefunction on coarser mesh
            //dividing j by the number of points the finer mesh has within one panel of the coarser mesh
            for (int m=0; m < N; m++){
                //rescaling quadrature rule to finer mesh
                //contribution of first shape fct.
                //std::cout << sqrt(pow(panels_coarse[j]->operator[](GaussQR.x(m)).x(),2)
                //                + pow(panels_coarse[j]->operator[](GaussQR.x(m)).y(),2))
                //<< std::endl;
                complex_t temp = (u_t_N[j] * cont_space.evaluateShapeFunction(1, GaussQR.x(m))
                                 +
                                 //contribution of second shape fct.
                                 u_t_N[(j + 1) % numpanels[i]] * cont_space.evaluateShapeFunction(0, GaussQR.x(m))
                                 -
                                 //exact solution
                                 u_t_dir(panels_coarse[j]->operator[](GaussQR.x(m)).x(), panels_coarse[j]->operator[](GaussQR.x(m)).y()));
                //squared norm multiplied by scaling factors
                res3 += (temp*(temp.real()-ii*temp.imag())) * GaussQR.w(m) * panels_coarse[j]->length();
            }
        }
        std::cout << "Test3: " << sqrt(res3) << std::endl;

        //double res4 = 0;
        //for (unsigned j=0; j < 4*numpanels[n_runs-1]; j++){
        //    int j_specific = j/4/pow(2,(n_runs-1-i));
        //    for (unsigned k=0; k < N; k++){
        //        double x =
        //                (GaussQR.x(k)
        //                 +(j%(4*int(pow(2,(n_runs-1-i))))))
        //                *1.0/double(4*pow(2,n_runs-1-i));
        //        complex_t temp = (u_t_N[j_specific]*cont_space.evaluateShapeFunction(1,x)
        //                           //*panels_fine[i]->Derivative_01(k*0.25).norm()
        //                           +
        //                           u_t_N[(j_specific+1)%numpanels[i]]*cont_space.evaluateShapeFunction(0,x)
        //                           //*panels_fine[i]->Derivative_01(k*0.25).norm()
        //                           -
        //                           u_t_dir(panels_fine[j]->operator[](GaussQR.x(k)).x(), panels_fine[j]->operator[](GaussQR.x(k)).y()));
        //        res4 += (temp*(temp.real()-ii*temp.imag())).real() * panels_fine[j]->length() * GaussQR.w(k);
        //    }
        //}
        //std::cout << "Test4: " << sqrt(res4) << std::endl;
        // write difference to computed solution in L^2 norm to file
        file_out.open(file_order, std::ios_base::app);
        //file_out << mesh.getPanels()[0]->length() << " " << sqrt(abs((sol - u_t_N).dot(M * (sol - u_t_N)))) << std::endl;
        file_out << mesh.getPanels()[0]->length() << " " << sqrt(res2) <<  " " << sqrt(res3) << std::endl;
        file_out.close();
#ifdef CMDL
        std::cout << "#######################################################" << std::endl;
        std::cout << "Computed Cauchy data on " << numpanels[i] << " panels_fine." << std::endl;
        std::cout << "Residual error of FEM-space interpolation coefficients:" << std::endl;
		std::cout << sqrt(abs((sol-u_t_N).dot(M*(sol-u_t_N)))) << std::endl;
        std::cout << "#######################################################" << std::endl;
        std::cout << std::endl;
#endif
    }
    return 0;
}
