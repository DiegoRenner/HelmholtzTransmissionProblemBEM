/**
 * \file sv_circle_arnoldi.cpp
 * \brief This target builds a script that computes the singular values
 * of the Galerkin BEM approximated BIO for the
 * second-kind direct BIEs of the Helmholtz
 * transmission problem. The Arnoldi algorithm from arpack is used to compute the
 * singular values. The scatterer is set to be a circle.
 * The results are written to file.
 * The script can be run as follows:
 *
 * <tt>
 * /path/to/sv_circle \<radius of circle\> \<refraction inside\>
 *      \<refraction outside\> \<initial wavenumber\> \<final wavenumber\>
 *      \<\#points to evaluate\> \<scan complex wavenumbers\> \<\#panels\>
 *      \<order of quadrature rule\> \<accuracy of arnoldi algorithm\>.
 * </tt>
 *
 * The resulting file will contain the value of \f$k\f$ in the first column.
 * The rest of the columns contain the singular values from 
 * smallest to largest for this \f$k\f$.
 * The user will be updated through the command line about the
 * progress of the algorithm
 * if \f$ \verb|-DCMDL| \f$ is set.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "singular_values_arnoldi.hpp"
#include "singular_values.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "continuous_space.hpp"

// define shorthand for complex data type and imaginary unit
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

int main(int argc, char** argv) {

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[2]);
    double c_o = atof(argv[3]);
    complex_t k_0 = atof(argv[4]);
    complex_t k_N = atof(argv[5]);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = atof(argv[6]);
    double h_x = (k_N-k_0).real()/double(n_points_x-1);
    bool scan_complex_wavenumbers = atoi(argv[7]);
    unsigned n_points_y;
    double h_y;
    if (scan_complex_wavenumbers) {
        n_points_y = n_points_x;
        h_y = (k_N-k_0).imag()/double(n_points_y-1);
    } else {
        n_points_y = 1;
        h_y = 0;
    }
    unsigned numpanels = atoi(argv[8]);
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries
    unsigned order = atoi(argv[9]);
    // define accuracy of arnoldi algorithm
    double acc = atof(argv[10]);

    // generate output filename with set parameters
    std::string base_order = "../data/file_SVs_circle_";
    std::string suffix = ".dat";
    std::string divider = "_";
    std::string file_SVs = base_order.append(argv[8])
                             + divider.append(argv[10]) + suffix;
    // clear existing file
    std::ofstream file_out;
    file_out.open(file_SVs, std::ofstream::out | std::ofstream::trunc);
    file_out.close();

    // Inform user of started computation.
	#ifdef CMDL
    std::cout << "Parameters are:" << std::endl;
    std::cout << "radius of circle [arg]: "  << eps << std::endl;
    std::cout << "refraction inside [arg]: "  << c_i << std::endl;
    std::cout << "refraction outside [arg]: "  << c_o << std::endl;
    std::cout << "initial wavenumber [arg]: "  << k_0 << std::endl;
    std::cout << "final wavenumber [arg]: "  << k_N << std::endl;
    std::cout << "#points to evaluate [arg]: "  << n_points_x << std::endl;
    std::cout << "scan complex wavenumbers [arg]: "  << scan_complex_wavenumbers << std::endl;
    std::cout << "real step size [computed]: "  << h_x << std::endl;
    std::cout << "imaginary step size [computed]: "  << h_y << std::endl;
    std::cout << "#panels [arg]: "  << numpanels << std::endl;
    std::cout << "order of quadrature rule [arg]: "  << order << std::endl;
    std::cout << "accuracy of arnoldi algorithm [arg]: "  << acc << std::endl;
    std::cout << "file name [computed]: "  << file_SVs << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Computing singular values of BIO." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
	#endif

    ContinuousSpace<1> cont_space;
    SolutionsOperator so(mesh, order, cont_space, cont_space);

    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            Eigen::MatrixXd res(2*numpanels,3);
            // define wavenumber for current loop
            complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

            // compute solutions operator
            Eigen::MatrixXcd T;
            so.gen_sol_op(k_temp, c_o, c_i, T);


            // set singular values to be computed, all
            // from arpack docs:
            // NEV     Integer.  (INPUT)
            // Number of eigenvalues of OP to be computed. 0 < NEV < N.
            // https://www.caam.rice.edu/software/ARPACK/UG/node136.html
            unsigned count = 2*numpanels-1;
            double list[count];
            for (unsigned i = 0; i < count; i++){
                list[i] = i;
            }

            // compute singular value
            res = arnoldi::sv(T,count,acc);

            // write singular values to file
            file_out.open(file_SVs, std::ios_base::app);
            file_out << k_temp.real() << " ";
            file_out << res.block(0, 0, count, 1).transpose() << std::endl;
            file_out.close();
			
			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			std::cout << "Singular values at " << k_temp << " computed." << std::endl;
			std::cout << "Smallest singular value is: "
				<< res.block(0, 0, 1, 1).transpose() << std::endl;
            std::cout << "#######################################################" << std::endl;
			std::cout << std::endl;
			#endif
        }
    }
    return 0;
}
