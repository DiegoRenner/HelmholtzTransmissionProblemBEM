/**
 * \file sv_circle.cpp
 * \brief This target builds a script that computes the singular values
 * of the Galerkin BEM approximated BIO for the
 * second-kind direct BIEs of the Helmholtz
 * transmission problem. The direct algorithm from Eigen is used to compute the
 * singular values.
 * The scatterer is set to be a circle.
 * The results are written to file.
 * The script can be run as follows:
 *
 * <tt>
 * /path/to/sv_circle \<radius of circle\> \<refraction inside\>
 *      \<refraction outside\> \<initial wavenumber\> \<final wavenumber\>
 *      \<\#panels\> \<order of quadrature rule\> \<outputfile\>.
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
#include "singular_values.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"

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
    unsigned n_points_x = atoi(argv[6]);
    unsigned n_points_y = atoi(argv[6]);
    unsigned numpanels = atoi(argv[7]);
    double h_x = (k_N-k_0).real()/double(n_points_x-1);
    double h_y = (k_N-k_0).imag()/double(n_points_x-1);
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries
    unsigned order = atoi(argv[8]);

    // clear existing file
    std::ofstream file_out;
    file_out.open(argv[9], std::ofstream::out | std::ofstream::trunc);
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
    //std::cout << "scan complex wavenumbers [arg]: "  << scan_complex_wavenumbers << std::endl;
    std::cout << "real step size [computed]: "  << h_x << std::endl;
    std::cout << "imaginary step size [computed]: "  << h_y << std::endl;
    std::cout << "#panels [arg]: "  << numpanels << std::endl;
    std::cout << "order of quadrature rule [arg]: "  << order << std::endl;
    //std::cout << "accuracy of arnoldi algorithm [arg]: "  << acc << std::endl;
    std::cout << "file name [arg]: "  << argv << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Computing singular values of BIO." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
    #endif

    SolutionsOperator so(mesh, order);

    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            Eigen::MatrixXd res(2*numpanels,3);
            // define wavenumber for current loop
            complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

            // compute solutions operator
            Eigen::MatrixXcd T = so.gen_sol_op(k_temp, c_o, c_i);

            // set singular values to be computed, all
            unsigned count = T.cols();
            double list[count];
            for (unsigned i = 0; i < count; i++){
                list[i] = i;
            }

            // compute singular value
            res = direct::sv(T,list,count);

            // write singular values to file
            file_out.open(argv[7], std::ios_base::app);
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
