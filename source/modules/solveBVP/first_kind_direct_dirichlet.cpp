#include <complex>
#include "first_kind_direct_dirichlet.hpp"
#include "parametrized_mesh.hpp"
#include "solvers.hpp"
#include "gen_sol.hpp"
#include "mass_matrix.hpp"
#include "discontinuous_space.hpp"

FirstKindDirectDirichlet::FirstKindDirectDirichlet(Json::Value job_config): ipt(Eigen::Vector2d(0.125,0.0)){
    if (job_config["domain"] != "circle"){
        if (job_config["radius"] != "None"){
            eps = job_config["radius"].asDouble();
        }
    boundary = ParametrizedCircularArc(Eigen::Vector2d(0.0,0.0), eps, 0, 2*M_PI);
    } else {
        boundary = ParametrizedCircularArc(Eigen::Vector2d(0.0,0.0), eps, 0, 2*M_PI);
    }
    if (job_config["wavenumber"] != "None"){
        k = job_config["wavenumber"].asDouble();
    }
    if (job_config["n_runs"] != "None"){
        n_runs = job_config["n_runs"].asUInt();
    }
    if (job_config["init_panels"] != "None"){
        init_panels = job_config["init_panels"].asUInt();
    }
    if (job_config["eval_point"] != "None"){
        ipt(0) = job_config["eval_point"][0].asDouble();
        ipt(1) = job_config["eval_point"][1].asDouble();
    }

    numpanels.resize(n_runs);
    numpanels[0] = init_panels;
    for (int i = 1; i < n_runs; i++) {
        numpanels[i] = 2 * numpanels[i - 1];
    }
}

void FirstKindDirectDirichlet::run() {
    typedef std::complex<double> complex_t;
    complex_t ii = complex_t(0,1.);
    auto fund_sol_dir = [&] (double x1, double x2){
        return sol::fund_sol_dir(k,x1,x2,ipt[0],ipt[1]);
    };
    auto fund_sol_neu = [&] (double x1, double x2){
        return sol::fund_sol_neu(k,x1,x2,ipt[0],ipt[1]);
    };
    // define FEM-spaces for result validation later on
    // Inform user of started computation.
    DiscontinuousSpace<0> discont_space;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Solving Dirichlet problem for increasing grid sizes." << std::endl;
    std::cout << "Using first-kind direct BIEs." << std::endl;
    std::cout << "Using lowest order FEM-spaces." << std::endl;
    std::cout << std::endl;
    // loop over #panels
    for (unsigned i = 0; i <= n_runs; i++) {
        // compute mesh
        ParametrizedMesh mesh(boundary.split(numpanels[i]));
        // compute Neumann data from Dirichlet data using 1st kind BIE
        Eigen::VectorXcd res = bvp::direct_first_kind::solve_dirichlet(
                mesh, fund_sol_dir, order, k);
        // compute interpolation coefficients of known solution Neumann data in FEM-space
        Eigen::VectorXcd res_known = discont_space.Interpolate_helmholtz(fund_sol_neu,mesh);
        // compute mass matrix for projection onto orthonormal basis functions
        Eigen::MatrixXcd M = mass_matrix::GalerkinMatrix(mesh,discont_space,discont_space,order);

        // setup mesh and QR for computing residuals
        PanelVector panels_coarse = mesh.getPanels();
        unsigned N = 20;
        QuadRule GaussQR = getGaussQR(N,0.,1.);

        // compute residual w.r.t. computed and exact solution
        complex_t res_val = 0.0;
        for (int j=0; j < numpanels[i]; j++){
            //compute index for shapefunction on coarser mesh
            //dividing j by the number of points the finer mesh has within one panel of the coarser mesh
            for (int m=0; m < N; m++){
                //rescaling quadrature rule to finer mesh
                //contribution of first shape fct.
                complex_t temp = (res[j] * discont_space.evaluateShapeFunction(0, GaussQR.x(m))
                                  -
                                  //exact solution
                                  fund_sol_neu(panels_coarse[j]->operator[](GaussQR.x(m)).x(), panels_coarse[j]->operator[](GaussQR.x(m)).y()));
                //squared norm multiplied by scaling factors
                res_val += (temp*(temp.real()-ii*temp.imag())) * GaussQR.w(m) * panels_coarse[j]->length();
            }
        }

        // compute residual w.r.t. best solution in FEM space and exact solution
        complex_t res_val1 = 0.0;
        for (int j=0; j < numpanels[i]; j++){
            //compute index for shapefunction on coarser mesh
            //dividing j by the number of points the finer mesh has within one panel of the coarser mesh
            for (int m=0; m < N; m++){
                //rescaling quadrature rule to finer mesh
                //contribution of first shape fct.
                complex_t temp = (res_known[j] * discont_space.evaluateShapeFunction(0, GaussQR.x(m))
                                  -
                                  //exact solution
                                  fund_sol_neu(panels_coarse[j]->operator[](GaussQR.x(m)).x(), panels_coarse[j]->operator[](GaussQR.x(m)).y()));
                //squared norm multiplied by scaling factors
                res_val1 += (temp*(temp.real()-ii*temp.imag())) * GaussQR.w(m) * panels_coarse[j]->length();
            }
        }
        // update user on residual error between computed and known FEM-space
        // interpolation coefficients that have been projected onto orthornomal basis
        // also output total error from found to exact solution and from known best solution
        // in FEM space to exact solution
        std::cout << "#######################################################" << std::endl;
        std::cout << res.transpose() << std::endl;
        std::cout << res_known.transpose() << std::endl;
        std::cout << "Computed Neumann data on " << numpanels[i] << " panels." << std::endl;
        std::cout << "Residual error of FEM-space interpolation coefficients:" << std::endl;
        std::cout << sqrt(((res-res_known)).dot(M*(res-res_known))).real() << std::endl;
        std::cout << "Residual error w.r.t. computed and exact solution" << std::endl;
        std::cout << sqrt(abs(res_val)) << std::endl;
        std::cout << "Residual error w.r.t. best solution in FEM space and exact solution" << std::endl;
        std::cout << sqrt(abs(res_val1)) << std::endl;

        std::cout << "#######################################################" << std::endl;
        std::cout << std::endl;
    }

}
