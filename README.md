# HelmholtzTransmissionProblemBEM
## Configuration and Dependencies
The library can be configured by running 
~~~
cmake CMakeLists.txt
~~~
in the base directory or
~~~
cmake ..
~~~
in a subdirectory for an out-of-source build.
This should automatically generate a target with which the Eigen library for matrix computations can be built by running
~~~
make Eigen
~~~
if it is not already available.
The [complex_bessel library](https://github.com/joeydumont/complex_bessel)
that is used for passing complex arguments to the Hankel and Bessel functions 
aswell as the [arpackpp library](https://github.com/m-reuter/arpackpp) which gives a <tt>c++</tt> interface to the [arpack library](https://github.com/opencollab/arpack-ng) are installed automatically and don't need to be built.
<tt>arpack</tt> and <tt>lapack</tt> need to be installed separately and can usually be done so with your distributions packagemanager.

For <tt>arch</tt> based distros:
~~~
sudo pacman -S arpack
sudo pacman -S lapack
~~~
For <tt>debian</tt> based distros:
~~~
sudo apt install libboost-all-dev
sudo apt install libarpack2-dev 
sudo apt install liblapack3-dev
~~~

To generate the documentation <tt>latex</tt> and <tt>doxygen</tt> have to be installed as well.
For <tt>arch</tt> based distros:
~~~
sudo pacman -S texlive-most
sudo pacman -S doxygen
~~~
For <tt>debian</tt> based distros:
~~~
sudo apt install texlive-full
sudo apt install doxygen
~~~
Running CMake also configures some examples of how to use the library as <tt>make</tt> targets.
These can then be built by running 

~~~
make <target_name>
~~~

The compiled binary can be found in the <tt>bin</tt> directory.

## Usage
We will show how the built targets are to be used.
We commonly refer to the wavenumber by k.

#### <tt>doxygen_HelmholtzTransmissionProblemBEM</tt>
This target generates a documentation of the library in the <tt>doxygen/generated_doc</tt> directory.
The documentation can be browsed using any common browser.

#### <tt>debugging_SVs</tt>
This target builds a script that computes the Eigenvalues of the BIO for Helmholtz Transmission Problem. The results are written to file. The script can be run as follows:
~~~
 /path/to/debugging_SVs
~~~
The output file will contain a section for each set mesh resolution and each of those sections will contain one section each for every BIO where all Eigenvalues for different wavenumbers will be listed in columns. The Eigenvalues are computed using the facts stated in Lemma 3.22. [TODO: find reference]

#### <tt>direct_v_arnoldi</tt>
This target builds a script that computes the singular values of the Galerkin BEM approximated BIO for the second-kind direct BIEs of the Helmholtz transmission problem, once using the Arnoldi algorithm and once using s direct solver. The scatterer is set to be a circle. The results are written to file. The script can be run as follows: 
~~~
/path/to/direct_v_arnoldi <radius of circle> <number of SVs to be computed> <accurracy of arnoldi algorithm>. 
~~~
 The script will generate four files: file_vals_eig_<number of SVs to be computed>_<accurracy of arnoldi algorithm>.dat, file_vals_arpp_<number of SVs to be computed>_<accurracy of arnoldi algorithm>.dat, file_timings_<number of SVs to be computed>_<accurracy of arnoldi algorithm>.dat, file_iter_<number of SVs to be computed>_<accurracy of arnoldi algorithm>.dat. These will contain the SVs computed using the direct solver, the SVs computed using the Arnoldi algorithm, the time taken by the direct solver and the Arnoldi algorithm, and the number of iterations the Arnoldi algorithm took to converge respectively.
#### <tt>direct_v_arnoldi_1st_der</tt>
This target builds a script that computes the first derivative of the singular values of the Galerkin BEM approximated BIO for the second-kind direct BIEs of the Helmholtz transmission problem, once using the Arnoldi algorithm and once using s direct solver. The scatterer is set to be a circle. The results are written to file. The script can be run as follows: 
~~~
 /path/to/direct_v_arnoldi_1st_der <radius of circle> <number of SV derivatives to be computed>
 <accurracy of arnoldi algorithm>.
~~~
 The script will generate four files:
 file_vals_eig_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_1stDer.dat,
 file_vals_arpp_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_1stDer.dat,
 file_timings_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_1stDer.dat,
 file_iter_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_1stDer.dat.
 These will contain the derivatives computed using the direct solver,
 the derivatives computed using the Arnoldi algorithm,
 the time taken by the direct solver and the Arnoldi algorithm,
 and the number of iterations the Arnoldi algorithm took to converge respectively.

#### <tt>direct_v_arnoldi_2nd_der</tt>
This target builds a script that computes the second derivative of the singular values of the Galerkin BEM approximated BIO for the second-kind direct BIEs of the Helmholtz transmission problem, once using the Arnoldi algorithm and once using s direct solver. The scatterer is set to be a circle. The results are written to file. The script can be run as follows: <tt> /path/to/direct_v_arnoldi_2nd_der <radius of circle> <number of SV derivatives to be computed> <accurracy of arnoldi algorithm>. </tt> The script will generate four files: file_vals_eig_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_2ndDer.dat, file_vals_arpp_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_2ndDer.dat, file_timings_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_2ndDer.dat, file_iter_<number of SV derivatives to be computed>_<accurracy of arnoldi algorithm>_2ndDer.dat. These will contain the derivatives computed using the direct solver, the derivatives computed using the Arnoldi algorithm, the time taken by the direct solver and the Arnoldi algorithm, and the number of iterations the Arnoldi algorithm took to converge respectively.

#### <tt>dirichlet_example</tt>
This target builds a script that computes the solution to a Dirichlet problem
using first kind direct BIEs.
No command line parameters are necessary.
Once built the script can be run as follows: 
~~~
/path/to/library/bin/dirichlet_example.
~~~
The user will be updated over the residual error in the euclidean norm of the computed FEM-space interpolation coefficients to the known FEM-space interpolation coefficients for the current number of panels through the command line.
 
#### <tt>neumann_example</tt>
This target builds a script that computes the solution to a Neumann problem
using first kind direct BIEs.
No command line parameters are necessary.
Once built the script can be run as follows: 
~~~
/path/to/library/bin/neumann_example.
~~~
The user will be updated over the residual error in the euclidean norm of the computed FEM-space interpolation coefficients to the known FEM-space interpolation coefficients for the current number of panels through the command line.

#### <tt>parabolic_approximation</tt>
This target builds a script that tries to find minimas in the smallest sinuglar values
of the Galerkin BEM approximated solutions operator for the second-kind direct BIEs of 
the Helmholtz Transmission problem.
The minimas are searched for using a parabolic approximation
based on evaluating the smallest singular values and their first
two derivatives.
The results are written to disk.
No command line arguments are necessary.
The script can be run as follows:
~~~
/path/to/library/bin/parabolic_approximation <outfile>.
~~~
In the file the first column contains the initial point used for the parabolic approximation.
The next three columns contain the function value and the first two derivatives at the initial point that were used to compute the parabolic approximation.
The user also will get updates on the current best approximation for a minima and the value of the first derivatie at this point through the command line if <tt>-DCMDL</tt> is set.

#### <tt>roots_brent_circle</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the Van Wijngaarden-Dekker-Brent method.
The scatterer is set to be a circle.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <radius of circle> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the direct
Eigen algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>roots_brent_circle_arnoldi</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the Van Wijngaarden-Dekker-Brent method.
The scatterer is set to be a circle.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <radius of circle> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>roots_brent_square</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the sedond-kind direct BIEs of the Helmholtz
transmission problem using the Van Wijngaarden-Dekker-Brent method.
The scatterer is set to be a square.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <half side length of square> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the direct
Eigen algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>roots_brent_square_arnoldi</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the sedond-kind direct BIEs of the Helmholtz
transmission problem using the Van Wijngaarden-Dekker-Brent method.
The scatterer is set to be a square.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <half side length of square> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>roots_mixed_circle_arnoldi</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using a precursor to the Algorithm described in Listing 1 of
 https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2022/2022-38.pdf.
The scatterer is set to be a circle.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <radius of circle> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>roots_mixed_square_arnoldi</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the sedond-kind direct BIEs of the Helmholtz
transmission problem using a precursor to the Algorithm described in Listing 1 of
https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2022/2022-38.pdf.
The scatterer is set to be a square.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <half side length of square> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>roots_newton_circle</tt>
This target builds a sript that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the Newton-Raphson method.
The scatterer is set to be a circle.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_newton_circle <radius of circle> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point,
the function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the direct
Eigen algorithm.
The user will be updated through the command line about the
progress of the algorithm.
if <tt>-DCMDL</tt> is set.

#### <tt>roots_newton_circle_arnoldi</tt>
This target builds a sript that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the Newton-Raphson method.
The scatterer is set to be a circle.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_newton_circle <radius of circle> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point,
the function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm.
if <tt>-DCMDL</tt> is set.

#### <tt>roots_newton_square</tt>
This target builds a sript that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the Newton-Raphson method.
The scatterer is set to be a square.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_newton_circle <side length of square> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point,
the function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the direct
Eigen algorithm.
The user will be updated through the command line about the
progress of the algorithm.
if <tt>-DCMDL</tt> is set.

#### <tt>roots_newton_square_arnoldi</tt>
This target builds a sript that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the Newton-Raphson method.
The scatterer is set to be a square.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_newton_circle <side length of square> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point,
the function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm.
if <tt>-DCMDL</tt> is set.

#### <tt>roots_seq_circle_arnoldi</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the algorithm described in Listing 1 of
https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2022/2022-38.pdf.
The scatterer is set to be a circle.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <radius of circle> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>roots_seq_square_arnoldi</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the sedond-kind direct BIEs of the Helmholtz
transmission problem using the algorithm described in Listing 1 of
https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2022/2022-38.pdf.
The scatterer is set to be a square.
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <half side length of square> 
    <refraction inside> <refraction outside> <initial wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the
interval used to compute the root in the first column.
Then in the next three columns will be the point, the
function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.
The singular values and their derivatives are computed using the Arnoldi algorithm.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>sv_circle</tt>
This target builds a script that computes the singular values
of the Galerkin BEM approximated BIO for the
second-kind direct BIEs of the Helmholtz
transmission problem. The direct algorithm from Eigen is used to compute the
sinuglar values.
The scatterer is set to be a circle.
The results are written to file.
The script can be run as follows:

~~~
/path/to/sv_circle <radius of circle> <refraction inside>
     <refraction outside> <initial wavenumber> <final wavenumber>
     <#panels> <order of quadrature rule> <outputfile>.
~~~

The resulting file will contain the value of <tt>k</tt> in the first column.
The rest of the columns contain the singular values from
smallest to largest for this <tt>k</tt>.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>sv_circle_arnoldi</tt>
This target builds a script that computes the singular values
of the Galerkin BEM approximated BIO for the
second-kind direct BIEs of the Helmholtz
transmission problem. The arnoldi algorithm from arpack is used to compute the
sinuglar values. The scatterer is set to be a circle.
The results are written to file.
The script can be run as follows:

~~~
/path/to/sv_circle <radius of circle> <refraction inside>
     <refraction outside> <initial wavenumber> <final wavenumber>
     <#points to evaluate> <scan complex wavenumbers> <#panels>
     <order of quadrature rule> <accuracy of arnoldi algorithm>.
~~~

The resulting file will contain the value of <tt>k</tt> in the first column.
The rest of the columns contain the singular values from
smallest to largest for this <tt>k</tt>.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>sv_derivative_full</tt>
This target builds a script that computes the singular values and
their first two derivatives of the Galerkin BEM
approximated BIO for the second-kind direct BIEs of the Helmholtz
transmission problem.
Minimas in the smallest singular value are determined as well
by use of the Newton-Raphson method.
The scatterer is set to be a circle.
The results are written to file.
The script can be run as follows:
~~~
/path/to/library/bin/sv_derivative_full <radius of circle> 
    <refraction inside> <refraction outside> <initial wavenumber>
    <#panels> <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the value of k in the first column.
Then the singular values and their first two derivatives at k will be listed from smallest to largest in the columns.
The singular values and their derivatives occupy three neighboring columns.
The final three columns will contain the value of the root, the value of the first derivative at the root and the number of iterations taken to find the root in the interval between the current and the next evaluation point.
If no root was found these three columns will contain <tt>NAN</tt>.
The user will be updated through the command line about the progress of the algorithm if <tt>-DCMDL</tt> is set.

#### <tt>sv_derivative_verification_circle</tt>
This target builds a script that verifies the derivatives of the singular
values of the Galerkin BEM approximated BIO for the
second-kind direct BIEs of the Helmholtz transmsission problem
using extrapolation.
The scatterer is set to be a circle.
The results are written to file.
The script can be run as follows:
~~~
/path/to/library/bin/sv_derivative_verification_circle 
    <radius of circle> <refraction inside> <refraction outside> 
    <initial wavenumber> <#panels> <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the value of k in the first column.
The second column will contain the value of the smallest singular value at this k.
Then the columns will contain the computed derivative, the extrapolated derivative, the computed second derivative and the extrapolated second derivative in this order.
The user will be updated through the command line about the progress of the algorithm if <tt>-DCMDL</tt> is set.

#### <tt>sv_derivative_verification_square</tt>
This target builds a script that verifies the derivatives of the singular
values and their derivatives of the Galerkin BEM BIO for the
second-kind direct BIEs of the Helmholtz transmsission problem
using extrapolation.
The scatterer is set to be a square.
The results are written to file.
The script can be run as follows:
~~~
/path/to/library/bin/sv_derivative_verification_circle 
    <half side length of square> <refraction inside> 
    <refraction outside> <initial wavenumber> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the value of k in the first column.
The second column will contain the value of the smallest singular value at this k.
Then the columns will contain the computed derivative, the extrapolated derivative, the computed second derivative and the extrapolated second derivative in this order.
The user will be updated through the command line about the progress of the algorithm if <tt>-DCMDL</tt> is set.

#### <tt>sv_square</tt>
This target builds a script that computes the singular values
of the Galerkin BEM approximated BIO for the
second-kind direct BIEs of the Helmholtz
transmission problem. The direct algorithm from Eigen is used to compute the
sinuglar values.
The scatterer is set to be a square.
The results are written to file.
The script can be run as follows:

~~~
/path/to/sv_square <half of side length of square> <refraction inside>
     <refraction outside> <initial wavenumber>
     <#panels> <order of quadrature rule> <outputfile>.
~~~

The resulting file will contain the value of <tt>k</tt> in the first column.
The rest of the columns contain the singular values from
smallest to largest for this <tt>k</tt>.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>sv_square_arnoldi</tt>
This target builds a script that computes the singular values
of the Galerkin BEM approximated BIO for the
second-kind direct BIEs of the Helmholtz
transmission problem. The arnoldi algorithm from arpack is used to compute the
sinuglar values. The scatterer is set to be a square.
The results are written to file.
The script can be run as follows:

~~~
/path/to/sv_circle <radius of circle> <refraction inside>
     <refraction outside> <initial wavenumber> <final wavenumber>
     <#points to evaluate> <scan complex wavenumbers> <#panels>
     <order of quadrature rule> <accuracy of arnoldi algorithm>.
~~~

The resulting file will contain the value of <tt>k</tt> in the first column.
The rest of the columns contain the singular values from
smallest to largest for this <tt>k</tt>.
The user will be updated through the command line about the
progress of the algorithm
if <tt>-DCMDL</tt> is set.

#### <tt>transmission_problem_verification</tt>
This target builds a script that computes solutions to
the analytically solvable case of the Helmholtz transmission
problem where the scatterer is a circle using second-kind direct
BIEs and Galerkin BEM.
The results are written to file.
The script can be run as follows:
~~~
/path/to/library/bin/transmission_problem_verification 
    <radius of circle> <#coeffs for series expansion of solution> 
    <refraction inside> <refraction outside> <initial wavenumber>
    <order of quadrature rule> <outputfile>.
~~~
This output file will contain two columns.
The first will contain the current panel size.
The second will contain the residual error in the euclidean norm of the computed FEM-space interpolation coefficients to the known FEM-space interpolation coefficients for the current number of panels.
The user will be updated through the command line about the progress of the algorithm if <tt>-DCMDL</tt> is set.

