# TransmissionScatteringProblemBEM
## Usage
The library can be configured by running <tt>cmake CMakeLists.txt</tt> in the base directory.
This should automatically install the Eigen library for matrix computations if this is not already available.
HelmholtzTransmissionProblemBEM also relies on the library

https://github.com/joeydumont/complex_bessel

for passing complex arguments to the Hankel and Bessel functions.
To generate the documentation \LaTeX and <tt>doxygen</tt> have to be installed.
Running CMake also configures <tt>make</tt> targets.
These can then be built by running 

<tt>make <target_name></tt>.

The compiled binary can be found in the <tt>bin</tt> directory.
We will show how the built targets are to be used.

#### <tt>doxygen_HelmholtzTransmissionProblemBEM</tt>
This target generates a documentation of the library in the <tt>doxygen/generated_doc</tt> directory.
The documentation can be browsed using any common browser.

#### <tt>dirichlet_example</tt>
This target builds a script that computes the solution to a Dirichlet problem
using first kind direct BIEs.
No command line parameters are necessary.
Once built the script can be run as follows: 
~~~
/path/to/library/bin/dirichlet_example.
~~~
The user will be updated over the residual error in the euclidean norm of the computed FEM-space interpolation coefficients to the known FEM-space interpolation coefficients for the current number of panels through the command line
 
#### <tt>neumann_example</tt>
This target builds a script that computes the solution to a Neumann problem
using first kind direct BIEs.
No command line parameters are necessary.
Once built the script can be run as follows: 
~~~
/path/to/library/bin/neumann_example.
~~~
The user will be updated over the residual error in the euclidean norm of the computed FEM-space interpolation coefficients to the known FEM-space interpolation coefficients for the current number of panels through the command line

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
/path/to/library/bin/parabolic_approximation.
~~~
In the file the first column contains the initial point used for the parabolic approximation.
The next three columns contain the function value and the first two derivatives at the initial point that were used to compute the parabolic approximation.
The user also will get updates on the current best approximation for a minima and the value of the first derivatie at this point through the command line.

#### <tt>roots_brent_circle</tt>
This target builds a script that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the sedond-kind direct BIEs of the Helmholtz
transmission problem using the Van Wijngaarden-Dekker-Brent method.
The scatterer is set to be a circle
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <radius of circle> 
    <refraction inside> <refraction outside> <wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the interval used to compute the root in the first column. 
Then in the next three columns will be the point, the function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.

#### <tt>roots_newton_circle</tt>
This target builds a sript that computes minimas in the smallest singular value of the
Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
transmission problem using the Newton-Raphson method.
The scatterer is set to be a circle
The results are written to disk.
The script can be run as follows:
~~~
/path/to/library/bin/roots_brent_circle <radius of circle> 
    <refraction inside> <refraction outside> <wavenumber> 
    <#grid points for root search> <#panels> 
    <order of quadrature rule> <outputfile>.
~~~
The resulting file will contain the left boundary of the interval used to compute the root in the first column. 
Then in the next three columns will be the point, the function value and the derivative at which the root was found.
The last column will contain the number of iterations used to find the root.
If no root was found the last four columns will be set to <tt>NAN</tt>.

#### <tt>sv_circle</tt>
This target builds a script that computes the singular values
of the Galerkin BEM approximated BIO for the
second-kind direct BIEs of the Helmholtz
transmission problem.
The scatterer is set to be a circle.
The results are written to file.
The script can be run as follows:
~~~
/path/to/library/bin/sv_circle <radius of circle> 
    <refraction inside> <refraction outside> <wavenumber>
    <#panels> <order of quadrature rule> <outputfile>
~~~
The resulting file will contain the value of $k$ in the first column.
The rest of the columns contain the singular values from smallest to largest for this $k$.

#### <tt>sv_square</tt>
This target builds a script that computes the singular values
of the Galerkin BEM approximated BIO for the
second-kind direct BIEs of the Helmholtz
transmission problem.
The scatterer is set to be a square.
The results are written to file.
The script can be run as follows:
~~~
/path/to/library/bin/sv_circle <half of side length of square> 
    <refraction inside> <refraction outside> <wavenumber>
    <#panels> <order of quadrature rule> <outputfile>
~~~
The resulting file will contain the value of $k$ in the first column.
The rest of the columns contain the singular values from smallest to largest for this $k$.


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
    <wavenumber> <#panels> <order of quadrature rule> <outputfile>
~~~
The resulting file will contain the value of $k$ in the first column.
The second column will contain the value of the smallest singular value at this $k$.
Then the columns will contain the computed derivative, the extrapolated derivative, the computed second derivative and the extrapolated second derivative in this order.

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
    <refraction outside> <wavenumber> <#panels> 
    <order of quadrature rule> <outputfile>
~~~
The resulting file will contain the value of $k$ in the first column.
The second column will contain the value of the smallest singular value at this $k$.
Then the columns will contain the computed derivative, the extrapolated derivative, the computed second derivative and the extrapolated second derivative in this order.

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
    <refraction inside> <refraction outside> <wavenumber>
    <#panels> <order of quadrature rule> <outputfile>
~~~
The resulting file will contain the value of $k$ in the first column.
Then the singular values and their first two derivatives at $k$ will be listed from smallest to largest in the columns.
The singular values and their derivatives occupy three neighboring columns.
The final three columns will contain the value of the root, the value of the first derivative at the root and the number of iterations taken to find the root in the interval between the current and the next evaluation point.
If no root was found these three columns will contain <tt>NAN</tt>.

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
    <refraction inside> <refraction outside> <wavenumber>
    <order of quadrature rule> <outputfile>
~~~
This output file will contain two columns.
The first will contain the current panel size.
The second will contain the residual error in the euclidean norm of the computed FEM-space interpolation coefficients to the known FEM-space interpolation coefficients for the current number of panels.
\end{document}
