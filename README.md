# AdaptiveQuad - An Adaptive Gaussian Quadrature Implementation

This is an implementation of a modified version of the classic Gander and Gautschi adaptive Gauss-Lobatto/Kronrod quadrature algorithm [1]. It has the same broad properties: Estimation of integrals with efficient use of integrand evaluations and performance which is responsive to the user's selected precision requirement. The main modifications are a change to subdividing intervals into three subintervals, rather than six, and an alteration to how the first round of function evaluations are used. In addition, some minor modifications are made to the stopping criteria. 

This implementation is designed to be more parsimonious with function evaluations than the original in many cases, making it advantageous for very simple problems in which no subdivisions of the original interval are required, and for integrands which are expensive to evaluate. To do this, it trades for some additional overhead which may make it slightly slower when evaluating a cheap integrand which requires a small number of interval subdivisions. 

## Usage

This implementation is header-only, so it is sufficient to `#include <AdaptiveQuad.h>`. Simple use requires only calling the primary function:

	double parameter1=5, parameter2=7;
	auto MyFunction=[=](double x)->double{
		return parameter1*sin(x) + x/parameter2;
	};
	double integral=AdaptiveQuad::integrate(MyFunction,4.8,11.2);
	std::cout << "Integral: " << integral << std::endl;
	//Integral: 6.73676

## Differences from the Original

The main difference between this implementation and the original `adaptlob` (a translation of which is included as part of the test suite) is the choice to subdivide intervals into three rather than six. The original rationale for dividing into six subintervals was to reuse every point in the original interval which was available from evaluating the four-point Gauss-Lobatto rule and the seven-point Kronrod extension (which itself uses all of the points used by the four-point rule). This seems highly advantageous as no points are 'wasted' when subdividing. However, whenever a subdivision occurs new points must be evaluated in the subintervals, and the division must eventually terminate, at which point no points in the final subintervals will be reused. The six-fold division requires 6 * 5 = 30 new points to be evaluated. In contrast, a three-fold subdivision reusing only the points for the base four-point rule may not reuse the three points used for the seven-point extension, but requires generating only 3 * 5 = 15 new points. This can be viewed as an optimistic bias which hopes that further subdivision will not be required, and so does less of it to start with. It can also be useful when the regions of the integral which are difficult are small and few in number, as subdivisions can be more readily focused on them for the same number of function evaluations. 

Another difference is specific to very simple problems which require no subdivision. In such cases, this implementation will make no further function evaluations after those required for the initial 13-point Gauss-Kronrod calculation, and it will return the result of that calculation, rather than one of the lower-order estimates. The reduction in the minimum number of function evaluations from 18 to 13 is accomplished by returning early from the main integration function, unlike `adaptlob`'s unconditionally calling `adaptlobstp`. Doing so makes it simple to also return the result of the 13-point calculation, since it is available and the termination conditions being satisfied implies that it is consistent with the lower order estimates, so one may as well take it in the hope that it will give extra precision 'for free'. When integrating easy functions this often works out to be the case. 

Structurally, this version of the algorithm keeps more logic in the top-level driver function, and uses heap-based storage for subintervals, rather than using recursion. This does impose some additional overhead when evaluating a small number of subdivisions with an integrand which is cheap to compute. Persons for whom the greatest possible speed is required in such a situation may wish to use a different, or more specialized implementation. 

## Tests

A simple benchmark tool is included in the `test` directory, which can be built with the included Makefile. It compares the performance of this implementation to the original and a few others, including a basic fifth-order Romberg method, the [GNU Scientific Library](https://www.gnu.org/software/gsl/)'s QAG implementation (if GSL is available when the benchmark is compiled), and the version of the Gander and Gautschi algorithm presented in the third edition of [Numerical Recipes](https://www.numerical.recipes) if available. (The interested user will have to supply their own copy of this latter implementation, as its license does not permit it to be redistributed.)

----

[1]: Gander, W and Gautschi, W, "Adaptive Quadrature Revisited", [https://doi.org/10.3929/ethz-a-006652954](https://doi.org/10.3929/ethz-a-006652954) (1998)