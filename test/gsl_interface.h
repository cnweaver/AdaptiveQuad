#ifndef GSL_INTERFACE_H

#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

namespace GSL{
	static std::string errMsg;
	static int gsl_errno;
	void error_handler(const char* reason, const char* file, int line, int gsl_errno){
		GSL::gsl_errno=gsl_errno;
		errMsg=file+std::string(":")+std::to_string(line)+": "+reason;
	}
	
	void reset_gsl_error(){ gsl_errno=0; }
	
	///\class
	///\brief Container for GSL workspace to be use with the integrators.
	class IntegrateWorkspace {
	private:
		
	public:
		gsl_integration_workspace* ws;
		IntegrateWorkspace(size_t limit) {
			ws=gsl_integration_workspace_alloc(limit);
		}
		~IntegrateWorkspace() {
			gsl_integration_workspace_free(ws);
		}
	};
	
	///\brief One dimensional integral using GSL.
	/// @param ws GSL integration workspace.
	/// @param f Function to integrate.
	/// @param a Lower integration limit.
	/// @param b Upper integration limit.
	/// @param acc Accuracy parameter.
	/// @param max_iter Maximum number of iterations to perform the integral.
	template<typename FunctionType>
	double integrate(IntegrateWorkspace& ws, FunctionType&& f, double a, double b, double acc=1e-7, unsigned int max_iter=10000){
		//precondition copied from gsl-1.16/integration/qag.c:122
		//assuming epsabs=0 and epsrel=acc
		if ((acc < 50 * GSL_DBL_EPSILON || acc < 0.5e-28))
			throw std::runtime_error("Specified tolerance too small");
		
		using FPtr=decltype(&f);
		double (*wrapper)(double,void*)=[](double x, void* params)->double{
			auto& f=*static_cast<FPtr>(params);
			return(f(x));
		};
		
		double result, error;
		gsl_function F;
		F.function = wrapper;
		F.params = &f;
		
		gsl_integration_qag(&F, a, b, 0, acc, max_iter, GSL_INTEG_GAUSS15, ws.ws, &result, &error);
		
		return(result);
	}
	
	///\brief One dimensional integral using GSL.
	/// @param f Function to integrate.
	/// @param a Lower integration limit.
	/// @param b Upper integration limit.
	/// @param acc Accuracy parameter.
	/// @param max_iter Maximum number of iterations to perform the integral.
	template<typename FunctionType>
	double integrate(FunctionType&& f, double a, double b, double acc=1e-7, unsigned int max_iter=10000, size_t memory_alloc=10000){
		IntegrateWorkspace ws(memory_alloc);
		return integrate(ws, f, a, b, acc, max_iter);
	}
} //namespace GSL
#endif //HAVE_GSL

#endif //GSL_INTERFACE_H