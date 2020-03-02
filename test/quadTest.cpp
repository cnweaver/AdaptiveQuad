#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "../AdaptiveQuad.h"
#include "gg.h"
#include "gsl_interface.h"
#include "romberg.h"

#ifdef HAVE_GSL
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_dilog.h>
#endif

//The NR implementation cannot be redistributed even for testing purposes. . . 
#if __has_include("nr_adapt.h")
namespace NR{ //inject into a namespace for neatness
#include "nr_adapt.h"
}
#define HAVE_NR 1
#endif

struct FuncBase;

using FuncPtr=std::unique_ptr<FuncBase>;

struct FuncBase{
	FuncBase():calls(0){}
	virtual ~FuncBase();
	
	double operator()(double x) const{
		calls++;
		return Evaluate(x);
	}
	
	virtual FuncPtr copy() const{
		std::cout << __PRETTY_FUNCTION__ << std::endl;
		return FuncPtr(new FuncBase);
	}
	
	virtual double Evaluate(double) const{ return 0; }
	
	mutable unsigned int calls;
	
	virtual double leftBoundary() const{ return 0; }
	virtual double rightBoundary() const{ return 0; }
	virtual double exactValue() const{ return 0; }
	virtual const char* name() const{ return "none"; }
};

FuncBase::~FuncBase(){}

template<typename FuncType>
struct FuncCopy : public FuncBase{
	virtual FuncPtr copy() const{
		return FuncPtr(new FuncType(static_cast<const FuncType&>(*this)));
	}
};

struct SplitTrig1 : public FuncCopy<SplitTrig1>{
	virtual double Evaluate(double x) const override{
		if(x<=2)
			return std::sin(x);
		else
			return std::cos(x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 4*atan(1); }
	virtual double exactValue() const override{ return -cos(2)+1+0-sin(2); }
	virtual const char* name() const override{ 
		return "if x<=2 then sin(x) else cos(x) fi, Interval:[0,pi]"; 
	}
};

struct SplitTrig2 : public FuncCopy<SplitTrig2>{
	virtual double Evaluate(double x) const override{
		if(x<=1)
			return std::sin(x);
		else
			return std::cos(x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 2; }
	virtual double exactValue() const override{ return -cos(1)+1+sin(2)-sin(1); }
	virtual const char* name() const override{ return "if x<=1 then sin(x) else cos(x) fi, Interval:[0,2]"; }
};

struct Exp : public FuncCopy<Exp>{
	virtual double Evaluate(double x) const override{
		return exp(x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{ return exp(1)-exp(0); }
	virtual const char* name() const override{ return "exp(x), Interval:[0,1]"; }
};

struct Step : public FuncCopy<Step>{
	virtual double Evaluate(double x) const override{
		return (x>=3?1:0);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 10; }
	virtual double exactValue() const override{ return 7; }
	virtual const char* name() const override{ return "if x>=3 then 1 else 0 fi, Interval:[0,10]"; }
};

struct Sqrt : public FuncCopy<Sqrt>{
	virtual double Evaluate(double x) const override{
		return sqrt(x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{ return 2./3; }
	virtual const char* name() const override{ return "sqrt(x), Interval:[0,1]"; }
};

struct CoshCos : public FuncCopy<CoshCos>{
	virtual double Evaluate(double x) const override{
		return (23./25)*cosh(x)+cos(x);
	}
	virtual double leftBoundary() const override{ return -1; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{ return (23./25)*(sinh(1)-sinh(-1))+sin(1)-sin(-1); }
	virtual const char* name() const override{ return "(23/25)*cosh(x)-cos(x), Interval:[-1,1]"; }
};

struct Rat1 : public FuncCopy<Rat1>{
	virtual double Evaluate(double x) const override{
		return 1./(x*x*x*x+x*x+0.9);
	}
	virtual double leftBoundary() const override{ return -1; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{
		auto antiDeriv=[](double x){
			std::complex<double> i(0,1);
			auto ir65=i*std::sqrt(65);
			auto r5m=std::sqrt(5.-ir65);
			auto r5p=std::sqrt(5.+ir65);
			return -sqrt(5./117)*imag((r5m*std::atan(x/(sqrt(0.1)*r5p)) -
			                           r5p*std::atan(x/(sqrt(0.1)*r5m))));
		};
		return antiDeriv(rightBoundary())-antiDeriv(leftBoundary());
	}
	virtual const char* name() const override{ return "1/(x^4+x^2+0.9), Interval:[-1,1]"; }
};

struct Sqrt3 : public FuncCopy<Sqrt3>{
	virtual double Evaluate(double x) const override{
		return sqrt(x*x*x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{ return 2./5; }
	virtual const char* name() const override{ return "sqrt(x^3), Interval:[0,1]"; }
};

struct SqrtInv : public FuncCopy<SqrtInv>{
	virtual double Evaluate(double x) const override{
		if(x==0)
			return 0;
		return 1./sqrt(x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{ return 2; }
	virtual const char* name() const override{ return "if x>0 then 1/sqrt(x) else 0 fi, Interval:[0,1]"; }
};

struct Rat2 : public FuncCopy<Rat2>{
	virtual double Evaluate(double x) const override{
		return 1./(1+x*x*x*x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{
		auto antiDeriv=[](double x)->double{
			double s=sqrt(2);
			return (-log(x*x-s*x+1)+log(x*x+s*x+1)-2*atan(1-s*x)+2*atan(1+s*x))/(4*s);
		};
		return antiDeriv(rightBoundary())-antiDeriv(leftBoundary());
	}
	virtual const char* name() const override{ return "1/(1+x^4), Interval:[0,1]"; }
};

struct SinInv : public FuncCopy<SinInv>{
	virtual double Evaluate(double x) const override{
		static const double pi=4*atan(1);
		return 2/(2+sin(10*pi*x));
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{ return 2/sqrt(3); }
	virtual const char* name() const override{ return "2/(2+sin(10*pi*x)), Interval:[0,1]"; }
};

struct Hyperb : public FuncCopy<Hyperb>{
	virtual double Evaluate(double x) const override{
		return 1/(1+x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{ return log(2); }
	virtual const char* name() const override{ return "1/(1+x), Interval:[0,1]"; }
};

struct ExpInv : public FuncCopy<ExpInv>{
	virtual double Evaluate(double x) const override{
		return 1./(1+exp(x));
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{
		auto antiDeriv=[](double x)->double{
			return x-log(1+exp(x));
		};
		return antiDeriv(rightBoundary())-antiDeriv(leftBoundary());
	}
	virtual const char* name() const override{ return "1/(1+exp(x)), Interval:[0,1]"; }
};

#ifdef HAVE_GSL //need GSL to compute the dilogarthim
struct XExpInv : public FuncCopy<XExpInv>{
	virtual double Evaluate(double x) const override{
		if(x==0)
			return 0;
		return x/(exp(x)-1);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{
		static const double pi=4*atan(1.);
		double Li2e=gsl_sf_dilog(exp(1));
		double l1me=std::real(std::log(std::complex<double>(1,0)-exp(1)));
		return (Li2e-1/2.+l1me)-(pi*pi/6);
	}
	virtual const char* name() const override{ return "if x>0 then x/(1+exp(x)) else x=0 fi, Interval:[0,1]"; }
};
#endif //HAVE_GSL

#ifdef HAVE_GSL //need GSL to compute the Sine integral
struct Sinc : public FuncCopy<Sinc>{
	virtual double Evaluate(double x) const override{
		static const double pi=4.*atan(1);
		return std::sin(100*pi*x)/(pi*x);
	}
	virtual double leftBoundary() const override{ return 0.1; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{
		static const double pi=4.*atan(1);
		return (gsl_sf_Si(100*pi*rightBoundary())-gsl_sf_Si(100*pi*leftBoundary()))/pi;
	}
	virtual const char* name() const override{ return "sin(100*pi*x)/(pi*x), Interval:[0.1,1]"; }
};
#endif //HAVE_GSL

struct Gauss : public FuncCopy<Gauss>{
	virtual double Evaluate(double x) const override{
		static const double pi=4.*atan(1);
		return sqrt(50)*exp(-50*pi*x*x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 10; }
	virtual double exactValue() const override{ return erf(sqrt(2000*atan(1)))/2; }
	virtual const char* name() const override{ return "sqrt(50)*exp(-50*pi*x^2), Interval:[0,1]"; }
};

struct ExpFall : public FuncCopy<ExpFall>{
	virtual double Evaluate(double x) const override{
		return 25*exp(-25*x);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 10; }
	virtual double exactValue() const override{ return 1-exp(-250); } //this will round to 1. . . 
	virtual const char* name() const override{ return "25*exp(-25*x), Interval:[0,10]"; }
};

struct Quad : public FuncCopy<Quad>{
	virtual double Evaluate(double x) const override{
		static const double pi=4.*atan(1);
		return (50/pi)*(2500*x*x+1);
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 10; }
	virtual double exactValue() const override{ return (50/(4.*atan(1)))*((2500./3)*1e3+10); }
	virtual const char* name() const override{ return "(50/pi)*(2500*x^2+1), Interval:[0,10]"; }
};

#ifdef HAVE_GSL //need GSL to compute the Sine integral
struct SincSq : public FuncCopy<SincSq>{
	virtual double Evaluate(double x) const override{
		static const double pi=4.*atan(1);
		double z=50*pi*x;
		double t=std::sin(z)/z;
		return 50*t*t;
	}
	virtual double leftBoundary() const override{ return 0.01; }
	virtual double rightBoundary() const override{ return 1; }
	virtual double exactValue() const override{
		static const double pi=4.*atan(1);
		auto antiderivative=[](double x){
			return (100*pi*x*gsl_sf_Si(100*pi*x)+cos(100*pi*x)-1)/(100*pi*pi*x);
		};
		return antiderivative(rightBoundary())-antiderivative(leftBoundary());
	}
	virtual const char* name() const override{ return "50*(sin(50*pi*x)/(50*pi*x))^2, Interval:[0.01,1]"; }
};
#endif //HAVE_GSL

struct TrigNest : public FuncCopy<TrigNest>{
	virtual double Evaluate(double x) const override{
		return cos(cos(x)+3*sin(x)+2*cos(2*x)+3*sin(2*x)+3*cos(3*x));
	}
	virtual double leftBoundary() const override{ return 0; }
	virtual double rightBoundary() const override{ return 4.*atan(1); }
	virtual double exactValue() const override{ return 0; throw std::runtime_error("Not implemented"); }
	virtual const char* name() const override{ return "cos(cos(x)+3sin(x)+2cos(2x)+3sin(2x)+3cos(3x), Interval:[0,pi]"; }
};

struct TestResult{
	double relativeError;
	double time;
	unsigned long evaluations;
	bool failed;
	
	TestResult():
	relativeError(std::numeric_limits<double>::quiet_NaN()),
	time(std::numeric_limits<double>::quiet_NaN()),
	evaluations(0),
	failed(false)
	{}
};

using ResultMap=std::map<std::string,TestResult>;

ResultMap testFunction(const FuncBase& func, double tolerance){
	ResultMap results;
	
	std::cout << "Integrating f(x)=" << func.name() << std::endl;
	double exact=func.exactValue();
	std::cout << " Exact:     " << exact << std::endl;
	
	const double a=func.leftBoundary();
	const double b=func.rightBoundary();
	
	try{
		auto f=func.copy();
		double i=GG::adaptlob(*f,a,b,tolerance);
		TestResult result;
		result.relativeError=(i-exact)/exact;
		result.evaluations=f->calls;
		results.emplace("adaptlob",result);
		std::cout << " adaptlob:  " << i << " with " << f->calls << " evaluations, relative error " << (i-exact)/exact << std::endl;
	}catch(std::exception& ex){
		TestResult result;
		result.failed=true;
		results.emplace("adaptlob",result);
		std::cout << " adaptlob:  " << ex.what() << std::endl;
	}
	
	try{
		auto f=func.copy();
		double i=Romberg::rombergIntegrate(*f,a,b,tolerance);
		TestResult result;
		result.relativeError=(i-exact)/exact;
		result.evaluations=f->calls;
		results.emplace("Romberg",result);
		std::cout << " Romberg:   " << i << " with " << f->calls << " evaluations, relative error " << result.relativeError << std::endl;
	}catch(std::exception& ex){
		TestResult result;
		result.failed=true;
		results.emplace("Romberg",result);
		std::cout << " Romberg:   " << ex.what() << std::endl;
	}
	
#ifdef HAVE_GSL
	try{
		auto f=func.copy();
		GSL::reset_gsl_error();
		double i=GSL::integrate(*f,a,b,tolerance);
		if(GSL::gsl_errno)
			throw std::runtime_error(GSL::errMsg);
		TestResult result;
		result.relativeError=(i-exact)/exact;
		result.evaluations=f->calls;
		results.emplace("GSL::QAG",result);
		std::cout << " GSL::QAG:  " << i << " with " << f->calls << " evaluations, relative error " << (i-exact)/exact << std::endl;
	}catch(std::exception& ex){
		TestResult result;
		result.failed=true;
		results.emplace("GSL::QAG",result);
		std::cout << " GSL::QAG:  " << ex.what() << std::endl;
	}
#endif //HAVE_GSL
	
#ifdef HAVE_NR
	try{
		auto f=func.copy();
		NR::Adapt adapt(tolerance);
		double i=adapt.integrate(*f,a,b);
		TestResult result;
		result.relativeError=(i-exact)/exact;
		result.evaluations=f->calls;
		results.emplace("NR::Adapt",result);
		std::cout << " NR::Adapt: " << i << " with " << f->calls << " evaluations, relative error " << (i-exact)/exact << std::endl;
	}catch(std::exception& ex){
		TestResult result;
		result.failed=true;
		results.emplace("NR::Adapt",result);
		std::cout << " NR::Adapt: " << ex.what() << std::endl;
	}
#endif //HAVE_NR
	
	try{
		auto f=func.copy();
		AdaptiveQuad::Options opt;
		opt.fa=(*f)(a);
		opt.fb=(*f)(b);
		double i=AdaptiveQuad::integrate(*f,a,b,tolerance,&opt);
		TestResult result;
		result.relativeError=(i-exact)/exact;
		result.evaluations=f->calls;
		results.emplace("AdaptQuad",result);
		std::cout << " AdaptQuad: " << i << " with " << f->calls << " evaluations, relative error " << (i-exact)/exact << std::endl;
	}catch(std::exception& ex){
		TestResult result;
		result.failed=true;
		results.emplace("AdaptQuad",result);
		std::cout << " AdaptQuad: " << ex.what() << std::endl;
	}
	
	return results;
}

struct PrecisionSummary{
	unsigned int successes, minorFailures, majorFailures;
	
	PrecisionSummary():successes(0),minorFailures(0),majorFailures(0){}
};

std::ostream& operator<<(std::ostream& os, const PrecisionSummary& sum){
	os << sum.successes << " success" << (sum.successes!=1?"es, ":", ");
	os << sum.minorFailures << " minor failure" << (sum.minorFailures!=1?"s, ":", ");
	os << sum.majorFailures << " major failure" << (sum.majorFailures!=1?"s, ":", ");
	return os;
}

template<typename T=double>
class RunningStat{
public:
	RunningStat():n(0){}
	void clear(){ n=0; }
	void insert(T x){
		n++;
		
		//Knuth TAOCP vol 2, 3rd edition, page 232
		if (n == 1){
			oldM = newM = x;
			oldS = 0.0;
		}
		else{
			newM = oldM + (x - oldM)/n;
			newS = oldS + (x - oldM)*(x - newM);
			//prepare for next iteration
			oldM = newM; 
			oldS = newS;
		}
	}
	uint64_t count() const{ return n; }
	T mean() const{ return (n ? newM : T(0.)); }
	T variance() const{ return ( n ? newS/(n-1) : T(0.) ); }
	T stddev() const{ return sqrt( variance() ); }
private:
	uint64_t n;
	T oldM, newM, oldS, newS;
};

std::map<std::string,double> timeFunction(const FuncBase& func, double tolerance, std::size_t iterations){
	std::map<std::string,double> results;
	
	std::cout << "Timing integrating f(x)=" << func.name() << std::endl;
	std::cout.setf(std::ios::fixed,std::ios::floatfield);
	std::cout.precision(0);
	
	const double a=func.leftBoundary();
	const double b=func.rightBoundary();
	
	try{
		auto f=func.copy();
		RunningStat<> stats;
		unsigned int rep=0;
		//run for at least 20 and not more than 1010 repetitions, and stop if
		//uncertainty falls below 1% or time taken exceeds 10 seconds
		for(; (rep<20 || stats.stddev()>stats.mean()/100) && rep<1010 
			&& (rep<20 || iterations*stats.mean()*stats.count()<1e10); rep++){
			std::chrono::high_resolution_clock::time_point t1, t2;
			t1 = std::chrono::high_resolution_clock::now();
			for(unsigned int j=0; j<iterations; j++)
				double i=GG::adaptlob(*f,a,b,tolerance);
			t2 = std::chrono::high_resolution_clock::now();
			double time=std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/double(iterations);
			if(rep>=10) //ignore early values, hopefully get system to steady state
				stats.insert(time);
		}
		std::cout << " adaptlob:  " << stats.mean() << '(' << stats.stddev() << ") ns with " 
			<< f->calls/(iterations*rep) << " evaluations" << std::endl;
		results.emplace("adaptlob",stats.mean());
	}catch(std::exception& ex){
		std::cout << " adaptlob:  " << ex.what() << std::endl;
		results.emplace("adaptlob",std::numeric_limits<double>::quiet_NaN());
	}
	
	try{
		auto f=func.copy();
		RunningStat<> stats;
		unsigned int rep=0;
		//run for at least 20 and not more than 1010 repetitions, and stop if
		//uncertainty falls below 1% or time taken exceeds 10 seconds
		for(; (rep<20 || stats.stddev()>stats.mean()/100) && rep<1010 
			&& (rep<20 || iterations*stats.mean()*stats.count()<1e10); rep++){
			std::chrono::high_resolution_clock::time_point t1, t2;
			t1 = std::chrono::high_resolution_clock::now();
			for(unsigned int j=0; j<iterations; j++)
				double i=Romberg::rombergIntegrate(*f,a,b,tolerance);
			t2 = std::chrono::high_resolution_clock::now();
			double time=std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/double(iterations);
			if(rep>=10) //ignore early values, hopefully get system to steady state
				stats.insert(time);
		}
		std::cout << " Romberg:   " << stats.mean() << '(' << stats.stddev() << ") ns with " 
		<< f->calls/(iterations*rep) << " evaluations" << std::endl;
		results.emplace("Romberg",stats.mean());
	}catch(std::exception& ex){
		std::cout << " Romberg:   " << ex.what() << std::endl;
		results.emplace("Romberg",std::numeric_limits<double>::quiet_NaN());
	}
	
#ifdef HAVE_GSL
	try{
		auto f=func.copy();
		RunningStat<> stats;
		unsigned int rep=0;
		//run for at least 20 and not more than 1010 repetitions, and stop if
		//uncertainty falls below 1% or time taken exceeds 10 seconds
		for(; (rep<20 || stats.stddev()>stats.mean()/100) && rep<1010 
			&& (rep<20 || iterations*stats.mean()*stats.count()<1e10); rep++){
			std::chrono::high_resolution_clock::time_point t1, t2;
			t1 = std::chrono::high_resolution_clock::now();
			for(unsigned int j=0; j<iterations; j++){
				GSL::reset_gsl_error();
				double i=GSL::integrate(*f,a,b,tolerance);
				if(GSL::gsl_errno)
					throw std::runtime_error(GSL::errMsg);
			}
			t2 = std::chrono::high_resolution_clock::now();
			double time=std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/double(iterations);
			if(rep>=10) //ignore early values, hopefully get system to steady state
				stats.insert(time);
		}
		std::cout << " GSL::QAG:  " << stats.mean() << '(' << stats.stddev() << ") ns with " 
		<< f->calls/(iterations*rep) << " evaluations" << std::endl;
		results.emplace("GSL::QAG",stats.mean());
	}catch(std::exception& ex){
		std::cout << " GSL::QAG:  " << ex.what() << std::endl;
		results.emplace("GSL::QAG",std::numeric_limits<double>::quiet_NaN());
	}
#endif //HAVE_GSL
	
#ifdef HAVE_NR
	try{
		auto f=func.copy();
		RunningStat<> stats;
		unsigned int rep=0;
		//run for at least 20 and not more than 1010 repetitions, and stop if
		//uncertainty falls below 1% or time taken exceeds 10 seconds
		for(; (rep<20 || stats.stddev()>stats.mean()/100) && rep<1010 
			&& (rep<20 || iterations*stats.mean()*stats.count()<1e10); rep++){
			std::chrono::high_resolution_clock::time_point t1, t2;
			t1 = std::chrono::high_resolution_clock::now();
			for(unsigned int j=0; j<iterations; j++){
				NR::Adapt adapt(tolerance);
				double i=adapt.integrate(*f,a,b);
			}
			t2 = std::chrono::high_resolution_clock::now();
			double time=std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/double(iterations);
			if(rep>=10) //ignore early values, hopefully get system to steady state
				stats.insert(time);
		}
		std::cout << " NR::Adapt: " << stats.mean() << '(' << stats.stddev() << ") ns with " 
		<< f->calls/(iterations*rep) << " evaluations" << std::endl;
		results.emplace("NR::Adapt",stats.mean());
	}catch(std::exception& ex){
		std::cout << " NR::Adapt: " << ex.what() << std::endl;
		results.emplace("NR::Adapt",std::numeric_limits<double>::quiet_NaN());
	}
#endif //HAVE_NR
	
	try{
		auto f=func.copy();
		RunningStat<> stats;
		unsigned int rep=0;
		//run for at least 20 and not more than 1010 repetitions, and stop if
		//uncertainty falls below 1% or time taken exceeds 10 seconds
		for(; (rep<20 || stats.stddev()>stats.mean()/100) && rep<1010 
			&& (rep<20 || iterations*stats.mean()*stats.count()<1e10); rep++){
			std::chrono::high_resolution_clock::time_point t1, t2;
			t1 = std::chrono::high_resolution_clock::now();
			for(unsigned int j=0; j<iterations; j++)
				double i=AdaptiveQuad::integrate(*f,a,b,tolerance);
			t2 = std::chrono::high_resolution_clock::now();
			double time=std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/double(iterations);
			if(rep>=10) //ignore early values, hopefully get system to steady state
				stats.insert(time);
		}
		std::cout << " AdaptQuad: " << stats.mean() << '(' << stats.stddev() << ") ns with " 
		<< f->calls/(iterations*rep) << " evaluations" << std::endl;
		results.emplace("AdaptQuad",stats.mean());
	}catch(std::exception& ex){
		std::cout << " AdaptQuad: " << ex.what() << std::endl;
		results.emplace("AdaptQuad",std::numeric_limits<double>::quiet_NaN());
	}
	
	std::cout.unsetf(std::ios::floatfield);
	std::cout.precision(16);
	return results;
}

int main(int argc, char* argv[]){
	double tolerance=1e-8;
	bool doTiming=false;
	std::size_t timingIterations=10000;
	std::cout.precision(16);
	
	FuncPtr st1(new SplitTrig1);
	std::map<std::string,FuncPtr> testFunctions;
	testFunctions.emplace("SplitTrig1",FuncPtr(new SplitTrig1));
	testFunctions.emplace("Exp",FuncPtr(new Exp));
	testFunctions.emplace("Step",FuncPtr(new Step));
	testFunctions.emplace("Sqrt",FuncPtr(new Sqrt));
	testFunctions.emplace("CoshCos",FuncPtr(new CoshCos));
	testFunctions.emplace("Rat1",FuncPtr(new Rat1));
	testFunctions.emplace("Sqrt3",FuncPtr(new Sqrt3));
	testFunctions.emplace("SqrtInv",FuncPtr(new SqrtInv));
	testFunctions.emplace("Rat2",FuncPtr(new Rat2));
	testFunctions.emplace("SinInv",FuncPtr(new SinInv));
	testFunctions.emplace("Hyperb",FuncPtr(new Hyperb));
	testFunctions.emplace("ExpInv",FuncPtr(new ExpInv));
	testFunctions.emplace("XExpInv",FuncPtr(new XExpInv));
	testFunctions.emplace("Sinc",FuncPtr(new Sinc));
	testFunctions.emplace("Gauss",FuncPtr(new Gauss));
	testFunctions.emplace("ExpFall",FuncPtr(new ExpFall));
	testFunctions.emplace("Quad",FuncPtr(new Quad));
	testFunctions.emplace("SincSq",FuncPtr(new SincSq));
	//testFunctions.emplace("TrigNest",FuncPtr(new TrigNest));
	
	std::set<std::string> functionsToTest;
	
	for(int i=1; i<argc; i++){
		std::string arg=argv[i];
		if(arg.empty())
			continue;
		if(arg=="--help" || arg=="-h" || arg=="help"){
			std::cout << "quadTest [--help] [--list-functions] [tolerance] [function ...]" << std::endl;
			return 0;
		}
		if(arg=="--list-functions" || arg=="-l"){
			unsigned long nameWidth=0;
			for(const auto& func : testFunctions)
				nameWidth=std::max(func.first.size(),nameWidth);
			for(const auto& func : testFunctions)
				std::cout << std::setw(nameWidth+2) << std::left << func.first+":" << func.second->name() << '\n';
			return 0;
		}
		if(arg=="--time" || arg=="-t"){
			doTiming=true;
			continue;
		}
		if(std::isdigit(arg.front())){
			try{
				tolerance=std::stod(arg);
			}catch(...){
				std::cerr << "Failed to interpret " << arg << " as a tolerance" << std::endl;
				return 1;
			}
			continue;
		}
		if(arg=="all"){
			for(const auto& func : testFunctions)
				functionsToTest.insert(func.first);
			continue;
		}
		if(testFunctions.find(arg)!=testFunctions.end()){
			functionsToTest.insert(arg);
			continue;
		}
		std::cerr << "quadTest: Unrecognized argument '" << arg << "'\n";
		return 1;
	}
	//if the user requested nothing give him/her _everything_
	if(functionsToTest.empty()){
		for(const auto& func : testFunctions)
			functionsToTest.insert(func.first);
	}
	
#ifdef HAVE_GSL
	//The GSL docs admonish us not to do error handling this way, and to instead
	//check the return codes of library functions. However, the integration 
	//functions do not return error codes for us to check, and will simply 
	//abort() the whole program if we don't set our own error handler. 
	gsl_set_error_handler(GSL::error_handler);
#endif
	
	std::cout << "Testing with tolerance " << tolerance << std::endl;
	std::map<std::string,PrecisionSummary> precSums;
	for(const auto& funcName : functionsToTest){
		const auto& func=testFunctions[funcName];
		std::cout << '(' << funcName << ") ";
		auto results=testFunction(*func,tolerance);
		for(const auto& result : results){
			PrecisionSummary& summary=precSums[result.first];
			double precFactor=tolerance/std::abs(result.second.relativeError);
			if(precFactor>1)
				summary.successes++;
			else if(precFactor>0.1)
				summary.minorFailures++;
			else
				summary.majorFailures++;
		}
	}
	
	std::cout << '\n';
	for(const auto& summary : precSums)
		std::cout << summary.first << ": " << summary.second << std::endl;
	
	std::cout << '\n';
	for(const auto& funcName : functionsToTest){
		const auto& func=testFunctions[funcName];
		std::cout << '(' << funcName << ") ";
		auto results=timeFunction(*func,tolerance,timingIterations);
	}
}
