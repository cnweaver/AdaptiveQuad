#include <cmath>

namespace GG{	

//A transliteration of Gander & Gautschi's original MATLAB code. 
//As-is this is suitable only to light-weight use (i.e.) testing as it is not 
//thread safe, and make little effort to properly communicate errors. 
namespace{ //hack to side-step the ODR
bool termination2;
}
	
template<typename FuncType, typename... VarArgs>
double adaptlobstp(FuncType&& f, double a, double b, double fa, double fb, double is, bool trace=false, VarArgs... varargin){
	double h=(b-a)/2, m=(a+b)/2;
	double alpha=sqrt(2./3), beta=1/sqrt(5);
	double mll=m-alpha*h, ml=m-beta*h, mr=m+beta*h, mrr=m+alpha*h;
	double x[5]={mll,ml,m,mr,mrr};
	double y[5];
	for(int i=0; i<5; i++)
		y[i]=f(x[i],varargin...);
	double fmll=y[0], fml=y[1], fm=y[2], fmr=y[3], fmrr=y[4];
	double i2=(h/6.)*(fa+fb+5.*(fml+fmr));
	double i1=(h/1470.)*(77.*(fa+fb)+432.*(fmll+fmrr)+625.*(fml+fmr)
	                     +672.*fm);
	if(is+(i1-i2)==is || mll<=a || b<=mrr){
		if((m<=a || b<=m) && !termination2){
			//Warning: Interval contains no more machine numbers. 
			//         Required tolerance may not be met.
			termination2=true;
		}
		if(trace)
			std::cout << a << ' ' << b-a << ' ' << i1 << '\n';
		return i1;
	}
	else{
		return adaptlobstp(f,a,mll,fa,fmll,is,trace,varargin...)+
		       adaptlobstp(f,mll,ml,fmll,fml,is,trace,varargin...)+
		       adaptlobstp(f,ml,m,fml,fm,is,trace,varargin...)+
		       adaptlobstp(f,m,mr,fm,fmr,is,trace,varargin...)+
		       adaptlobstp(f,mr,mrr,fmr,fmrr,is,trace,varargin...)+
		       adaptlobstp(f,mrr,b,fmrr,fb,is,trace,varargin...);
	}
}

template<typename FuncType, typename... VarArgs>
double adaptlob(FuncType&& f, double a, double b, double tol=0, bool trace=false, VarArgs... varargin){
	static const std::size_t nargin=sizeof...(VarArgs);
	termination2=false;
	if(tol<std::numeric_limits<double>::epsilon())
		tol=std::numeric_limits<double>::epsilon();
	double m=(a+b)/2, h=(b-a)/2;
	double alpha=sqrt(2./3), beta=1/sqrt(5);
	double x1=0.942882415695480;
	double x2=0.641853342345781;
	double x3=0.236383199662150;
	double x[13]={a,m-x1*h,m-alpha*h,m-x2*h,m-beta*h,m-x3*h,m,m+x3*h,
	              m+beta*h,m+x2*h,m+alpha*h,m+x1*h,b};
	double y[13];
	for(int i=0; i<13; i++)
		y[i]=f(x[i],varargin...);
	double fa=y[0], fb=y[12];
	double i2=(h/6.)*(y[0]+y[12]+5.*(y[4]+y[8]));
	double i1=(h/1470.)*(77.*(y[0]+y[12])+432.*(y[2]+y[10])+
	                     625.*(y[4]+y[8])+672.*y[6]);
	double is=h*(0.0158271919734802*(y[0]+y[12])+0.0942738402188500
	             *(y[1]+y[11])+0.1550719873365850*(y[2]+y[10])+
	             0.1888215739601820*(y[3]+y[9])+0.1997734052268590
	             *(y[4]+y[8])+0.224926465333340*(y[5]+y[7])
	             +0.242611071901408*y[6]);
	int s=(is>=0?1:-1);
	double erri1=std::abs(i1-is);
	double erri2=std::abs(i2-is);
	double R=erri1/erri2;
	if(R>0. && R<1.) tol/=R;
	is=s*std::abs(is)*tol/std::numeric_limits<double>::epsilon();
	if(is==0) is=b-a;
	return adaptlobstp(f,a,b,fa,fb,is,trace,varargin...);
}
}