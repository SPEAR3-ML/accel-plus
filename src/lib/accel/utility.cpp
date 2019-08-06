#include <math.h>
#include <float.h>
#include "nr3.h"
#include "fourier.h"
#include "mins.h"
#include "gamma.h"
#include "incgammabeta.h"
#include "fitab.h"
#include "interp_1d.h"
#include "interp_linear.h"
#include "interp_2d.h"
#include "utility.h"
#include <time.h>
#include <cmath>
//using namespace std;
using std::cout;
using std::endl;
using std::cin;
using std::isnan;	
/* NAFF *********************************************/
double fprojection(double nv, int N, double *X)
{
	double ry,iy;
	ry=iy=0;
	
	for(int i=0;i<N;i++)
	{
		double rz, iz;
		rz = cos(2*PI*nv*i);
		iz =-sin(2*PI*nv*i);
		
		double W;
		W = pow(sin(PI*i/N/1.0),4);
		
		ry += rz*X[i]*W;
		iy += iz*X[i]*W;
		
	}
	return -sqrt(ry*ry+iy*iy);	

}
struct Func {
	int N;
	double *X;
	Func(int n, double *nX) : N(n), X(nX){};
	Doub operator()(const Doub x)
	{ return fprojection(x, N, X); }	
};

double naff(int N, double *data,double xlow, double xhi)
{
	VecDoub vd(N, data);
	realft(vd,1);
	
	xlow = MAX(xlow, 0.9/N);
	xhi = MIN(0.5, xhi);
	unsigned int ilow = (int(xlow*N*2)>>1)<<1;
	unsigned int ihi = (int(xhi*N*2)>>1)<<1;
	
	double fm,im;
	im = ilow;
	fm = double(SQR(vd[im])+SQR(vd[im+1]));
	for(int i=ilow+2;i<ihi;i+=2)
	{
		double fnew = SQR(vd[i])+SQR(vd[i+1]);
		if(fnew>fm)
		{	im = i; fm = fnew;}
	}
	xlow = MAX(xlow, (im-2)*1.0/N/2.0);
	xhi = MIN(xhi, (im+2)*1.0/N/2.0);
	
	
	//Brent ginst;
	Golden ginst;
	Func func(N, data);
	//ginst.bracket(0,0.5,func);
	/*
	xlow = MAX(xlow, 0.9/N);
	xhi = MIN(0.5, xhi);
	double nv1,nv2,nv3,fm;
	unsigned int Ndiv = 20;
	nv1 = xlow; nv3 = xhi; nv2 = (nv1+nv3)/2.0;
	fm = func(nv1);
	for(int i=0;i<Ndiv;i++)
	{
		nv2 = nv1+(xhi-xlow)*i/Ndiv*1.0;
	}
	*/
	{
		ginst.ax = xlow; ginst.cx=xhi;  ginst.bx = (ginst.ax+ginst.cx)/2.0;
		double nv1 = ginst.minimize(func);
		//cout<<ginst.xmin<<"\t"<<ginst.fmin<<endl;
		return nv1;
	}
		
}


/*double naff(int N, double *data,double xlow=0.0, double xhi=0.5)
{
	Brent ginst;
	//Golden ginst;
	Func func(N, data);
	//ginst.bracket(0,0.5,func);
	if (xhi-xlow>0.32){
		ginst.ax = MAX(xlow, 0.9/N); ginst.cx=(ginst.ax*0.4+xhi*0.6);  ginst.bx = (ginst.ax+ginst.cx)/2.0;
		double nv1 = ginst.minimize(func);
		//cout<<ginst.xmin<<"\t"<<ginst.fmin<<endl;
		double f1 = ginst.fmin;
		
		ginst.ax = (ginst.ax*.6+xhi*.4)/2.0; ginst.cx=MIN(0.5, xhi);  ginst.bx = (ginst.ax+ginst.cx)/2.0;
		double nv2 = ginst.minimize(func);
		double f2 = ginst.fmin;
		//cout<<ginst.xmin<<"\t"<<ginst.fmin<<endl;
		
		return (f1<f2)?nv1:nv2;	
	}
	else
	{
		ginst.ax = MAX(xlow, 0.9/N); ginst.cx=MIN(0.5, xhi);  ginst.bx = (ginst.ax+ginst.cx)/2.0;
		double nv1 = ginst.minimize(func);
		//cout<<ginst.xmin<<"\t"<<ginst.fmin<<endl;
		return nv1;
	}
		
}
*/

/*End of NAFF *********************************************/

/* Linear Fit *********************************************/

void linearfit(CLinfit &ca,int N, double *xx,double  *yy, double *sig)
{
	VecDoub vx(N,xx), vy(N,yy); 
	if(sig==NULL)
	{
		Fitab fab(vx,vy);
		ca.a = fab.a;
		ca.b = fab.b;
		ca.siga = fab.siga;
		ca.sigb = fab.sigb;
		ca.q = fab.q;
		ca.chi2 = fab.chi2;
	}else
	{
		VecDoub vsig(N,sig);
		Fitab fab(vx,vy,vsig);
		ca.a = fab.a;
		ca.b = fab.b;
		ca.siga = fab.siga;
		ca.sigb = fab.sigb;
		ca.q = fab.q;
		ca.chi2 = fab.chi2;
	}
	
}

/*End of Linear Fit *********************************************/

/*interpolation *********************************************/
//linear interpolation, make sure xl and yl are in ascending order
double interp1(double ytab[],double xl[], int nx, double x)
{
	VecDoub xv(nx), yv(nx);
	for(int i=0;i<nx;i++)
		{
			xv[i] = xl[i];
			yv[i] = ytab[i];
		}
		Linear_interp func(xv,yv);
		return func.interp(x);
}
double interp2(double *ztab[], double xl[],int nx, double yl[],int ny,double x, double y)
{
	MatDoub zm(nx,ny);
	VecDoub xv(nx), yv(ny);
	for(int i=0;i<ny;i++) yv[i] = yl[i];
	for(int j=0;j<nx;j++) xv[j] = xl[j];
	
	for(int i=0;i<ny;i++)
		for(int j=0;j<nx;j++)
			zm[j][i] = ztab[i][j]; 
	Bilin_interp func(xv,yv,zm);
	return func.interp(x,y);
	
}

/*End of interpolation *********************************************/

/*Rand number ****************************************************/


double grandn(double gmean,double gstd)
{
	const double epsilon = DBL_EPSILON;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
	   return z1*gstd+gmean ;

	double u1, u2;
	do
	 {
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(TWOPI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(TWOPI * u2);
	return z0*gstd+gmean ;
	
}
/*
double grandn(double gmean,double gstd)
{
	//NormalGenerator g(gmean,gstd);
	//return g();
	return 0;	
}*/	


/* End of Random number ************************************************/

/*Generation of distriubtion*********************************************/

//double* gen_distr_6D(int N, double bx, double emitx, double by, double emity, double sigz, double sigdpp)
double* gen_distr_6D(int N, double sigx, double sigxp, double sigy, double sigyp, double sigz, double sigdpp)
{
	if(N<=0) return NULL;
	
	double *r = new double[6*N];
/*	double sigx, sigxp, sigy, sigyp;
	sigx = sqrt(bx*emitx);
	sigxp = sqrt(emitx/bx);
	sigy = sqrt(by*emity);
	sigyp = sqrt(emity/by);
*/	
	for(int i=0;i<N;i++)
	{
		double *rn = r+6*i;
		rn[0] = grandn(0,sigx);
		rn[1] = grandn(0,sigxp);
		rn[2] = grandn(0,sigy);
		rn[3] = grandn(0,sigyp);
		rn[4] = grandn(0,sigz);
		rn[5] = grandn(0,sigdpp);
	}
	
	return r;
} 

/*End of Generation of distriubtion*********************************************/

/*String operations *********************************************/

/*vector<string> str_split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}*/

int str_split2(std::string& s, const std::string& delimiter, std::string *tokenlist, int ntoken_max)
{
	int cnt = 0;
	std::string token;
	
	while(token != s){
	  token = s.substr(0,s.find_first_of(delimiter));
	  s = s.substr(s.find_first_of(delimiter) + 1);
	  //printf("%s ",token.c_str());
	  tokenlist[cnt++] = token;
	}
	//std::cout << s << std::endl;
	return cnt;
}	

bool str_icompare(const std::string& a, const std::string& b)
{
    unsigned int sz = a.size();
    if (b.size() != sz)
        return false;
    for (unsigned int i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}
	
/*End of String operations *********************************************/

/*Time operations*********************************************/
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
	//localtime_s(&tstruct, &now);
	strftime(buf, sizeof(buf), "UTC %Y-%m-%d %X", &tstruct);
	
	//asctime_s(buf, sizeof(buf), &tstruct);
    return buf;
}

/*End of time operations *********************************************/

/*Beam statistics**************************/
int calc_beam_stat( double *r, int N, double ar[],double sigr[],double *sigrr=NULL)
{
	for (int i = 0; i < 6; i++) { ar[i] = 0; sigr[i] = 0; }

	int nlive = 0;
	for (int j = 0; j < N; j++)
	{
		double *pr = r + j * 6;
		if (isnan(pr[0])) continue;
		else nlive++;

		for (int i = 0; i < 6; i++) { ar[i] += pr[i]; }
	}
	for (int i = 0; i < 6; i++) { ar[i] /= nlive;  }

	for (int j = 0; j < N; j++)
	{
		double *pr = r + j * 6;
		if (isnan(pr[0])) continue;
		for (int i = 0; i < 6; i++) 
		{
			sigr[i] += SQR(pr[i] - ar[i]);
		}
	}
	for (int i = 0; i < 6; i++) { sigr[i] /= nlive;  sigr[i] = sqrt(sigr[i]); }

	if (sigrr != NULL)
	{
		sigrr[0] = 0; //xx'
		sigrr[1] = 0; //yy'
		sigrr[2] = 0; //xdelta
		sigrr[3] = 0; //x'delta
		sigrr[4] = 0; //yz
		sigrr[5] = 0; //y'z
		for (int j = 0; j < N; j++)
		{
			double *pr = r + j * 6;
			if (isnan(pr[0])) continue;

			sigrr[0] += (pr[0] - ar[0])*(pr[1] - ar[1]);  //xx'
			sigrr[1] += (pr[2] - ar[2])*(pr[3] - ar[3]);  //yy'
			sigrr[2] += (pr[0] - ar[0])*(pr[5] - ar[5]);  //xdelta
			sigrr[3] += (pr[1] - ar[1])*(pr[5] - ar[5]);  //x'delta
			sigrr[4] += (pr[2] - ar[2])*(pr[4] - ar[4]);  //yz
			sigrr[5] += (pr[3] - ar[3])*(pr[4] - ar[4]);  //y'z
		}
		for (int i = 0; i < 6; i++) { sigrr[i] /= nlive;   }

	}
	return nlive;

}
/*End of Beam statistics**************************/