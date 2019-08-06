#ifndef __ACCEL_H
#define __ACCEL_H
/***
Accel header file
Created by Xiaobiao Huang, SLAC, Dec 2011
***/ 
//#define OPENMP_BY_ELEMENT
#define OPENMP_BY_LINE
#ifdef		OPENMP_BY_ELEMENT
#undef OPENMP_BY_LINE
#endif
#ifdef		OPENMP_BY_LINE
#undef OPENMP_BY_ELEMENT
#endif


#include <iostream>
#include <string>
#include <map>
#include <cmath>
//#include "utility.h"

#define PI  3.1415926535897932384626
#define TWOPI  (2*PI)
#define MAX(A,B) ((A)>(B))?(A):(B)
#define MIN(A,B) ((A)<(B))?(A):(B)
#define SQR(X) ((X)*(X))

#define C0  	2.99792458e8 
#define DEFAULTE0  3000
		
//using namespace std;
using std::cout;
using std::endl;
using std::cin;
using std::string;
using std::istream;
using std::ostream;
using std::map;
using std::ofstream;
using std::ifstream;
using std::isnan;

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656
#define MAX_MPOLE   10 

#define CPPNAN (std::numeric_limits<double>::quiet_NaN())

class ALine;

class AccelElement {
	protected:
		double len;
		string famname;
		string type;
	public:
		AccelElement(void){len = 0; famname="undefined"; type="undefined";};
		AccelElement(double len);
		AccelElement(double len,  const string& famname, const string& type);
		
		string& gettype(){return type; };
		string& getfamname(){return famname; };
		double getlength(){return len;};
		virtual void serialize(ostream& os)=0;
		virtual AccelElement* create(istream& is,string& famname)=0;
		static AccelElement* deserialize(istream& is);
		virtual void pass(double *r, int N)=0;
		virtual void setfield(void *p, const string pname);
		virtual void getfield(void *p, const string pname);
		
		friend class ALine;
		static map<string, AccelElement*> maptype;
		static void initializemap();
};

class ADrift : public AccelElement {
	public:
		ADrift(){len=0; famname="undefined";type="Drift";};
		ADrift(double len,const string& fam="undefined") : AccelElement(len){famname=fam;type="Drift";};
		ADrift& operator=(ADrift&);
		void serialize(ostream& os);
		ADrift* create(istream& is, string& famname);
		void pass(double *r, int N);
};

class AQuad : public AccelElement {
	protected:
		double K;
	public:
		AQuad(){len=0; K=0; famname="undefined";type="Quad";};
		AQuad(double len, double K,const string& fam="undefined") : AccelElement(len)
			{this->K = K; famname=fam;type="Quad";};
		AQuad& operator=(AQuad& a){len = a.len; K=a.K; famname=a.famname; type=a.type; return *this;};
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		AQuad* create(istream& is, string& famname);
		void pass(double *r, int N);
};
class ACorrector : public AccelElement {
	double kickx,kicky;
	public:
		ACorrector(){len=0; kickx=kicky=0;famname="undefined";type="Corrector";};
		ACorrector(double len,double kx, double ky, const string& fam="undefined") : AccelElement(len){kickx=kx; kicky=ky;famname=fam;type="Corrector";};
		ACorrector& operator=(ACorrector& a){len=a.len; kickx=a.kickx; kicky=a.kicky;famname=a.famname; type=a.type; return *this;};
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		ACorrector* create(istream& is, string& famname);
		void pass(double *r, int N);
};

class AAperture : public AccelElement {
	double xlim[2];
	double ylim[2];
	public:
		AAperture(){len=0; xlim[0]=ylim[0]=-1.0;xlim[1]=ylim[1]=1.0; famname="undefined";type="Aperture";};
		AAperture(double len,double *xl, double* yl, const string& fam="undefined") : AccelElement(len){xlim[0]=xl[0];ylim[0]=yl[0]; xlim[1]=xl[1]; ylim[1]=yl[1];famname=fam;type="Aperture";};
		AAperture& operator=(AAperture& a){len=a.len; xlim[0]=a.xlim[0]; xlim[1]=a.xlim[1]; ylim[0]=a.ylim[0]; ylim[1]=a.ylim[1]; famname=a.famname; type=a.type; return *this;};
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		AAperture* create(istream& is, string& famname);
		void pass(double *r, int N);
};

class AMonitor : public AccelElement {
	private:
		unsigned int flag_print; //0, do nothing; 1, write coordinates; 2, average and moment matrix; 3, histogram
		ofstream* pfs;
		int  cnt_pass;
	public:
		double ratio_print; //ratio of particles whose coordinates to be printed out
		int interval_turn;
		string outfile;
		
		AMonitor(){len=0; flag_print=0; pfs=NULL; famname="undefined";type="Monitor"; ratio_print=1.0; interval_turn=1; cnt_pass=0;
			outfile=string("");};
		AMonitor(const string& fname, const string& fam="undefined") : AccelElement(0){flag_print=0; this->len = 0;pfs=NULL; famname=fam;type="Monitor"; ratio_print=1.0; interval_turn=1;cnt_pass=0;outfile=fname;};
		AMonitor(int iturn, double ratio,const string& fname, const string& fam="undefined") : AccelElement(0){flag_print=0; this->len = 0;pfs=NULL; famname=fam;type="Monitor"; ratio_print=ratio; interval_turn=iturn;cnt_pass=0;outfile=fname;};
		AMonitor(int iturn, double ratio,unsigned int flag, const string& fname, const string& fam="undefined") : AccelElement(0){flag_print=flag; this->len = 0;pfs=NULL; famname=fam;type="Monitor"; ratio_print=ratio; interval_turn=iturn;cnt_pass=0;outfile=fname;};
		//~AMonitor(){if(pfs!=NULL) ;/*(*pfs).close();*/ };
		AMonitor& operator=(AMonitor& a){len=a.len; famname=a.famname; pfs=a.pfs; type=a.type; flag_print=a.flag_print;
		ratio_print=a.ratio_print; interval_turn=a.interval_turn;cnt_pass=a.cnt_pass;outfile=a.outfile;
		return *this;};
		unsigned int set_print_on(unsigned int flag) {flag_print = flag; return flag_print;};
		void serialize(ostream& os);
		//void setstream(ofstream& f) {fout = f;};
		void setpassnumber(int num=0) {cnt_pass = num;};
		
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		AMonitor* create(istream& is, string& famname);
		void pass(double *r, int N);
};

class ABend : public AccelElement {
	protected:
		double irho;
		double K;
		double theta;
		double e1,e2;
	public:
		ABend(){len=0; K=0; famname="undefined";type="Bend";};
		ABend(double len, double theta, double e1, double e2, double K=0,const string& fam="undefined") : AccelElement(len)
			{this->theta = theta; this->e1=e1; this->e2=e2; this->K=K;
				irho = theta/len;
				famname=fam;type="Bend";
			};
		ABend& operator=(ABend& a){len = a.len; K=a.K;irho=a.irho;theta=a.theta;e1=a.e1;e2=a.e2; famname=a.famname; type=a.type; return *this;};
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		ABend* create(istream& is, string& famname);
		void pass(double *r, int N);
};

 
class ABendMPole : public AccelElement {
	private:
		unsigned int flag_rad_on; //0 for no radiation, 1 for radiation damping, 2 for radiation damping and quantum excitation
		
	protected:
		double irho;
		double theta;
		double e1,e2;
		double polyB[MAX_MPOLE];
		int    nslice;
		double energy;
		int    MaxNPole;
	public:
		ABendMPole(){len=0; irho=0; theta=0;e1=e2=0;nslice=10;energy=DEFAULTE0;
			famname="undefined";type="BendMPole"; flag_rad_on = 0; for(int i=0;i<MAX_MPOLE;i++) this->polyB[i]=0; MaxNPole=MAX_MPOLE;};
		ABendMPole(double len, double theta, double e1, double e2, double *polyB=NULL,int MaxNPole=MAX_MPOLE,const string& fam="undefined", 
			double en=DEFAULTE0) : AccelElement(len)
			{this->theta = theta; this->e1=e1; this->e2=e2; this->MaxNPole=MaxNPole;
				irho = theta/len;
				nslice = 10;
				energy = en;
				flag_rad_on = 0;
				
				if(polyB==NULL)
				{	for(int i=0;i<MAX_MPOLE;i++) this->polyB[i]=0; 
				}
				else
				{	for(int i=0;i<MaxNPole;i++) this->polyB[i]=polyB[i];
					for(int i=MaxNPole;i<MAX_MPOLE;i++) this->polyB[i]=0;
				}
				famname=fam;type="BendMPole";
			};
		
		double getenergy(){return energy;};
		void   setenergy(double en) {this->energy = en;};	
		unsigned int set_rad_onoff(unsigned int flag) 
		{ if(flag<=2 && flag>=0) flag_rad_on = flag;
			else flag_rad_on = 0; //off
			return flag_rad_on; };
		ABendMPole& operator=(const ABendMPole& a){len = a.len; nslice=a.nslice; irho=a.irho;theta=a.theta;e1=a.e1;e2=a.e2; famname=a.famname; type=a.type;
			flag_rad_on = a.flag_rad_on; MaxNPole=a.MaxNPole;
			for(int i=0;i<MAX_MPOLE;i++) this->polyB[i]=a.polyB[i]; 
			return *this;};
		//void setfield(void *p, const char* pname);
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		void setfield(void *p,const string pname, unsigned int index);
		void getfield(void *p,const string pname, unsigned int index);
		ABendMPole* create(istream& is, string& famname);
		void pass(double *r, int N);
};

class AStrMPole : public AccelElement {
	protected:
		double polyA[MAX_MPOLE];
		double polyB[MAX_MPOLE];
		int    nslice;
		int    MaxNPole;
	public:
		AStrMPole(){nslice=10; MaxNPole=MAX_MPOLE; for(int i=0;i<MAX_MPOLE;i++) this->polyB[i]=0; 
			for(int i=0;i<MAX_MPOLE;i++) this->polyA[i]=0; famname="undefined";type="StrMPole"; };
		AStrMPole(double len, double *polyB=NULL, double *polyA=NULL,int MaxNPole=MAX_MPOLE, const string& fam="undefined") : AccelElement(len)
			{
				nslice = 10;
				this->MaxNPole = MaxNPole;
				if(polyB==NULL)
				{	for(int i=0;i<MAX_MPOLE;i++) this->polyB[i]=0; 
				}
				else
				{	for(int i=0;i<MaxNPole;i++) this->polyB[i]=polyB[i];
					for(int i=MaxNPole;i<MAX_MPOLE;i++) this->polyB[i]=0;
					}
					
				if(polyA==NULL)
				{	for(int i=0;i<MAX_MPOLE;i++) this->polyA[i]=0; 
				}
				else
				{ for(int i=0;i<MAX_MPOLE;i++) this->polyA[i]=polyA[i];
					for(int i=MaxNPole;i<MAX_MPOLE;i++) this->polyA[i]=0;	}
				famname=fam;type="StrMPole";
					
			};
		
		AStrMPole& operator=(AStrMPole& a){len = a.len; nslice=a.nslice;  famname=a.famname; type=a.type; MaxNPole=a.MaxNPole;
			for(int i=0;i<MAX_MPOLE;i++) 
			{this->polyB[i]=a.polyB[i]; this->polyA[i]=a.polyA[i];}
			return *this; 
			};	
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		void setfield(void *p,const string pname, unsigned int index);
		void getfield(void *p,const string pname, unsigned int index);
		
		AStrMPole* create(istream& is, string& famname);
		void pass(double *r, int N);
};

class ACavity : public AccelElement {
	protected:
		double voltage;  //MV
		double freq;    //Hz
		double phi_s;  //radian 
		double energy;  //MeV
		bool  flag_onoff;
	public:
		ACavity() {len=0; famname="undefined";type="Cavity";voltage=0; freq=476.314e6; phi_s=0;energy=DEFAULTE0; flag_onoff=true;};
		ACavity(double len, double volt, double freq, double phis, double en=DEFAULTE0, const string& fam="undefined") 
		{this->len=len; famname=fam; type="Cavity";this->voltage=volt; this->freq=freq; this->phi_s=phis;this->energy = en;flag_onoff=true;};
		
		ACavity& operator=(ACavity& a){len = a.len; famname=a.famname; type=a.type;
			voltage = a.voltage; freq = a.freq; phi_s = a.phi_s; energy = a.energy;
			flag_onoff=a.flag_onoff;
			return *this;
			};	
		double getenergy(){return energy;};
		double getvoltage(){return voltage;};
		void   setenergy(double en) {this->energy = en;};	
		bool set_cavity_onoff(bool flag){bool tmp = flag_onoff; flag_onoff = flag; return tmp;};
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		ACavity* create(istream& is, string& famname);
		void pass(double *r, int N);
};

//Crab cavity was added Dec. 1, 2016
class ACrabCav : public AccelElement {
	protected:
		double devoltx;  //MV, horizontal deflecting voltage
		double devolty;  //MV, vertical deflecting voltage
		double freq;    //Hz
		double phi_s;  //radian 
		double energy;  //MeV
		
		double sigphi; //rad, phase noise sigma
		double sigdVV; //, dV/V noise sigma
		bool  flag_onoff;
		//NormalGenerator gn;
	public:
		ACrabCav() {len=0; famname="undefined";type="CrabCav"; freq=476.314e6*6; phi_s=0;energy=DEFAULTE0;
		devoltx = 0; devolty = 0; sigphi = 0; sigdVV = 0; flag_onoff = 0; 
		};
		ACrabCav(double len, double devoltx, double devolty, double freq, double phis, double en=DEFAULTE0, double sigphi=0, double sigdVV=0, const string& fam="undefined") 
		{this->len=len; famname=fam; type="CrabCav";this->devoltx=devoltx; 
			this->devolty=devolty;this->freq=freq; this->phi_s=phis;this->energy = en;
			this->sigphi=sigphi; this->sigdVV=sigdVV; flag_onoff = 0;
		};
		
		ACrabCav& operator=(ACrabCav& a){len = a.len; famname=a.famname; type=a.type; devoltx = a.devoltx; 
			devolty = a.devolty; freq = a.freq; phi_s = a.phi_s; energy = a.energy; a.sigphi=sigphi; a.sigdVV=sigdVV;
			flag_onoff=a.flag_onoff;
			return *this;
			};	
		double getenergy(){return energy;};
		void   setenergy(double en) {this->energy = en;};
		bool set_crab_onoff(bool flag){bool tmp = flag_onoff; flag_onoff = flag; return tmp;};
		double getdevoltx(){return devoltx;};
		double getdevolty(){return devolty;};
		void   setnoisesigma(double sigphi=0, double sigdVV=0) {this->sigphi=sigphi; this->sigdVV=sigdVV;};	
			
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		ACrabCav* create(istream& is, string& famname);
		void pass(double *r, int N);
};

//RF multipole was added Dec. 1, 2016
class ARFMult : public AccelElement {
	protected:
		double freq;    //Hz
		double phi_s;  //radian 
		double energy;  //MeV
		
		double polyA[MAX_MPOLE];  //Integrated multipole, even if len is nonzero!
		double polyB[MAX_MPOLE];  //Integrated multipole, even if len is nonzero!
		int MaxNPole;
		double sigphi; //rad, phase noise sigma
		double sigdVV; //, dV/V noise sigma
		
		//NormalGenerator gn;
	public:
		ARFMult() {len=0; famname="undefined";type="RFMult"; freq=476.314e6*6; phi_s=0;energy=DEFAULTE0;
			for(int i=0;i<MAX_MPOLE;i++) { this->polyB[i]=0;this->polyA[i]=0; }
			 sigphi=0; sigdVV=0;MaxNPole=0;};
		ARFMult(double len, double freq, double phis, double en=DEFAULTE0, double *polyB=NULL, double *polyA=NULL,int MaxNPole=0, const string& fam="undefined") 
		{this->len=len; famname=fam; type="RFMult";this->freq=freq; this->phi_s=phis;this->energy = en;
			for(int i=0;i<MaxNPole;i++) 
			{
				this->polyB[i]= (polyB==NULL)?0:polyB[i];
				this->polyA[i]= (polyA==NULL)?0:polyA[i];
			}
			this->sigphi=0; this->sigdVV=0;
			this->MaxNPole=MaxNPole;
			
			};
		ARFMult& operator=(ARFMult& a){len = a.len; famname=a.famname; type=a.type;  
			freq = a.freq; phi_s = a.phi_s; energy = a.energy; a.sigphi=sigphi; a.sigdVV=sigdVV;
			MaxNPole = a.MaxNPole;
			for(int i=0;i<MaxNPole;i++) 
			{this->polyB[i]=a.polyB[i]; this->polyA[i]=a.polyA[i];}
			return *this;
			};	
		double getenergy(){return energy;};
		void   setenergy(double en) {this->energy = en;};
		
		void   setnoisesigma(double sigphi=0, double sigdVV=0) {this->sigphi=sigphi; this->sigdVV=sigdVV;};	
			
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		void setfield(void *p,const string pname, unsigned int index);
		void getfield(void *p,const string pname, unsigned int index);
		ARFMult* create(istream& is, string& famname);
		void pass(double *r, int N);
};

#define MAX_WIGTAB_X 101
#define MAX_WIGTAB_Y 21
class AWigTable : public AccelElement{
	protected:
		int    nslice;
		unsigned int mx, my;
		double xl[MAX_WIGTAB_X], yl[MAX_WIGTAB_X];
		double *xkick[MAX_WIGTAB_Y];
		double *ykick[MAX_WIGTAB_Y];
		
	public:
		
		AWigTable() {len=0; famname="undefined";type="WigTable";nslice=1; mx=0; my=0; };
		AWigTable(double len, int nslice, unsigned int mx, unsigned int my, const string& fam="undefined");
		~AWigTable() { for(int i=0;i<my;i++) { delete []xkick[i];  delete []ykick[i]; }};
		
		AWigTable& operator=(AWigTable& a);
		
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		AWigTable* create(istream& is, string& famname);
		void pass(double *r, int N);
			
};

#define IDEALWIG_TYPE_H  0   //horizontal wiggler, vertical field, regular
#define IDEALWIG_TYPE_V  1   //vertical wiggler, horizontal field, generate vertical dispersion
class AIdealWig : public AccelElement{
	private:
		unsigned int flag_rad_on;
		unsigned int flag_wigtype; 
		
	protected:
		int    nslice;
		double lambda,h;  //h = 1/rho
		double energy;
				
	public:
		AIdealWig() {len=0; famname="undefined";type="IdealWig";nslice=1; lambda=0.10; h=0; flag_rad_on = 0;energy=DEFAULTE0;flag_wigtype=IDEALWIG_TYPE_H; };
		AIdealWig(double len, int nslice, double lambda, double h, unsigned int wig_type=IDEALWIG_TYPE_H,const string& fam="undefined");
//		~AIdealWig() { };
		
		AIdealWig& operator=(AIdealWig& a);
		
		double getenergy(){return energy;};
		void   setenergy(double en) {this->energy = en;};	
		unsigned int set_rad_onoff(unsigned int flag) 
		{
			if(flag<=2 && flag>=0) flag_rad_on = flag;
			else flag_rad_on = 0; //off
			return flag_rad_on; 
		};
		
		void serialize(ostream& os);
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		AIdealWig* create(istream& is, string& famname);
		void pass(double *r, int N);
			
};

class ATransWake : public AccelElement{
	protected:
		double energy;
		double Ne;  //charge in nC
		double Wx, Wy;  //Wake amplitude in V/pC/m
		double fRx, fRy; //frequency in Hz
		double Qx, Qy;   //Q-value
		bool flag_onoff;

	public:
		ATransWake() { len = 0; famname = "undefined"; type = "TransWake"; Ne = 0.0; Wx = 0; Wy = 0; fRx = 1.0e10; fRy = 1.0e10; Qx = 1; Qy = 1; energy = DEFAULTE0; flag_onoff = 0; };
		ATransWake(double len, double Ne, double Wx, double fRx, double Qx, double Wy, double fRy, double Qy, double en=DEFAULTE0, const string& fam="undefined"){
			this->len = len; this->Ne = Ne; this->Wx=Wx; this->fRx=fRx; this->Qx=Qx; this->energy=en;
			this->Wy=Wy; this->fRy=fRy; this->Qy=Qy;
			famname=fam;type="TransWake"; flag_onoff = 0;
		};
	
		ATransWake& operator=(ATransWake& a){len = a.len; famname=a.famname; type=a.type;  
			Ne = a.Ne; Wx=a.Wx; fRx=a.fRx; Qx=a.Qx; Wy=a.Wy; fRy=a.fRy; Qy=a.Qy; flag_onoff = a.flag_onoff;
			return *this;
			};	
		bool set_wake_onoff(bool flag) { bool tmp = flag_onoff; flag_onoff = flag; return tmp; };
		void serialize(ostream& os);
		double getenergy(){return energy;};
		void   setenergy(double en) {this->energy = en;};	
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		ATransWake* create(istream& is, string& famname);
		void pass(double *r, int N);
			
};

class AMatrix : public AccelElement{
	protected:
		double mat[6][6];
		

	public:
		AMatrix() { len = 0; famname = "undefined"; type = "Matrix"; for(int i=0;i<6;i++) 
				{ for(int j=0;j<6;j++) mat[i][j] = 0; 
					mat[i][i] = 1; 
				} 
			};
		AMatrix(double len, double mat[][6], const string& fam="undefined"){
			this->len = len; 
			for(int i=0;i<6;i++) 
				{ for(int j=0;j<6;j++) this->mat[i][j] = mat[i][j]; 
				} 
			famname=fam;type="Matrix"; 
		};
	
		AMatrix& operator=(AMatrix& a){len = a.len; famname=a.famname; type=a.type;  
			for(int i=0;i<6;i++) 
				{ for(int j=0;j<6;j++) this->mat[i][j] = a.mat[i][j]; 
				} 
			return *this;
			};	
		void serialize(ostream& os);
		AMatrix* create(istream& is, string& famname);
		void setfield(void *p, int i, int j);
		void getfield(void *p, int i, int j);
		void pass(double *r, int N);
			
};

class ASyncRad : public AccelElement{
	protected:
		double energy;  //MeV
		double emitx, emity, sigdpp;  //equilibrium emittance in meter*rad for x and y planes, and sigma_\delta
		double U0;  //energy loss per turn in MeV
		double Jx, Jy, Jz; //damping partition, Jz=4-Jx-Jy
		double betax, betay;   //beta function in x an y at the location
		bool flag_damping_onoff, flag_excitation_onoff;

	public:
		ASyncRad() { len = 0; famname = "undefined"; type = "SyncRad"; 
			energy = DEFAULTE0; 
			emitx = 1e-9; emity=10e-12; sigdpp=0.001;
			U0 = 0.881; Jx=1.35; Jy=1; Jz=4-Jx-Jy;
			betax = 1; betay=1;
			flag_damping_onoff = 0; //off
			flag_excitation_onoff = 0; //off
			};
		ASyncRad(double len, double emitx, double emity, double sigdpp, double U0, double Jx, double Jy, double betax, double betay, double en=DEFAULTE0, const string& fam="undefined"){
			this->len = len; this->emitx = emitx; this->emity = emity; this->sigdpp = sigdpp; 
			this->U0 = U0; this->Jx = Jx;this->Jy = Jy; this->betax = betax;this->betay = betay; this->energy = en;
			flag_damping_onoff = 0; //off
			flag_excitation_onoff = 0; //off
			Jz=4-Jx-Jy;
			famname=fam;type="SyncRad"; 
		};
	
		ASyncRad& operator=(ASyncRad& a){len = a.len; famname=a.famname; type=a.type;  
			emitx = a.emitx; emity = a.emity; sigdpp = a.sigdpp; 
			U0 = a.U0; Jx = a.Jx;Jy = a.Jy; betax = a.betax;betay = a.betay; energy = a.energy;
			Jz = a.Jz;  flag_damping_onoff=a.flag_damping_onoff;  flag_excitation_onoff=a.flag_excitation_onoff; 
			return *this;
			};	
			
		void serialize(ostream& os);
		ASyncRad* create(istream& is, string& famname);
		bool set_rad_onoff(bool flag) { flag_damping_onoff = flag; flag_excitation_onoff = flag; return 0; };
		bool set_rad_onoff(bool flag_damp,bool flag_QE) { flag_damping_onoff = flag_damp; flag_excitation_onoff = flag_QE; return 0; };
		void setfield(void *p, const string pname);
		void getfield(void *p, const string pname);
		void pass(double *r, int N);
			
};


struct Twiss;
class Mat;      // matrix class
class Vec;      // vector class

#define FLAG_RAD_ONOFF_OFF       0
#define FLAG_RAD_ONOFF_RAD       1
#define FLAG_RAD_ONOFF_RADQE     2
#define FLAG_RAD_ONOFF_LUMP      3 //Use lump element SyncRad, both damping and excitation on
#define FLAG_RAD_ONOFF_LUMP_DAMP 4 //Use lump element SyncRad, damping on, excitation off
#define FLAG_RAD_ONOFF_LUMP_QE   5 //Use lump element SyncRad, damping off, excitation on
#define FLAG_CAVITY_OFF false
#define FLAG_CAVITY_ON true
#define MAX_ELEMENT 100000
#define DELTA_DIFF_TRACK (1.0E-8)
class ALine {
	AccelElement *pElement[MAX_ELEMENT];
	unsigned int nElement;
	
	public:
	ALine(){ nElement = 0;};
	ALine(ALine &line);
	~ALine(){for(int i=0;i<nElement;i++) delete pElement[i];};
	
	void copyline(ALine &line);		
	int addelement(AccelElement*);
	int deleteelement(unsigned int);
	int insertelement(unsigned int i, AccelElement* pnew);
	AccelElement* getelement(unsigned int i);
	unsigned int findelements(unsigned int *list, string type, string famname, unsigned int index=0); //index = 1 for the 1st element, index=0 for all elements
	unsigned int findelements(unsigned int *list, string type);
	unsigned int resetfirstelement(string type, string famname, unsigned int index=1); //move the matched element to the beginning of the ring
	unsigned int resetfirstelement(unsigned int eindex); //move the element at eindex to the beginning of the ring
	
	void serialize(ostream& os);
	static ALine* deserialize(istream& is);
	unsigned int getnumber(){return nElement;};
	double getlength();
	double getlength(unsigned int ed);
	double getlength(unsigned int st, unsigned int ed);
	
	int setenergy(double en);
	double setrfphase(double phis);
	double setrfphase();
	int set_rad_onoff(unsigned int flag);
	int set_wake_onoff(bool flag);
	int set_cavity_onoff(bool flag);
	int set_crab_onoff(bool flag);
	int set_monitor(unsigned int flag, unsigned int index=0); //index=0 for all monitors, index>0 for one monitor
	int getfield(void *fvalue, unsigned nfam, unsigned *lfam,string pname,unsigned index=0);
	int setfield(void *fvalue, unsigned nfam, unsigned *lfam,string pname,unsigned index=0);
	int set_monitor_pass(int num=0);
	
	void pass(double *r, int N,unsigned st, unsigned ed);
	void pass(double *r, int N);
	void pass(double *r, int N, int nturn);
	double trackDA(double xmax, double xmin, double ymax, int nturn, unsigned int nray, int nstep, double *rx, double *ry,bool flag_negweight=false);
	double trackDiffusion(int nturn, unsigned int np, double *x0, double *y0, double *nux, double *nuy, double *dif,double xlow, double xhi, double ylow, double yhi);
	double trackMA(double dppmax, double dppmin, int nturn, unsigned int nstep, double *dpp_plus, double *dpp_minus,unsigned nfam, unsigned *lfam);
	//bool setfield(string type, string famname, unsigned int index, void* p);  //index=ALL for all matched elements

	Vec& getco4(Vec& rc, double dp=0,double delta=DELTA_DIFF_TRACK);
	Vec& getco6(Vec& rc, double delta=DELTA_DIFF_TRACK);
	Mat& getm66(Mat& m66,double dp=0,double delta=DELTA_DIFF_TRACK);
	Mat& getm66(Mat& m66,unsigned int ed,double dp=0,double delta=DELTA_DIFF_TRACK);
	Mat& getm66(Mat& m66,unsigned int st, unsigned int ed,double dp=0,double delta=DELTA_DIFF_TRACK);

	void calctune(double &tunex, double &tuney, double dp=0,double delta=DELTA_DIFF_TRACK);
	void calctwiss(Twiss& tw,double dp=0,double delta=DELTA_DIFF_TRACK);;
	void calctwiss(Twiss *tw,unsigned int indx[], int n,double dp=0,double delta=DELTA_DIFF_TRACK);
	void calcchrom(double &chromx, double &chromy,double dp=0,double delta=DELTA_DIFF_TRACK);
	void calcchrom2(double &chromx, double &chromy,double &chromx2, double &chromy2);
	void correcttune(double &tunex, double &tuney, unsigned int nfam1, unsigned int nfam2,unsigned int *lfam1, unsigned int *lfam2);
	void correctchrom(double &chromx, double &chromy, unsigned int nfam1, unsigned int nfam2,unsigned int *lfam1, unsigned int *lfam2);
	void calctuneresp(double resp[][2], unsigned nknob, unsigned *nfam, unsigned **lfam);
	void calcchromresp(double resp[][2], unsigned nknob, unsigned *nfam, unsigned **lfam);
	
	void tracktune(double &tunex, double &tuney, unsigned int nturn,const double x0, const double y0,double xlow=0.0, double xhi=0.5,double ylow=0.0, double yhi=0.5);
	//void trackdtune(double &tunex, double &tuney, double &dtunex, double &dtuney,unsigned int nturn,const double x0, const double y0,double xlow=0.0, double xhi=0.5,double ylow=0.0, double yhi=0.5);
	void trackdnudexy(double *dnudexy, unsigned int nturn=512,const double xm=0.004, const double ym=0.003);
	
	static unsigned int flag_msg;
	unsigned int set_flag_msg(unsigned int flag) {flag_msg |= flag; return flag_msg;};
	unsigned int clear_flag_msg(unsigned int flag) {flag_msg &= ~flag; return flag_msg;};	
	
};
#define FLAG_MSG_DA 0x1
#define FLAG_MSG_MA (0x1<<1)
#define FLAG_MSG_DF (0x1<<2)
#define FLAG_MSG_PASS (0x1<<3)
#define FLAG_MSG_CHROM_CORR (0x1<<4)

class ARing : public ALine {
	protected:
		
};

struct Twiss {
	double s;
	double psix, psiy;
	double betx, alfx;
	double bety, alfy;
	
	void serialize(ostream& os);
};

#endif 
