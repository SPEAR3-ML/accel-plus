#ifndef __UTILITY_H
#define __UTILITY_H

//#include <random>
#include <stddef.h>
#include <string>
#include <cstdlib>


#define PI  3.1415926535897932384626
#define TWOPI  (2*PI)
#define MAX(A,B) ((A)>(B))?(A):(B)
#define MIN(A,B) ((A)<(B))?(A):(B)
#define SQR(X) ((X)*(X))
//using namespace std;

double naff(int N, double *data,double xlow=0, double xhi=0.5);

struct CLinfit {
	double a, b, siga, sigb, chi2, q;
};
void linearfit(CLinfit &ca,int N, double *xx,double  *yy, double *sig=NULL);
double interp1(double ytab[],double xl[], int nx, double x);
double interp2(double *ztab[], double xl[],int nx, double yl[],int ny,double x, double y);

double grandn(double gmean=0,double gstd=1.0);
//double* gen_distr_6D(int N, double bx, double emitx, double by, double emity, double sigz, double sigdpp);
double* gen_distr_6D(int N, double sigx, double sigxp, double sigy, double sigyp, double sigz, double sigdpp);

//std::string str_trim(const std::string& str, const std::string& whitespace = " \t");

// trim from left
inline std::string& str_trim_left(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& str_trim_right(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& str_trim_both(std::string& s, const char* t = " \t\n\r\f\v")
{
    return str_trim_left(str_trim_right(s, t), t);
}

inline int str2int(std::string s) //convert string to integer
{
    return atoi(s.c_str());
}

//vector<string> str_split(const string& str, const string& delim);
//make sure tokenlist has enough space for the tokens
int str_split2(std::string &s, const std::string& delimiter, std::string *tokenlist, int ntoken_max);
bool str_icompare(const std::string& a, const std::string& b); //case insensitive comparison
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime();
int calc_beam_stat( double *r, int N, double ar[], double sigr[], double *sigrr);

#endif
