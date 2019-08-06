#ifndef __ATPASS_H
#define __ATPASS_H

#define SQR(X) ((X)*(X))

void drift(double *r, double L);
void quad6 (double *r, double L, double K);
void bend6(double* r, double L, double b_angle, double grd, double ByError);
void edge(double* r, double inv_rho, double edge_angle);

void bndthinkick(double* r, double* B, double L, double irho, int max_order);
void bndthinkickrad(double* r, double* B, double L, double irho, double E0,int max_order, unsigned int flag_rad=1);
void ATbendhxdrift6(double* r, double L,double h);
void fastdrift(double* r, double NormL);
void strthinkick(double* r, double* A, double* B, double L, int max_order);

#endif
