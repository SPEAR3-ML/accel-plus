#include <math.h>
#include "atpass.h"
#include "utility.h"
#include <iostream>  //for debugging

inline void drift(double *r6, double L)
{
	double p_norm = 1.0/(1+r6[5]);
		r6[0] += r6[1]*L;
		r6[2] += r6[3]*L;
		r6[4] -= SQR(p_norm)*L*(SQR(r6[1]+SQR(r6[2])))/2.0;
		
}

void quad6 (double *r, double L, double K)
{	/* K - is the quadrupole strength defined as
	   (e/Eo)(dBz/dx) [1/m^2] 
	   another notation: g0 [DESY paper]
	*/
	
	double p_norm = 1/(1+r[5]);
	double x, xpr, y ,ypr, g, t ,lt;
	double M12,M21,M34,M43,MVD,MHD;  /* non-0 elements of transfer matrix */
		
	if(K==0) /* Track as a drift */
	{	drift(r,L);
		return;
	}   
	
	x   = r[0];
	xpr = r[1]*p_norm;
	y   = r[2];
	ypr = r[3]*p_norm;
	
	g  = fabs(K)/(1+r[5]);
	t  = sqrt(g);
	lt = L*t;
				
   if(K>0)
		{	/* Horizontal  */
				MHD = cos(lt); 		
				M12 = sin(lt)/t;
				M21 = -M12*g;		
   			/* Vertical */
				MVD = cosh(lt);		
				M34 = sinh(lt)/t;
				M43 = M34*g;			
		}
	else 
		{	/* Horizontal  */
				MHD = cosh(lt);		
				M12 = sinh(lt)/t;
				M21 = M12*g;			
			/* Vertical */
				MVD = cos(lt); 		
				M34 = sin(lt)/t;
				M43 = -M34*g;			
		}			
	

	/* M transformation compute change in tau first */
	r[0]=  MHD*x + M12*xpr;
	r[1]= (M21*x + MHD*xpr)/p_norm; 
    r[2]=  MVD*y + M34*ypr;
	r[3]= (M43*y + MVD*ypr)/p_norm;  
  
    /* no change in r[4] (delta) */
	
	r[4]-= g*(x*x*(L-MHD*M12)-y*y*(L-MVD*M34))/4;
	r[4]-= (xpr*xpr*(L+MHD*M12)+ypr*ypr*(L+MVD*M34))/4;
	r[4]-= (x*xpr*M12*M21 + y*ypr*M34*M43)/2;
}


/**********  BEND *************************/
void bend6(double* r, double L, double b_angle, double grd, double ByError)
{	
	double M12,M21,M34,M43,MVD,MHD;  /*  non-0 elements of transfer matrix */
	double x, xpr,y,ypr,delta;
	
	
	
	double sqrtG1, sqrtG2, arg1, arg2;
	double Kx = b_angle/L; /* Curvature of the design trajectory */
	double p_norm = 1/(1+r[5]);
	double G1 = (Kx*Kx+grd)*p_norm;
	double G2 = -grd*p_norm;

	
	/* Horizontal Transverse Matrix */
	if(G1==0)		/* Special case Kx^2 + grd = 0 */
 
		{	MHD = 1;				
			M12 = L;
			M21 = 0;	
		}
	else		
		{	if(G1 > 0)
				{	sqrtG1 = sqrt(G1);
					arg1 = L*sqrtG1;
					MHD = cos(arg1);				
					M12 = sin(arg1)/sqrtG1;
					M21 = -sin(arg1)*sqrtG1;	
				}
			else
				{	sqrtG1 = sqrt(-G1);
					arg1 = L*sqrtG1;
					MHD = cosh(arg1); 				
					M12 = sinh(arg1)/sqrtG1;
					M21 = sinh(arg1)*sqrtG1;	
				}
		}

	

	/*  Vertical Transverse Matrix */
	
	if(G2==0) /*  No gradient - vertical motion is a drift  */
	
		{	MVD = 1;				
			M34 = L;
			M43 = 0;	
		}
	else	
		{	if(G2 > 0)	/* Vertical focusing */
				{	sqrtG2 = sqrt(G2);
					arg2 = L*sqrtG2;
					MVD = cos(arg2);;				
					M34 = sin(arg2)/sqrtG2;
					M43 = -sin(arg2)*sqrtG2;	
				}
			else		/*  Vertical defocusing	*/
				{	sqrtG2 = sqrt(-G2);
					arg2 = L*sqrtG2;
					MVD = cosh(arg2); 				
					M34 = sinh(arg2)/sqrtG2;
					M43 = sinh(arg2)*sqrtG2;	
				}
		}

	x   = r[0];
	xpr = r[1]*p_norm;
	y   = r[2];
	ypr = r[3]*p_norm;
	delta = r[5]; 

	r[0]=  MHD*x + M12*xpr ;
	r[1]= (M21*x + MHD*xpr)/p_norm; 
	
	if(G1==0)	
		{	r[0]+= (delta*p_norm-ByError)*L*L*Kx/2;
			r[1]+= (delta*p_norm-ByError)*L*Kx/p_norm ;
		}
	else
		{	if(G1>0)
				{	r[0]+= (delta*p_norm-ByError)*(1-cos(arg1))*Kx/G1;
					r[1]+= (delta*p_norm-ByError)*sin(arg1)*Kx/(sqrtG1*p_norm) ;
				}	
			else
				{	r[0]+= (delta*p_norm-ByError)*(1-cosh(arg1))*Kx/G1;
					r[1]+= (delta*p_norm-ByError)*sinh(arg1)*Kx/(sqrtG1*p_norm) ;
				}
		}

	r[2]=  MVD*y + M34*ypr;
	r[3]= (M43*y + MVD*ypr)/p_norm ;
		
	r[4]-= xpr*xpr*(L+MHD*M12)/4;
   	
	if (G1==0) {
	    /* Do nothing */   
	} else {
    	r[4]-= (L-MHD*M12)*(x*x*G1+(delta*p_norm-ByError)*(delta*p_norm-ByError)*Kx*Kx/G1 
										-2*x*Kx*(delta*p_norm-ByError))/4;

    	r[4]-= M12*M21*( x*xpr - xpr*(delta*p_norm-ByError)*Kx/G1)/2;
    	
    	r[4]-= Kx*x*M12  +   xpr*(1-MHD)*Kx/G1   +   (delta*p_norm-ByError)*(L-M12)*Kx*Kx/G1;
    }



	r[4]-= ((L-MVD*M34)*y*y*G2 + ypr*ypr*(L+MVD*M34))/4;
   	
	r[4]-= M34*M43*y*ypr/2;	

}

/***********************/
void edge(double* r, double inv_rho, double edge_angle)
{	/* Edge focusing in dipoles with hard-edge field */
    double psi = inv_rho*tan(edge_angle);

    r[1]+=r[0]*psi;
	r[3]-=r[2]*psi;
}

/*************************************/
void bndthinkick(double* r, double* B, double L, double irho, int max_order)
/*****************************************************************************
(1) The vector potential is expanded up to 4th order of x and y. 
(2) Coefficients in PolynomB higher than 4th order is treated as if they are on straight geometry.
(3) The Hamiltonian is H2 = - h x delta - (1+h x)As/Brho-B0 x/Brho      
*/
{ int i;
	double ReSum = 0; /*B[max_order];*/
 	double ImSum = 0; /*A[max_order];*/

	double ReSumTemp;
	double K1,K2,K3,h;
 
	K1 = B[1]; 
	if (max_order>=2)
	 K2=B[2]; 
	else
	 K2=0;
	 
	if (max_order>=3)
	 K3=B[3];
	else
	 K3=0;
	 h=irho;

  ReSum = B[max_order];
	for(i=max_order-1;i>=0;i--)
		{	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
			ImSum = ImSum*r[0] +  ReSum*r[2] ;
			ReSum = ReSumTemp;
		}
	
	r[1] -=  L*(-h*r[5] + ReSum + h*(h*r[0]+K1*(r[0]*r[0]-0.5*r[2]*r[2])+K2*(r[0]*r[0]*r[0]-4.0/3.0*r[0]*r[2]*r[2]))    ); 
	r[3] +=  L*(ImSum+h*(K1*r[0]*r[2]+4.0/3.0*K2*r[0]*r[0]*r[2]+(h/6.0*K1-K2/3.0)*r[2]*r[2]*r[2])) ;
	r[4] -=  L*h*r[0]; /* pathlength */

}

double B2perp(double bx, double by, double irho, 
                            double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|) , where e is a unit vector in the direction of velocity  */
    
{	double v_norm2;
	v_norm2 = 1/(SQR(1+x*irho)+ SQR(xpr) + SQR(ypr));

	/* components of the  velocity vector
	   double ex, ey, ez;
	    ex = xpr; 
	    ey = ypr; 
	    ez = (1+x*irho);
	*/
  	
	return((SQR(by*(1+x*irho)) + SQR(bx*(1+x*irho)) + SQR(bx*ypr - by*xpr) )*v_norm2) ;

} 

void bndthinkickrad(double* r, double* B, double L, double irho, double E0,int max_order, unsigned int flag_rad)
/*****************************************************************************
(1) The vector potential is expanded up to 4th order of x and y. 
(2) Coefficients in PolynomB higher than 4th order is treated as if they are on straight geometry.
*/
{ int i;
	double ReSum = 0; /*B[max_order];*/
 	double ImSum = 0; /*A[max_order];*/

	double ReSumTemp;
	double K1,K2,K3,h;
	double x ,xpr, y, ypr, p_norm,dp_0, B2P;

#define PI  3.1415926535897932384626
#define TWOPI  (2*PI) 
#define CGAMMA 	8.846056192e-05 
	
	double CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e9);	/* [m]/[GeV^3] M.Sands (4.1) */
	double uc_E0 = 0.665*E0*E0/1.0e6/2.99792458E8*fabs(irho)*1000; /*=uc/E0, with uc being the nominal critical photon energy   */
	double delta_dpp,mean_dpp, sigma_dpp;
	
	K1 = B[1]; 
	if (max_order>=2)
	 K2=B[2]; 
	else
	 K2=0;
	 
	if (max_order>=3)
	 K3=B[3];
	else
	 K3=0;
	 h=irho;

  ReSum = B[max_order];
	for(i=max_order-1;i>=0;i--)
		{	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
			ImSum = ImSum*r[0] +  ReSum*r[2] ;
			ReSum = ReSumTemp;
		}
	
		/* calculate angles from momentums 	*/
	p_norm = 1/(1+r[5]);
	x   = r[0];
	xpr = r[1]*p_norm;
	y   = r[2];
	ypr = r[3]*p_norm;

	/* see Iselin Part. Accel. 1985  */
	ImSum += h*(K1*h-K2)*y*y*y/6.0;
	ReSum += -K1*h*y*y/2.0 + h*(K1*h-K2)*x*y*y/2.0;

	B2P = B2perp(ImSum, ReSum +irho, irho, x , xpr, y ,ypr);
	//cout<<B2P<<endl;

	dp_0 = r[5];
	if(flag_rad==1)
	{
		r[5] = r[5] - CRAD*SQR(1+r[5])*B2P*(1 + x*irho + (SQR(xpr)+SQR(ypr))/2 )*L;
	}
	else if(flag_rad==2)
	{
		
		mean_dpp = CRAD*SQR(1+r[5])*B2P*(1 + x*irho + (SQR(xpr)+SQR(ypr))/2 )*L;
		if (mean_dpp>=0)
			sigma_dpp = sqrt(55*sqrt(3)/72.0*uc_E0*mean_dpp);
		else
			sigma_dpp = sqrt(-55*sqrt(3)/72.0*uc_E0*mean_dpp);
		
		double gr = grandn();
		delta_dpp = mean_dpp+gr*sigma_dpp;
		//cout<<uc_E0<<" "<<mean_dpp<<" "<<sigma_dpp<<" "<<delta_dpp<<" "<<gr<<" "<<r[5]<<endl;
		r[5] = r[5] - delta_dpp;
	}

	/* recalculate momentums from angles after losing energy for radiation 	*/
	p_norm = 1/(1+r[5]);
	r[1] = xpr/p_norm;
	r[3] = ypr/p_norm;

	
	r[1] -=  L*(-h*r[5] + ReSum + h*(h*r[0]+K1*(r[0]*r[0]-0.5*r[2]*r[2])+K2*(r[0]*r[0]*r[0]-4.0/3.0*r[0]*r[2]*r[2]))    ); 
	r[3] +=  L*(ImSum+h*(K1*r[0]*r[2]+4.0/3.0*K2*r[0]*r[0]*r[2]+(h/6.0*K1-K2/3.0)*r[2]*r[2]*r[2])) ;
	r[4] -=  L*h*r[0]; /* pathlength */
	
}


/* the pseudo-drift element described by Hamiltonian H1 = (1+hx) (px^2+py^2)/2(1+delta),     */
void ATbendhxdrift6(double* r, double L,double h)
{
	double hs = h*L;
	double i1pd = 1.0/(1+r[5]);
	double x=r[0],px=r[1],y=r[2],py=r[3];

	r[0] += (1+h*x)*px*i1pd*L+1/4.*hs*L*(px*px-py*py)*i1pd*i1pd; /* (1.0/h+x)*((1.0+hs*px*i1pd/2.)*(1.0+hs*px*i1pd/2.)-(hs*py*i1pd/2.)*(hs*py*i1pd/2.))-1./h;*/
	r[1] -= hs*(px*px+py*py)*i1pd/2.0;
	
	r[2]+= (1.0+h*x)*i1pd*py*L*(1.+px*hs/2.0);
	r[4]-= (1.0+h*x)*i1pd*i1pd*L/2.0*(px*px+py*py);

}

/*************************************************  */
void strthinkick(double* r, double* A, double* B, double L, int max_order)
/***************************************************************************** 
Calculate and apply a multipole kick to a 6-dimentional
phase space vector in a straight element ( quadrupole)

IMPORTANT !!!
The reference coordinate system is straight but the field expansion may still
contain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
A[0], B[0] - C,C++ notation


   Note: in the US convention the transverse multipole field is written as:

                         max_order+1
                           ----
                           \                       n-1
	   (B + iB  )/ B rho  =  >   (ia  + b ) (x + iy)
         y    x            /       n    n
	                        ----
                          n=1
	is a polynomial in (x,y) with the highest order = MaxOrder
	

	Using different index notation 
   
                         max_order
                           ----
                           \                       n
	   (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
         y    x            /       n    n
	                       ----
                          n=0

	A,B: i=0 ... max_order
   [0] - dipole, [1] - quadrupole, [2] - sextupole ...
   units for A,B[i] = 1/[m]^(i+1)
	Coeficients are stroed in the PolynomA, PolynomB field of the element
	structure in MATLAB

	A[i] (C++,C) =  PolynomA(i+1) (MATLAB) 
	B[i] (C++,C) =  PolynomB(i+1) (MATLAB) 
	i = 0 .. MaxOrder

******************************************************************************/
{  int i;
	double ReSum = B[max_order];
 	double ImSum = A[max_order];
	double ReSumTemp;
    	for(i=max_order-1;i>=0;i--)
        {   ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
            ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
            ReSum = ReSumTemp;
        }

    r[1] -=  L*ReSum;
    r[3] +=  L*ImSum;
}

void fastdrift(double* r, double NormL)

/*   NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
     in the loop if momentum deviation (delta) does not change
     such as in 4-th order symplectic integrator w/o radiation
*/

{   double dx = NormL*r[1];
    double dy = NormL*r[3];
    r[0]+= dx;
    r[2]+= dy;
    r[4]-= NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[5]));
}

