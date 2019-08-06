#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <limits>
#include <assert.h>
#include "atpass.h"
#include "matvec.h"
#include "accel.h"
#include "utility.h"
#include <omp.h>

/***
Accel lattice functions
Created by Xiaobiao Huang, SLAC, Dec 2011
***/ 
/*******************************************************/
/* ALINE   *********************************************/
/*******************************************************/
ALine::ALine(ALine& line)
{
	nElement =0;
	copyline(line);
}
void ALine::copyline(ALine& line)
{
	string etype;
	//cout<<line.nElement<<endl;
	for(int i=0;i<line.nElement;i++)
	{
		etype = line.pElement[i]->gettype();
		//cout<<i<<"\t"<<etype<<endl;
		
		if( !etype.compare("Drift"))
		{
			ADrift *a = new ADrift;
			*a = *(ADrift*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("Quad"))
		{
			AQuad *a = new AQuad;
			*a = *(AQuad*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("Bend"))
		{
			ABend *a = new ABend;
			*a = *(ABend*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("BendMPole"))
		{
			ABendMPole *a = new ABendMPole;
			*a = *(ABendMPole*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("StrMPole"))
		{
			AStrMPole *a = new AStrMPole;
			*a = *(AStrMPole*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("Cavity"))
		{
			ACavity *a = new ACavity;
			*a = *(ACavity*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("WigTable"))
		{
			AWigTable *a = new AWigTable;
			*a = *(AWigTable*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("IdealWig"))
		{
			AIdealWig *a = new AIdealWig;
			*a = *(AIdealWig*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}	
		else if( !etype.compare("Aperture"))
		{
			AAperture *a = new AAperture;
			*a = *(AAperture*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("Corrector"))
		{
			ACorrector *a = new ACorrector;
			*a = *(ACorrector*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("Monitor"))
		{
			AMonitor *a = new AMonitor;
			*a = *(AMonitor*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("CrabCav"))
		{
			ACrabCav *a = new ACrabCav;
			*a = *(ACrabCav*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("RFMult"))
		{
			ARFMult *a = new ARFMult;
			*a = *(ARFMult*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}
		else if( !etype.compare("TransWake"))
		{
			ATransWake *a = new ATransWake;
			*a = *(ATransWake*)line.pElement[i];
			this->pElement[i] = (AccelElement*)a; 
			}	
		else
			cout<<"unknown element type: "<<line.pElement[i]->gettype()<<endl;
		
		this->nElement ++;
		
	}
		
}
double ALine::getlength(unsigned int st, unsigned int ed)
{
	assert(ed<=nElement);
	assert(st<=ed);
	
	double len = 0.;
	for(int i=st;i<ed;i++)
		len += pElement[i]->getlength();
	return len;
}
double ALine::getlength(unsigned int ed)
{
	assert(ed<=nElement);
	return getlength(0,ed);
}	
double ALine::getlength()
{
	return getlength(0,nElement);
}		
int ALine::addelement(AccelElement* pnew)
{
	pElement[nElement] = pnew;
	nElement += 1;
	return nElement;
}
int ALine::deleteelement(unsigned int i)
{
	if( i<nElement) 
	{
		for(int j=i;j<nElement-1;j++)
		{	pElement[j] = pElement[j+1];}
		nElement--;
	}
	else
		cout<<"warning(ALine::deleteelement): index out of range "<<endl;
	
	return nElement;
}
int ALine::insertelement(unsigned int i, AccelElement* pnew)
{
	if(nElement+1 > MAX_ELEMENT)
	{	cout<<"error: number of elements over flow"<<endl;
		return MAX_ELEMENT;
	}
		
	if( i<=nElement) 
	{
		for(int j=nElement;j>i;j--)
		{	pElement[j] = pElement[j-1];}
		pElement[i] = pnew;
		nElement++;
	}
	else
		cout<<"warning(ALine::insertelement): index out of range "<<endl;
	
	return nElement;
}
AccelElement* ALine::getelement(unsigned int i)
{
	if( i<=nElement) 
		return pElement[i];
	else
	{
		cout<<"warning(ALine::insertelement): index out of range "<<endl;
		return NULL;
	}
}	
void ALine::serialize(ostream& os)
{
	for(int i=0;i<nElement; i++)
	{
		pElement[i]->serialize(os);
	}
	
}

ALine* ALine::deserialize(istream& is)
{
	ALine* al = new ALine;
	
	al->nElement = 0;
	while(!is.eof())
	{
		AccelElement* np = AccelElement::deserialize(is);
			if(np!=NULL)
				al->pElement[al->nElement++] = np;
			else
				break;	
			//np->serialize(cout);
			
	}
	return al;
	
}
int ALine::setenergy(double en)
{
	int cnt=0;
	for(int i=0;i<nElement; i++)
	{
		if (pElement[i]->gettype()=="Cavity" ) 
		{	ACavity* p = (ACavity*)pElement[i];
			p->setenergy(en);
			cnt++;
		}
		if (pElement[i]->gettype()=="CrabCav" ) 
		{	ACrabCav* p = (ACrabCav*)pElement[i];
			p->setenergy(en);
			cnt++;
		}
		if (pElement[i]->gettype()=="RFMult" ) 
		{	ARFMult* p = (ARFMult*)pElement[i];
			p->setenergy(en);
			cnt++;
		}
		else if ( pElement[i]->gettype()=="BendMPole")
		{	ABendMPole* p = (ABendMPole*)pElement[i];
			p->setenergy(en);
			cnt++;
		}
		else if ( pElement[i]->gettype()=="IdealWig")
		{	AIdealWig* p = (AIdealWig*)pElement[i];
			p->setenergy(en);
			cnt++;
		} 
		else if ( pElement[i]->gettype()=="TransWake")
		{	ATransWake* p = (ATransWake*)pElement[i];
			p->setenergy(en);
			cnt++;
		}
	}
	return cnt;
}	
int ALine::set_rad_onoff(unsigned int flag)
{
	int cnt=0;
	for(int i=0;i<nElement; i++)
	{
		if ( pElement[i]->gettype()=="BendMPole")
		{	ABendMPole* p = (ABendMPole*)pElement[i];
			p->set_rad_onoff(flag);
			cnt++;
		}
		if ( pElement[i]->gettype()=="IdealWig")
		{	AIdealWig* p = (AIdealWig*)pElement[i];
			p->set_rad_onoff(flag);
			cnt++;
		}
		if ( pElement[i]->gettype()=="SyncRad")
		{	
			//cout<<"Found SyncRad:  "<<flag<<" "<<FLAG_RAD_ONOFF_LUMP<<" "<< (flag==FLAG_RAD_ONOFF_LUMP)<<endl;
			ASyncRad* p = (ASyncRad*)pElement[i];
			
			if(flag==FLAG_RAD_ONOFF_LUMP)
			{	
				//cout<<"Set Rad flag:  "<<endl;
				p->set_rad_onoff(1);
			}
			else if(flag==FLAG_RAD_ONOFF_LUMP_DAMP)
				p->set_rad_onoff(1,0);
			else if(flag==FLAG_RAD_ONOFF_LUMP_QE)
				p->set_rad_onoff(0,1);	
			else
				p->set_rad_onoff(FLAG_RAD_ONOFF_OFF); //do not use SyncRad elements
			cnt++;
		}
	}
	return cnt;
}

int ALine::set_monitor_pass(int num)
{
	int cnt=0;
	for(int i=0;i<nElement; i++)
	{
		if ( pElement[i]->gettype()=="Monitor")
		{	AMonitor* p = (AMonitor*)pElement[i];
			p->setpassnumber(num);
			cnt++;
		}
	}
	return cnt;
}


//parallel pass
#ifdef OPENMP_BY_LINE

void ALine::pass(double *r, int N,unsigned st, unsigned ed)
{
	#pragma omp parallel //num_threads(4)
	{	
	int this_thread = omp_get_thread_num();
	int num_threads = omp_get_num_threads();

	int my_start = (this_thread  ) * N / num_threads;
 	int my_end   = (this_thread+1) * N / num_threads;
 	int block_size = my_end-my_start;
 	//cout<<"thread: "<<this_thread<<" out of "<<num_threads<<" "<<my_start<<endl;
 	
	for(int i=st;i<ed; i++){
		pElement[i]->pass(r+my_start*6,block_size);
		}
	}
}

void ALine::pass(double *r, int N)
{
	if (pElement[0]->gettype() == "Monitor")
	{
		AMonitor* pm = (AMonitor*)pElement[0];
		pm->pass(r, N);
		pass(r, N, 1, nElement);
	}
	else
		pass(r, N, 0, nElement);
		
	/*
#pragma omp parallel //num_threads(4)
	{	
	int this_thread = omp_get_thread_num();
	int num_threads = omp_get_num_threads();

	int my_start = (this_thread  ) * N / num_threads;
 	int my_end   = (this_thread+1) * N / num_threads;
 	int block_size = my_end-my_start;
 	//cout<<"thread: "<<this_thread<<" out of "<<num_threads<<" "<<my_start<<endl;
 	
	for(int i=0;i<nElement; i++){
		pElement[i]->pass(r+my_start*6,block_size);
		}
	}*/
}

void ALine::pass(double *r, int N, int nturn)
{
	for (int i = 0; i<nturn; i++)
		pass(r, N);
}
/*
void ALine::pass(double *r, int N, int nturn)
{
#pragma omp parallel //num_threads(4)
	{
		int this_thread = omp_get_thread_num();
		int num_threads = omp_get_num_threads();

		int my_start = (this_thread)* N / num_threads;
		int my_end = (this_thread + 1) * N / num_threads;
		int block_size = my_end - my_start;

		for (int i = 0; i < nturn; i++)
			//pass(r + my_start * 6, block_size);
		{
			for (int j = 0; j<nElement; j++) {
				pElement[j]->pass(r + my_start * 6, block_size);
			}
		}

	}
}*/
#else

void ALine::pass(double *r, int N,unsigned st, unsigned ed)
{
	for(int i=st;i<ed; i++)
	{
		pElement[i]->pass(r,N);
	}
}
void ALine::pass(double *r, int N)
{
	for(int i=0;i<nElement; i++)
	{
		pElement[i]->pass(r,N);
	}

}
void ALine::pass(double *r, int N, int nturn)
{
	for(int i=0;i<nturn; i++)
		pass(r,N);	
}
#endif

int ALine::set_monitor(unsigned int flag, unsigned int index) //index=0 for all monitors, index>0 for one monitor, e.g., index=1 for the 1st monitor
{
	int cnt = 0;
	//cout<<"test"<<endl;
	for(int i=0;i<nElement; i++)
	{
		if(pElement[i]->gettype()=="Monitor")
		{AMonitor* p = (AMonitor*)pElement[i];
			//cout<<p->gettype()<<endl;
			cnt++;
			if(index==0)
				p->set_print_on(flag);
			else if (index==(cnt))
			{	p->set_print_on(flag); break; }
		}
	}
	return 0;
}
#define MAX_NUM_RFCAVITY 50
int ALine::set_cavity_onoff(bool flag)
{
	unsigned int nrf, rflist[MAX_NUM_RFCAVITY];
	//nrf = findelements(rflist,string("Cavity"),string("RF"),1);
	nrf = findelements(rflist, string("Cavity"));
	for(int i=0;i<nrf; i++){
			ACavity *pc = (ACavity *)getelement(rflist[i]);
			pc->set_cavity_onoff(flag);
		}
	return nrf;
}
#define MAX_NUM_CRABCAVITY 50
int ALine::set_crab_onoff(bool flag)
{
	unsigned int ncrab, crablist[MAX_NUM_CRABCAVITY];
	ncrab = findelements(crablist, string("CrabCav"));
	for(int i=0;i<ncrab; i++){
			ACrabCav *pc = (ACrabCav *)getelement(crablist[i]);
			pc->set_crab_onoff(flag);
		}
	return ncrab;
}

#define MAX_NUM_WAKE 50
int ALine::set_wake_onoff(bool flag)
{
	unsigned int nwake, wakelist[MAX_NUM_WAKE];
	nwake = findelements(wakelist, string("TransWake"));
	bool fprev;
	for (int i = 0; i<nwake; i++) {
		ATransWake *pc = (ATransWake *)getelement(wakelist[i]);
		fprev = pc->set_wake_onoff(flag);
	}
	return nwake;
}

double ALine::setrfphase()
{
		double r[]={0.0,0.0,0,0.0,0,0};
		set_cavity_onoff(0);
		set_crab_onoff(0);
		set_rad_onoff(1);
		
		pass(r,1); //1 turn
		double delta = -r[5];
		set_rad_onoff(0);
		//cout<<"energy loss/energy:  "<<delta<<endl;
		
		unsigned int nrf, rflist[MAX_NUM_RFCAVITY];
		nrf = findelements(rflist,string("Cavity"));
		if (!nrf) {
			cerr << "cavity not found" << endl;
			return 0;
		}
		
		double tot_volt=0, energy;
		for(int i=0;i<nrf; i++){
			ACavity *pc = (ACavity *)getelement(rflist[i]);
			tot_volt += pc->getvoltage();
			energy = pc->getenergy();
		}
	
		double phis;
		phis =  asin(delta*energy/tot_volt);
		//cout<<energy<<"\t"<<tot_volt<<"\t"<<delta<<endl;
		
		for(int i=0;i<nrf; i++){
			ACavity *pc = (ACavity *)getelement(rflist[i]);
			pc->setfield((void*)&phis, string("Phis"));			
		}
		return phis;

}
double ALine::setrfphase(double phis)
{
		unsigned int nrf, rflist[MAX_NUM_RFCAVITY];
		//nrf = findelements(rflist,string("Cavity"),string("RF"),1);
		nrf = findelements(rflist,string("Cavity"));
		if(!nrf) {
			cerr << "cavity not found" << endl;
			return 0;
		}
		
		for(int i=0;i<nrf; i++){
			ACavity *pc = (ACavity *)getelement(rflist[i]);
			pc->setfield((void*)&phis, string("Phis"));			
		}
		return phis;

}

unsigned int ALine::findelements(unsigned int *list, string type, string famname, unsigned int index)  //index = 1 for the 1st element, index=0 for all elements
{
	unsigned int cnt = 0;
		
	for(unsigned i=0;i<nElement; i++)
	{
		if((pElement[i]->gettype()==type)  && (pElement[i]->famname == famname ) )
		{
			cnt++;
			if(index==0)
				list[cnt-1] = i;
			else if (index==cnt)
			{list[cnt-1] = i;
				cnt = 1;
				break;
			}	
		}
	}
	return cnt; //number of element in list
}

unsigned int ALine::findelements(unsigned int *list, string type) 
{
	unsigned int cnt = 0;
		
	for(unsigned i=0;i<nElement; i++)
	{
		//cout<<pElement[i]->gettype()<<"=?"<<type<<endl;
		if((pElement[i]->gettype()==type)  )
		{
			cnt++;
			list[cnt-1] = i;
		}
	}
	return cnt; //number of element in list
}

unsigned int ALine::resetfirstelement(unsigned int eindex) //move the element at eindex to the beginning of the ring
{
	if(eindex==0)
		return 0;  //nothing need to be done
		
	//AccelElement **pTemp;
	AccelElement *pTemp[MAX_ELEMENT];
	int i,j;
	if(eindex<=nElement/2)
	{
		//pTemp= new (*AccelElement)[eindex];
		for(i=0;i<eindex; i++)
			pTemp[i] = pElement[i];
		for(i=0;i<nElement-eindex; i++)
			pElement[i] = pElement[i+eindex];
		for(j=0;i<nElement; i++,j++)	
			pElement[i] = pTemp[j];
	}
	else
	{
		//pTemp= new (*AccelElement)[nElement - eindex];
		for(i=eindex,j=0;i<nElement; i++,j++)
			pTemp[j] = pElement[i];
		for(i=nElement-1,j=eindex-1;j>=0; i--,j--)
			pElement[i] = pElement[j];
		for( i=0;i<nElement-eindex; i++)	
			pElement[i] = pTemp[i];
		
	}
	
	//delete []pTemp;	
	return nElement-eindex;
	
}

// return the index of the previous first element in the new order
unsigned int ALine::resetfirstelement(string type, string famname, unsigned int index) //move the matched element to the beginning of the ring
{
	unsigned int ne, eindex;
	if(index==0)
		cerr<<"resetfirstelement: element index cannot be zero"<<endl;
	ne = findelements(&eindex, type, famname,index);
	if(ne==0)
	{	
		cout<<"resetfirstelement: element not found (type ="<<type<<", famname ="<<famname<<", index ="<<index<<")"<<endl;
		return 0; //element not found
	}
	return resetfirstelement(eindex);
}

double ALine::trackDA(double xmax, double xmin, double ymax, int nturn, unsigned int nray, int nstep, double *rx, double *ry,bool flag_negweight)
{ //make sure rx and ry are double vector with size nray
	
	assert(xmax>0);
	assert(ymax>0);
	assert(xmin>0);
	
	if(nray<3)
		nray = 3;
	if(nray>51)
		nray = 51;
	if(nstep<2)
		nstep = 2;
			
	double deltax = xmax - xmin;
	if (deltax<0)
		{ cerr<<"xmax must be greater than xmin"<<endl; return -1;}
		
	double *xl,*yl;
	xl = new double[nray*nstep];
	yl = new double[nray*nstep];
	double *losst = new double[nray*nstep];
	double *r = new double[nray*nstep*6];
	for(int i=0;i<nray*nstep*6;i++) r[i] = 0;		
	
	for(int i=0;i<nray;i++)
		for(int j=0;j<nstep;j++)
			{
				xl[i*nstep+j] = (xmin+deltax*j/(nstep-1))*cos(i*PI/(nray-1));
				yl[i*nstep+j] = (xmin+deltax*j/(nstep-1))*sin(i*PI/(nray-1))*ymax/xmax;
				r[(i*nstep+j)*6+0] = xl[i*nstep+j];
				r[(i*nstep+j)*6+2] = yl[i*nstep+j];
				r[(i*nstep+j)*6+4] = 1.0e-5;
				losst [i*nstep+j] = 0;
			}
	
  //for loss over turn plot
	if(flag_msg & FLAG_MSG_DA)
	{	
		for (int t = 0; t<nturn; t++)
		{
			pass(r, nray*nstep);
			//cout<<r[0]<<" "<<r[1]<<"\t"<<r[2]<<" \t"<<r[3]<<"\t "<<r[4]<<"\t "<<r[5]<<endl;
			for (int i = 0; i<nray; i++)
				for (int j = 0; j<nstep; j++)
				{
					if (isnan(r[(i*nstep + j) * 6]) && !losst[i*nstep + j])
					{
						losst[i*nstep + j] = t;
					}

				}
		}
		cout<<"$$$DA-Loss "<<endl;
		for(int i=0;i<nray;i++)
				for(int j=0;j<nstep;j++)
				{
					cout<<xl[i*nstep+j]<<"\t"<<yl[i*nstep+j]<<"\t"<<losst[i*nstep+j]<<endl;
				}
	}
	else //no need to get loss turn number
	{
		pass(r, nray*nstep,nturn);

	}
	
	for(int i=0;i<nray;i++)
	{		rx[i] = xl[i*nstep+nstep-1];
			ry[i] = yl[i*nstep+nstep-1];

			for(int j=0;j<nstep;j++)
			{
				if(isnan(r[(i*nstep+j)*6]))
					{
						rx[i] = xl[i*nstep+j];
						ry[i] = yl[i*nstep+j];
						//cout<<i<<"\t"<<j<<endl;
						break;
					}
			}
		}
		
	double Sxy = 0;
	for(int i=1;i<nray;i++)
	{
		double rd1, rd2, dphi;
		rd1 = sqrt(rx[i-1]*rx[i-1]+ry[i-1]*ry[i-1]);
		rd2 = sqrt(rx[i]*rx[i]+ry[i]*ry[i]);
		
		//dphi = PI/(nray-1); //this is incorrect because the y-direction is squeezed
		dphi = atan2(ry[i],rx[i]) - atan2(ry[i-1],rx[i-1]);
		
	//	Sxy += .5*rd1*rd2*sin(dphi);
		if (rx[i]<=0)
			Sxy += .5*rd1*rd2*sin(dphi);
		else
		{
			if(flag_negweight) Sxy += .25*rd1*rd2*sin(dphi);  //reduce weight for positive side
			else Sxy += .5*rd1*rd2*sin(dphi);
		}
			
	}
	delete []losst;
	delete []r;
	delete []xl;
	delete []yl;
	
	return Sxy;
}
double ALine::trackDiffusion(int nturn, unsigned int np, double *x0, double *y0, double *nux, double *nuy, double *dif,double xlow, double xhi, double ylow, double yhi)
{
	double r[6];

	unsigned int nt=1;
	while(nt<nturn) 	nt <<=1;
	double mDiffu = -4*log(nt*1.0);
	
	double *xl = new double[nt];
	double *yl = new double[nt];	
	if(flag_msg & FLAG_MSG_DF)
		cout<<"Diffusion \t"<<nturn<<" turns\t"<<np<<endl;
		
	double tdif = 0;
	for(int ip=0;ip<np;ip++)
	{
		double nux0,nux1,nuy0,nuy1;
		
		r[1]=r[3]=r[4]=r[5] = 0.0;
		r[0]= x0[ip];
		r[2]= y0[ip];
		
		xl[0] = r[0];
		yl[0] = r[2];
		for(int j=1;j<nt;j++)
		{	
			pass(r, 1);
			xl[j] = r[0];
			yl[j] = r[2];
		}
		if(isnan(r[0]))
		{
			nux[ip] = CPPNAN;
			nuy[ip] = CPPNAN;
			dif[ip] = CPPNAN;
			nux0 = nuy0 = nux1 = nuy1 = CPPNAN;
			tdif += -3.0;
		}
		else
		{
			
			nux0 = naff(nt>>1, xl,xlow, xhi);
			nuy0 = naff(nt>>1, yl,ylow, yhi);
			nux1 = naff(nt>>1, xl+(nt>>1),xlow, xhi);
			nuy1 = naff(nt>>1, yl+(nt>>1),ylow, yhi);
		
			nux[ip] = nux0;
			nuy[ip] = nuy0;
			dif[ip] = MAX(mDiffu,log(sqrt( SQR(nux1-nux0)+SQR(nuy1-nuy0) )*2.0/nt));
			tdif += dif[ip];
		}	
		if(flag_msg & FLAG_MSG_DF)
		{	
			cout << xl[0]<<"\t"<< yl[0]<<"\t"<< nux[ip]<<"\t"<< nuy[ip]<<"\t"<< log(SQR(nux1-nux0))<<"\t"<< log(SQR(nuy1-nuy0))<<"\t"<< dif[ip]<<"\n";
		}

		
	}
	
	delete []xl;
	delete []yl;
	
	return tdif;
}

double ALine::trackMA(double dppmax, double dppmin, int nturn, unsigned int nstep, double *dpp_plus, double *dpp_minus,unsigned nfam, unsigned *lfam)
{
	assert(dppmax>0);
	assert(dppmin>0);
	
	double *r = new double[6*(nstep+1)*2];
	
	double stepsize = (dppmax-dppmin)/nstep;
	double dppavg=0;
	for(int i=0;i<nfam;i++)
	{
		int ind1st = resetfirstelement(lfam[i]);
		//cout<<lfam[i]<<"\t"<<pElement[0]->famname<<endl;
		for(int j=0;j<6*(nstep+1)*2;j++)
			r[j] = 0.0;
			
		for(int j=0;j<(nstep+1);j++)	
		{
			r[j*12+5] = dppmin + stepsize*j;
			r[j*12+6+5] = -(dppmin + stepsize*j);
		}
		
		for(int t=0;t<nturn;t++)
		{	
			pass(r, (nstep+1)*2);
			/*cout<<t<<endl;
			for(int j=1;j<(nstep+1);j++)
				cout<<r[j*12+0]<<"\t"<<r[j*12+5]<<"\t"<<r[j*12+6]<<"\t"<<r[j*12+6+5]<<endl;
			cout<<endl;
			*/
		}
		
		resetfirstelement(ind1st);
		
		dpp_plus[i] = dppmin;
		dpp_minus[i] = -dppmin;
		for(int j=1;j<(nstep+1);j++)	
		{
			//cout<<r[j*12]<<"\t"<<r[j*12+6]<<endl;
			
			if(!isnan(r[j*12]))
				dpp_plus[i] = dppmin + stepsize*j;
				
			if(!isnan(r[j*12+6]))
				dpp_minus[i] = -(dppmin + stepsize*j);
		}
		
		dppavg += fmin(dpp_plus[i], -dpp_minus[i]);
		//cout << i <<"\t"<<dppavg<<endl;
	}
	dppavg /= nfam;
	
	delete []r;
	
	return dppavg;
}

Vec& ALine::getco4(Vec& rc, double dp,double delta)
{
	Mat m66(6,6);
	set_rad_onoff(0);
	set_cavity_onoff(0);
	
	getm66(m66,dp,delta);
	Mat m44(4,4);
	m44 = m66.submatrix(m44,1,4,1,4);
	
	double r[6];
	for(int i=0;i<6;i++)r[i] = 0;
	r[5] = dp;
	pass(r,1);
	
	Vec v4(4);
	for(int i=0;i<4;i++) v4[i+1] = r[i];
	
	Mat I44(4,4,MAT_IDENTITY);
	rc = inv(I44-m44)*v4;
	
	/*#define MIN_ORBIT_PRECISION  1.0E-6
	int cnt=0;
	while(rc.L2_norm()>MIN_ORBIT_PRECISION && cnt++<3)
	{
		for(int i=0;i<4;i++) r[i] = rc[i+1];
		r[4]=0; r[5]=0;
		pass(r,1);	
		
	}
	*/
	return rc;
}
Vec& ALine::getco6(Vec& rc, double delta)
{
	Mat m66(6,6);
	set_rad_onoff(1);
	set_cavity_onoff(1);
	
	getm66(m66,0.,delta);
	
	double r[6];
	for(int i=0;i<6;i++)r[i] = 0;
	pass(r,1);
	
	Vec v6(6);
	for(int i=0;i<6;i++) v6[i+1] = r[i];
	
	Mat I66(6,6,MAT_IDENTITY);
	rc = inv(I66-m66)*v6;
	return rc;
}
	
Mat& ALine::getm66(Mat& m66,double dpp,double delta)
{
	return getm66(m66,0,nElement,dpp,delta);
}	

Mat& ALine::getm66(Mat& m66,unsigned int ed, double dpp,double delta)
{
	return getm66(m66,0,ed,dpp,delta);	
}	
Mat& ALine::getm66(Mat& m66,unsigned int st, unsigned int ed,double dpp,double delta)
{
	if(st<0 || st>nElement)
		cerr<<"element index of range"<<endl;
	if(ed<0 || ed>nElement)
		cerr<<"element index of range"<<endl;
	if(ed<st)
		cerr<<"start > end"<<endl;
		
	double r0[72];
	for(int i=0;i<72; i++)  r0[i]=0.0;
	r0[0] = r0[7] = r0[14] = r0[21]	= r0[28] = r0[35] = delta;
	r0[36+0] = r0[36+7] = r0[36+14] = r0[36+21]	= r0[36+28] = r0[36+35] = -delta;
	for(int i=0;i<12; i++)  r0[i*6+5] +=dpp;
	
	for(int i=st;i<ed; i++)
	{
		pElement[i]->pass(r0,12);
	}
	//cout<< r0[0]<<endl;
	//Mat m66(6,6);
	for(int i=0;i<6;i++) //row
		for (int j=0;j<6;j++) //col
		{
			m66[i+1][j+1] = (r0[j*6+i]-r0[j*6+i+36])/2.0/delta;
		}
	return m66;
}	
int ALine::getfield(void *fvalue, unsigned nfam, unsigned *lfam,string pname,unsigned index) //index, position of the field in an array
{
	unsigned nfbyte;
	/*if(pname=="Length" || pname=="K"  || pname=="PolynomA" || pname=="PolynomB" || pname=="PolynomA" 
		|| pname=="BendingAngle" || pname=="EntranceAngle" || pname=="ExitAngle" || pname=="energy"
		|| pname=="KickX" || pname=="KickY" || pname=="Voltage" || pname=="Frequency" || pname=="Phis"		)
		{
			nfbyte = sizeof(double);
		}
		*/
		if (pname=="NSlice")
			nfbyte = sizeof(int);
		else
		{
			nfbyte = sizeof(double);
			//cerr<<"type not supported by ALine->getfield:  "<<pname<<endl;
			//return -1;
		}
	char* pfi = (char*) fvalue;		
	for(int i=0;i<nfam;i++)
	{
		if(pname=="PolynomA" || pname=="PolynomB")
		{
			if(pElement[lfam[i]]->gettype()=="BendMPole")
			{
				ABendMPole* p= (ABendMPole*)pElement[lfam[i]];
				p->getfield((void*) pfi,pname,index);		
			}
			else if(pElement[lfam[i]]->gettype()=="StrMPole")
			{
				AStrMPole* p= (AStrMPole*)pElement[lfam[i]];
				p->getfield((void*) pfi,pname,index);	
			}
			else if(pElement[lfam[i]]->gettype()=="ARFMult")
			{
				ARFMult* p= (ARFMult*)pElement[lfam[i]];
				p->getfield((void*) pfi,pname,index);	
			}
			else
			{
				cerr<<"field "<<pname<<" not present in element type "<<pElement[lfam[i]]->gettype()<<endl;
				return -1;
			}
		}
		else
			pElement[lfam[i]]->getfield((void*) pfi,pname);		
		
		pfi += nfbyte;
	}
	
	return 0;
}
int ALine::setfield(void *fvalue, unsigned nfam, unsigned *lfam,string pname,unsigned index)
{
	unsigned nfbyte;
	/*if(pname=="Length" || pname=="K"  || pname=="PolynomA" || pname=="PolynomB" || pname=="PolynomA" 
		|| pname=="BendingAngle" || pname=="EntranceAngle" || pname=="ExitAngle" || pname=="energy"
		|| pname=="KickX" || pname=="KickY" || pname=="Voltage" || pname=="Frequency" || pname=="Phis"	
			) 
		{
			nfbyte = sizeof(double);
		}
		else */
		if (pname=="NSlice")
			nfbyte = sizeof(int);
		else
		{
			nfbyte = sizeof(double);
			//cerr<<"type not supported by ALine->setfield:  "<<pname<<endl;
			//return -1;
		}
	char* pfi = (char*) fvalue;		
	for(int i=0;i<nfam;i++)
	{
		if(pname=="PolynomA" || pname=="PolynomB")
		{
			if(pElement[lfam[i]]->gettype()=="BendMPole")
			{
				ABendMPole* p= (ABendMPole*)pElement[lfam[i]];
				p->setfield((void*) pfi,pname,index);		
			}
			else if(pElement[lfam[i]]->gettype()=="StrMPole")
			{
				AStrMPole* p= (AStrMPole*)pElement[lfam[i]];
				p->setfield((void*) pfi,pname,index);	
			}
			else if(pElement[lfam[i]]->gettype()=="ARFMult")
			{
				ARFMult* p= (ARFMult*)pElement[lfam[i]];
				p->setfield((void*) pfi,pname,index);	
			}
			else
			{
				cerr<<"field "<<pname<<" not present in element type "<<pElement[lfam[i]]->gettype()<<endl;
				return -1;
			}
		}
		else
			pElement[lfam[i]]->setfield((void*) pfi,pname);		
		
		pfi += nfbyte;	
	}
	return 0;
}	
void ALine::calctune(double &tunex, double &tuney,double dp, double delta)
{
	Mat m66(6,6);
	getm66(m66,dp,delta);
	
	tunex = acos((m66[1][1]+m66[2][2])*0.5)/2.0/PI;
	tuney = acos((m66[3][3]+m66[4][4])*0.5)/2.0/PI;
//	double betx, bety;
//	betx = m66[1][2]/sin(Phix);
//	bety = m66[3][4]/sin(Phiy);
	
}
void ALine::calctwiss(Twiss *tw,unsigned int indx[], int n,double dp, double delta)
{
	Mat m66(6,6);
	getm66(m66,dp,delta);
	Twiss tw0;
	tw0.s = 0;
	tw0.psix = acos((m66[1][1]+m66[2][2])*0.5);
	tw0.psiy = acos((m66[3][3]+m66[4][4])*0.5);

	tw0.betx = m66[1][2]/sin(tw0.psix);
	tw0.bety = m66[3][4]/sin(tw0.psiy);
	tw0.alfx = (m66[1][1] - m66[2][2])/sin(tw0.psix)/2.0;
	tw0.alfy = (m66[3][3] - m66[4][4])/sin(tw0.psiy)/2.0;
	tw0.psix = 0;
	tw0.psiy = 0;
	//
	Mat T01(6,6); //matrix from point 0 to indx
	unsigned int id = indx[0];
	getm66(T01,0,id,dp,delta);
	//cout<<T01<<endl;
	double len01 = getlength(0,id);
	for(int i=0;i<n;i++)
	{
		Twiss tw1; 
		tw1.s = tw0.s+len01;
		tw1.betx = 1.0/tw0.betx*(SQR(T01[1][1]*tw0.betx-T01[1][2]*tw0.alfx)+SQR(T01[1][2]));
		tw1.bety = 1.0/tw0.bety*(SQR(T01[3][3]*tw0.bety-T01[3][4]*tw0.alfy)+SQR(T01[3][4]));
		tw1.alfx = -1.0/tw0.betx*((T01[1][1]*tw0.betx-T01[1][2]*tw0.alfx)*(T01[2][1]*tw0.betx-T01[2][2]*tw0.alfx)+T01[1][2]*T01[2][2]);
		tw1.alfy = -1.0/tw0.bety*((T01[3][3]*tw0.bety-T01[3][4]*tw0.alfy)*(T01[4][3]*tw0.bety-T01[4][4]*tw0.alfy)+T01[3][4]*T01[4][4]);
		tw1.psix = tw0.psix+atan2(T01[1][2],T01[1][1]*tw0.betx-T01[1][2]*tw0.alfx); 
		tw1.psiy = tw0.psiy+atan2(T01[3][4],T01[3][3]*tw0.bety-T01[3][4]*tw0.alfy); 
		
		//tw0.serialize(cout);
		//tw1.serialize(cout);
		tw[i] = tw1;
		tw0 = tw1;
		                
		unsigned int id2;
		if(i<n-1)
		{ id2 = indx[i+1];
			getm66(T01,id,id2,dp,delta);
			len01 = getlength(id,id2);
			id = id2;		
		}
		
	}
			
}	
void ALine::calctwiss(Twiss& tw,double dp, double delta)
{
	Mat m66(6,6);
	getm66(m66,dp,delta);
	
	tw.s = 0;
	tw.psix = acos((m66[1][1]+m66[2][2])*0.5);
	tw.psiy = acos((m66[3][3]+m66[4][4])*0.5);

	tw.betx = m66[1][2]/sin(tw.psix);
	tw.bety = m66[3][4]/sin(tw.psiy);
	tw.alfx = (m66[1][1] - m66[2][2])/sin(tw.psix)/2.0;
	tw.alfy = (m66[3][3] - m66[4][4])/sin(tw.psiy)/2.0;
	
}	

void ALine::calcchrom(double &chromx, double &chromy,double dp,double delta)
{
	double tunex0,tuney0,tunex,tuney;
  double dpp = 1.0e-5;
  
	Mat m66(6,6);
	
	getm66(m66,dp-dpp,delta);	
	tunex0 = acos((m66[1][1]+m66[2][2])*0.5)/2.0/PI;
	tuney0 = acos((m66[3][3]+m66[4][4])*0.5)/2.0/PI;
	
	getm66(m66,dp+dpp,delta);
	tunex = acos((m66[1][1]+m66[2][2])*0.5)/2.0/PI;
	tuney = acos((m66[3][3]+m66[4][4])*0.5)/2.0/PI;
	
	chromx = (tunex-tunex0)/dpp/2.0;
	chromy = (tuney-tuney0)/dpp/2.0;

}	
void ALine::calcchrom2(double &chromx, double &chromy,double &chromx2, double &chromy2)
{
	double chrx0,chry0;
  const int np = 11; 
   
  double *cxl = new double[np];
	double *cyl = new double[np];
	double *dpl = new double[np];
	
	for(int i=0;i<np;i++)
	{
		double dpp = 1.0e-3*(i-6.0);
		calcchrom(chrx0, chry0,dpp);
		cxl[i] = chrx0;
		cyl[i] = chry0;
		dpl[i] = dpp;
	}
	 
	 CLinfit ca; //double a, b, siga, sigb, chi2, q;
	
   linearfit(ca, np, dpl, cxl);
	 chromx2 = ca.b/2.0; //chromx2
	 
	 linearfit(ca, np, dpl, cyl);
	 chromy2 = ca.b/2.0; //chromy2
	 
	calcchrom(chromx, chromy, 0);
	
	delete []cxl;
	delete []cyl;
	delete []dpl;
}	
void ALine::tracktune(double &tunex, double &tuney, unsigned int nturn,const double x0, const double y0,double xlow, double xhi,double ylow, double yhi)
{
	double r0[6];
	for(int i=0;i<6; i++)  r0[i]=0.0;
	r0[0] = x0;
	r0[2] = y0;
	
	unsigned int nt=1;
	while(nt<nturn) 	nt <<=1;
	
	double *xl = new double[nt];
	double *yl = new double[nt];
	
	xl[0] = r0[0];
	yl[0] = r0[2];
	for(int i=1;i<nt;i++)
	{
		pass(r0,1);
		xl[i] = r0[0];
		yl[i] = r0[2];
	}
	tunex = naff(nt, xl,xlow, xhi);
	tuney = naff(nt, yl,ylow, yhi);
	
	delete []xl;
	delete []yl;
	
}
void ALine::calctuneresp(double resp[][2], unsigned nknob, unsigned *nfam, unsigned **lfam)
{
	double tunex0, tuney0;
	calctune(tunex0,tuney0);
	
	double delta = 1.0e-5;
	for(int i=0;i<nknob;i++)
	{
		double *K1old = new double[nfam[i]];
		for(int j=0;j<nfam[i];j++)
			{
				double K1new;
				unsigned index = lfam[i][j];
				//first get original values, then set value with delta added
				if(pElement[index]->gettype() == "Quad")
				{	 
					pElement[index]->getfield(K1old+j,"K");
					K1new = K1old[j]+delta;
					pElement[index]->setfield((void*) &K1new,"K");
				}
				else if(pElement[index]->gettype() == "StrMPole")
				{
					AStrMPole *p = (AStrMPole *)pElement[index];
					p->getfield(K1old+j,"PolynomB",1);
					K1new = K1old[j]+delta;
					p->setfield((void*) &K1new,"PolynomB",1);
					
				}
				
			}	//end innoer for 1
			double tunex, tuney;
			calctune(tunex,tuney);
			resp[i][0] = (tunex-tunex0)/delta;
			resp[i][1] = (tuney-tuney0)/delta;
			
			for(int j=0;j<nfam[i];j++)
			{
				unsigned index = lfam[i][j];
				//set to original values
				if(pElement[index]->gettype() == "Quad")
				{	 
					pElement[index]->setfield((void*) &K1old[j],"K");
				}
				else if(pElement[index]->gettype() == "StrMPole")
				{
					AStrMPole *p = (AStrMPole *)pElement[index];
					
					p->setfield((void*) &K1old[j],"PolynomB",1);
					
				}
				
			}	//end innoer for 2
				
		delete []K1old;
	}	//end outer for
}
void ALine::correcttune(double &tunex, double &tuney, unsigned int nfam1, unsigned int nfam2,unsigned int *lfam1, unsigned int *lfam2)
{
	double tunex0, tuney0;
	calctune(tunex0,tuney0);
	
	double resp[2][2];
	unsigned *lfam[2];
	unsigned nfam[2];
	nfam[0] = nfam1;
	nfam[1] = nfam2;
	lfam[0] = lfam1;
	lfam[1] = lfam2;
	
	calctuneresp(resp, 2, nfam, lfam);
	
	double detr = resp[0][0]*resp[1][1] - resp[0][1]*resp[1][0];
	double delta1,delta2;
	delta1 = -( resp[1][1]*(tunex0-tunex) - resp[1][0]*(tuney0-tuney))/detr;
	delta2 = -(-resp[0][1]*(tunex0-tunex) + resp[0][0]*(tuney0-tuney))/detr;
	
	
	double *K1fam1,*K1fam2;
	K1fam1 = new double[nfam1];
	K1fam2 = new double[nfam2];
	
	string fam1type,fam2type;
	fam1type = pElement[lfam1[0]]->gettype();
	fam2type = pElement[lfam2[0]]->gettype();
	
	if(fam1type=="Quad")
	{
			getfield((void*)K1fam1, nfam1, lfam1,"K");
		}
	else if(fam1type=="StrMPole") 	
		getfield((void*)K1fam1, nfam1, lfam1,"PolynomB",1);
	else
		{ cerr<<"element type " <<fam1type<<" not valid for tune correction (only Quad or StrMPole are)"<<endl;  
			delete []K1fam1;
			delete []K1fam2; return; }
	if(fam2type=="Quad")
		getfield((void*)K1fam2, nfam2, lfam2,"K");
	else if(fam2type=="StrMPole") 	
		getfield((void*)K1fam2, nfam2, lfam2,"PolynomB",1);
	else
		{ cerr<<"element type " <<fam1type<<" not valid for tune correction (only Quad or StrMPole are)"<<endl;
			delete []K1fam1;
			delete []K1fam2; return; }
	
	for(int i=0;i<nfam1;i++)
		K1fam1[i] += delta1;
	for(int i=0;i<nfam2;i++)
		K1fam2[i] += delta2;	
		
	if(fam1type=="Quad")
		setfield((void*)K1fam1, nfam1, lfam1,"K");
	else if(fam1type=="StrMPole") 	
		setfield((void*)K1fam1, nfam1, lfam1,"PolynomB",1);
	if(fam2type=="Quad")
		setfield((void*)K1fam2, nfam2, lfam2,"K");
	else if(fam2type=="StrMPole") 	
		setfield((void*)K1fam2, nfam2, lfam2,"PolynomB",1);
	
	calctune(tunex,tuney);	
	
	delete []K1fam1;
	delete []K1fam2;
}
void ALine::calcchromresp(double resp[][2], unsigned nknob, unsigned *nfam, unsigned **lfam)
{
	double chromx0, chromy0;
	calcchrom(chromx0,chromy0);
	
	double delta = 1.0e-5;
	for(int i=0;i<nknob;i++)
	{
		double *K2old = new double[nfam[i]];
		for(int j=0;j<nfam[i];j++)
			{
				double K2new;
				unsigned index = lfam[i][j];
				//first get original values, then set value with delta added
				if(pElement[index]->gettype() == "StrMPole")
				{
					AStrMPole *p = (AStrMPole *)pElement[index];
					p->getfield(K2old+j,"PolynomB",2);
					K2new = K2old[j]+delta;
					p->setfield((void*) &K2new,"PolynomB",2);
					
				}
				
			}	//end innoer for 1
			double chromx, chromy;
			calcchrom(chromx,chromy);
			resp[i][0] = (chromx-chromx0)/delta;
			resp[i][1] = (chromy-chromy0)/delta;
			
			for(int j=0;j<nfam[i];j++)
			{
				unsigned index = lfam[i][j];
				//set to original values
				if(pElement[index]->gettype() == "StrMPole")
				{
					AStrMPole *p = (AStrMPole *)pElement[index];
					
					p->setfield((void*) &K2old[j],"PolynomB",2);
					
				}
				
			}	//end innoer for 2
			
		delete []K2old;
	}	//end outer for
}
void ALine::correctchrom(double &chromx, double &chromy, unsigned int nfam1, unsigned int nfam2,unsigned int *lfam1, unsigned int *lfam2)
{
	double chromx0, chromy0;
	calcchrom(chromx0,chromy0);
	
	double resp[2][2];
	unsigned *lfam[2];
	unsigned nfam[2];
	nfam[0] = nfam1;
	nfam[1] = nfam2;
	lfam[0] = lfam1;
	lfam[1] = lfam2;
	
	calcchromresp(resp, 2, nfam, lfam);
	
	double detr = resp[0][0]*resp[1][1] - resp[0][1]*resp[1][0];
	double delta1,delta2;
	delta1 = -( resp[1][1]*(chromx0-chromx) - resp[1][0]*(chromy0-chromy))/detr;
	delta2 = -(-resp[0][1]*(chromx0-chromx) + resp[0][0]*(chromy0-chromy))/detr;
	
	
	double *K2fam1,*K2fam2;
	K2fam1 = new double[nfam1];
	K2fam2 = new double[nfam2];
	
	string fam1type,fam2type;
	fam1type = pElement[lfam1[0]]->gettype();
	fam2type = pElement[lfam2[0]]->gettype();
	
	if(fam1type=="StrMPole") 	
		getfield((void*)K2fam1, nfam1, lfam1,"PolynomB",2);
	else
		{ cerr<<"element type " <<fam1type<<" not valid for chrom correction (only StrMPole )"<<endl;  
			delete []K2fam1;
			delete []K2fam2; return; }
	if(fam2type=="StrMPole") 	
		getfield((void*)K2fam2, nfam2, lfam2,"PolynomB",2);
	else
		{ cerr<<"element type " <<fam1type<<" not valid for chrom correction (only StrMPole )"<<endl;
			delete []K2fam1;
			delete []K2fam2; return; }
	
	for(int i=0;i<nfam1;i++)
		K2fam1[i] += delta1;
	for(int i=0;i<nfam2;i++)
		K2fam2[i] += delta2;	
	
	
	if(flag_msg & FLAG_MSG_CHROM_CORR)
	{	cout<<"Chrom corr"<<endl;
		cout<<K2fam1[0]<<"\t"<<K2fam2[0]<<"\t"<<endl;
	}
			
	if(fam1type=="StrMPole") 	
		setfield((void*)K2fam1, nfam1, lfam1,"PolynomB",2);
	if(fam2type=="StrMPole") 	
		setfield((void*)K2fam2, nfam2, lfam2,"PolynomB",2);
	
	calcchrom(chromx,chromy);	
	
	delete []K2fam1;
	delete []K2fam2;
}


void  ALine::trackdnudexy(double *dnudexy, unsigned int nturn,const double xm, const double ym)
{
	 double tunex, tuney;
	 Twiss tw;
	 calctwiss(tw);
	 tunex = tw.psix/2/PI;
	 tuney = tw.psiy/2/PI;
	 //cout<<tunex<<" "<<tuney<<endl;
	 
	 const unsigned int np = 6;
	 double xtunel[np], ytunel[np], de[np];
	 
	 double r0[6];
	 
	 unsigned int nt=1;
	 while(nt<nturn) 	nt <<=1;
	 		
	 double *xl = new double[nt];
	 double *yl = new double[nt];
	 CLinfit ca; //double a, b, siga, sigb, chi2, q;
	
	 const double TINY = 1.0E-5;

	 //cout<<xm<<" "<<ym<<" "<<nturn<<endl;
	 for(int k=0;k<np;k++)
	 {	
	 	 for(int i=0;i<6; i++)  r0[i]=0.0;
			r0[0] = TINY + k*1.0/(np-1)*xm;
			r0[2] = TINY + r0[0]*.02;
			de[k] = r0[0]*r0[0]/tw.betx;
			
			for(int i=0;i<nt;i++)
			{
				xl[i] = r0[0];
				yl[i] = r0[2];
				pass(r0,1);
				//cout<<xl[i]<<"\t"<<yl[i]<<endl;
				
			}
			
			xtunel[k] = naff(nt, xl,MAX(0,tunex-0.1), MIN(0.5,tunex+0.1));
			ytunel[k] = naff(nt, yl,MAX(0,tuney-0.1), MIN(0.5,tuney+0.1));
			cout<<de[k]<<"\t"<<xtunel[k]<<"\t"<<ytunel[k]<<endl;
		}
   linearfit(ca, np, de, xtunel);
	 dnudexy[0] = ca.b; //dnux/dex
	 
	 linearfit(ca, np, de, ytunel);
	 dnudexy[1] = ca.b; //dnuy/dex   
	 
	 
	 for(int k=0;k<np;k++)
	 {
	 	  for(int i=0;i<6; i++)  r0[i]=0.0;	
			r0[2] = TINY + k*1.0/(np-1)*ym;
		  r0[0] = TINY + r0[2]*.02;
			de[k] = r0[2]*r0[2]/tw.bety;
			
			//cout<<"%*******"<<endl;
			for(int i=0;i<nt;i++)
			{
				xl[i] = r0[0];
				yl[i] = r0[2];
				pass(r0,1);	
				//cout<<	xl[i]<<"\t "<<	yl[i]<<"\t "<<endl;
			}
			
			xtunel[k] = naff(nt, xl,MAX(0,tunex-0.1), MIN(0.5,tunex+0.1));
			ytunel[k] = naff(nt, yl,MAX(0,tuney-0.1), MIN(0.5,tuney+0.1));
			cout<<de[k]<<"\t"<<xtunel[k]<<"\t"<<ytunel[k]<<endl;
		}
   linearfit(ca, np, de, xtunel);
	 dnudexy[2] = ca.b; //dnux/dey
	 
	 linearfit(ca, np, de, ytunel);
	 dnudexy[3] = ca.b; //dnuy/dey   
   
	delete []xl;
	delete []yl;
}

unsigned int ALine::flag_msg=0;

void Twiss::serialize(ostream& os)
{
		os<<s<<"\t"<<betx<<"\t"<<bety<<"\t"<<alfx<<"\t"<<alfy<<"\t"<<psix<<"\t"<<psiy<<"\t"<<endl;
}
