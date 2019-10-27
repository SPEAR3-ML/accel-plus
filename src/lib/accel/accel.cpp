#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <limits>
#include <assert.h>
#include <omp.h>
#include "atpass.h"
#include "matvec.h"
#include "accel.h"
#include "utility.h"

/***
Accel elements
Created by Xiaobiao Huang, SLAC, Dec 2011,
Last modified Jan. 2017
All rights reserved. 
***/ 

AccelElement::AccelElement(double len)
{ this-> len = len;
	this-> famname = "";
	this-> type = "";
}
void AccelElement::initializemap()
{
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Drift", new ADrift));
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Quad", new AQuad));
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Corrector", new ACorrector));	
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Aperture", new AAperture));	
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Monitor", new AMonitor));	
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Bend", new ABend));
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("BendMPole", new ABendMPole));	
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("StrMPole", new AStrMPole));		
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Cavity", new ACavity));
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("WigTable", new AWigTable));					
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("IdealWig", new AIdealWig));
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("CrabCav", new ACrabCav));	
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("RFMult", new ARFMult));	
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("TransWake", new ATransWake));
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("Matrix", new AMatrix));	
	AccelElement::maptype.insert(std::pair<string, AccelElement*>("SyncRad", new ASyncRad));
}

AccelElement::AccelElement(double len, const string & famname, const string & type)
{ this-> len = len;
	this-> famname = famname;
	this-> type = type;
}
map<string, AccelElement*> AccelElement::maptype;

void AccelElement::setfield(void *p, const string pname)
{
	if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
	
}
void AccelElement::getfield(void *p, const string pname)
{
	if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
}	
AccelElement* AccelElement::deserialize(istream& is)
{
	string str,famname,type;
	if (is.eof())
	{	cout<<"end of file"<<endl; return NULL; }
	
	try {
	is>>famname;
	//std::string skip; // dummy
	//std::getline(std::getline(is, skip, '"'), famname, '"');
	//cout << famname<<"\t";
	if (is.eof())
	{	//cout<<"end of file"<<endl; 
		return NULL; }
	
	if(famname=="RETURN")
	{	//cout<<"end of LINE"<<endl; 
		return NULL; }
	if(famname[0]=='!')
	{	//cout<<"skip a line"<<endl; 
		getline(is,str); }
	
	if(famname.at(famname.length()-1)==':')
	{
		famname.erase(famname.length()-1,1);
	}
	else
	{
		getline(is, str,':');
		//cout<< endl;
		//cout<<": not found\n";
		if (is.eof())
		{	
			cout<<"end of file"<<endl; return NULL; 
		}
	}
	/*
	size_t pos = famname.find(':');
		cout<<pos<<' '<<famname.length();
	if(pos>=famname.length())
	{	getline(is, str,':');
		//cout<< endl;
		//cout<<": not found\n";
		if (is.eof())
	{	cout<<"end of file"<<endl; return NULL; }
	}
	else if(pos==famname.length()-1)
		famname.erase(pos,1);
	else
	*/

	is>>type;
	//cout<<type<<endl;
	if (is.eof())
	{	cout<<"end of file"<<endl; return NULL; }
	
	//cout << type;
	if (maptype.count(type)==0)
		throw string("type not found");
	return (AccelElement*) maptype[type]->create(is, famname);
}
catch (string& str)  //(exception &e)
{	cout<<"exception caught: "<<str<<endl;
	return NULL;
}

}

/**** DRIFT  **************************************/
ADrift& ADrift::operator=(ADrift& a)
{
			len = a.len;
			famname = a.famname;
			type = a.type;
			return *this;
}
void ADrift::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<endl;
}
ADrift* ADrift::create(istream& is, string& famname)
{
	string str;
	//is>>famname;
	getline(is, str);
	
	std::stringstream(str)>>len;

	return new ADrift(len,famname);
	//return new ADrift;
}
void ADrift::pass(double *r, int N)
{
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
	{
		if(isnan(*(r+6*i)))
			continue;
			
		double *r6 = r+6*i;
		fastdrift(r6,len/(1.0+r6[5]));
		
		/*r6[0] += r6[1]*this->len;
		r6[2] += r6[3]*this->len;
		r6[4] -= SQR(1.0/(1+r6[5]))*this->len*(SQR(r6[1]+SQR(r6[2])))/2.0;
		*/
	}		
}

/**** QUAD  **************************************/
void AQuad::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<K<<endl;
}
AQuad* AQuad::create(istream& is, string& famname)
{
	string str;
	//is>>famname;
	getline(is, str);
	std::stringstream(str)>>len>>K;
	
	return new AQuad(len,K,famname);
}
void AQuad::pass(double *r, int N)
{
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
	{
		if(isnan(*(r+6*i)))
			continue;
		
		double *r6;	
		r6 = r+6*i;
		quad6(r6, this->len, this->K);
	
	}		
}
void AQuad::setfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("K"))
		{ K = *(double*)p; }
		
	}
void AQuad::getfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("K"))
		{ *(double*)p = K ; }
		
	}

/**** Corrector  **************************************/
void ACorrector::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<kickx<<"\t"<<kicky<<endl;
}	
ACorrector* ACorrector::create(istream& is, string& famname)
{
	string str;
	//is>>famname;
	getline(is, str);
	std::stringstream(str)>>len>>kickx>>kicky;
	
	return new ACorrector(len,kickx,kicky,famname);
}		
void ACorrector::setfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("KickX"))
		{ kickx = *(double*)p; }
		else if (!pname.compare("KickY"))
		{ kicky = *(double*)p; }
		
	}
void ACorrector::getfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname ; }
		else if (!pname.compare("KickX"))
		{ *(double*)p = kickx ; }
		else if (!pname.compare("KickY"))
		{ *(double*)p = kicky ; }
		
	}
	
void ACorrector::pass(double *r, int N)
{
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;

			double *r6;
			double pn;				
			r6 = r+6*i;
			pn = (1+r6[5]);
			
			if(len!=0)fastdrift(r6, 0.5*len/pn);
			r6[1] += kickx/pn;
			r6[3] += kicky/pn;
			if(len!=0)fastdrift(r6, 0.5*len/pn);
		}		
	
}	

/**** Aperture  **************************************/
void AAperture::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<xlim[0]<<"\t"<<xlim[1]<<"\t"<<ylim[0]<<"\t"<<ylim[1]<<endl;
}	
AAperture* AAperture::create(istream& is, string& famname)
{
	string str;
	//is>>famname;
	getline(is, str);
	std::stringstream(str)>>len>>xlim[0]>>xlim[1]>>ylim[0]>>ylim[1];
	
	return new AAperture(len,xlim,ylim,famname);
}	
void AAperture::setfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("XLim"))
		{ double *a = (double*)p; 
			xlim[0] = a[0]; xlim[1] = a[1];}
		else if (!pname.compare("YLim"))
		{ double *a = (double*)p; 
			ylim[0] = a[0]; ylim[1] = a[1];}
		
	}
void AAperture::getfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{  *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{  *(string*)p = famname; }
		else if (!pname.compare("XLim"))
		{ double *a = (double*)p; 
			a[0] = xlim[0] ; a[1] = xlim[1];}
		else if (!pname.compare("YLim"))
		{ double *a = (double*)p; 
			a[0] = ylim[0]; a[1] = ylim[1];}
		
	}
	
void AAperture::pass(double *r, int N)
{
	if(len==0){
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif
		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;
				
			double *r6;
			r6 = r+6*i;
			if(r6[0]<=xlim[0] || r6[0]>=xlim[1] || r6[2]<=ylim[0] || r6[2]>=ylim[1])
			{
				r6[4] = (r6[0]<=xlim[0])*1 + (r6[0]>=xlim[1])*2 + (r6[2]<=ylim[0])*3 +  (r6[2]>=ylim[1])*4; //
				r6[0] = CPPNAN;				
			}
		}		
	}
	else
	{
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;
				
			double *r6;
			r6 = r+6*i;
			if(r6[0]<=xlim[0] || r6[0]>=xlim[1] || r6[2]<=ylim[0] || r6[2]>=ylim[1])
			{
				r6[4] = (r6[0]<=xlim[0])*1 + (r6[0]>=xlim[1])*2 + (r6[2]<=ylim[0])*3 +  (r6[2]>=ylim[1])*4; //
				r6[0] = CPPNAN;
				continue;
			}
						
			fastdrift(r6, len/(1+r6[5]));
			
			if(r6[0]<=xlim[0] || r6[0]>=xlim[1] || r6[2]<=ylim[0] || r6[2]>=ylim[1])
			{
				r6[4] = (r6[0]<=xlim[0])*1 + (r6[0]>=xlim[1])*2 + (r6[2]<=ylim[0])*3 +  (r6[2]>=ylim[1])*4; //
				r6[0] = CPPNAN;
			}
			
		}		
	}
}	

/**** Monitor  **************************************/
void AMonitor::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<interval_turn<<"\t"<<ratio_print<<"\t"<<flag_print<<endl;
}		
AMonitor* AMonitor::create(istream& is, string& famname)
{
	string str;
	getline(is, str);
	std::stringstream(str)>>len>>interval_turn>>ratio_print;//>>flag_print;
	//cout<<str<<endl;
	//cout<<len<<" "<<interval_turn<<" "<<ratio_print<<endl;
	if(len!=0){len=0; cerr<<"Monitor length should always be 0"<<endl; }
	if(interval_turn<=1) interval_turn=1;
	if(ratio_print<0 ) ratio_print=0.1;
	if(ratio_print>1.0) ratio_print=1.0;
	
	return new AMonitor(interval_turn,ratio_print,string(""),famname);
}	

void AMonitor::setfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("Interval"))
		{ interval_turn = *(int*)p; }
		else if (!pname.compare("Ratio"))
		{ ratio_print = *(double*)p; }
		else if (!pname.compare("PrintFlag"))
		{ flag_print = *(unsigned int*)p; }
		
	}
void AMonitor::getfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{  *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{  *(string*)p = famname; }
		else if (!pname.compare("Interval"))
		{ *(int*)p = interval_turn ;}
		else if (!pname.compare("Ratio"))
		{ *(double*)p = ratio_print ;}	
		else if (!pname.compare("PrintFlag"))
		{ *(unsigned int*)p = flag_print ;}	
		
	}
void AMonitor::pass(double *r, int N)
{
	assert(len==0);
	
	cnt_pass++;
	if( (interval_turn>1) && (((cnt_pass-1)%interval_turn)!=0)) return;
	if(!flag_print) return;
	
	int iskip = int(1./ratio_print);
	if(iskip<1)iskip=1;
	if(iskip>N)iskip=N;
	
	if(flag_print ==1) 
	{ //print coordinate
		for (int i = 0; i<N; i += iskip) 
		{
			double *rp = r + 6 * i;
			cout << cnt_pass-1 << "\t" << rp[0] << "\t" << rp[1] << "\t" << rp[2] << "\t" << rp[3] << "\t" << rp[4] << "\t" << rp[5] << endl;
		}
	}
	else if (flag_print == 2 && N>=2)
	{ //print statistics, no off-diagonal elements
		double ar[6], sigr[6];
		int nlive = calc_beam_stat(r, N, ar, sigr, NULL);
		cout << std::scientific;
		cout << cnt_pass - 1 << "\t" <<nlive<< "\t" << ar[0] << "\t" << ar[1] << "\t" << ar[2] << "\t" << ar[3] << "\t" << ar[4] << "\t" << ar[5] ;
		cout << "\t" << sigr[0] << "\t" << sigr[1] << "\t" << sigr[2] << "\t" << sigr[3] << "\t" << sigr[4] << "\t" << sigr[5] << endl;
	}
	else if (flag_print == 3 && N>=2)
	{ //print statistics, with off-diagonal elements
		double ar[6], sigr[6], sigrr[6];
		int nlive = calc_beam_stat(r, N, ar, sigr, sigrr);
		cout << std::scientific;
		cout << cnt_pass - 1 << "\t" << nlive << "\t" << ar[0] << "\t" << ar[1] << "\t" << ar[2] << "\t" << ar[3] << "\t" << ar[4] << "\t" << ar[5] ;
		cout <<  "\t" << sigr[0] << "\t" << sigr[1] << "\t" << sigr[2] << "\t" << sigr[3] << "\t" << sigr[4] << "\t" << sigr[5] ;
		cout <<  "\t" << sigrr[0] << "\t" << sigrr[1] << "\t" << sigrr[2] << "\t" << sigrr[3] << "\t" << sigrr[4] << "\t" << sigrr[5] << endl;
	}
			
}						
/**** BEND  **************************************/
void ABend::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<theta<<"\t"<<e1<<"\t"<<e2<<"\t"<<K<<endl;
}
ABend* ABend::create(istream& is, string& famname)
{
	string str;
	//is>>famname;
	getline(is, str);
	std::stringstream(str)>>len>>theta>>e1>>e2>>K;

	return new ABend(len,theta,e1,e2,K,famname);
}
void ABend::setfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("BendingAngle"))
		{ theta = *(double*)p; irho=theta/len;}
		else if (!pname.compare("K"))
		{ K = *(double*)p; }
		else if (!pname.compare("EntranceAngle"))
		{ e1 = *(double*)p; }
		else if (!pname.compare("ExitAngle"))
		{ e2 = *(double*)p; }
		
	}
void ABend::getfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname ; }
		else if (!pname.compare("BendingAngle"))
		{ *(double*)p = theta ; }
		else if (!pname.compare("K"))
		{ *(double*)p = K ; }
		else if (!pname.compare("EntranceAngle"))
		{ *(double*)p = e1 ; }
		else if (!pname.compare("ExitAngle"))
		{ *(double*)p = e2 ; }
		
	}

void ABend::pass(double *r, int N)
{
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
	{
		if(isnan(*(r+6*i)))
			continue;
			
		double *r6;
		r6 = r+6*i;
		//if(useFringe1)
		//     edge_fringe(r6, ba/le, entrance_angle, fint1, gap);
		//else
		edge(r6, irho, e1);		 
				  
		bend6(r6, len, theta, K,0);	    
		//bend6(r6, le, ba, grd, bye);
			    
		 //if(useFringe2)
		 //   edge_fringe(r6, ba/le, exit_angle,fint2,gap);
		 //else
		 edge(r6, irho, e2);
		
	}		
}

/**** BEND Symplectic Integrator  **************************************/
void ABendMPole::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<theta<<"\t"<<e1<<"\t"<<e2<<"\t"<<MaxNPole<<endl;
	os<<"\t\t\t";
	for(int i=0;i<MaxNPole;i++) os<<polyB[i]<<"\t";
	os<<endl;
}
ABendMPole* ABendMPole::create(istream& is, string& famname)
{
	string str;
	double polyB[MAX_MPOLE];
	//is>>famname;
	getline(is, str);
	std::stringstream(str)>>len>>theta>>e1>>e2>>MaxNPole;
	//getline(is, str);
	//cout<<str<<endl;
	//for(int i=0;i<MAX_MPOLE;i++) std::stringstream(str)>>polyB[i];
	for(int i=0;i<MaxNPole;i++) is>>polyB[i];
	getline(is, str);
	return new ABendMPole(len,theta,e1,e2,polyB,MaxNPole,famname);
}
void ABendMPole::setfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("BendingAngle"))
		{ theta = *(double*)p; irho=theta/len;}
		else if (!pname.compare("EntranceAngle"))
		{ e1 = *(double*)p; }
		else if (!pname.compare("ExitAngle"))
		{ e2 = *(double*)p; }
		else if (!pname.compare("energy"))
		{ energy = *(double*)p; }
		else if (!pname.compare("NSlice"))
		{ nslice = *(int*)p; }
		else if (!pname.compare("MaxNPole"))
		{ MaxNPole = *(int*)p; }
		else if (!pname.compare("PolynomB"))
		{ for(int i=0;i<MaxNPole;i++) polyB[i] = *( ((double*)p)+i ); }
		else
		{ cout<<"ABendMPole::setfield() --> no matched field "<<pname<<endl; }
		
	}
void ABendMPole::getfield(void *p, const string pname) 
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname ; }
		else if (!pname.compare("BendingAngle"))
		{ *(double*)p = theta ; }
		else if (!pname.compare("EntranceAngle"))
		{ *(double*)p = e1 ; }
		else if (!pname.compare("ExitAngle"))
		{ *(double*)p = e2 ; }
		else if (!pname.compare("energy"))
		{ *(double*)p = energy ; }
		else if (!pname.compare("NSlice"))
		{ *(int*)p = nslice ; }
		else if (!pname.compare("MaxNPole"))
		{ *(int*)p = MaxNPole ; }
		else if (!pname.compare("PolynomB"))
		{ for(int i=0;i<MaxNPole;i++) *( ((double*)p)+i ) = polyB[i]; }
		else
		{ cout<<"ABendMPole::getfield() --> no matched field "<<pname<<endl; }
		
	}
void ABendMPole::setfield(void *p,const string pname, unsigned int index)
	{
		if (!pname.compare("PolynomB"))
		{ if(index<MaxNPole)
				polyB[index] = *( (double*)p); 
			else
				cout<<"error: index="<<index<<">MaxNPole="<<MaxNPole<<endl;
		}
		//else if (!pname.compare("PolynomA"))
		//{ polyA[index] = *( (double*)p); }
		else
		{ cout<<"ABendMPole::setfield() --> no matched field "<<pname<<endl; }
	}
void ABendMPole::getfield(void *p,const string pname, unsigned int index)
	{
		if (!pname.compare("PolynomB"))
		{  
			if(index<MaxNPole)
				*( (double*)p) = polyB[index] ; 
			else
				cout<<"error: index="<<index<<">MaxNPole="<<MaxNPole<<endl;
		}
		//else if (!pname.compare("PolynomA"))
		//{ *( (double*)p) = polyA[index] ; }
		else
		{ cout<<"ABendMPole::getfield() with index --> no matched field "<<pname<<endl; }
	}
void ABendMPole::pass(double *r, int N)
{
	double SL,L1, L2,K1,K2;
	SL = len/nslice;
	L1 = SL*DRIFT1;
	L2 = SL*DRIFT2;
	K1 = SL*KICK1;
	K2 = SL*KICK2;
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
	{
		if(isnan(*(r+6*i)))
			continue;
			
		double *r6;
		r6 = r+6*i;
		edge(r6, irho, e1);	
		/* integrator */
		for(int m=0; m < nslice; m++) /* Loop over slices*/			
			{		
				if(flag_rad_on<=1){ //fourth order integrator, w/ or w/o radiation damping
					ATbendhxdrift6(r6,L1,irho);
					if(flag_rad_on)
						bndthinkickrad(r6, polyB, K1, irho, energy, MaxNPole-1,flag_rad_on);
					else
          	bndthinkick(r6, polyB, K1, irho, MaxNPole-1);

					ATbendhxdrift6(r6,L2,irho);
					if(flag_rad_on)
						bndthinkickrad(r6, polyB, K2, irho, energy, MaxNPole-1,flag_rad_on);
					else
          	bndthinkick(r6, polyB, K2, irho, MaxNPole-1);
					ATbendhxdrift6(r6,L2,irho);

					if(flag_rad_on)
						bndthinkickrad(r6, polyB, K1, irho, energy, MaxNPole-1,flag_rad_on);
					else
						bndthinkick(r6, polyB,  K1, irho, MaxNPole-1);
					ATbendhxdrift6(r6,L1,irho);
				}
				else if(flag_rad_on==2) {//second order integrator, w/ radiation damping and quantum excitation
					ATbendhxdrift6(r6,0.5*SL,irho);
					double *r = r6;
					//cout<<r[0]<<" "<<r[1]<<" "<<r[2]<<" "<<r[3]<<" "<<r[4]<<" "<<r[5]<<endl;
					bndthinkickrad(r6, polyB, SL, irho, energy, MaxNPole-1,flag_rad_on);
					//cout<<r[0]<<" "<<r[1]<<" "<<r[2]<<" "<<r[3]<<" "<<r[4]<<" "<<r[5]<<endl;
					ATbendhxdrift6(r6,0.5*SL,irho);
				}
			}  
		edge(r6, irho, e2);
		
		}
}

/**** Straight Symplectic Integrator  **************************************/
void AStrMPole::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<MaxNPole<<endl;
	os<<"\t\t\t";
	for(int i=0;i<MaxNPole;i++) os<<polyB[i]<<"\t";
	os<<endl;
	os<<"\t\t\t";
	for(int i=0;i<MaxNPole;i++) os<<polyA[i]<<"\t";
	os<<endl;
}
AStrMPole* AStrMPole::create(istream& is, string& famname)
{
	string str;
	double polyB[MAX_MPOLE],polyA[MAX_MPOLE];
	//is>>famname;
	getline(is, str);
	//cout<<str;
	std::stringstream(str)>>len>>MaxNPole;
	//getline(is, str);
	//for(int i=0;i<MAX_MPOLE;i++) std::stringstream(str)>>polyB[i];
	//getline(is, str);
	//for(int i=0;i<MAX_MPOLE;i++) std::stringstream(str)>>polyA[i];
	for(int i=0;i<MaxNPole;i++) is>>polyB[i];
	getline(is, str);
	for(int i=0;i<MaxNPole;i++) is>>polyA[i];
	getline(is, str);

	return new AStrMPole(len,polyB,polyA,MaxNPole,famname);
}
void AStrMPole::setfield(void *p, const string pname)
	{
		if (!pname.compare("NSlice"))
		{ nslice = *(int*)p; }
		else if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("MaxNPole"))
		{ MaxNPole = *(int*)p; }
		else if (!pname.compare("PolynomB"))
		{ for(int i=0;i<MAX_MPOLE;i++) polyB[i] = *( ((double*)p)+i ); }
		else if (!pname.compare("PolynomA"))
		{ for(int i=0;i<MAX_MPOLE;i++) polyA[i] = *( ((double*)p)+i ); }
		else
		{ cout<<"AStrMPole::setfield() --> no matched field "<<pname<<endl; }
		
	}
void AStrMPole::getfield(void *p, const string pname)
	{
		if (!pname.compare("NSlice"))
		{ *(int*)p = nslice ; }
		else if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname ; }
		else if (!pname.compare("MaxNPole"))
		{ *(int*)p = MaxNPole ; }
		else if (!pname.compare("PolynomB"))
		{ for(int i=0;i<MaxNPole;i++) *( ((double*)p)+i ) = polyB[i]; }
		else if (!pname.compare("PolynomA"))
		{ for(int i=0;i<MaxNPole;i++) *( ((double*)p)+i ) = polyA[i] ; }
		else
		{ cout<<"AStrMPole::getfield() --> no matched field "<<pname<<endl; }
		
	}
void AStrMPole::setfield(void *p,const string pname, unsigned int index)
	{
		if(index>MaxNPole){
				cout<<"error: index="<<index<<">MaxNPole="<<MaxNPole<<endl;
				return;
		}
				
		if (!pname.compare("PolynomB"))
			polyB[index] = *( (double*)p);
		else if (!pname.compare("PolynomA"))
			polyA[index] = *( (double*)p);
		
	}
void AStrMPole::getfield(void *p,const string pname, unsigned int index)
	{
		if(index>MaxNPole)
		{		cout<<"error: index="<<index<<">MaxNPole="<<MaxNPole<<endl;
			return;
		}
		if (!pname.compare("PolynomB"))
		{  *( (double*)p) = polyB[index] ; }
		else if (!pname.compare("PolynomA"))
		{ *( (double*)p) = polyA[index] ; }
		else
		{ cout<<"AStrMPole::getfield() with index --> no matched field "<<pname<<endl; }
	}
	
void AStrMPole::pass(double *r, int N)
{
	
	double SL,L1, L2,K1,K2;
	SL = len/nslice;
	L1 = SL*DRIFT1;
	L2 = SL*DRIFT2;
	K1 = SL*KICK1;
	K2 = SL*KICK2;
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
	{
		if(isnan(*(r+6*i)))
			continue;
			
		double *r6;
		r6 = r+6*i;
		
		double norm, NormL1, NormL2;

		if(len==0)
		{
			strthinkick(r6, polyA, polyB, 1.0, MAX_MPOLE-1);
		}
		else
		{
			/* integrator */
			for(int m=0; m < nslice; m++) /* Loop over slices*/			
				{		
						norm = 1/(1+r6[5]);
						NormL1 = L1*norm;
						NormL2 = L2*norm;
						fastdrift(r6, NormL1);
	    			strthinkick(r6, polyA, polyB,  K1, MAX_MPOLE-1);
	    			fastdrift(r6, NormL2);
	    			strthinkick(r6, polyA, polyB, K2, MAX_MPOLE-1);
	    			fastdrift(r6, NormL2);
	    			strthinkick(r6, polyA, polyB,  K1, MAX_MPOLE-1);
	    			fastdrift(r6, NormL1);	
				}  
			}
		
		}
}
/*********************************/
//Cavity
/*********************************/
void ACavity::serialize(ostream& os)
{
	std::streamsize ss = os.precision();
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<voltage<<"\t";
	os.precision(std::streamsize(16));
	os<<freq<<"\t"<<phi_s<<"\t";
	os.precision(ss);
	os<<energy<<endl;
}
ACavity* ACavity::create(istream& is, string& famname)
{
	string str;
	double volt, freq, phis,en;
	
	getline(is, str);
	std::stringstream(str)>>len>>volt>>freq>>phis>>en;
	
	return new ACavity(len,volt,freq,phis,en, famname);
}
void ACavity::setfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("Voltage"))
		{ voltage = *(double*)p; }
		else if (!pname.compare("Frequency"))
		{ freq = *(double*)p; }
		else if (!pname.compare("Phis"))
		{ phi_s = *(double*)p; }	
		else if (!pname.compare("Energy"))
		{ energy = *(double*)p; }	
		
	}
void ACavity::getfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("Voltage"))
		{ *(double*)p = voltage ; }
		else if (!pname.compare("Frequency"))
		{ *(double*)p = freq ; }
		else if (!pname.compare("Phis"))
		{ *(double*)p = phi_s ; }	
		else if (!pname.compare("Energy"))
		{ *(double*)p = energy ; }	
		
	}

void ACavity::pass(double *r, int N)
{
	double nv = voltage/energy;
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;
			//cout << "Hello from thread " << omp_get_thread_num() << endl;		
			double *r6;
			r6=r+i*6;
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
			if(flag_onoff) r6[5] += nv*sin(TWOPI*freq*r6[4]/C0+phi_s);	
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));	
		}		
	
}

/*********************************/
//Crab Cavity
/*********************************/
void ACrabCav::serialize(ostream& os)
{
	std::streamsize ss = os.precision();
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<devoltx<<"\t"<<devolty<<"\t";
	os.precision(16);
	os<<freq<<"\t"<<phi_s<<"\t";
	os.precision(ss);
	os<<energy<<"\t"<<sigphi<<"\t"<<sigdVV<<endl;
}
ACrabCav* ACrabCav::create(istream& is, string& famname)
{
	string str;
	double dvoltx, dvolty, freq, phis,en;
	
	getline(is, str);
	std::stringstream(str)>>len>>dvoltx>>dvolty>>freq>>phis>>en>>sigphi>>sigdVV;
	
	return new ACrabCav(len,dvoltx,dvolty,freq,phis,en, sigphi, sigdVV,famname);
}
void ACrabCav::setfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("DeVoltX"))
		{ devoltx = *(double*)p; }
		else if (!pname.compare("DeVoltY"))
		{ devolty = *(double*)p; }
		else if (!pname.compare("Frequency"))
		{ freq = *(double*)p; }
		else if (!pname.compare("Phis"))
		{ phi_s = *(double*)p; }	
		else if (!pname.compare("Energy"))
		{ energy = *(double*)p; }
		else if (!pname.compare("SigPhi"))
		{ sigphi = *(double*)p; }	
		else if (!pname.compare("SigdVV"))
		{ sigdVV = *(double*)p; }
		
	}
void ACrabCav::getfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("DeVoltX"))
		{ *(double*)p = devoltx ; }
		else if (!pname.compare("DeVoltY"))
		{ *(double*)p = devolty ; }
		else if (!pname.compare("Frequency"))
		{ *(double*)p = freq ; }
		else if (!pname.compare("Phis"))
		{ *(double*)p = phi_s ; }	
		else if (!pname.compare("Energy"))
		{ *(double*)p = energy ; }
		else if (!pname.compare("SigPhi"))
		{ *(double*)p = sigphi ; }
		else if (!pname.compare("SigdVV"))
		{ *(double*)p = sigdVV ; }
	}

void ACrabCav::pass(double *r, int N)
{
	double nvx = devoltx/energy;
	double nvy = devolty/energy;
	double k=TWOPI*freq/C0;
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;	
			
			double *r6;
			r6=r+i*6;
			double nvxa, nvya, omega_t;
			
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
			
			if (flag_onoff) {
				if (sigphi > 0)
				{
					omega_t = -k*r6[4] + phi_s + sigphi*grandn();
				}
				else
					omega_t = -k*r6[4] + phi_s;
				if(sigdVV > 0)
				{
					nvxa = nvx*(1.0 + sigdVV*grandn());
					nvya = nvy*(1.0 + sigdVV*grandn());
				}
				else
				{
					nvxa = nvx;
					nvya = nvy;
				}
				//cout<<nvya<<" "<<sin(omega_t)<<endl;

				r6[1] += nvxa*sin(omega_t);
				r6[3] += nvya*sin(omega_t);
				r6[5] += -(nvya*r6[2] + nvxa*r6[0])*k*cos(omega_t) / (1 + r6[5]);
			}
			if(len!=0)fastdrift(r6, 0.5*len/(1+r6[5]));	
		}		
	
}

/*********************************/
//RF Multipole
/*********************************/
void ARFMult::serialize(ostream& os)
{
	std::streamsize ss = os.precision();
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t";
	os.precision(16);
	os<<freq<<"\t"<<phi_s<<"\t";
	os.precision(ss);
	os<<energy<<"\t"<<sigphi<<"\t"<<sigdVV<<"\t"<<MaxNPole<<endl;
	os<<"\t\t\t";
	for(int i=0;i<MaxNPole;i++) os<<polyB[i]<<"\t";
	os<<endl;
	os<<"\t\t\t";
	for(int i=0;i<MaxNPole;i++) os<<polyA[i]<<"\t";
	os<<endl;
}
ARFMult* ARFMult::create(istream& is, string& famname)
{
	string str;
	double freq, phis,en;
	double polyB[MAX_MPOLE],polyA[MAX_MPOLE];
	
	getline(is, str);
	std::stringstream(str)>>len>>freq>>phis>>en>>sigphi>>sigdVV>>MaxNPole;
	for(int i=0;i<MaxNPole;i++) is>>polyB[i];
	getline(is, str);
	for(int i=0;i<MaxNPole;i++) is>>polyA[i];
	getline(is, str);

	return new ARFMult(len,freq,phis,en, polyB, polyA,MaxNPole,famname);
}
void ARFMult::setfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("Frequency"))
		{ freq = *(double*)p; }
		else if (!pname.compare("Phis"))
		{ phi_s = *(double*)p; }	
		else if (!pname.compare("Energy"))
		{ energy = *(double*)p; }
		else if (!pname.compare("SigPhi"))
		{ sigphi = *(double*)p; }	
		else if (!pname.compare("SigdVV"))
		{ sigdVV = *(double*)p; }
		else if (!pname.compare("MaxNPole"))
		{ MaxNPole = *(int*)p; }
		else if (!pname.compare("PolynomB"))
		{ for(int i=0;i<MaxNPole;i++) polyB[i] = *( ((double*)p)+i ); }
		else if (!pname.compare("PolynomA"))
		{ for(int i=0;i<MaxNPole;i++) polyA[i] = *( ((double*)p)+i ); }
		
	}
void ARFMult::getfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("Frequency"))
		{ *(double*)p = freq ; }
		else if (!pname.compare("Phis"))
		{ *(double*)p = phi_s ; }	
		else if (!pname.compare("Energy"))
		{ *(double*)p = energy ; }
		else if (!pname.compare("SigPhi"))
		{ *(double*)p = sigphi ; }
		else if (!pname.compare("SigdVV"))
		{ *(double*)p = sigdVV ; }
		else if (!pname.compare("MaxNPole"))
		{ *(int*)p = MaxNPole ; }
		else if (!pname.compare("PolynomB"))
		{ for(int i=0;i<MaxNPole;i++) *( ((double*)p)+i ) = polyB[i]; }
		else if (!pname.compare("PolynomA"))
		{ for(int i=0;i<MaxNPole;i++) *( ((double*)p)+i ) = polyA[i] ; }
	}
void ARFMult::setfield(void *p,const string pname, unsigned int index)
	{
		if(index>MaxNPole){
				cout<<"error: index="<<index<<">MaxNPole="<<MaxNPole<<endl;
				return;
		}
		if (!pname.compare("PolynomB"))
			polyB[index] = *( (double*)p);
		else if (!pname.compare("PolynomA"))
			polyA[index] = *( (double*)p);
		
	}
void ARFMult::getfield(void *p,const string pname, unsigned int index)
	{
		if(index>MaxNPole)
		{		cout<<"error: index="<<index<<">MaxNPole="<<MaxNPole<<endl;
			return;
		}
		if (!pname.compare("PolynomB"))
		{  *( (double*)p) = polyB[index] ; }
		else if (!pname.compare("PolynomA"))
		{ *( (double*)p) = polyA[index] ; }
		else
		{ cout<<"ARFMult::getfield() with index --> no matched field "<<pname<<endl; }
	}

void ARFMult::pass(double *r, int N)
{
	double k=TWOPI*freq/C0;
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;
			
			double omega_t;
			double rho, phis; //here phis=atan(y/x) as angle in cylindrical coordinate
			double AsL, ArhoL, AphiL;
			double nrho,nphi,tdArhoL; //temporary
				
			double *r6;
			r6=r+i*6;
			
			if(len!=0)
				fastdrift(r6, 0.5*len/(1+r6[5]));
			
			if(sigphi>0)
				omega_t = -k*r6[4]+phi_s+sigphi*grandn();
			else
				omega_t = -k*r6[4]+phi_s;
			
			rho = sqrt(SQR(r6[0])+SQR(r6[2]));
			phis = atan2(r6[2],r6[0]);
			AsL=0;
			ArhoL=0;
			AphiL=0;
			nrho = 1.0;
			for(int n=0;n<MaxNPole;n++) 
			{
				nphi = (n+1)*phis;
				tdArhoL = ( polyB[n]*cos(nphi)+polyA[n]*sin(nphi));
				ArhoL += nrho*tdArhoL;
				AphiL += nrho*(-polyB[n]*sin(nphi)+polyA[n]*cos(nphi));
				nrho *= rho;
				AsL += 1.0/(n+1)*nrho*tdArhoL;
			}
			//cout<<ArhoL<<" "<<AphiL<<" "<<AsL<<" "<<phis<<" "<<cos(omega_t)<<endl;
			
			
			r6[1] += (ArhoL*cos(phis)-AphiL*sin(phis))*cos(omega_t);
			r6[3] += (ArhoL*sin(phis)+AphiL*cos(phis))*cos(omega_t);
			r6[5] += k*AsL*sin(omega_t)/(1.+r6[5]);
			
			
			if(len!=0)
				fastdrift(r6, 0.5*len/(1+r6[5]));
		}		
	
}

/********************************/
//WigTable
/********************************/
AWigTable::AWigTable(double len, int nslice, unsigned int mx, unsigned int my, const string& fam) 
{
		this->len=len; famname=fam; type="WigTable";this->nslice=nslice; this->mx=mx; this->my=my; 
			for(int i=0;i<my;i++) 
			{
				xkick[i] = new double[mx]; ykick[i] = new double[mx];
				for(int j=0;j<mx;j++)
					{
						xkick[i][j] = 0.0;
						ykick[i][j] = 0.0;
					}
				}
}
			
AWigTable& AWigTable::operator=(AWigTable& a)
{
	len = a.len; famname=a.famname; type=a.type;
	nslice = a.nslice; 
			if(my>0)
			{ for(int i=0;i<my;i++) { delete []xkick[i];  delete []ykick[i]; }};
			
			mx = a.mx; my = a.my; 
			for(int i=0;i<mx;i++) 
				xl[i] = a.xl[i];
			for(int i=0;i<my;i++) 
				yl[i] = a.yl[i];	
			
			for(int i=0;i<my;i++) 
			{
				xkick[i] = new double[mx]; ykick[i] = new double[mx];
				for(int j=0;j<mx;j++)
					{
						xkick[i][j] = a.xkick[i][j];
						ykick[i][j] = a.ykick[i][j];
					}
				}
			
			return *this;
}

void AWigTable::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<nslice<<"\t"<<mx<<"\t"<<my<<endl;
	for(int i=0;i<mx;i++)
		os<<xl[i]<<"\t";
	os<<endl;
	for(int i=0;i<my;i++)
		os<<yl[i]<<"\t";
	os<<endl;
	for(int i=0;i<my;i++)
	{
		for(int j=0;j<mx;j++)
			os<<xkick[i][j]<<"\t";
		os<<endl;
	}
	os<<endl;
			
	for(int i=0;i<my;i++)
	{
		for(int j=0;j<mx;j++)
			os<<ykick[i][j]<<"\t";
		os<<endl;
	}
	
}	
void AWigTable::setfield(void *p, const string pname)
{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("nslice"))
		{ nslice = *(int*)p; }
}	
void AWigTable::getfield(void *p, const string pname)
{
	if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("nslice"))
		{ *(int*)p = nslice ; }
		else if (!pname.compare("mx"))
		{ *(unsigned int*)p = mx ; }
		else if (!pname.compare("my"))
		{ *(unsigned int*)p = my ; }
}	
AWigTable* AWigTable::create(istream& is, string& famname)
{
	string str;
	std::stringstream sstr;
	
	getline(is, str);
	std::stringstream(str)>>len>>nslice>>mx>>my;
	
	AWigTable* pa = new AWigTable(len,nslice,mx,my, famname);
	
	getline(is, str);
	sstr<<str;
	for(int i=0;i<mx;i++)
		sstr>>pa->xl[i];
		
	getline(is, str);
	sstr<<str;
	for(int i=0;i<my;i++)
		sstr>>pa->yl[i];
			

	for(int i=0;i<my;i++)
	{
		getline(is, str);
		sstr<<str;
		for(int j=0;j<mx;j++)
			sstr>>pa->xkick[i][j];
	}
	getline(is, str); //an empty line
	for(int i=0;i<my;i++)
	{
		getline(is, str);
		sstr<<str;
		for(int j=0;j<mx;j++)
			sstr>>pa->ykick[i][j];
	}

	return pa;
	
}

void AWigTable::pass(double *r, int N)
{

	if(len==0)
	{
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;
			
			double *r6=r+i*6;
			double x = r6[0];
			double y = r6[2];
			
			r6[1] += interp2(xkick, xl,mx, yl,my, x, y);
			r6[3] += interp2(ykick, xl,mx, yl,my, x, y);
		}		
	}
	else //with finite length
	{
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;	
			double *r6=r+i*6;
			
			for(int j=0;j<nslice;j++)
			{
				fastdrift(r6, 0.5*len/(1+r6[5])/nslice);
				double x = r6[0];
				double y = r6[2];
				r6[1] += interp2(xkick, xl,mx, yl,my, x, y)/nslice;
				r6[3] += interp2(ykick, xl,mx, yl,my, x, y)/nslice;
				fastdrift(r6, 0.5*len/(1+r6[5])/nslice);	
			}
		}		
		
	}				
		
}
/**************************************/
// AIdealWig
/**************************************/
AIdealWig::AIdealWig(double len, int nslice, double lambda, double h, unsigned int wig_type,const string& fam) 
{
		this->len=len; famname=fam; type="IdealWig";
		this->nslice=nslice; this->lambda=lambda; this->h=h; 
		energy=DEFAULTE0;	
		flag_rad_on = 0;
		flag_wigtype = wig_type;
	
}
			
AIdealWig& AIdealWig::operator=(AIdealWig& a)
{
	len = a.len; famname=a.famname; type=a.type;
	nslice = a.nslice; 
	lambda = a.lambda; 
	h = a.h;		
	energy=a.energy;
	flag_rad_on = a.flag_rad_on;
	flag_wigtype = a.flag_wigtype;
			
			return *this;
}

void AIdealWig::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<nslice<<"\t"<<lambda<<"\t"<<h<<flag_wigtype<<endl;
	
}	
void AIdealWig::setfield(void *p, const string pname)
{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("nslice"))
		{ nslice = *(int*)p; }
		else if (!pname.compare("lambda"))
		{ lambda = *(double*)p; }
		else if (!pname.compare("h"))
		{ h = *(double*)p; }
		else if (!pname.compare("energy"))
		{ energy = *(double*)p; }
		else if (!pname.compare("wigtype"))
		{ flag_wigtype = *(unsigned int*)p; }
		
}	
void AIdealWig::getfield(void *p, const string pname)
{
	if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("nslice"))
		{ *(int*)p = nslice ; }
		else if (!pname.compare("lambda"))
		{ *(double*)p = lambda ; }
		else if (!pname.compare("h"))
		{ *(double*)p = h ; }
		else if (!pname.compare("energy"))
		{ *(double*)p = energy ; }
		else if (!pname.compare("wigtype"))
		{ *(unsigned int*)p = flag_wigtype ; }
}	
AIdealWig* AIdealWig::create(istream& is, string& famname)
{
	string str;
	std::stringstream sstr;
	
	getline(is, str);
	std::stringstream(str)>>len>>nslice>>lambda>>h>>flag_wigtype;
	
	AIdealWig* pa = new AIdealWig(len,nslice,lambda,h, flag_wigtype,famname);
	
	return pa;
	
}

void AIdealWig::pass(double *r, int N)
{
	double k = 2.0*PI/lambda;
	
	if(len==0)
	{
		//nothing to be done
	}
	else 
	{
					#define PI  3.1415926535897932384626
					#define TWOPI  (2*PI) 
					#define CGAMMA 	8.846056192e-05 
		
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

		for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;	
			
			double *r6;
			double x, y,xp,yp;
			r6=r+i*6;
			
			for(int j=0;j<nslice;j++)
			{
				fastdrift(r6, 0.5*len/(1+r6[5])/nslice);
				y = r6[2];
				x = r6[0];

				if(flag_rad_on)
				{
					double CRAD = CGAMMA*energy*energy*energy/(TWOPI*1e9);	/* [m]/[GeV^3] M.Sands (4.1) */
					double dp = r6[5];
					
					yp = r6[3]/(1+dp);
					xp = r6[1]/(1+dp);
					r6[5] = r6[5] - CRAD*SQR(1+r6[5])*h*h*len/nslice/2.0;
					r6[1] = xp*(1+r6[5]);
					r6[3] = yp*(1+r6[5]);
				}
								
				if(flag_wigtype==IDEALWIG_TYPE_H){				
					double y3 = y*y*y;
					double k2 = k*k;
					r6[3] -= h*h*(0.5*y+1.0/3.0*y3*k2+1.0/15.0*y3*y*y*k2*k2)*len/nslice;
				}
				else if(flag_wigtype==IDEALWIG_TYPE_V){
					double x3 = x*x*x;
					double k2 = k*k;
					r6[1] -= h*h*(0.5*x+1.0/3.0*x3*k2+1.0/15.0*x3*x*x*k2*k2)*len/nslice;
				}
				
				if(flag_rad_on)
				{
					double CRAD = CGAMMA*energy*energy*energy/(TWOPI*1e9);	/* [m]/[GeV^3] M.Sands (4.1) */
					double dp = r6[5];
					
					yp = r6[3]/(1+dp);
					xp = r6[1]/(1+dp);
					r6[5] = r6[5] - CRAD*SQR(1+r6[5])*h*h*len/nslice/2.0;
					r6[1] = xp*(1+r6[5]);
					r6[3] = yp*(1+r6[5]);
				}
				
				fastdrift(r6, 0.5*len/(1+r6[5])/nslice);	
			}
		}		
		
	}				
		
}

/*********************************/
//TransWake: Resonant transverse wake model
/*********************************/
void ATransWake::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<energy<<"\t"<<Ne<<"\t"<<Wx<<"\t"<<fRx<<"\t"<<Qx<<"\t"<<Wy<<"\t"<<fRy<<"\t"<<Qy<<endl;
}
ATransWake* ATransWake::create(istream& is, string& famname)
{
	string str;
	double en;
	getline(is, str);
	std::stringstream(str)>>len>>en>>Ne>>Wx>>fRx>>Qx>>Wy>>fRy>>Qy;
	
	return new ATransWake(len,Ne,Wx,fRx,Qx, Wy, fRy, Qy, en, famname);
}
void ATransWake::setfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("Energy"))
		{ energy = *(double*)p; }
		else if (!pname.compare("Charge"))
		{ Ne = *(double*)p; }	
		else if (!pname.compare("WakeX"))
		{ Wx = *(double*)p; }
		else if (!pname.compare("FreqX"))
		{ fRx = *(double*)p; }
		else if (!pname.compare("Qx"))
		{ Qx = *(double*)p; }
		else if (!pname.compare("WakeY"))
		{ Wy = *(double*)p; }
		else if (!pname.compare("FreqY"))
		{ fRy = *(double*)p; }
		else if (!pname.compare("Qy"))
		{ Qy = *(double*)p; }	
		
	}
void ATransWake::getfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("Energy"))
		{ *(double*)p = energy ; }
		else if (!pname.compare("Charge"))
		{ *(double*)p = Ne ; }	
		else if (!pname.compare("WakeX"))
		{ *(double*)p = Wx ; }
		else if (!pname.compare("FreqX"))
		{ *(double*)p = fRx ; }
		else if (!pname.compare("Qx"))
		{ *(double*)p = Qx ; }
		else if (!pname.compare("WakeY"))
		{ *(double*)p = Wy ; }
		else if (!pname.compare("FreqY"))
		{ *(double*)p = fRy ; }
		else if (!pname.compare("Qy"))
		{ *(double*)p = Qy ; }
	}

void ATransWake::pass(double *r, int N)
{
	double q = Ne/N*1000;  //charge per macro-particle, converted to pC
	double kRx, kRy,ax,ay,kbarRx,kbarRy;
	
	kRx = 2*PI*fRx*1.0e9/C0; //1/m, k=omega/c
  kRy = 2*PI*fRy*1.0e9/C0; //1/m,
    
  ax = kRx/2/Qx;
  ay = kRy/2/Qy;
  kbarRx = kRx*sqrt(1.0-0.25/SQR(Qx));
  kbarRy = kRy*sqrt(1.0-0.25/SQR(Qy));
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;	
			double *r6;
			r6=r+i*6;
			double deltaxp, deltayp;
			
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
			
			if(flag_onoff){
				deltaxp=0;
				deltayp=0;
			
				for(int j=0;j<N;j++)	
				{
					if(isnan(*(r+6*j)))
						continue;
					if (i == j) continue;
					double *rj6;
				
					rj6=r+j*6;
					double zij = r6[4]-rj6[4];
					if(zij>0) continue;
				
					deltaxp = deltaxp -q*Wx*exp(ax*zij)*sin(kbarRx*zij)*rj6[0]/energy/1.0e6;
					deltayp = deltayp -q*Wy*exp(ay*zij)*sin(kbarRy*zij)*rj6[2]/energy/1.0e6;
				}
				//cout << q << " " << Wy << " " << kRy << " " << ay << " " <<kbarRy<<" "<< deltayp << endl;
				r6[1] += deltaxp*(1+r6[5]);
				r6[3] += deltayp*(1+r6[5]);			
			}
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
		}
		
		
	
}

/*********************************/
//Matrix: matrix element
/*********************************/
void AMatrix::serialize(ostream& os)
{
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<endl;
	for(int i=0;i<6;i++)
	{	
		os<<"\t";
		for(int j=0;j<6;j++)
			os<<mat[i][j]<<"\t";
		os<<endl;
	}
}
AMatrix* AMatrix::create(istream& is, string& famname)
{
	string str;
	double mat[6][6];
	getline(is, str);
	std::stringstream(str)>>len;
	for(int i=0;i<6;i++)
	{
		getline(is, str);	
		//cout<<str<<endl;
		std::stringstream(str)>>mat[i][0]>>mat[i][1]>>mat[i][2]>>mat[i][3]>>mat[i][4]>>mat[i][5];
	}
	
	return new AMatrix(len,mat, famname);
}

void AMatrix::setfield(void *p, int i, int j)
{
	this->mat[i][j] = *(double*)p;
}
		
void AMatrix::getfield(void *p, int i, int j)
{
	*(double*)p = this->mat[i][j];
}

void AMatrix::pass(double *r, int N)
{
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;	
			double *r6;
			r6=r+i*6;
			
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
			
			//Matrix multiplication
			double tr[6];
			for(int j=0;j<6;j++)	
				{
					tr[j] = mat[j][0]*r6[0]+mat[j][1]*r6[1]+mat[j][2]*r6[2]+mat[j][3]*r6[3]+mat[j][4]*r6[4]
					+mat[j][5]*r6[5];
					
				}
			for(int j=0;j<6;j++)	
				r6[j] = tr[j];
					
			
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
		}
	
}

/*********************************/
//SyncRad: Synchrotron radiation lump element
/*********************************/
void ASyncRad::serialize(ostream& os)
{
	//this->len = len; this->emitx = emitx; this->emity = emity; this->sigdpp = sigdpp; 
	//		this->U0 = U0; this->Jx = Jx;this->Jy = Jy; this->betax = betax;this->betay = betay; this->energy = en;
			
	os<<famname<<"  : \t" <<type<<"\t"<<len<<"\t"<<energy<<"\t"<<U0<<"\t"<<Jx<<"\t"<<Jy<<"\t"<<endl;
	os<<"\t"<<emitx<<"\t"<<emity<<"\t"<<sigdpp<<"\t"<<betax<<"\t"<<betay<<"\t"<<endl;
}
ASyncRad* ASyncRad::create(istream& is, string& famname)
{
	string str;
	double en,U0, emitx,emity, sigdpp,Jx,Jy,betax,betay;
	getline(is, str);
	std::stringstream(str)>>len>>en>>U0>>Jx>>Jy;
	getline(is, str);
	std::stringstream(str)>>emitx>>emity>>sigdpp>>betax>>betay;
	
	return new ASyncRad(len,emitx,emity,sigdpp,U0,Jx,Jy,betax,betay, en, famname);
}
void ASyncRad::setfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ len = *(double*)p; }
		else if (!pname.compare("FamName"))
		{ famname = *(string*)p; }
		else if (!pname.compare("Energy"))
		{ energy = *(double*)p; }
		else if (!pname.compare("U0"))
		{ U0 = *(double*)p; }	
		else if (!pname.compare("EmitX"))
		{ emitx = *(double*)p; }
		else if (!pname.compare("EmitY"))
		{ emity = *(double*)p; }
		else if (!pname.compare("SigDelta"))
		{ sigdpp = *(double*)p; }
		else if (!pname.compare("Jx"))
		{ Jx = *(double*)p; }
		else if (!pname.compare("Jy"))
		{ Jy = *(double*)p; }
		else if (!pname.compare("BetaX"))
		{ betax = *(double*)p; }
		else if (!pname.compare("BetaY"))
		{ betay = *(double*)p; }		
		
	}
void ASyncRad::getfield(void *p, const string pname)
	{
		if (!pname.compare("Length"))
		{ *(double*)p = len ; }
		else if (!pname.compare("FamName"))
		{ *(string*)p = famname; }
		else if (!pname.compare("Energy"))
		{ *(double*)p = energy ; }
		else if (!pname.compare("U0"))
		{ *(double*)p = U0 ; }	
		else if (!pname.compare("EmitX"))
		{ *(double*)p = emitx ; }
		else if (!pname.compare("EmitY"))
		{ *(double*)p = emity ; }
		else if (!pname.compare("SigDelta"))
		{ *(double*)p = sigdpp ; }
		else if (!pname.compare("Jx"))
		{ *(double*)p = Jx; }
		else if (!pname.compare("Jy"))
		{ *(double*)p = Jy; }
		else if (!pname.compare("BetaX"))
		{ *(double*)p = betax; }
		else if (!pname.compare("BetaY"))
		{ *(double*)p = betay; }
	}

void ASyncRad::pass(double *r, int N)
{
	double alfx, alfy, alfz;
	double sx, sy,sz;
	double sigx, sigy;
	alfx = U0/2./energy*Jx;
	alfy = U0/2./energy*Jy;
	alfz = U0/2./energy*Jz;
	sx = sqrt(2.*alfx); 
	sy = sqrt(2.*alfy);
	sz = sqrt(2.*alfz);
	
	sigx = sqrt(emitx*betax);
	sigy = sqrt(emity*betay);
	
#ifdef		OPENMP_BY_ELEMENT
#pragma omp parallel for
#endif

	for(int i=0;i<N;i++)
		{
			if(isnan(*(r+6*i)))
				continue;	
			double *r6;
			r6=r+i*6;
			
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
			
			//radiation damping
			if(flag_damping_onoff)
			{
				r6[0] *= (1-alfx);
				r6[2] *= (1-alfy);
				r6[5] *= (1-alfz);
			}
			//cout<<flag_damping_onoff<<"\t"<<flag_excitation_onoff<<"\t"<<endl;
			//quantum excitation
			if(flag_excitation_onoff)
			{
				r6[0] += sigx*sx*grandn();
				r6[2] += sigy*sy*grandn();
				r6[5] += sigdpp*sz*grandn();
			}
					
			
			if(len!=0) fastdrift(r6, 0.5*len/(1+r6[5]));
		}
	
}
