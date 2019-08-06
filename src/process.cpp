#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <ctime>
#include "matvec.h"
#include "accel.h"
#include "utility.h"
#include <omp.h>
#include <algorithm> 
#include <functional>
#include <cctype>
#include <vector>

//using namespace std;

#define CHARCOMMENT ("#%")
#define STRINGTOUPPER(s)  (std::transform(s.begin(), s.end(), s.begin(), std::ptr_fun<int, int>(std::toupper)))
#define STR_LATTICE (string("LATTICE"))
#define STR_OPTICS (string("OPTICS"))
#define STR_SURVEY (string("SURVEY"))
#define STR_DAMA (string("DAMA"))
#define STR_TRACK (string("TRACK"))
#define BLOCK_END (string("END"))

int lookforblock(ifstream& fin, string &smatch,string key[],int ns=1);
int proc_lattice(ifstream& fin, ALine **pLine);
int proc_optics(ifstream& fin, ALine **pLine);
int proc_track(ifstream& fin, ALine **pLine);
unsigned int sel_elements(unsigned int* plist,ALine *pLine,int nElem, string& sline);
void str_remove_comments(string &sline);
bool apply_line_setpt(const string &filename, ALine *pLine);
int proc_DAMA(ifstream& fin, ALine *pLine);
int get_type_famname_str(string &sline, string &type, string &famname);
unsigned int get_elem_list(unsigned int* plist, ALine *pLine, int nElem, string& type, string &famname, unsigned index, string str_indices=string(""));
string& mgetline(ifstream& fin, string &sline);
void print_m66_to_line(Mat& m66, string &str_linestyle);

int main(int argc, char* argv[]){

srand (time(NULL));
AccelElement::initializemap();

//initially set maximum thread to 1
omp_set_num_threads(1);
ALine *aLine = NULL;
ifstream fin(argv[1]);

string key_blocks[]={STR_LATTICE,STR_OPTICS,STR_SURVEY,STR_DAMA,STR_TRACK};
int n_keys = 5;
string sline;
int sw_indx = lookforblock(fin,sline, key_blocks,n_keys);
if(sw_indx>0)
{	
	cout<<"The first block must be LATTICE."<<endl;
	return 1;
}
else if(sw_indx<0)
{
	cout<<"No LATTICE block is found."<<endl;
	return 1;
}

cout << "#" << currentDateTime() << endl;
proc_lattice(fin, &aLine);
cout<<"Total elements: "<<aLine->getnumber()<<endl;
while(!fin.eof())
{
	string smatch;
	int sw_indx = lookforblock(fin, smatch, key_blocks,n_keys);
	if(sw_indx<0) break; //end of file or stop/return
	cout<<key_blocks[sw_indx]<<endl; 
	
	switch(sw_indx)
	{
		case 0: //lattice
			proc_lattice(fin, &aLine);
			//cout<<sw_indx<<endl;
			break;
		case 1:  //optics
			//cout<<sw_indx<<endl;
			
			proc_optics(fin, &aLine);
			break;
		case 2: //survey
			//cout<<sw_indx<<endl;
			break;
		case 3:  //DAMA
			//cout<<sw_indx<<endl;
			proc_DAMA(fin, aLine);
			break;
		case 4: //Track	
			//cout<<sw_indx<<endl;
			cout << "#" << currentDateTime() << endl;
			proc_track(fin, &aLine);
			cout << "#" << currentDateTime() << endl;
			break;
				
	}

}

fin.close();
delete aLine;
return 0;
}

void str_remove_comments(string &sline)
{
	size_t pos = sline.find_first_of(CHARCOMMENT);
	if(pos!=string::npos)
		sline.erase(pos);
			
}
#define SYMBOL_LINE_CONTD ('\\')
string& mgetline(ifstream& fin, string &sline)
{
	getline(fin, sline);
	//return sline;
	
	str_remove_comments(sline);
	str_trim_right(sline);
	while(!fin.eof() && (sline.back()==SYMBOL_LINE_CONTD) && sline.length()>0)
	{
		sline.back() = ' ';
		string nl;
		getline(fin,nl);
		str_remove_comments(nl);
		str_trim_right(nl);	
		sline = sline+nl;
	}
	return sline;
	
}
int lookforblock(ifstream& fin, string& smatch,string key[],int ns)
{
	while(!fin.eof())
	{
		string sline;
		mgetline(fin, sline);
		
		str_remove_comments(sline);
		str_trim_both(sline);
		smatch = sline;
		
		STRINGTOUPPER(sline);
		for(int i=0;i<ns;i++)
			if(sline.compare(0,key[i].length(),key[i])==0)
			{
					//cout<<sline<<endl;
					return i;
			}
			else if(str_icompare(sline,"STOP") || str_icompare(sline,"RETURN"))
			{ return -1;}
	}
	smatch = string("");
	return -1;
	
}

#define STR_LAT_FILENAME (string("FILE_NAME"))
#define STR_LAT_RESETFIRST (string("RESET_FIRST"))
#define STR_LAT_PRINTSUMMARY (string("PRINT_LAT_SUMMARY"))
#define STR_LAT_PRINT_ELEM_VALUE (string("PRINT_ELEM_VALUE"))
#define STR_LAT_EXPORT_ELEMENT (string("EXPORT_ELEMENT"))
#define STR_LAT_CORRECT_TUNE (string("CORRECT_TUNE"))
#define STR_LAT_CORRECT_CHROM (string("CORRECT_CHROM"))
#define STR_LAT_MAX_THREADS (string("SET_MAX_NUM_THREADS"))
#define STR_LAT_MODIFY_ELEM (string("MODIFY_ELEMENT"))
#define STR_LAT_APPLY_SETPT (string("APPLY_SETPT"))
#define STR_LAT_SET_FLAGS (string("SET_FLAGS"))
int proc_lattice(ifstream& fin, ALine *(*pLine))
{
	string key_blocks[]={STR_LAT_FILENAME,STR_LAT_RESETFIRST,STR_LAT_PRINT_ELEM_VALUE,STR_LAT_PRINTSUMMARY,STR_LAT_EXPORT_ELEMENT,
		STR_LAT_CORRECT_TUNE, STR_LAT_CORRECT_CHROM,STR_LAT_MAX_THREADS,STR_LAT_MODIFY_ELEM,
		STR_LAT_SET_FLAGS,STR_LAT_APPLY_SETPT,BLOCK_END};
	int n_key=12;
	
	string sline;
	while(!fin.eof())
	{
		int sw_indx=lookforblock(fin, sline,key_blocks,n_key);
		//cout<<sline<<" : "<<sw_indx<<endl;
		if(sw_indx!=0 && (*pLine)==NULL)
		{
			cout<<sw_indx<<" Lattice file name needs to be given first"<<endl;
			return -1;
		}
		if(sw_indx==n_key-1) break; //end of block 
		
		switch(sw_indx)
		{
			case 0: //filename
			{	
				if((*pLine)!=NULL)
					{
						delete (*pLine);
					}
				
				str_remove_comments(sline);
				replace(sline.begin(), sline.end(), ',', ' ');
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0, pos + 1);
				str_trim_left(sline);
				pos = sline.find_first_of((" \t"));
				string filename = sline.substr(0,pos);
				//str_trim_both(filename);
				//cout<<filename<<endl;
				sline.erase(0, pos + 1);
				str_trim_left(sline);
				
				ifstream fin_lat;
				try {
					fin_lat.open(filename.c_str());
					(*pLine)=ALine::deserialize(fin_lat);
				}
				catch (const ifstream::failure& e)
				{ cout << "Exception opening/reading file:"<<filename<<endl;}

				

				string keyws[] = { "TYPE"};
				string type = "";

				int cnt = 0;
				while (sline.length()>0 && cnt++<1)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0, pos);
					sline.erase(0, pos + 1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if (valtype == "TYPE")
					{
						type = sline.substr(0, pos);
						if (str_icompare(type, "ERING")) //electron storage ring
						{
							double phis = (*pLine)->setrfphase();
							cout << "RF cavities phase set to: phi_s=" << phis << " rad" << endl;
							
							(*pLine)->set_rad_onoff(0);
							cout << "radiation off" << endl;
						}
					}
					sline.erase(0, pos + 1);
					str_trim_both(sline);
				}
				
				int ncav = (*pLine)->set_cavity_onoff(0);
				if(ncav>0)
					cout<<ncav<<" cavity(ies) turned off"<<endl;
				int nwake = (*pLine)->set_wake_onoff(0);
				if(nwake>0)
						cout << nwake << " wake element(s) turned off" << endl;
				int ncrab = (*pLine)->set_crab_onoff(false);
				if(ncrab>0)
					cout<<ncrab<<" crab cavity(ies) turned off"<<endl;
			}
				
				break;
				
			case 1: //reset first element
				{
				str_remove_comments(sline);
					/*
					size_t pos = sline.find_first_of(("=:"));
					string s_elem = sline.substr(pos+1);
					str_trim_both(s_elem);
					pos = s_elem.find_first_of(("."));
					size_t pos2 = s_elem.find_first_of((" \t\n\r\f\v"));
					
					string type = s_elem.substr(0,pos);
					string name = s_elem.substr(pos+1,pos2-1-pos);
					
				  int eid = str2int(s_elem.substr(pos2+1));
				  if(eid==0)eid=1; //default to 1
				  //cout<<type<<" "<<name<<" "<<eid<<endl;
				  (*pLine)->resetfirstelement(type,name,eid);
				  */
				  //cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0, pos + 1);
				str_trim_both(sline);
				replace(sline.begin(), sline.end(), ',', ' ');

				string keyws[] = { "TYPE","FAMNAME","INDEX" };
				string type = "", famname = ""; 
				unsigned int index = 1;

				int cnt = 0;
				while (sline.length()>0 && cnt++<3)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0, pos);
					sline.erase(0, pos + 1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if (valtype == "TYPE")
					{
						type = sline.substr(0, pos);
					}
					else if (valtype == "FAMNAME")
					{
						famname = sline.substr(0, pos);
					}
					else if (valtype == "INDEX")
						std::stringstream(sline) >> index;

					sline.erase(0, pos + 1);
					str_trim_both(sline);
				}
				if(!type.empty() && !famname.empty())
				{
					cout << "First element moved to: type="<<type << ", famname=" << famname << ", index= " << index << endl;
					(*pLine)->resetfirstelement(type,famname,index);
				}
				else if (type.empty() && !famname.empty())
				{
					unsigned int nElem = (*pLine)->getnumber();
					for (unsigned i = 0; i < nElem; i++)
					{
						AccelElement* pc = (AccelElement*)(*pLine)->getelement(i);
						if(pc->getfamname()==famname)
						{
							type = pc->gettype();
							cout << "First element moved to: type=" << type << ", famname=" << famname << ", index= " << 1 << endl;
							(*pLine)->resetfirstelement(type, famname, 1);
						}
							
						
					}
				}



			
				}
				break;
			
			case 2: //print element value
				{
					str_remove_comments(sline);
					
					size_t pos = sline.find_first_of(("=:"));
					//cout<<sline<<" "<<pos<<endl;
					sline.erase(0,pos+1);
					//cout<<sline<<endl;
					replace(sline.begin(),sline.end(),',',' ');
					str_trim_both(sline);
					
					string type, famname;
					unsigned int nlist;
					int nElem = (*pLine)->getnumber();
					//cout<<nElem<<endl;
					unsigned int *plist = new unsigned int[nElem];
					//cout<<sline<<endl;
					
					int flag = get_type_famname_str(sline, type,famname);
					//cout<<sline<<endl;
					//cout<<type<<" "<<famname<<" "<<flag<<endl;
					if(flag<1)
						break;
					else if(flag==1 || flag==2)
					{
						nlist = (*pLine)->findelements(plist, type);
					}
					else if(flag==11)
					{
						nlist = (*pLine)->findelements(plist, type, famname);
					}
					cout<<"$$$Print-Element_Values***********************"<<endl;
					cout<<"#"<<type<<"\t"<<famname<<"\t"<<sline<<endl;
					for(int i=0;i<nlist;i++)
					{
						if(type=="StrMPole")
						{
							AStrMPole* pE=(AStrMPole*) (*pLine)->getelement(plist[i]);
							
							size_t pos = sline.find("Index");
							if(pos!=string::npos)
							{
								cout<<plist[i]<<"\t";
							}
							pos = sline.find("Length");
							if(pos!=string::npos)
							{
								double len;
								pE->getfield((void*)&len,"Length");
									cout<<len<<"\t";
							}
							pos = sline.find("PolyB[1]");
							if(pos!=string::npos)
							{
								double K1;
								pE->getfield((void*)&K1,"PolynomB",1);
									cout<<K1<<"\t";
							}
							pos = sline.find("PolyB[2]");
							if(pos!=string::npos)
							{
								double K2;
								pE->getfield((void*)&K2,"PolynomB",2);
									cout<<K2<<"\t";
							}
							pos = sline.find("PolyB[3]");
							if(pos!=string::npos)
							{
								double K3;
								pE->getfield((void*)&K3,"PolynomB",3);
									cout<<K3<<"\t";
							}
							
							pos = sline.find("PolyA[1]");
							if(pos!=string::npos)
							{
								double K1;
								pE->getfield((void*)&K1,"PolynomA",1);
									cout<<K1<<"\t";
							}
							pos = sline.find("PolyA[2]");
							if(pos!=string::npos)
							{
								double K2;
								pE->getfield((void*)&K2,"PolynomA",2);
									cout<<K2<<"\t";
							}
							pos = sline.find("PolyA[3]");
							if(pos!=string::npos)
							{
								double K3;
								pE->getfield((void*)&K3,"PolynomA",3);
									cout<<K3<<"\t";
							}
								
						}
						else if(type=="Quad")
						{
							AQuad* pE=(AQuad*) (*pLine)->getelement(plist[i]);
							
							size_t pos = sline.find("Index");
							if(pos!=string::npos)
							{
								cout<<plist[i]<<"\t";
							}
							pos = sline.find("Length");
							if(pos!=string::npos)
							{
								double len;
								pE->getfield((void*)&len,"Length");
									cout<<len<<"\t";
							}
							pos = sline.find("K");
							if(pos!=string::npos)
							{
								double K1;
								pE->getfield((void*)&K1,"K");
									cout<<K1<<"\t";
							}
						}
						cout<<endl;
					}
					cout<<"******************************************"<<endl;
					
					delete []plist;
				}
			
				break;	
			case 3: //print ring summary
				{
				cout << "$$$Lattice-Summary***********************" << endl;
					cout << "Total " << (*pLine)->getnumber() << " elements, total length = " << (*pLine)->getlength() << " m" << endl;
					cout << "First element is " << (*pLine)->getelement(0)->gettype() << "."
						<< (*pLine)->getelement(0)->getfamname() << endl;

					(*pLine)->set_rad_onoff(0);
					double tunex, tuney;
					double chromx, chromy;

					(*pLine)->calctune(tunex, tuney);
					cout << "Tunes: " << tunex << "\t" << tuney << "\t" << endl;

					(*pLine)->calcchrom(chromx, chromy);
					cout << "Chromaticity: " << chromx << "\t" << chromy << "\t" << endl;
					//cout<<sline<<endl;

					Vec co(6);
					(*pLine)->getco6(co);
					cout<<"Closed orbit 6D "<<endl;
					//cout<<co<<endl;
					for(int i=0;i<6;i++) cout<<co[i+1]<<"\t";
					cout<<endl;
					
					Twiss tw;
					(*pLine)->calctwiss(tw);
					cout << "betax\tbetay\talphax\talphay" << endl;
					cout << tw.betx << " \t" << tw.bety << " \t" << tw.alfx << " \t" << tw.alfy << endl;
					cout << "******************************************" << endl;
				}
				break;
			case 4: //Export element
			{
				str_remove_comments(sline);

				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0, pos + 1);
				replace(sline.begin(), sline.end(), ',', ' ');
				str_trim_both(sline);

				string type, famname;
				unsigned int nlist;
				unsigned int nElem = (*pLine)->getnumber();

				unsigned int *plist = new unsigned int[nElem];

				int flag = get_type_famname_str(sline, type, famname);
				//cout<<sline<<endl;
				//cout<<type<<" "<<famname<<" "<<flag<<endl;
				if (flag<1)
					break;
				else if (flag == 1 || flag == 2)
				{
					nlist = (*pLine)->findelements(plist, type);
				}
				else if (flag == 11)
				{
					nlist = (*pLine)->findelements(plist, type, famname);
				}
				cout << "$$$Export-Elements***********************" << endl;
				//cout << "#" << sline << endl;
				if(str_icompare(type,"RING") || str_icompare(type,"LINE"))
				{
					for (unsigned int i = 0; i<nElem; i++)
						(*pLine)->getelement(i)->serialize(cout);
				}
				else{
					for (int i = 0; i<nlist; i++)
						(*pLine)->getelement(plist[i])->serialize(cout);
					cout << "******************************************" << endl;
				}
				delete []plist;
			}
			break;
			case 5: //Correct the tunes
			{
				str_remove_comments(sline);
				cout << "$$$Correct Tunes***********************" << endl;

				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0, pos + 1);
				replace(sline.begin(), sline.end(), ',', ' ');
				str_trim_both(sline);

				string type, fam1, fam2;
				pos = sline.find_first_of(("="));
				string keyw = sline.substr(0, pos);
				sline.erase(0, pos + 1);
				str_trim_left(sline);
				pos = sline.find_first_of(" \t");
				str_trim_both(keyw);
				if (str_icompare(keyw, "type"))
				{
					type = sline.substr(0, pos);
					str_trim_both(type);
					cout << type << endl;
				}
				else
				{
					cout << "error, specify type first" << endl;
					break;
				}
				sline.erase(0, pos + 1);

				pos = sline.find_first_of(("="));
				keyw = sline.substr(0, pos);
				sline.erase(0, pos + 1);
				str_trim_left(sline);
				
				str_trim_both(keyw);
				if (str_icompare(keyw, "famname"))
				{
					pos = sline.find_first_of("(");
					size_t pos2 = sline.find_first_of(")");
					string str_fam = sline.substr(pos+1, pos2 - pos - 1);
					sline.erase(0, pos2 + 1);
					
					str_trim_both(str_fam);
					std::stringstream(str_fam) >> fam1 >> fam2;
					/* //if not compile with -std=c++0x flag
					pos = str_fam.find_first_of(" \t");
					fam1 = str_fam.substr(0,pos);
					str_fam.erase(0,pos+1);
					str_trim_both(str_fam);
					fam2 = str_fam;
					*/
					cout << fam1 << " " << fam2 << endl;
				}
				else
				{
					cout << "error, specify famname secondly" << endl;
					break;
				}
				
				pos = sline.find_first_of("[");
				size_t pos2 = sline.find_first_of("]");
				string str_target = sline.substr(pos+1, pos2 - pos - 1);
				double nux, nuy, nux_target, nuy_target;
				std::stringstream(str_target) >> nux_target >> nuy_target;
				cout << "Target tunes: " << nux_target << "\t" << nuy_target << "\t" << endl;
				(*pLine)->set_rad_onoff(0);
				(*pLine)->calctune(nux, nuy);
				cout << "Initial tunes: " << nux << "\t" << nuy << "\t" << endl;

				unsigned int nfam1, nfam2;
				int nElem = (*pLine)->getnumber();

				unsigned int *plist1 = new unsigned int[nElem];
				unsigned int *plist2 = new unsigned int[nElem];

				nfam1 = (*pLine)->findelements(plist1, type, fam1);
				nfam2 = (*pLine)->findelements(plist2, type, fam2);
				if (nfam1 < 1 || nfam2 < 1)
				{
					cout << "error: no elements found " << nfam1 << "\t" << nfam2 << endl;
					delete[]plist1;
					delete[]plist2;
					break;
				}
				else
					cout << nfam1 << "\t" << nfam2 << endl;

				nux = nux_target;
				nuy = nuy_target;

				(*pLine)->correcttune(nux, nuy, nfam1, nfam2, plist1, plist2);
				cout << "corrected tunes: " << nux << "\t" << nuy << "\t" << endl;
				STRINGTOUPPER(sline);
				pos = sline.find("REPEAT");
				if (pos != string::npos)
				{
					nux = nux_target;
					nuy = nuy_target;
					(*pLine)->correcttune(nux, nuy, nfam1, nfam2, plist1, plist2);
					cout << "corrected tunes: " << nux << "\t" << nuy << "\t" << endl;
				}

				cout << "*********************************" << endl;

				delete []plist1;
				delete []plist2;
			}
			break;
			
			case 6: //Correct the chromaticities
			{
				str_remove_comments(sline);
				cout << "$$$Correct chromaticities***********************" << endl;

				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0, pos + 1);
				replace(sline.begin(), sline.end(), ',', ' ');
				str_trim_both(sline);

				string type, fam1, fam2;
				pos = sline.find_first_of(("="));
				string keyw = sline.substr(0, pos);
				sline.erase(0, pos + 1);
				str_trim_left(sline);
				pos = sline.find_first_of(" \t");
				str_trim_both(keyw);
				if (str_icompare(keyw, "type"))
				{
					type = sline.substr(0, pos);
					str_trim_both(type);
					cout << type << endl;
				}
				else
				{
					cout << "error, specify type first" << endl;
					break;
				}
				sline.erase(0, pos + 1);

				pos = sline.find_first_of(("="));
				keyw = sline.substr(0, pos);
				sline.erase(0, pos + 1);
				str_trim_left(sline);

				str_trim_both(keyw);
				if (str_icompare(keyw, "famname"))
				{
					pos = sline.find_first_of("(");
					size_t pos2 = sline.find_first_of(")");
					string str_fam = sline.substr(pos+1, pos2 - pos - 1);
					sline.erase(0, pos2 + 1);
					
					str_trim_both(str_fam);
					std::stringstream(str_fam) >> fam1 >> fam2;
					
					/* //if not compile with -std=c++0x flag
					pos = str_fam.find_first_of(" \t");
					fam1 = str_fam.substr(0,pos);
					str_fam.erase(0,pos+1);
					str_trim_both(str_fam);
					fam2 = str_fam;
					*/
					cout << fam1 << " " << fam2 << endl;
				}
				else
				{
					cout << "error, specify famname secondly" << endl;
					break;
				}

				pos = sline.find_first_of("[");
				size_t pos2 = sline.find_first_of("]");
				string str_target = sline.substr(pos + 1, pos2 - pos - 1);
				double chromx, chromy, chromx_target, chromy_target;
				std::stringstream(str_target) >> chromx_target >> chromy_target;
				cout << "Target chromaticities: " << chromx_target << "\t" << chromy_target << "\t" << endl;
				(*pLine)->set_rad_onoff(0);
				(*pLine)->set_cavity_onoff(0);

				(*pLine)->calcchrom(chromx, chromy);
				cout << "Initial chromaticities: " << chromx << "\t" << chromy << "\t" << endl;

				unsigned int nfam1, nfam2;
				int nElem = (*pLine)->getnumber();

				unsigned int *plist1 = new unsigned int[nElem];
				unsigned int *plist2 = new unsigned int[nElem];

				nfam1 = (*pLine)->findelements(plist1, type, fam1);
				nfam2 = (*pLine)->findelements(plist2, type, fam2);
				if (nfam1 < 1 || nfam2 < 1)
				{
					cout << "error: no elements found " << nfam1 << "\t" << nfam2 << endl;
					delete[]plist1;
					delete[]plist2;
					break;
				}
				else
					cout << nfam1 << "\t" << nfam2 << endl;

				chromx = chromx_target;
				chromy = chromy_target;
				(*pLine)->correctchrom(chromx, chromy, nfam1, nfam2, plist1, plist2);
				cout << "corrected chromaticities: " << chromx << "\t" << chromy << "\t" << endl;
				STRINGTOUPPER(sline);
				pos = sline.find("REPEAT");
				if (pos != string::npos)
				{
					chromx = chromx_target;
					chromy = chromy_target;
					(*pLine)->correctchrom(chromx, chromy, nfam1, nfam2, plist1, plist2);
					cout << "corrected chromaticities: " << chromx << "\t" << chromy << "\t" << endl;
				}
				cout << "*********************************" << endl;

				delete[]plist1;
				delete[]plist2;
			}
			break;
			case 7: //set maximum number of openmp threads (which overrides the environment variable)
				{
					str_remove_comments(sline);
	
					size_t pos = sline.find_first_of(("=:"));
					sline.erase(0, pos + 1);
					int num;
					std::stringstream(sline)>>num;
					if(num<1) num=1;
					omp_set_num_threads(num);
					cout<<"Maximum thread number set to "<<num<<endl;
				}
				break;
			case 8: //modify elements
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"TYPE","FAMNAME","INDEX","FIELD","VALUE","DELTA"};
				string type="", famname="", field;
				unsigned int index=0;
				double value;
				int flag_value=0;
				
				//type=Cavity, famname=RF, index=1, field=Voltage, value=3.0
				int cnt=0;
				while(sline.length()>0 && cnt++<6)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="TYPE")
					{
						type = sline.substr(0,pos);						
					}	
					else if(valtype=="FAMNAME")
					{
						famname = sline.substr(0,pos);						
					}	
					else if(valtype=="INDEX")
						std::stringstream(sline)>>index;
					else if(valtype=="VALUE")
					{
							std::stringstream(sline)>>value;
							flag_value = 1;
					}
					else if(valtype=="DELTA")
					{
							std::stringstream(sline)>>value;
							flag_value = 2;
					}
					else if (valtype == "FIELD")
					{
						field = sline.substr(0, pos);
						//cout<<str_printfinal<<endl;
					}
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				unsigned int nfam;
				int nElem = (*pLine)->getnumber();

				unsigned int *plist = new unsigned int[nElem];

				nfam = (*pLine)->findelements(plist, type, famname,index);
				if(nfam>0 && flag_value>0){
					if(field=="NSlice")
					{
						int *pns = new int[nfam];
						
						if(flag_value==2) //delta
						{
							(*pLine)->getfield((void*)pns, nfam, plist,field,index);
							for(int i=0;i<nfam;i++) pns[i] += (int) value;
						}
						else
						{ for(int i=0;i<nfam;i++) pns[i] = (int) value; }
						
						(*pLine)->setfield((void*)pns, nfam, plist,field,index);
						delete []pns;
					}
					else
					{
						double *pns = new double[nfam];
						if(flag_value==2) //delta
						{
							(*pLine)->getfield((void*)pns, nfam, plist,field,index);
							for(int i=0;i<nfam;i++) pns[i] += value;
						}
						else
						{ for(int i=0;i<nfam;i++) pns[i] = value; }
						
						(*pLine)->setfield((void*)pns, nfam, plist,field,index);
						delete []pns;
					}
				}
				else{
					if(nfam==0) cout<<"No element is found: "<<type<<"\t"<<famname<<"\t"<<index<<endl;
					if(flag_value==0) cout<<"No value is given"<<endl;
				}
			
				delete []plist;
			}	
			break;
			case 9: //set flags
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"RADIATION","CAVITY","CRABCAV","WAKE"};
				
				int cnt=0;
				while(sline.length()>0 && cnt++<4)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					
					string str_flag = sline.substr(0,pos);
					if(valtype=="RADIATION")
					{
						if(str_icompare(str_flag,"off"))
						{		
							(*pLine)->set_rad_onoff(FLAG_RAD_ONOFF_OFF);
							cout<<"radiation turned off"<<endl;
						}
						else if (str_icompare(str_flag,"rad") || str_icompare(str_flag,"on"))
						{
									(*pLine)->set_rad_onoff(FLAG_RAD_ONOFF_RAD);		
									cout<<"radiation turned on without quantum excitation"<<endl;
						}
						else if (str_icompare(str_flag,"radQE"))
						{
								(*pLine)->set_rad_onoff(FLAG_RAD_ONOFF_RADQE);			
								cout<<"radiation turned on with quantum excitation"<<endl;
						}
					}
					else if(valtype=="CAVITY")
					{
						if(str_icompare(str_flag,"off"))
						{	
							int nl = (*pLine)->set_cavity_onoff(false);
							if(nl>0)
								cout<<nl<<" cavity(ies) turned off"<<endl;
						}
						else if (str_icompare(str_flag,"on"))
						{
							int nl = (*pLine)->set_cavity_onoff(true);
							if(nl>0)
								cout<<nl<<" cavity(ies) turned on"<<endl;
						}
					}
					else if(valtype=="CRABCAV")
					{
						if(str_icompare(str_flag,"off"))
						{		
							int nl = (*pLine)->set_crab_onoff(false);
							if(nl>0)
								cout<<nl<<" crab cavity(ies) turned off"<<endl;
						}
						else if (str_icompare(str_flag,"on"))
						{
							int nl = (*pLine)->set_crab_onoff(true);
							if(nl>0)
								cout<<nl<<" crab cavity(ies) turned on"<<endl;
						}
					}
					else if(valtype=="WAKE")
					{
						if(str_icompare(str_flag,"off"))
						{		
							int nl = (*pLine)->set_wake_onoff(false);
							if(nl>0)
								cout<<nl<<" wake element(s) turned off"<<endl;
						}
						else if (str_icompare(str_flag,"on"))
						{
							int nl = (*pLine)->set_wake_onoff(true);
							if(nl>0)
								cout<<nl<<" wake element(s) turned on"<<endl;
						}
					}	
					
					sline.erase(0,pos+1);
					str_trim_both(sline);	
				} //while sline
			}
			break;
			case 10: //apply set points, same it does in DAMA block,
			{
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				string filename = sline.substr(pos + 1);
				str_trim_both(filename);

				pos = filename.find_first_of(("="));
				if (pos != string::npos)filename.erase(0, pos + 1);
				str_trim_both(filename);
				int flag = apply_line_setpt(filename, (*pLine));

			}
			break;
		}
	}
	
	return 0;
}
//select elements for print out
unsigned int sel_elements(unsigned int* plist,ALine *pLine,int nElem, string& sline)
{
	str_remove_comments(sline);
	size_t pos = sline.find_first_of(("=:"));
	string sel = sline.substr(pos+1);
	str_trim_both(sel);
	
	unsigned int nlist;
				
	string tokenlist[5];
	int  ntoken_max=5;
				
	int ntoken = str_split2(sel, string(" \t"), tokenlist, ntoken_max);
	//cout<<sel<<endl;
	//cout<<ntoken<<" tokens"<<endl;
	for(int i=0;i<ntoken;i++)
		cout<<tokenlist[i]<<" "<<endl;
								
	if(ntoken==0) //first element only
	{
		plist[0] = 0;
		nlist = 1;
	}
	else if(ntoken==1) 
	{
		string type;
		type = tokenlist[0];
		if(str_icompare(tokenlist[0],"ALL"))
		{
		//all elements
			for(int i=0;i<nElem;i++)
				plist[i] = i;
			nlist = nElem;
		}
		else 
		{
			int sw_indx = str2int(type);
			if(sw_indx>0 && sw_indx<nElem) //by index number
				{ plist[0] = sw_indx-1; nlist=1;}
			else if(sw_indx!=0)
				{	
					nlist = -1;//invalid range
					//cout<<"index="<<sw_indx<<" out of range "<<endl;
					//break;
				}
			else //by type
				nlist = pLine->findelements(plist, type);
		}
	}
	else if(ntoken==2) //type and famname
	{
		string type, famname;
		//unsigned int st, ed;
		type = tokenlist[0];
		famname = tokenlist[1];
		nlist = pLine->findelements(plist, type, famname);
	}
	return nlist;
}

unsigned int get_elem_list(unsigned int* plist, ALine *pLine, int nElem, string& type, string &famname, unsigned index, string str_indices)
{
	unsigned int nfam;
	size_t pos;

	if (type.length() == 0)
	{
		plist[0] = index;
		nfam = 1;
	}
	else if (str_icompare(type, "ALL"))
	{
		for (int i = 0; i<nElem; i++) plist[i] = i;
		nfam = nElem;
	}
	else if (famname.length() == 0)
	{
		nfam = pLine->findelements(plist, type);
		if (index>0) //one element
		{
			plist[0] = plist[index - 1];
			nfam = 1;
		}
	}
	else if (!str_indices.empty())
	{
		nfam = pLine->findelements(plist, type, famname, index);
		int i = 0;
		while (!str_indices.empty())
		{
			int tind;
			pos = str_indices.find_first_of(" \t");
			if (pos != string::npos)
			{
				std::stringstream(str_indices.substr(0, pos)) >> tind;
				if (tind > 0 && tind < nfam)
					plist[i++] = plist[tind - 1];
				str_indices.erase(0, pos + 1);
				str_trim_left(str_indices);
			}
			else
			{
				break;
			}
		}
		nfam = i;
	}
	else
		nfam = pLine->findelements(plist, type, famname, index);

	return nfam;
}
#define STR_OPTICS_PRINT_OPTICS (string("PRINT_OPTICS"))
#define STR_OPTICS_PRINT_MATRIX (string("PRINT_MATRIX"))
#define STR_OPTICS_PRINT_TUNE (string("PRINT_TUNE"))
#define STR_OPTICS_PRINT_CHROM (string("PRINT_CHROM"))
#define STR_OPTICS_PRINT_CHROM2 (string("PRINT_2ND_CHROM"))
int proc_optics(ifstream& fin, ALine *(*pLine))
{
	string key_blocks[]={STR_OPTICS_PRINT_OPTICS,STR_OPTICS_PRINT_MATRIX,STR_OPTICS_PRINT_TUNE,
		STR_OPTICS_PRINT_CHROM,STR_OPTICS_PRINT_CHROM2,BLOCK_END};
	int n_key=6;
	string sline;
	while(!fin.eof())
	{
		int sw_indx=lookforblock(fin, sline,key_blocks,n_key);
		cout<<sline<<" : "<<sw_indx<<endl;
		if(sw_indx==n_key-1) break; //end of block 
		
		switch(sw_indx)
		{
			case 0: //STR_OPTICS_PRINT
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"TYPE","FAMNAME","INDEX","DPP","DELTA"};
				string type="", famname=""; //, str_linestyle="6by6", str_option="";
				unsigned int index=0;
				string str_indices = "";
				double dpp=0, delta=DELTA_DIFF_TRACK; //1.0e-8;
				
				int cnt=0;
				while(sline.length()>0 && cnt++<5)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="TYPE")
					{
						type = sline.substr(0,pos);						
					}	
					else if(valtype=="FAMNAME")
					{
						famname = sline.substr(0,pos);						
					}
					else if (valtype == "INDEX")
					{
						if (sline.front() == '[')
						{
							size_t pos2 = sline.find_first_of(("]"));
							str_indices = sline.substr(1, pos2 - 1);
							str_trim_both(str_indices);
							str_indices.append(" ");
						}
						else
							std::stringstream(sline) >> index;
					}
					else if(valtype=="DPP")
						std::stringstream(sline)>>dpp;		
					else if(valtype=="DELTA")
						std::stringstream(sline)>>delta;		
					
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				//cout<<type<<" "<<famname<<" "<<index<<endl;
				
				unsigned int nElem = (*pLine)->getnumber();
				unsigned int *plist = new unsigned int[nElem];
				unsigned int nlist;
				
				/*if(type.length()==0)
				{
					plist[0] = index;
					nlist = 1;
				}
				else if(str_icompare(type,"ALL"))
				{
					for(int i=0;i<nElem;i++) plist[i] = i;
					nlist = nElem;
				}
				else if(famname.length()==0)
				{
					nlist = (*pLine)->findelements(plist, type);
					if(index>0) //one element
					{
						plist[0] = plist[index-1];
						nlist = 1;
					}
				}
				else
					nlist = (*pLine)->findelements(plist, type, famname,index);
				*/
				nlist = get_elem_list(plist, *pLine, nElem, type, famname, index, str_indices);

				/*unsigned int nlist = sel_elements(plist,(*pLine),nElem, sline);
				if(nlist<0)
					{
						cout<<"index="<<sw_indx<<" out of range "<<endl;
						break;
					}
					*/
			
				Twiss *twlist = new Twiss[nlist];
				cout<<"nlist="<<nlist<<endl;
				(*pLine)->calctwiss(twlist,plist, nlist,dpp,delta);
				cout<<"$$$Optics-Twiss, dp="<<dpp<<"*************************** "<<nlist<<" points"<<endl;
				cout<<"index\ts (m)\tbetx (m)\tbety (m)\talfx\talfy\tpsix\tpsiy\n";
				for(int i=0;i<nlist;i++)
				{
					cout<<plist[i]<<"\t"
					  <<twlist[i].s<<"\t"<<twlist[i].betx<<"\t"<<twlist[i].bety<<"\t"
						<<twlist[i].alfx<<"\t"<<twlist[i].alfy<<"\t"<<twlist[i].psix<<"\t"<<twlist[i].psiy<<"\t"<<endl;
					
				}
				
				delete []twlist;
				delete []plist;
				
			}
			break;
			case 1: //print matrix
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"TYPE","FAMNAME","INDEX","DPP","DELTA","LINE_STYLE","PRINT_OPTION"};
				string type="", famname="", str_linestyle="6by6", str_option="";
				unsigned int index=0;
				string str_indices = "";
				double dpp=0, delta=1.0e-8;
				
				int cnt=0;
				while(sline.length()>0 && cnt++<3)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="TYPE")
					{
						type = sline.substr(0,pos);						
					}	
					else if(valtype=="FAMNAME")
					{
						famname = sline.substr(0,pos);						
					}
					else if(valtype=="LINE_STYLE")
					{
						str_linestyle = sline.substr(0,pos);						
					}
					else if(valtype=="PRINT_OPTION")
					{
						str_option = sline.substr(0,pos);						
					}	
					else if (valtype == "INDEX")
					{
						if (sline.front() == '[')
						{
							size_t pos2 = sline.find_first_of(("]"));
							str_indices = sline.substr(1, pos2 - 1);
							str_trim_both(str_indices);
							str_indices.append(" ");
						}
						else
							std::stringstream(sline) >> index;
					}
					else if(valtype=="DPP")
						std::stringstream(sline)>>dpp;		
					else if(valtype=="DELTA")
						std::stringstream(sline)>>delta;		
					
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				//cout<<type<<" "<<famname<<" "<<index<<endl;

				if(str_icompare(type,"RING")||str_icompare(type,"ONETURN"))
				{
					Mat m66(6,6);
					m66=(*pLine)->getm66(m66,dpp,delta);
					cout<<"$$$OPTICS-MATRIX*******************ONE-TURN"<<endl;
					cout<<"type=RING at "<<(*pLine)->getelement(0)->gettype()<<", famname="<<
						(*pLine)->getelement(0)->getfamname()<<", ring index="<<0<<", spos="<<(*pLine)->getlength()<<endl;
					//cout<<m66<<endl;
					print_m66_to_line(m66, str_linestyle);
					cout<<"****************************"<<endl;
				}
				else{
					int nElem = (*pLine)->getnumber();
					unsigned int *plist = new unsigned int[nElem];
					unsigned int nlist =
						get_elem_list(plist, *pLine, nElem, type, famname, index, str_indices);

					if(nlist>0)
						cout<<"$$$OPTICS-MATRIX, dpp="<<dpp<<"******************* "<< nlist <<" points"<<endl;
					else
						cout<<"#no element found for matrix print out, type="<<type<<", famname="<<famname<<",index"<<index<<endl;
					//cout<<plist[0]<<" "<<dpp<<" "<<delta<<endl;
					for(int i=0;i<nlist;i++)
					{
						Mat m66(6,6);
						
						m66=(*pLine)->getm66(m66,plist[i], dpp,delta);
						//cout<<plist[i]<<endl;
						double spos =  (*pLine)->getlength(plist[i]);
						cout<<"type="<<type<<", famname="<<famname<<", ring index="<<plist[i]<<", spos="<<spos<<endl;
						//cout<<m66<<endl;
						print_m66_to_line(m66, str_linestyle);
						
					}
					delete []plist;
					if(nlist>0)
						cout<<"**********************"<<endl;
				}
				
			}
			break;
			case 2: //PRINT_TUNE
			case 3: //PRINT_CHROM
			case 4: //PRINT_CHROM2
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"DP","DELTA"};
				double dpp=0, delta=DELTA_DIFF_TRACK; //1.0e-8;
				
				int cnt=0;
				while(sline.length()>0 && cnt++<2)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="DPP")
						std::stringstream(sline)>>dpp;		
					else if(valtype=="DELTA")
						std::stringstream(sline)>>delta;		
					
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				if(sw_indx==2)
				{
					(*pLine)->set_rad_onoff(0);
					double tunex, tuney;
					(*pLine)->calctune(tunex, tuney,dpp,delta);
					cout << "Tunes: [" << tunex << "\t" << tuney <<"], dpp="<<dpp<< endl;
				}
				else if(sw_indx==3)
				{
					(*pLine)->set_rad_onoff(0);
					double chromx, chromy;

					(*pLine)->calcchrom(chromx, chromy,dpp,delta);
					cout << "Chromaticity: [" << chromx << "\t" << chromy << "], dpp=" <<dpp<< endl;
				}
				else //if(sw_indx==4)
				{
					(*pLine)->set_rad_onoff(0);
					double chromx, chromy,chromx2, chromy2;

					(*pLine)->calcchrom2(chromx, chromy,chromx2, chromy2);
					cout << "Chromaticity: [" << chromx << "\t" << chromy << "]" << endl;
					cout << "2nd-order Chromaticity: [" << chromx2 << "\t" << chromy2 << "], " << endl;
				}
				
				
			}
			break;

		} //SWITCH
	} //while feof
	return 0;
}

void print_m66_to_line(Mat& m66, string &str_linestyle)
{
	if(str_icompare(str_linestyle,"oneline"))
	{
		for(int i=1;i<=6;i++)
			for(int j=1;j<=6;j++)
				cout<<m66[i][j]<<"\t";
		cout<<endl;
	}
	else if (str_icompare(str_linestyle,"4by4"))
	{
		for(int i=1;i<=4;i++)
		{	
			for(int j=1;j<=4;j++)
				cout<<m66[i][j]<<"\t";
			cout<<endl;
		}
	}
	else if (str_icompare(str_linestyle,"6by6"))
	{
		for(int i=1;i<=6;i++)
		{	
			for(int j=1;j<=6;j++)
				cout<<m66[i][j]<<"\t";
			cout<<endl;
		}
	}

			
}
//generate distribution
int gen_distri(double **r, string& sline, ALine **pLine=NULL)
{
	cout<<sline<<endl;
	str_remove_comments(sline);
	//cout<<sline<<endl;
	size_t pos = sline.find_first_of(("=:"));
	string sel = sline.substr(pos+1);
	str_trim_both(sel);
	
	pos = sel.find_first_of((" \t"));
	
	if(pos==string::npos) //one particle with all zeros
	{
		double *r6 = new double[6];
		for(int i = 0;i<6;i++) r6[i]=0;
		(*r)= r6;
		//cout<<"one particle: "<<r6[0]<<endl;
		
		int np = 1;
		return np;
	}
	else 
	{
		string option = sel.substr(0,pos);
		//cout<<"Option="<<option<<endl;
		if(str_icompare(option,"ZERO")) //N particle with all zeros
		{
			size_t pos2 = sel.find_first_of(("="));
			string str_num = sel.substr(pos2+1);
			
			//cout<<"number"<<str_num<<endl;
			
			int np; // = str2int(str_num);
			std::stringstream(str_num)>>np; 
			if(np>0) 
			{
				(*r)=new double[6*np];
				for(int i = 0;i<6*np;i++) (*r)[i]=0;
				return np;
			}
			else return 0;
		
		}
		else if(str_icompare(option,"LOAD")) //load file
		{
			string filename;
			sel.erase(0, pos + 1);
			size_t pos = sel.find_first_of(", \t\n");
			filename = sel.substr(0,pos);
			str_trim_both(filename);
			cout<<filename<<endl;
			
			std::ifstream fin;
			fin.open(filename.c_str());

			string str_format = "PLAIN";
			sel.erase(0, pos + 1);
			pos = sel.find("=");
			if (pos != string::npos)
			{
				string valtype = sel.substr(0, pos);
				str_trim_both(valtype);
				if (str_icompare(valtype, "FORMAT"))
				{
					//size_t pos = sel.find_first_of(", \t\n");
					str_format = sel.substr(pos + 1);
					str_trim_both(str_format);
					STRINGTOUPPER(str_format);
				}

			}
			if (str_icompare(str_format, "ACCEL-EXPORT"))
			{
				while (fin.good())
				{
					string sl;
					mgetline(fin, sl);
					if (sl.compare(0, 3, "$$$")) continue;
					else
					{
						if (sl.find("Full export")!=string::npos)
						{
							mgetline(fin, sl); //read the Np=??, NPrint=?? line
							cout << sl << endl;
							break;
						}
					}
				}
			}


			int np = 0;
			std::vector<double*> pr;
			while(fin.good())
			{
				string sl;
				mgetline(fin, sl);
				if (!sl.compare(0, 5, "*****"))break;
				//cout<<sl<<endl;
				str_trim_both(sl);
				if(sl.length()>11 && sl.compare(0,1,"%"))
				{
					np++;
					double *r6 = new double[6];
					if (str_icompare(str_format, "ACCEL-EXPORT"))
					{
						int tmp;
						cout << sl << endl;
						std::stringstream(sl) >>tmp>> r6[0] >> r6[1] >> r6[2] >> r6[3] >> r6[4] >> r6[5];
					}
					else
					{
						std::stringstream(sl) >> r6[0] >> r6[1] >> r6[2] >> r6[3] >> r6[4] >> r6[5];
					}
					pr.push_back(r6);
					//cout<<pr[np-1][0]<<endl;
					//cout<<pr.size()<<endl;
				}
				
			}
			fin.close();
			//cout<<"# of particles :"<<np<<endl;
			if(np>0)
			{
				(*r) = new double[6*np];
				for(int i=0;i<np;i++)
				{
					double *r6 = (*r)+6*i;
					for(int j=0;j<6;j++) r6[j] = pr[i][j];
					delete []pr[i];
				}
			}
			
			return np;
			
		} //if(str_icompare
		else if(str_icompare(option,"Coordinate"))
		{
			size_t bl = sel.find('(');
			size_t br = sel.find(')');
			string str_r = sel.substr(bl+1,br-bl-1);
			std::replace(str_r.begin(),str_r.end(),',',' ');
			//cout<<str_r<<endl;
			
			std::vector<double*> pr;
			int np = 1;
			
			double *r6 = new double[6];
			std::stringstream(str_r)>>r6[0]>>r6[1]>>r6[2]>>r6[3]>>r6[4]>>r6[5];
			pr.push_back(r6);
					
			pos=str_r.find_first_of(";");
			while(pos!=string::npos)
			{
				//cout<<np<<" "<<str_r<<endl;
				np++;
				str_r.erase(0,pos+1);
				double *r6 = new double[6];
				std::stringstream(str_r)>>r6[0]>>r6[1]>>r6[2]>>r6[3]>>r6[4]>>r6[5];
				pr.push_back(r6);
				pos=str_r.find_first_of(";");
			}
			//cout<<"# of particles :"<<np<<endl;
			if(np>0)
			{
				(*r) = new double[6*np];
				for(int i=0;i<np;i++)
				{
					double *r6 = (*r)+6*i;
					for(int j=0;j<6;j++) r6[j] = pr[i][j];
					delete []pr[i];
				}
			}
			
			return np;
			
		}
		else if(str_icompare(option,"Gaussian"))
		{
			sel.erase(0,pos+1);
			
			int np=1;
			pos = sel.find("NPar");
			if(pos!=string::npos)
			{
				string str_num = sel.substr(pos+4);
				size_t pos2=str_num.find_first_of("=");
				str_num.erase(0,pos2+1);
				std::stringstream(str_num)>>np;
			}
			
			double sigx=0, sigxp=0;
			pos = sel.find("emit_x");
			if(pos!=string::npos)
			{
				string str_num = sel.substr(pos+6);
				size_t pos2=str_num.find_first_of("=");
				str_num.erase(0,pos2+1);
				double emit_x, betx;
				std::stringstream(str_num)>>emit_x;

				pos = sel.find("betx");
				if (pos != string::npos)
				{
					pos2 = sel.find_first_of("=",pos+4);
					str_num.erase(0, pos2 + 1);
					std::stringstream(str_num) >> betx;
					if(betx<=0)
					{
						betx = 1.0;
						cout << "#betx=1.0 m is assumed" << endl;
					}
				}
				else
				{
					Twiss tw;
					(*pLine)->calctwiss(tw);
					if (std::isnan(tw.betx))
					{
						betx = 1.0;
						cout << "#betx=1.0 m is assumed" << endl;
					}
					else
						betx = tw.betx;
				}

				sigx = sqrt(betx*emit_x);
				sigxp = sqrt(emit_x/betx);
				
			}
			else
			{ //look for sigx and sigxp
				pos = sel.find("sigx");
				if(pos!=string::npos)
				{
					string str_num = sel.substr(pos+4);
					size_t pos2=str_num.find_first_of("=");
					str_num.erase(0,pos2+1);
					std::stringstream(str_num)>>sigx;
					
				}
				pos = sel.find("sigxp");
				if(pos!=string::npos)
				{
					string str_num = sel.substr(pos+5);
					size_t pos2=str_num.find_first_of("=");
					str_num.erase(0,pos2+1);
					std::stringstream(str_num)>>sigxp;					
				}
				
			}
			
			double sigy=0, sigyp=0;
			pos = sel.find("emit_y");
			if(pos!=string::npos)
			{
				string str_num = sel.substr(pos+6);
				size_t pos2=str_num.find_first_of("=");
				str_num.erase(0,pos2+1);
				double emit_y, bety;
				std::stringstream(str_num) >> emit_y;

				pos = sel.find("bety");
				if (pos != string::npos)
				{
					pos2 = sel.find_first_of("=", pos + 4);
					str_num.erase(0, pos2 + 1);
					std::stringstream(str_num) >> bety;
					if (bety <= 0)
					{
						bety = 1.0;
						cout << "#bety=1.0 m is assumed" << endl;
					}
				}
				else
				{
					Twiss tw;
					(*pLine)->calctwiss(tw);
					if (std::isnan(tw.bety))
					{
						bety = 1.0;
						cout << "#betx=1.0 m is assumed" << endl;
					}
					else
						bety = tw.bety;
				}
				sigy = sqrt(bety*emit_y);
				sigyp = sqrt(emit_y/bety);
			}
			else
			{ //look for sigy and sigyp
				pos = sel.find("sigy");
				if(pos!=string::npos)
				{
					string str_num = sel.substr(pos+4);
					size_t pos2=str_num.find_first_of("=");
					str_num.erase(0,pos2+1);
					std::stringstream(str_num)>>sigy;
					
				}
				pos = sel.find("sigyp");
				if(pos!=string::npos)
				{
					string str_num = sel.substr(pos+5);
					size_t pos2=str_num.find_first_of("=");
					str_num.erase(0,pos2+1);
					std::stringstream(str_num)>>sigyp;	
				}
				
			}
			
			double sigz=0, sigdpp=0;
			pos = sel.find("sigz");
			if(pos!=string::npos)
			{
				string str_num = sel.substr(pos+4);
				size_t pos2=str_num.find_first_of("=");
				str_num.erase(0,pos2+1);
				std::stringstream(str_num)>>sigz;
			
			}
			pos = sel.find("sigdpp");
			if(pos!=string::npos)
			{
				string str_num = sel.substr(pos+6);
				size_t pos2=str_num.find_first_of("=");
				str_num.erase(0,pos2+1);
				std::stringstream(str_num)>>sigdpp;	
			}
			
			(*r) = gen_distr_6D(np, sigx,sigxp,sigy,sigyp,sigz,sigdpp);
			return np;
		}
		else if(str_icompare(option,"Linear"))	
		{
			sel.erase(0,pos+1);
			str_trim_left(sel);
			pos = sel.find_first_of("=");
			string str_coord = sel.substr(0,pos);
			str_trim_both(str_coord);
			
			sel.erase(0,pos+1);
			str_trim_left(sel);
			pos = sel.find_first_of(" \t");
			replace(sel.begin(),sel.end(),':',' ');
			double st, ed, intval;
			std::stringstream(sel)>>st>>intval>>ed;	
			//cout<<st<<" "<<intval<<" "<<ed<<endl;
			if(intval==0)
			{
				cout<<"interval cannot be zero"<<endl;
				return 0;
			}
			int np=0;
			double val = st;
			while((val<=ed && intval>0) || (val>=ed && intval<0))
			{
				val += intval;
				np++;
			}		
			(*r) = new double[np*6];
			for(int i=0;i<np*6;i++)(*r)[i] = 0.0;
			
			int j=0;
			if(str_icompare(str_coord,"X")) j=0;
			else if(str_icompare(str_coord,"XP")) j=1;
			else if(str_icompare(str_coord,"Y")) j=2;
			else if(str_icompare(str_coord,"YP")) j=3;
			else if(str_icompare(str_coord,"Z")) j=4;
			else if(str_icompare(str_coord,"DELTA")) j=5;
			
			for(int i=0;i<np;i++)(*r)[i*6+j] = st+intval*i;
			
			return np;
		}
	}
	
	
	return 0;
	
}

#define STR_TRACK_DISTRIBUTION (string("DISTRIBUTION"))
#define STR_TRACK_MONITOR (string("MONITOR"))
#define STR_TRACK_DO_TRACK (string("DO_TRACK"))
#define STR_TRACK_TRANSPORT (string("TRANSPORT"))
int proc_track(ifstream& fin, ALine **pLine)
{
	string key_blocks[]={STR_TRACK_DISTRIBUTION,STR_TRACK_MONITOR, STR_TRACK_DO_TRACK,STR_TRACK_TRANSPORT,BLOCK_END};
	int n_key=5;
	string sline;
	double *r=NULL;
	int np = 0;
	int Nturn = 0;
	while(!fin.eof())
	{
		int sw_indx=lookforblock(fin, sline,key_blocks,n_key);
		//cout<<sline<<" : "<<sw_indx<<endl;
		if(sw_indx==n_key-1) break; //end of block 
		str_remove_comments(sline);
		switch(sw_indx)
		{
			case 0: //set up distribution
			{
				delete []r;
				np = gen_distri(&r, sline, pLine);
				STRINGTOUPPER(sline);
				
				size_t pos;
				if(np>0)
				{
					pos = sline.find("SHIFT");
					if(pos!=string::npos){
						pos = sline.find_first_of("(",pos+5);
						size_t pos2 = sline.find_first_of(")",pos+1);
						string str_sh = sline.substr(pos+1, pos2-pos-1);
						replace(str_sh.begin(),str_sh.end(),',',' ');
						//cout<<str_sh<<endl;
						double dlt[6];
						std::stringstream(str_sh)>>dlt[0]>>dlt[1]>>dlt[2]>>dlt[3]>>dlt[4]>>dlt[5];
						for(int i=0;i<np;i++)
						{	for(int j=0;j<6;j++)
								r[j+6*i] += dlt[j];
						}
					}
				}
				pos = sline.find("PRINT");
				if(pos!=string::npos){
					cout << "$$$Init-Distribution*************" << endl;
					for(int i=0;i<np;i++)
					{	for(int j=0;j<6;j++)
							cout<<r[j+6*i]<<" ";
						cout<<endl;
					}
					cout << "***********************" << endl;
				}
			}
			break;
			case 1: //set up monitor
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0, pos + 1);
				str_trim_both(sline);
				replace(sline.begin(), sline.end(), ',', ' ');
				
				AMonitor *am=new AMonitor;
				string str_amfam = "FirstMonitor";
				am->setfield((void*)&str_amfam, "FamName");

				//88888888888888
				string keyws[] = { "INTERVAL","FRACTION","PRINTFLAG" };
				int interval=1;
				double fraction = 1.0;
				string str_flag="COORDINATE";

				int cnt = 0;
				while (sline.length()>0 && cnt++<3)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0, pos);
					sline.erase(0, pos + 1);
					str_trim_both(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if (valtype == "INTERVAL")
						std::stringstream(sline) >> interval;
					else if (valtype == "FRACTION")
						std::stringstream(sline) >> fraction;
					else if (valtype == "PRINTFLAG") 
						str_flag = sline.substr(0, pos);
					sline.erase(0, pos + 1);
					str_trim_both(sline);
				}

				if (interval<1) interval = 1;
				am->setfield((void*)&interval, string("Interval"));

				if (fraction<0.01) fraction = 0.01;
				if (fraction>1.0)  fraction = 1.0;
				//cout<<fraction<<endl;
				am->setfield((void*)&fraction, string("Ratio"));

				//am->serialize(cout);				
				(*pLine)->insertelement(0,am);
				if(str_icompare(str_flag,"COORDINATE"))
					(*pLine)->set_monitor(1,1);//flag, index
				else if (str_icompare(str_flag, "STAT1"))
					(*pLine)->set_monitor(2, 1);//flag, index
				else if (str_icompare(str_flag, "STAT2"))
					(*pLine)->set_monitor(3, 1);//flag, index
			}
			break;
			case 2: //DO tracking
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"RADIATION","RFCAVITY","NTURN","PRINTFINAL"};
				
				string str_rad="", str_rf="",str_printfinal="NO";
				
				int cnt=0;
				while(sline.length()>0 && cnt++<4)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="RADIATION")
					{
						str_rad = sline.substr(0,pos);						
					}	
					else if(valtype=="RFCAVITY")
						str_rf = sline.substr(0,pos);
					else if(valtype=="NTURN")
						std::stringstream(sline)>>Nturn;
					else if (valtype == "PRINTFINAL")
					{
						str_printfinal = sline.substr(0, pos);
						//cout<<str_printfinal<<endl;
					}
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				
				unsigned flag_rad;
				if(str_icompare(str_rad,"OFF"))  {
					flag_rad = FLAG_RAD_ONOFF_OFF;
					(*pLine)->set_rad_onoff(flag_rad);
				}
				else if(str_icompare(str_rad,"RAD")) {
					flag_rad = FLAG_RAD_ONOFF_RAD;
					(*pLine)->set_rad_onoff(flag_rad);
				}
				else if(str_icompare(str_rad,"RADQE")) {
					flag_rad = FLAG_RAD_ONOFF_RADQE;
					(*pLine)->set_rad_onoff(flag_rad);
				}
				else if(str_icompare(str_rad,"LUMP")) {
					flag_rad = FLAG_RAD_ONOFF_LUMP;
					(*pLine)->set_rad_onoff(flag_rad);
				}
				else if(str_icompare(str_rad,"LUMPDAMP")) {
					flag_rad = FLAG_RAD_ONOFF_LUMP_DAMP;
					(*pLine)->set_rad_onoff(flag_rad);
				}
				else if(str_icompare(str_rad,"LUMPQE")) {
					flag_rad = FLAG_RAD_ONOFF_LUMP_QE;
					(*pLine)->set_rad_onoff(flag_rad);
				}
								
				if(str_icompare(str_rf,"OFF"))
				{
					(*pLine)->set_cavity_onoff(FLAG_CAVITY_OFF);
				}
				else if(str_icompare(str_rf,"ON"))
				{
					(*pLine)->set_cavity_onoff(FLAG_CAVITY_ON);
				}
				
				
				//DO TRACKING
				if(np>0 && Nturn>0)
				{
					AMonitor* pm = (AMonitor*)(*pLine)->getelement(0);
					unsigned int flag_print;
					pm->getfield((void*)&flag_print,"PrintFlag");

					if (flag_print == 1)
					{
						cout << "$$$Track-Data: Coordinate**************************" << endl;
						double fraction;
						pm->setfield((void*)&fraction, string("Ratio"));
						int iskip = int(1. / fraction);
						if (iskip<1)iskip = 1;
						if (iskip>np)iskip = np;

						cout << "#" << str_rad << " " << str_rf << " " << Nturn << endl;
						cout << "Np=" << np <<", Nprint="<<np/iskip<< endl;
					}
					else if (flag_print == 2)
					{
						cout << "$$$Track-Data: Stat1******************************" << endl;
						cout << "#" << str_rad << " " << str_rf << " " << Nturn << endl;
						cout << "Np=" << np << endl;
					}

					else if (flag_print == 3)
					{
						cout << "$$$Track-Data: Stat2******************************" << endl;
						cout << "#" << str_rad << " " << str_rf << " " << Nturn << endl;
						cout << "Np=" << np << endl;
					}

					
					
					//clock_t start=clock();
					//double duration;
					
					//delete []r;
					//r = gen_distr_6D(np, 9.5, 10.0E-9, 5, 10.0E-12, 0.006, 0.001);
					
					/*for(int i=0;i<np;i++)
					{	for(int j=0;j<6;j++)
							cout<<r[j+6*i]<<" ";
						cout<<endl;
					}
					 */	
					//cout<<"Nturn="<<Nturn<<endl;
					//cout<<"#";		
					(*pLine)->set_monitor_pass(0);
					/*int nwake = (*pLine)->set_wake_onoff(true);
					if (nwake > 0)
						cout << "#turn on " << nwake << " wake elements" << endl;
						*/

					(*pLine)->pass(r,np,Nturn);
					//(*pLine)->set_wake_onoff(false); //turn off wake after tracking
					if(str_icompare(str_printfinal, "Full"))
					{
						cout << "*************************************" << endl;
						cout << "$$$Track-Data: Full export ************************" << endl;
						cout << "Np=" << np << ", Nprint=" << np << endl;
						for (int i = 0; i<np; i++)
						{
							cout << Nturn << "\t";
							for (int j = 0; j<6; j++)
								cout << r[j + 6 * i] << "\t";
							cout << endl;
						}
					}
					else if (str_icompare(str_printfinal, "Yes") || str_icompare(str_printfinal, "same"))
					{
						//pm->pass(r, np);
						if (flag_print == 1)
						{ //print coordinate
							double ratio_print;
							pm->getfield((void*)&ratio_print, "Ratio");
							int iskip = int(1. / ratio_print);
							if (iskip<1)iskip = 1;
							if (iskip>np)iskip = np;
							for (int i = 0; i<np; i += iskip)
							{
								double *rp = r + 6 * i;
								cout << Nturn << "\t" << rp[0] << "\t" << rp[1] << "\t" << rp[2] << "\t" << rp[3] << "\t" << rp[4] << "\t" << rp[5] << endl;
							}
						}
						else if (flag_print == 2 && np >= 2)
						{ //print statistics, no off-diagonal elements
							double ar[6], sigr[6];
							int nlive = calc_beam_stat(r, np, ar, sigr, NULL);
							cout << std::scientific;
							cout << Nturn << "\t" << nlive << "\t" << ar[0] << "\t" << ar[1] << "\t" << ar[2] << "\t" << ar[3] << "\t" << ar[4] << "\t" << ar[5] ;
							cout << "\t" << sigr[0] << "\t" << sigr[1] << "\t" << sigr[2] << "\t" << sigr[3] << "\t" << sigr[4] << "\t" << sigr[5] << endl;
						}
						else if (flag_print == 3 && np >= 2)
						{ //print statistics, with off-diagonal elements
							double ar[6], sigr[6], sigrr[6];
							int nlive = calc_beam_stat(r, np, ar, sigr, sigrr);
							cout << std::scientific;
							cout << Nturn << "\t" << nlive << "\t" << ar[0] << "\t" << ar[1] << "\t" << ar[2] << "\t" << ar[3] << "\t" << ar[4] << "\t" << ar[5] ;
							cout << "\t" << sigr[0] << "\t" << sigr[1] << "\t" << sigr[2] << "\t" << sigr[3] << "\t" << sigr[4] << "\t" << sigr[5] ;
							cout << "\t" << sigrr[0] << "\t" << sigrr[1] << "\t" << sigrr[2] << "\t" << sigrr[3] << "\t" << sigrr[4] << "\t" << sigrr[5] << endl;
						}
					}
					
					//duration = ( clock() - start )/(double) CLOCKS_PER_SEC;
					cout<<"*************************************"<<endl;
				
				}
			}
			break;
			case 3: //TRANSPORT to locations in the ring
			{	
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"TYPE","FAMNAME","INDEX"};
				string type="", famname=""; //, str_linestyle="6by6", str_option="";
				unsigned int index=0;
				string str_indices = ""; //in case a list of indices is given in square bracket, such as [1,4,5]
				
				int cnt=0;
				while(sline.length()>0 && cnt++<3)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_left(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="TYPE")
					{
						type = sline.substr(0,pos);						
					}	
					else if(valtype=="FAMNAME")
					{
						famname = sline.substr(0,pos);						
					}
					else if(valtype=="INDEX")
					{
						if (sline.front() == '[')
						{
							size_t pos2= sline.find_first_of(("]"));
							str_indices = sline.substr(1, pos2-1);
							str_trim_both(str_indices);
							str_indices.append(" "); 
						}
						else
							std::stringstream(sline) >> index;
					}
						
										
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				if(str_indices.empty())
					cout<<type<<" "<<famname<<" "<<index<<endl;
				else
					cout << type << " " << famname << " " << str_indices << endl;
				
				unsigned int nElem = (*pLine)->getnumber();
				unsigned int *plist = new unsigned int[nElem];
				unsigned int nlist;
				
				
				nlist = get_elem_list(plist, *pLine, nElem, type, famname, index, str_indices);
				//
				double *nr = new double[np*6];
				for(int i=0;i<np*6;i++) nr[i] = r[i];
				
			  //cout<<nlist<<endl;
				if(nlist>0)
				{
					cout << "$$$Track-Transport: ************************" << endl;
					(*pLine)->pass(nr,np,0,plist[0]);
					double spos =  (*pLine)->getlength(plist[0]);
					cout<<"type="<<type<<", famname="<<famname<<", ring index="<<plist[0]<<", spos="<<spos<<endl;
					for (int i = 0; i<np; i++)
						{
							cout << plist[0] << "\t";
							for (int j = 0; j<6; j++)
								cout << r[j + 6 * i] << "\t";
							cout << endl;
						}
						cout << "#------end of ring index=" <<plist[0]<< endl;
						
					for(unsigned int k=1;k<nlist; k++) //loop over watch points
					{
						spos =  (*pLine)->getlength(plist[k]);
						cout<<"type="<<type<<", famname="<<famname<<", ring index="<<plist[k]<<", spos="<<spos<<endl;
						(*pLine)->pass(nr,np,plist[k-1],plist[k]);
						
						for (int i=0; i<np; i++)
						{
							cout << plist[k] << "\t";
							for (int j = 0; j<6; j++)
								cout << r[j + 6 * i] << "\t";
							cout << endl;
						}
						cout << "#------end of ring index=" <<plist[k]<< endl;
					}
					cout << "******************************" << endl;
				}
				else
					cout<<"No watch point found"<<endl;
				
				delete []nr;
				}
				break;
							
		} //end_switch
		
	} //end_while
	
	
	
	delete []r;
	return 0;
}

int get_type_famname_str(string &sline, string &type, string &famname)
{
			string keyw;		
			int flag= -1;
			
			size_t pos = sline.find_first_of(("="));
			keyw=sline.substr(0,pos);
			sline.erase(0,pos+1);
			str_trim_left(sline);
			pos = sline.find_first_of(" \t");
			str_trim_both(keyw);
			if(str_icompare(keyw,"type"))
			{
				type = sline.substr(0,pos);
				str_trim_both(type);
				flag=1;
			}
			else if(str_icompare(keyw,"famname"))
			{
				famname = sline.substr(0,pos);
				str_trim_both(famname);
				flag=10;
			}
			sline.erase(0,pos+1);
			
			pos = sline.find_first_of(("="));
			keyw=sline.substr(0,pos);
			sline.erase(0,pos+1);
			str_trim_left(sline);
			pos = sline.find_first_of(" \t");
			str_trim_both(keyw);
			if(str_icompare(keyw,"type"))
			{
				type = sline.substr(0,pos);
				str_trim_both(type);
				flag+=1;
			}
			else if(str_icompare(keyw,"famname"))
			{
				famname = sline.substr(0,pos);
				str_trim_both(famname);
				flag+=10;
			}
			sline.erase(0,pos+1);
			//cout<<type<<" "<<famname<<endl;
			//cout<<sline<<endl;
		return flag; //if -1, not found type or famname, if 1 or 2, only found type, if 11, found both
}
#define STR_DAMA_APPLY_SETPT (string("APPLY_SETPT"))
#define STR_DAMA_TRACK_DA (string("TRACK_DA"))
#define STR_DAMA_TRACK_MA (string("TRACK_MA"))
#define STR_DAMA_TRACK_DNUDJ (string("TRACK_DNUDJ"))
bool apply_line_setpt(const string &filename, ALine *pLine)
{
	std::ifstream dfin;
	dfin.open(filename.c_str());
	cout<<filename<<endl;
	string sline;
	bool flag_check_first=false;
	while(dfin.good())
	{
		mgetline(dfin,sline);
		str_remove_comments(sline);
		const string str_first = "$FirstElement";
		
		if(str_icompare(sline.substr(0,str_first.length()),str_first))
		{
			size_t pos=sline.find_first_of(":");
			sline.erase(0,pos+1);
			replace(sline.begin(),sline.end(),',',' ');
			str_trim_left(sline);
				
			string type, famname;
			get_type_famname_str(sline, type, famname);
			
			//cout<<type<<" "<<famname<<endl;
			int nElem = pLine->getnumber();
			for(int i=0;i<nElem;i++)
			{
				AccelElement* pElem = pLine->getelement(i);
				//cout<<i<<": "<<pElem->gettype()<<" "<<pElem->getfamname()<<endl;
				if(pElem->gettype()==type && pElem->getfamname()==famname)
				{
					flag_check_first = true; //correct
					break;
				}
				else if(pElem->getlength()>0)
					break;
				else if(pElem->gettype()!="Monitor")
					break;
					
				//if(i>3)break;
			}
			//cout<<"first element check: "<<flag_check_first<<endl;
			break;
		}	
		
	} //end while, checking first element
	if(!flag_check_first)
	{
		cout<<"#first element check: "<<flag_check_first<<endl;
		return false;
	}
	else cout<<"#first element matches"<<endl;
	
	while(dfin.good())
	{
		mgetline(dfin,sline);
		str_remove_comments(sline);
		const string str_first = "$SetValue";
		
		if(!str_icompare(sline.substr(0,str_first.length()),str_first))
			continue;
		size_t pos=sline.find_first_of(":");
			sline.erase(0,pos+1);
			replace(sline.begin(),sline.end(),',',' ');
			str_trim_left(sline);
		
		string type, famname;
		int flag = get_type_famname_str(sline, type, famname);
		
		int nlist;
		unsigned int nElem = (pLine)->getnumber();
						//cout<<nElem<<endl;
		unsigned int *plist = new unsigned int[nElem];
						//cout<<sline<<endl;
		//cout<<type<<" "<<famname<<" "<<flag<<endl;
			
		if(flag<1)
			return false;
		else if(flag==1 || flag==2)
		{
			nlist = (pLine)->findelements(plist, type);
		}
		else if(flag==11)
		{
			nlist = (pLine)->findelements(plist, type, famname);
		}
		double *pval = new double[nlist];
		for(int i=0;i<nlist;i++)
		{
			dfin>>pval[i];
		}
		if(type=="StrMPole")
		{
			size_t pos2=sline.find("Poly");
			string valtype = sline.substr(pos2,5);
			sline.erase(0,pos2+5);
			pos2 = sline.find_first_of('[');
			size_t pos3=sline.find(']');
			string str_num = sline.substr(pos2+1,pos3-pos2);
			int valindex = str2int(str_num);
			if(valtype=="PolyB") valtype = "PolynomB";
			if(valtype=="PolyA") valtype = "PolynomA";
			
			//cout<<valtype<<" "<<valindex<<endl;
			for(int i=0;i<nlist;i++)
			{
				AStrMPole* pE= (AStrMPole*)pLine->getelement(plist[i]);
				double tmpval = pval[i];
				pE->setfield((void*) &tmpval,valtype,valindex);
			}
		}
		
		
		delete []plist;	
		delete []pval;
	}
	
	//dfin.close();
	return false;
}
int proc_DAMA(ifstream& fin, ALine *pLine)
{
	string key_blocks[]={STR_DAMA_APPLY_SETPT,STR_DAMA_TRACK_DA, STR_DAMA_TRACK_MA,STR_DAMA_TRACK_DNUDJ, BLOCK_END};
	int n_key=5;
	string sline;
	
	while(!fin.eof())
	{
		
		int sw_indx=lookforblock(fin, sline,key_blocks,n_key);
		str_remove_comments(sline);
		//cout<<sline<<" : "<<sw_indx<<endl;
		if(sw_indx==n_key-1) break; //end of block 
		
		switch(sw_indx)
		{
			case 0: //apply set points,
			{
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				string filename = sline.substr(pos+1);
				str_trim_both(filename);
				
				pos = filename.find_first_of(("="));
				if(pos!=string::npos)filename.erase(0,pos+1);
				str_trim_both(filename);
				int flag = apply_line_setpt(filename, pLine);	
								
			}
			break;
			case 1: //track DA
			{
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"XMIN","XMAX","YMAX","NRAY","NSTEP","NTURN","PRINT","WEIGHTMINUSSIDE"};
				int nray, nstep, nturn;
				double xmin, xmax, ymax;
				bool flag_weightminusside = false;
				
				int cnt=0;
				while(sline.length()>0 && cnt++<8)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_both(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="XMIN")
						std::stringstream(sline)>>xmin;
					else if(valtype=="XMAX")
						std::stringstream(sline)>>xmax;
					else if(valtype=="YMAX")
						std::stringstream(sline)>>ymax;
					else if(valtype=="NRAY")
						std::stringstream(sline)>>nray;
					else if(valtype=="NSTEP")
						std::stringstream(sline)>>nstep;
					else if(valtype=="NTURN")
						std::stringstream(sline)>>nturn;
					else if(valtype=="PRINT")
					{
						string str_print;
						if(pos==string::npos)str_print = sline;
						else str_print = sline.substr(0,pos);
						if(str_icompare(str_print,"loss_turn"))
							pLine->set_flag_msg(FLAG_MSG_DA);
					}
					else if(valtype=="WEIGHTMINUSSIDE")
					{
						string str_weight;
						if(pos==string::npos)str_weight = sline;
						else str_weight = sline.substr(0,pos);
						cout<<str_weight<<endl;
						if(str_icompare(str_weight,"true"))
							flag_weightminusside=true;
					}
					
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				double rx[51], ry[51];
				if(nray<2 || nray>51) cout<<"Error: number of rays: "<<nray<<endl;
				
				pLine->set_monitor(0,1);
				cout<<"#"<<xmin<<" "<<xmax<<" "<<ymax<<" "<<nray<<" "<<nstep<<" "<<nturn<<" "<<flag_weightminusside<<endl;
				
				pLine->set_rad_onoff(FLAG_RAD_ONOFF_RAD);
				pLine->set_cavity_onoff(FLAG_CAVITY_ON);
				
				cout<<"#"<<currentDateTime()<<endl;
				double Sxy = pLine->trackDA(xmax, xmin, ymax,	nturn, nray, nstep, rx, ry,flag_weightminusside);
				if(flag_weightminusside)
					cout<<"$$$DA-profile "<<nturn<<" turns, 3/4 DA area="<<Sxy*1.0e6<<" mm^2"<<endl;
				else
					cout<<"$$$DA-profile "<<nturn<<" turns, DA area="<<Sxy*1.0e6<<" mm^2"<<endl;
				for(int i=0;i<nray;i++)
					cout<<rx[i]<<" \t"<<ry[i]<<" \t"<<endl;
				cout<<"#"<<currentDateTime()<<endl;
				cout<<"******************************"<<endl;
				
				
			}
			break;
			
			case 2: //track MA
			{
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				//type=Quad,famname=QFC,NMONI=7,dppmin=0.01, dppmax=0.035,NSTEP=25,NTURN=2000
				string type, famname;
				int nlist;
				int nElem = pLine->getnumber();
					//cout<<nElem<<endl;
				unsigned int *plist = new unsigned int[nElem];
					//cout<<sline<<endl;
					
				int flag = get_type_famname_str(sline, type,famname);
					//cout<<sline<<endl;
					//cout<<type<<" "<<famname<<" "<<flag<<endl;
				if(flag<1)
					break;
				else if(flag==1 || flag==2)
				{
					nlist = pLine->findelements(plist, type);
				}
				else if(flag==11)
				{
					nlist = pLine->findelements(plist, type, famname);
				}
				
				string keyws[]={"DPPMIN","DPPMAX","NSTEP","NTURN","NMONI"};
				unsigned int nmoni;
				int nstep, nturn;
				double dppmin, dppmax;
				
				int cnt=0;
				while(sline.length()>0 && cnt++<5)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_both(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="DPPMIN")
						std::stringstream(sline)>>dppmin;
					else if(valtype=="DPPMAX")
						std::stringstream(sline)>>dppmax;
					else if(valtype=="NMONI") //number of locations
						std::stringstream(sline)>>nmoni;
					else if(valtype=="NSTEP")
						std::stringstream(sline)>>nstep;
					else if(valtype=="NTURN")
						std::stringstream(sline)>>nturn;
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				cout<<"#"<<nmoni<<" "<<nstep<<" "<<nturn<<" "<<dppmin<<" "<<dppmax<<endl;
				if(nmoni>nlist) nmoni = nlist;
				
				double dppavg;
				double *dpp_plus = new double[nmoni];
				double *dpp_minus = new double[nmoni];
				
				pLine->set_rad_onoff(FLAG_RAD_ONOFF_RAD);
				pLine->set_cavity_onoff(FLAG_CAVITY_ON);
				cout<<"#"<<currentDateTime()<<endl;
				//cout<<nlist<<endl;
				dppavg = pLine->trackMA(dppmax, dppmin, nturn, nstep, dpp_plus, dpp_minus,nmoni, plist);
				
				if(1)
				{
					cout<<"$$$MA "<<nturn<<" turns, average MA="<<dppavg<<endl;
					for(int i=0;i<nmoni;i++)
						cout<<dpp_plus[i]<<" \t"<<dpp_minus[i]<<" \t"<<endl;
					cout<<"#"<<currentDateTime()<<endl;
					cout<<"******************************"<<endl;	
				}
				delete []plist;
				delete []dpp_plus;
				delete []dpp_minus;
			}
			break;
			
			case 3: //track dNudJ
			{
				//cout<<sline<<endl;
				size_t pos = sline.find_first_of(("=:"));
				sline.erase(0,pos+1);
				str_trim_both(sline);
				replace(sline.begin(),sline.end(),',',' ');
				
				string keyws[]={"XMAX","YMAX","NTURN"};
				
				unsigned int nturn;
				double xmax, ymax;
				
				int cnt=0;
				while(sline.length()>0 && cnt++<3)
				{
					pos = sline.find_first_of(("="));
					string valtype = sline.substr(0,pos);
					sline.erase(0,pos+1);
					str_trim_both(sline);
					pos = sline.find_first_of((" \t"));
					str_trim_both(valtype);
					STRINGTOUPPER(valtype);
					if(valtype=="XMAX")
						std::stringstream(sline)>>xmax;
					else if(valtype=="YMAX")
						std::stringstream(sline)>>ymax;
					else if(valtype=="NTURN")
						std::stringstream(sline)>>nturn;
					sline.erase(0,pos+1);
					str_trim_both(sline);					
				}
				cout<<"#"<<xmax<<" "<<ymax<<" "<<nturn<<endl;
				double dnudexy[4];
				pLine->set_monitor(0,0);
				pLine->set_rad_onoff(0); //turn radiation off
				pLine->set_cavity_onoff(0); //turn RF off
				pLine->trackdnudexy(dnudexy,nturn,xmax,ymax);
				double norm;
				norm = sqrt( SQR(dnudexy[0]) + SQR(dnudexy[3]) + SQR((dnudexy[1]+dnudexy[2])/2.0) );
				if(1)
				{
					cout<<"$$$dtuning coefficients"<<endl;
					cout<<dnudexy[0]<<" \t"<<dnudexy[1]<<" \t"<<dnudexy[2]<<" \t"<<dnudexy[3]<<" \t"<<endl;
					cout<<"*************************"<<endl;
				}
				
			}
			break;
		} //switch
	} //while
	return 0;
}

