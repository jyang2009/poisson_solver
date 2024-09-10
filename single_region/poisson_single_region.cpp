#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <stdlib.h>


using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> Mat;
typedef Eigen::Triplet<double> Trip;

int main()
{

int solve_whole(VectorXd& x0,VectorXd& E_field,VectorXd& rho,VectorXd& concentration, double delta_Ev,int n_ele,double* mu,int n_region, double* l,double* h,int* dim,int* dim_start,int* n_defect,int* defect_start,double* E_CBM,double* E_VBM,double* volume,double* epsilon,double* DOS_deno,double* E_MP,int n_MP,int* charge,double* E_defected,double* E_perfect,int* delta_n,int* num_site,double* delta_F_vib,int* DOS_start,double* DOS,double* energy_DOS,int* ENCUT,double temp, int* group_site,double c_dopant1,double c_dopant2);
int generate_mu(double* mu, int n_ele, int n_region, double pressure, double temp);

//input parameters
double temp=500;
//double pressure=1e-5;
double delta_Ev=0.5;


double mu_O_0=0;
double E_DFT_O2=0.096775438-.98542082E+01;
double epsilon0=8.8541878176e-12/(1.602176565e-19);   // F/m
double k_b=8.6173324E-5;

//the unit here should be change to per unit formula
double c_dopant1=0;
double c_dopant2=1e-5;

int array_size=200;
char info[array_size];

//read element list
int n_ele;
ifstream ele_list;
ele_list.open("element_list");
ele_list>>n_ele;
cout<<"n_ele="<<n_ele<<endl;

//read region list
ifstream region_list;
int n_region;
region_list.open("region_list");
region_list>>info;
region_list>>n_region;
region_list.ignore(1,'\n');
region_list.getline(info,array_size);
double l[n_region];
int dim[n_region];
double h[n_region];
for(int i=0;i<n_region;i++)
{
	region_list>>l[i];
	region_list>>dim[i];
	l[i]*=1e-9;
	h[i]=l[i]/(dim[i]+1);
}
cout<<"n_region="<<n_region;
cout<<"l=";
for(int i=0;i<n_region;i++)cout<<l[i]<<"  ";
cout<<endl;
cout<<"h=";
for(int i=0;i<n_region;i++)cout<<h[i]<<"  ";
cout<<endl;
//for each region, read in band, defect, DOSCAR
stringstream num_str;
string num_name;
string name="band";
string filename;
ifstream band_input;
int ENCUT[n_region];
double E_VBM[n_region];
double E_CBM[n_region];
int n_defect[n_region];
int n_MP=5;
double E_MP[n_region*n_MP];
double DOS_deno[n_region];
double volume[n_region];
double dielectric_constant[n_region];

for(int i=0;i<n_region;i++)
{
	num_name.clear();
    num_str.str("x");
    num_str<<i;
    num_name= num_str.str();
	filename=name+num_name;
	band_input.open(filename.c_str());
	band_input>>info;band_input>>ENCUT[i];
	band_input>>info;band_input>>E_CBM[i];
	band_input>>info;band_input>>E_VBM[i];
	band_input>>info;band_input>>n_defect[i];
	band_input>>info;band_input>>dielectric_constant[i];
	band_input>>info;
	for(int j=0;j<n_MP;j++){band_input>>E_MP[i*n_MP+j];E_MP[i*n_MP+j]/=dielectric_constant[i];}
	band_input>>info;band_input>>DOS_deno[i];
	band_input>>info;band_input>>volume[i];
	band_input.close();
}
/*
cout<<"E_CBM=";
for(int i=0;i<n_region;i++)cout<<E_CBM[i]<<"   ";
cout<<endl;
cout<<"volume=";
for(int i=0;i<n_region;i++)cout<<volume[i]<<"   ";
cout<<endl;
*/

 name="DOSCAR";
 ifstream DOSCAR_input;
 int num_DOSCAR_ele=0;
 for(int i=0;i<n_region;i++)num_DOSCAR_ele+=ENCUT[i];
 int DOS_start[n_region];
 DOS_start[0]=0;
 cout<<"DOS_start=";
 for(int i=1;i<n_region;i++)
 {
	 DOS_start[i]=DOS_start[i-1]+ENCUT[i-1];
	 cout<<DOS_start[i]<<"   ";
 }
 cout<<endl;
 double energy_DOS[num_DOSCAR_ele];
 double DOS[num_DOSCAR_ele];
 
 for(int i=0;i<n_region;i++)
{
	num_name.clear();
    num_str.str("x");
    num_str<<i;
    num_name= num_str.str();
	filename=name+num_name;
	DOSCAR_input.open(filename.c_str());
	for(int j=0;j<ENCUT[i];j++)
	{
	 DOSCAR_input>>energy_DOS[(DOS_start[i]+j)];
	 DOSCAR_input>>DOS[(DOS_start[i]+j)];
	 DOSCAR_input.ignore(1,'\n');
	}
	 DOSCAR_input.close();	
}


name="defect";
ifstream defect_input;
int num_defect_ele=0;
for(int i=0;i<n_region;i++)
{
	num_defect_ele+=n_defect[i];
}
int defect_start[n_region];
defect_start[0]=0;
for(int i=1;i<n_region;i++)
{
	defect_start[i]=defect_start[i-1]+n_defect[i-1];
}
int n;
int charge[num_defect_ele];
int delta_n[num_defect_ele*n_ele];
double E_defected[num_defect_ele];
int num_site[num_defect_ele];
double E_perfect[num_defect_ele];
int group_site[num_defect_ele];

for(int i=0;i<n_region;i++)
{
	num_name.clear();
    num_str.str("x");
    num_str<<i;
    num_name= num_str.str();
	filename=name+num_name;
	cout<<"reading..."<<filename<<endl;
	defect_input.open(filename.c_str());
	defect_input.getline(info,array_size);
	for(int j=0;j<n_defect[i];j++)
	{
		defect_input>>n;
		//defect_input>>info;
		defect_input>>charge[defect_start[i]+j];
		for(int k=0;k<n_ele;k++)defect_input>>delta_n[(defect_start[i]+j)*n_ele+k];
		defect_input>>E_defected[defect_start[i]+j];
		defect_input>>E_perfect[defect_start[i]+j];
		defect_input>>num_site[defect_start[i]+j];
		defect_input>>group_site[defect_start[i]+j];
		
	}
	defect_input.close();
}
double delta_F_vib[num_defect_ele];
for(int k=0;k<num_defect_ele;k++)
{
	delta_F_vib[k]=0;
}


stringstream pressure_str;
string pressure_name;


ofstream x0_output;
ofstream rho_output;
ofstream concentration_output;
ofstream field_output;
int dim_tot=0;
for(int i=0;i<n_region;i++)dim_tot+=dim[i];

int dim_start[n_region];
dim_start[0]=0;
for(int i=1;i<n_region;i++)dim_start[i]=dim_start[i-1]+dim[i-1];
VectorXd x0(dim_tot);
VectorXd E_field(dim_tot+n_region);
VectorXd rho(dim_tot);

int concentration_tot=0;
for(int i=0;i<n_region;i++)concentration_tot+=dim[i]*(n_defect[i]+2);
VectorXd concentration(concentration_tot);

double pressure=0.2;
//double pressure[50];
//for(int k=-50;k<0;k++)pressure[k+50]=pow(10,k);
double mu[n_ele*n_region];


    generate_mu(mu,n_ele,n_region,pressure,temp);

    for(int i=0;i<dim_tot;i++)x0(i)=-2.;
    
    //Here initialization of x0 doesn't matter, it's done in solve_whole
    solve_whole(x0, E_field, rho,concentration, delta_Ev,n_ele, mu, n_region,  l, h,dim, dim_start, n_defect, defect_start, E_CBM,E_VBM, volume, dielectric_constant,DOS_deno, E_MP, n_MP, charge, E_defected, E_perfect,delta_n,num_site,delta_F_vib,DOS_start,DOS,energy_DOS,ENCUT,temp,group_site,c_dopant1,c_dopant2);
	x0_output.open("x0_output");
	x0_output<<x0;
	x0_output.close();
	
    concentration_output.open("concentration_output");
	concentration_output<<concentration;
	concentration_output.close();
	
    field_output.open("field_output");
	field_output<<E_field;
	field_output.close();
	
	
	
	
	rho_output.open("rho_output");
	rho_output<<rho;
	rho_output.close();
	
	


return 0;

}

//solve for each region
//x0,E_field,rho has dim_tot dimension
//concentration has dim[i]*(n_defect[i]+2) dimension
int solve_whole(VectorXd& x0,VectorXd& E_field,VectorXd& rho,VectorXd& concentration, double delta_Ev,int n_ele,double* mu,int n_region, double* l,double* h,int* dim,int* dim_start,int* n_defect,int* defect_start,double* E_CBM,double* E_VBM,double* volume,double* epsilon,double* DOS_deno,double* E_MP,int n_MP,int* charge,double* E_defected,double* E_perfect,int* delta_n,int* num_site,double* delta_F_vib,int* DOS_start,double* DOS,double* energy_DOS,int* ENCUT,double temp, int* group_site,double c_dopant1,double c_dopant2)
{
	int solve_system(VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f, int* group_site);
	int find_Ef(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_CBM, int n_defect,int n_ele, int* charge, int* delta_n, double* E_DFT, int* num_site,double* E_MP, double T,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double* Ef, int* group_site);
	int generate_rho(double* Q,VectorXd& concentration,VectorXd& rho,VectorXd& E_field,VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f, int* group_site);
	int generate_region(int i_region,int n_ele,int n_MP,int* n_defect, int* defect_start, int* DOS_start, double* E_MP_region,double* E_MP,double* E_defect_region,double* E_defected,double* E_perfect_region, double* E_perfect, int* charge, int* charge_region,int* num_site_region, int* num_site, double* delta_F_vib_region, double* delta_F_vib, int* delta_n, int* delta_n_region,int* ENCUT,double* DOS_region, double* DOS,double* energy_DOS, double* energy_DOS_region,double* mu_region, double* mu, int* group_site, int* group_site_region);
	int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int* group_site);
	int find_Ef_external(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_CBM, int n_defect,int n_ele, int* charge, int* delta_n, double* E_DFT, int* num_site,double* E_MP, double T,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double* Ef, int* group_site,double Q_dopant);

	//c_dopant1 and c_dopant2 are the dopant concentrations on the two sides, in ZrO2 the doped cation, in water the fictitious A+ or B- ion
	//if c_dopant is zero, do not specify dopant concentration, do as iso-chemical potential case
	//if c_dopant is not zero, do iso-dopant concentration test
	
	ofstream output_solve_whole;
	
	int max_step=150;
	int flag_solve;
	
	int dim_tot=0;
	for(int i=0;i<n_region;i++)dim_tot+=dim[i];
	int dim_max=0;
	for(int i=0;i<n_region;i++){if(dim[i]>dim_max)dim_max=dim[i];}
	int n_defect_max=0;
	for(int i=0;i<n_region;i++){if(n_defect[i]>n_defect_max)n_defect_max=n_defect[i];}
	int ENCUT_max=0;
	for(int i=0;i<n_region;i++){if(ENCUT[i]>ENCUT_max)ENCUT_max=ENCUT[i];}
	VectorXd x0_solved(dim_tot);
	//outputs
	VectorXd x0_region(dim_max); for(int i=0;i<dim_max;i++)x0_region(i)=0;
	VectorXd E_field_region(dim_max+1);
	VectorXd rho_region(dim_max);
	VectorXd concentration_region(dim_max*(2+n_defect_max));
	
	int concentration_start[n_region];
	concentration_start[0]=0;
	for(int i=1;i<n_region;i++)
	{
		concentration_start[i]=concentration_start[i-1]+((2+n_defect[i-1])*dim[i-1]);
	}
	
	//inputs
	double E_MP_region[n_MP]; 
	double E_defect_region[n_defect_max];
	double E_perfect_region[n_defect_max];
	int charge_region[n_defect_max];
	int delta_n_region[n_defect_max*n_ele];
	int num_site_region[n_defect_max];
	int group_site_region[n_defect_max];
	double delta_F_vib_region[n_defect_max];
	double DOS_region[ENCUT_max];
	double energy_DOS_region[ENCUT_max];
	double mu_region[n_ele];
	double mu_l,mu_r;
	double c_dopant_calc;
	int i_region;
	i_region=0;	
	generate_region(i_region,n_ele,n_MP,n_defect,defect_start,DOS_start,E_MP_region,E_MP,E_defect_region,E_defected,E_perfect_region,E_perfect,charge,charge_region,num_site_region,num_site,delta_F_vib_region,delta_F_vib,delta_n,delta_n_region,ENCUT,DOS_region,DOS,energy_DOS, energy_DOS_region, mu_region, mu,group_site,group_site_region);
	//solve for bulk neutrality for the first and last region
    cout<<"region1 bulk"<<endl;
    double E_f1;
    cout<<"E_VBM1="<<E_VBM[i_region]<<endl;
    cout<<"E_CBM1="<<E_CBM[i_region]<<endl;
    cout<<"charge_region";
    cout<<"n_defect="<<n_defect[i_region]<<endl;
    for(int i=0;i<n_defect[i_region];i++)cout<<charge_region[i]<<"    ";
    cout<<endl;

	if(c_dopant1!=0)
	{
		//find chemical potential for dopant that achieves the given concentration
			mu_l=-30.;
			mu_r=-0.1;
			double * c_defect = new double[n_defect[i_region]];
			double * G_defect = new double[n_defect[i_region]];
			while(mu_r-mu_l>0.0001)
			{
				mu_region[n_ele-1]=(mu_r+mu_l)/2.;

				find_Ef(energy_DOS_region, DOS_region, ENCUT[i_region],DOS_deno[i_region], E_CBM[i_region], n_defect[i_region],n_ele, charge_region, delta_n_region, E_defect_region,num_site_region,E_MP_region,temp,mu_region,E_perfect_region,E_VBM[i_region],delta_F_vib_region,&E_f1,group_site_region);
				//cout<<"Ef="<<E_f<<endl;
				defect_concentration(n_ele,n_defect[i_region],charge_region,delta_n_region,E_defect_region,E_perfect_region,num_site_region,E_MP_region,E_f1,temp,mu_region,E_VBM[i_region],delta_F_vib_region,c_defect,G_defect,group_site_region);
				c_dopant_calc=0;
				for(int d=0;d<n_defect[i_region];d++){if(delta_n_region[d*n_ele+n_ele-1]==1){c_dopant_calc+=c_defect[d];}}
				if(c_dopant_calc<c_dopant1){mu_l=(mu_r+mu_l)/2.;}
				else{mu_r=(mu_r+mu_l)/2.;}
				cout<<"c_dopant_calc="<<c_dopant_calc<<endl;
			}
			mu[i_region*n_ele+n_ele-1]=(mu_l+mu_r)/2.;
			mu[(i_region+1)*n_ele+n_ele-1]=(mu_l+mu_r)/2.;
			delete[] c_defect;
			delete[] G_defect;
	}
	else 
	{
		//used chemical potential read from file
		find_Ef(energy_DOS_region, DOS_region, ENCUT[i_region],DOS_deno[i_region], E_CBM[i_region], n_defect[i_region],n_ele, charge_region, delta_n_region, E_defect_region,num_site_region,E_MP_region,temp,mu_region,E_perfect_region,E_VBM[i_region],delta_F_vib_region,&E_f1,group_site_region);
	}
	
	cout<<"E_f1="<<E_f1<<endl;
	ofstream band_output;
	band_output.open("band_output");
	band_output<<E_f1<<"       "<<E_CBM[0]-E_VBM[0]<<endl;
 	double beta=delta_Ev;	
	//in this version, solve for one region
	double E_f[n_region];
	E_f[0]=E_f1;
	double delta_V_left=0.0;
	double delta_V_right=beta;
	double Q[n_region];
	double Q_tot;
	double delta_V_test;
	int num_step=0;
	double E1;
	double E2;
	double V_left[n_region];
	double V_right[n_region];
	
		//write output
	
	    V_left[0]=0;
	    V_right[0]=beta;
		
	for(int i=0;i<n_region;i++)
	{
		i_region=i;
		generate_region(i_region,n_ele,n_MP,n_defect,defect_start,DOS_start,E_MP_region,E_MP,E_defect_region,E_defected,E_perfect_region,E_perfect,charge,charge_region,num_site_region,num_site,delta_F_vib_region,delta_F_vib,delta_n,delta_n_region,ENCUT,DOS_region,DOS,energy_DOS, energy_DOS_region, mu_region, mu,group_site,group_site_region);
		flag_solve=solve_system(x0_region,n_ele,max_step,dim[i_region],h[i_region],V_left[0],V_right[0],volume[i_region],epsilon[i_region],energy_DOS_region,DOS_region,ENCUT[i_region],E_CBM[i_region], n_defect[i_region],charge_region, delta_n_region, E_defect_region,num_site_region,E_MP_region,temp,mu_region,E_perfect_region,E_VBM[i_region],delta_F_vib_region,DOS_deno[i_region],E_f1,group_site_region);        
                        cout<<"solved!!"<<endl;
		for(int j=0;j<dim[i_region];j++)
		{
			x0_solved(dim_start[i_region]+j) = x0_region(j);
		}
		generate_rho(&Q[i_region], concentration_region,rho_region,E_field_region,x0_region,n_ele,max_step,dim[i_region],h[i_region],V_left[i_region],V_right[i_region],volume[i_region],epsilon[i_region],energy_DOS_region,DOS_region,ENCUT[i_region],E_CBM[i_region], n_defect[i_region],charge_region, delta_n_region, E_defect_region,num_site_region,E_MP_region,temp,mu_region,E_perfect_region,E_VBM[i_region],delta_F_vib_region,DOS_deno[i_region],E_f[i_region],group_site_region);
		
		for(int k=0;k<(dim[i_region]+1);k++)E_field(dim_start[i_region]+i_region+k)=E_field_region(k);
		for(int k=0;k<dim[i_region];k++)rho(dim_start[i_region]+k)=rho_region(k);
		for(int k=0;k<dim[i_region];k++)
		{   concentration(concentration_start[i_region]+k)=concentration_region(k);
			concentration(concentration_start[i_region]+dim[i_region]+k)=concentration_region(dim[i_region]+k);
			for(int j=0;j<n_defect[i_region];j++)
			concentration(concentration_start[i_region]+(j+2)*dim[i_region]+k)=concentration_region((j+2)*dim[i_region]+k);
		}
	}
	
			
	
	for(int i=0;i<dim_tot;i++)x0(i)=x0_solved(i);

		
	return 0;
}


int generate_region(int i_region,int n_ele,int n_MP,int* n_defect, int* defect_start, int* DOS_start, double* E_MP_region,double* E_MP,double* E_defect_region,double* E_defected,double* E_perfect_region, double* E_perfect, int* charge, int* charge_region,int* num_site_region, int* num_site, double* delta_F_vib_region, double* delta_F_vib, int* delta_n, int* delta_n_region,int* ENCUT,double* DOS_region, double* DOS,double* energy_DOS, double* energy_DOS_region,double* mu_region, double* mu,int* group_site, int* group_site_region)
{
		//give values for each region
	for(int i=0;i<n_MP;i++)E_MP_region[i]=E_MP[i_region*n_MP+i];
		//defects
	for(int i=0;i<n_defect[i_region];i++)E_defect_region[i]=E_defected[defect_start[i_region]+i];
	for(int i=0;i<n_defect[i_region];i++)E_perfect_region[i]=E_perfect[defect_start[i_region]+i];
	for(int i=0;i<n_defect[i_region];i++)charge_region[i]=charge[defect_start[i_region]+i];
	for(int i=0;i<n_defect[i_region];i++)num_site_region[i]=num_site[defect_start[i_region]+i];
	for(int i=0;i<n_defect[i_region];i++)delta_F_vib_region[i]=delta_F_vib[defect_start[i_region]+i];
	for(int i=0;i<n_defect[i_region];i++)
		{
			for(int k=0;k<n_ele;k++)
			delta_n_region[i*n_ele+k]=delta_n[(defect_start[i_region]+i)*n_ele+k];
	    }
	    //DOS
	for(int i=0;i<ENCUT[i_region];i++)DOS_region[i]=DOS[DOS_start[i_region]+i];
	for(int i=0;i<ENCUT[i_region];i++)energy_DOS_region[i]=energy_DOS[DOS_start[i_region]+i];   
		//mu
	for(int i=0;i<n_ele;i++)mu_region[i]=mu[i_region*n_ele+i];
	for(int i=0;i<n_defect[i_region];i++)group_site_region[i]=group_site[defect_start[i_region]+i];
}

//read mu[n_ele][n_region], calculate minimized mu
int generate_mu(double* mu, int n_ele, int n_region, double pressure, double temp)
{

	double k_b=8.61733238e-5;
	double mu_calc;
	int delta_n_mu[n_ele];   //read chemical species for potential calculations
	int n;
	int n_phase;
	int phase[n_ele*n_region];
	double E_phase_read;
	double E_phase;
	int array_size=100;
	char info[array_size];
	double pressure_dependence;
	ifstream mu_input;
	stringstream pressure_str;
	string pressure_name;
	string filename;
	string name="mu";
	pressure_str.str("x");
	pressure_name= pressure_str.str();
	mu_input.open("mu0");
	mu_input.close();
	for(int i=0;i<n_region;i++)
	{
		for(int j=0;j<n_ele;j++)
		{
			mu[i*n_ele+j]=100.0;
			phase[i*n_ele+j]=-1;
		}
		for(int j=0;j<n_ele;j++)
		{
			pressure_name.clear();
			pressure_str.str("x");
			pressure_str<<j;
			pressure_name= pressure_str.str();
			filename=name+pressure_name;
			
			mu_input.open(filename.c_str());
			mu_input>>n_phase;mu_input.ignore(1,'\n');
			cout<<"reading..."<<filename<<"   n_phase="<<n_phase<<endl;
			mu_input.getline(info,array_size);
			//cout<<info<<endl;
			for(int l=0;l<n_phase;l++)
		 {
		    mu_input>>n;
		    mu_input>>info;
		    for(int k=0;k<n_ele;k++)
		   {
			 mu_input>>delta_n_mu[k];
		   }
		    for(int k=0;k<n_region;k++)
		    {
		       mu_input>>E_phase_read;
		       if(k==i)E_phase=E_phase_read;
		    }
		     mu_input>>pressure_dependence;
		     cout<<"For element"<<j<<", phase"<<l<<"  ";
		     for(int k=0;k<n_ele;k++)cout<<delta_n_mu[k]<<"   ";
		     cout<<"E_phase="<<E_phase<<"   ";
		     cout<<"pressure="<<pressure_dependence<<endl;
		     mu_calc=E_phase;
		     for(int k=0;k<j;k++)
		     {
				if(k!=j)mu_calc-=delta_n_mu[k]*mu[i*n_ele+k];   
			 }
			 mu_calc+=pressure_dependence*k_b*temp*log(pressure);
			 mu_calc/=delta_n_mu[j];
		     if(mu_calc<=mu[i*n_ele+j]){phase[i*n_ele+j]=n;mu[i*n_ele+j]=mu_calc;}
		 }
			
			if(i==0||i==1){mu[i*n_ele+0]=(-0.547861+k_b*temp*log(pressure)-9.8580501+1.22)/2.;}
			
			mu_input.close();				
		}
			for(int j=0;j<n_ele;j++)
	    {
		    cout<<"calc_mu  region"<<i<<" element"<<j<<"  n_phase"<<"    "<<phase[i*n_ele+j]<<"    mu="<<mu[i*n_ele+j]<<endl;
	    }
		
	}
    
}

int el_concentration(double* energy_DOS, double* DOS, int ENCUT, int DOS_deno,double E_VBM, double E_CBM, double E_f, double T,double* c_e)
{
    double c=0;
    double k_b=8.61733238e-5;
    for(int k=1;k<ENCUT;k++)
    {
        if(energy_DOS[k]>E_CBM)
        {
            //c=c+DOS[k]*(energy_DOS[k]-energy_DOS[k-1])/(1+exp((energy_DOS[k]-E_f)/(k_b*T)));
            c=c+DOS[k]*(15./ENCUT)/(1+exp((energy_DOS[k]-(E_f+E_VBM))/(k_b*T)));
        }
    }
    *c_e=c/DOS_deno;
return 1;
}

int h_concentration(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_VBM, double E_CBM,double E_f, double T,double* c_h)
{
  double c=0;
    double k_b=8.61733238e-5;
    for(int k=1;k<ENCUT;k++)
    {
        if(energy_DOS[k]<E_VBM)
        {
            //c=c+DOS[k]*(energy_DOS[k]-energy_DOS[k-1])/(1+exp((-energy_DOS[k]+E_f)/(k_b*T)));
            c=c+DOS[k]*(15./ENCUT)/(1+exp((-energy_DOS[k]+(E_f+E_VBM))/(k_b*T)));
        }
    }
    *c_h=c/DOS_deno;
return 1;
}


int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int* group_site)
{

    double k_b=8.61733238e-5;
    int q;
    double deno;
    //calculate formation energy
    for(int k=0;k<n_defect;k++)
    {
		G_defect[k]=E_DFT[k]-E_perfect[k]+delta_F_vib[k]+charge[k]*(E_f+E_VBM);
		for(int j=0;j<n_ele;j++)
		{
			G_defect[k]-=delta_n[k*n_ele+j]*mu[j];
		}
        q=charge[k];
        q=abs(q);
        //cout<<"q="<<q<<endl;
        G_defect[k]+=E_MP[q];
        
        //c[k]=num_site[k]*exp(-G_defect[k]/(k_b*T));
        //cout<<"G_defect[k]="<<G_defect[k]<<"   c[k]="<<c[k]<<endl;
    }
    
        for(int k=0;k<n_defect;k++)
    {
        deno=1;
        for(int j=0;j<n_defect;j++)
        {
			if((group_site[k]==group_site[j])&&(j!=k))
			{
                deno=deno+exp(-G_defect[j]/(k_b*T));
			}
        }

        c[k]=num_site[k]*exp(-G_defect[k]/(k_b*T))/deno;
    }
    
    
return 1;
}




int rho(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_CBM, int n_defect, int n_ele,int* charge, int* delta_n,double* E_DFT, int* num_site,double* E_MP,double E_f, double T,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double* rho, int* group_site)
{
    int el_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int h_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int*);
    double c_rho[n_defect];
    double G_rho[n_defect];
    double c_e;
    double c_h;
    double charge_density=0;
    el_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f,T,&c_e);
    charge_density-=c_e;
    h_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f,T,&c_h);
    charge_density+=c_h;
    defect_concentration(n_ele,n_defect,charge,delta_n,E_DFT,E_perfect,num_site,E_MP,E_f,T,mu,E_VBM,delta_F_vib,c_rho,G_rho, group_site);
    for(int k=0;k<n_defect;k++)
    {
        charge_density+=charge[k]*c_rho[k];
    }
    *rho=charge_density;
    return 1;
}

int find_Ef(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_CBM, int n_defect,int n_ele, int* charge, int* delta_n, double* E_DFT, int* num_site,double* E_MP, double T,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double* Ef, int* group_site)
{
    int rho(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_CBM, int n_defect, int n_ele,int* charge, int* delta_n,double* E_DFT, int* num_site,double* E_MP,double E_f, double T,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double* rho, int* group_site);
    int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int* group_site);
    double E_f;
    double prec=1e-4;
    int num_iter=0;
    double rho1,rho2,rho3;
    double E_l=0;
    double E_r=E_CBM-E_VBM;
    double c_rho[n_defect];
    double G_rho[n_defect];
    ofstream iter_result;
    iter_result.open("iter_result.txt");
    E_f=(E_l+E_r)/2;
    int flag=1;   //flag=1 for negative formation energy exist
    while((flag==1)&&(E_l<(E_CBM-E_VBM)))
    {
		//check for E_l, avoid negative formation energy
		defect_concentration(n_ele,n_defect,charge,delta_n,E_DFT,E_perfect,num_site,E_MP,E_l,T,mu,E_VBM,delta_F_vib,c_rho,G_rho, group_site);
		for(int k=0;k<n_defect;k++)
		{
			if(G_rho[k]<0){E_l+=prec;break;}
		    if(k==(n_defect-1))flag=0;
		}
    }
    if(flag==1)cout<<"E_l doesn't exist!!"<<endl;
    for(int k=0;k<n_defect;k++)
    {
		if(G_rho[k]<0)cout<<"defect"<<k<<" is negative!!"<<endl;
	}
    flag=1;
    while((flag==1)&&(E_r>0))
    {
		//check for E_r, avoid negative formation energy
		defect_concentration(n_ele,n_defect,charge,delta_n,E_DFT,E_perfect,num_site,E_MP,E_r,T,mu,E_VBM,delta_F_vib,c_rho,G_rho, group_site);
		for(int k=0;k<n_defect;k++)
		{
			if(G_rho[k]<0){E_r-=prec;break;}
		    if(k==(n_defect-1))flag=0;
		}
    }
    if(flag==1)cout<<"E_r doesn't exist!!"<<endl;
        for(int k=0;k<n_defect;k++)
    {
		if(G_rho[k]<0)cout<<"defect"<<k<<" is negative!!"<<endl;
	}
    //rho1(E_l)
    rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_l,T,mu,E_perfect,E_VBM,delta_F_vib,&rho1,group_site);
    //rho2(E_r)
    rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_r,T,mu,E_perfect,E_VBM,delta_F_vib,&rho2, group_site);
    //cout<<"rho1="<<rho1<<"      "<<"rho2="<<rho2<<endl;
    
    if(rho1*rho2>0)
    {
        *Ef=E_l;
        return 0;
    }

    while(abs(rho1)>1e-50)
    {
        E_f=(E_l+E_r)/2;
        //rho1(E_f)
        rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_f,T,mu,E_perfect,E_VBM,delta_F_vib,&rho1,group_site);
        //rho2(E_r)
        rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_r,T,mu,E_perfect,E_VBM,delta_F_vib,&rho2,group_site);
        //rho3(E_l)
        rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_l,T,mu,E_perfect,E_VBM,delta_F_vib,&rho3,group_site);
        iter_result<<num_iter<<setprecision(10) <<"  E_l= "<<E_l<<"  E_r="<<E_r<<"  E_f="<<E_f<<"     rho="<<rho1<<"  rho_l="<<rho3<<"   rho_r="<<rho2<<endl;
        if(rho1*rho2<0)E_l=E_f;
        else E_r=E_f;

        num_iter+=1;
        if((abs(E_r-E_l)<1.0E-15)){E_f=(E_l+E_r)/2;break;}
    }
	//cout<<"final rho="<<rho3<<endl;
    *Ef=E_f;
    return 1;
}

int find_Ef_external(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_CBM, int n_defect,int n_ele, int* charge, int* delta_n, double* E_DFT, int* num_site,double* E_MP, double T,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double* Ef, int* group_site,double Q_dopant)
{
    int rho(double* energy_DOS, double* DOS, int ENCUT,int DOS_deno, double E_CBM, int n_defect, int n_ele,int* charge, int* delta_n,double* E_DFT, int* num_site,double* E_MP,double E_f, double T,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double* rho, int* group_site);
    int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int* group_site);
    double E_f;
    double prec=1e-4;
    int num_iter=0;
    double rho1,rho2,rho3;
    double E_l=0;
    double E_r=E_CBM-E_VBM;
    double c_rho[n_defect];
    double G_rho[n_defect];
    ofstream iter_result;
    iter_result.open("iter_result.txt");
    E_f=(E_l+E_r)/2;
    int flag=1;   //flag=1 for negative formation energy exist
    while((flag==1)&&(E_l<(E_CBM-E_VBM)))
    {
		//check for E_l, avoid negative formation energy
		defect_concentration(n_ele,n_defect,charge,delta_n,E_DFT,E_perfect,num_site,E_MP,E_l,T,mu,E_VBM,delta_F_vib,c_rho,G_rho, group_site);
		for(int k=0;k<n_defect;k++)
		{
			if(G_rho[k]<0){E_l+=prec;break;}
		    if(k==(n_defect-1))flag=0;
		}
    }
    if(flag==1){cout<<"E_l doesn't exist!!"<<endl;E_l=0;}
    for(int k=0;k<n_defect;k++)
    {
		if(G_rho[k]<0)cout<<"defect"<<k<<" is negative!!"<<endl;
	}
    flag=1;
    while((flag==1)&&(E_r>0))
    {
		//check for E_r, avoid negative formation energy
		defect_concentration(n_ele,n_defect,charge,delta_n,E_DFT,E_perfect,num_site,E_MP,E_r,T,mu,E_VBM,delta_F_vib,c_rho,G_rho, group_site);
		for(int k=0;k<n_defect;k++)
		{
			if(G_rho[k]<0){E_r-=prec;break;}
		    if(k==(n_defect-1))flag=0;
		}
    }
    if(flag==1){cout<<"E_r doesn't exist!!"<<endl;E_r=E_CBM-E_VBM;}
        for(int k=0;k<n_defect;k++)
    {
		if(G_rho[k]<0)cout<<"defect"<<k<<" is negative!!"<<endl;
	}
    //rho1(E_l)
    rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect-1,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_l,T,mu,E_perfect,E_VBM,delta_F_vib,&rho1,group_site);
    //rho2(E_r)
    rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect-1,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_r,T,mu,E_perfect,E_VBM,delta_F_vib,&rho2, group_site);
    //cout<<"rho1="<<rho1<<"      "<<"rho2="<<rho2<<endl;
    rho1=rho1+Q_dopant;
    rho2=rho2+Q_dopant;
    if(rho1*rho2>0)
    {
        *Ef=E_l;
        return 0;
    }

    while(abs(rho1)>1e-50)
    {
        E_f=(E_l+E_r)/2;
        //rho1(E_f)
        rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect-1,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_f,T,mu,E_perfect,E_VBM,delta_F_vib,&rho1,group_site);
        //rho2(E_r)
        rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect-1,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_r,T,mu,E_perfect,E_VBM,delta_F_vib,&rho2,group_site);
        //rho3(E_l)
        rho(energy_DOS,DOS,ENCUT,DOS_deno,E_CBM,n_defect-1,n_ele,charge,delta_n,E_DFT,num_site,E_MP,E_l,T,mu,E_perfect,E_VBM,delta_F_vib,&rho3,group_site);
        iter_result<<num_iter<<setprecision(10) <<"  E_l= "<<E_l<<"  E_r="<<E_r<<"  E_f="<<E_f<<"     rho="<<rho1<<"  rho_l="<<rho3<<"   rho_r="<<rho2<<endl;
        rho1=rho1+Q_dopant;
        rho2=rho2+Q_dopant;
        rho3=rho3+Q_dopant;
        if(rho1*rho2<0)E_l=E_f;
        else E_r=E_f;

        num_iter+=1;
        if((abs(E_r-E_l)<1.0E-15)){E_f=(E_l+E_r)/2;break;}
    }
	//cout<<"final rho="<<rho3<<endl;
    *Ef=E_f;
    return 1;
}

//x0 is an input
//x0 has dimension dim*(3+n_defect)
int solve_system(VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f,int* group_site)
{

	int generate_Fx(VectorXd& F_x,VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f, int* group_site);
	int generate_Jacobi(Mat& jacobi,VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f, int* group_site);	
	ofstream solve_system;
	
	cout<<"in solve_system..."<<endl;
	/*
	cout<<"volume="<<volume<<endl;
	cout<<"energy_DOS[0]="<<energy_DOS[0]<<endl;
	cout<<"delta_n=";
	for(int i=0;i<(n_defect*n_ele);i++)cout<<delta_n[i]<<"   ";cout<<endl;
    */
    int dim_J=dim;
	cout<<"dim="<<dim<<endl;
	VectorXd x0_solve(dim_J);
	for(int i=0;i<dim;i++)x0_solve(i)=x0(i);
	//initial guess : linear
	for(int i=0;i<dim;i++)x0_solve(i)=V_left+(V_right-V_left)*i/(dim+1);
    VectorXd x1(dim_J);
    VectorXd F_x0(dim_J);
    VectorXd F_x1(dim_J);
	
	Mat Jacobi_sparse(dim_J,dim_J);
	SimplicialLDLT<SparseMatrix<double> > solver;
    VectorXd y(dim_J);
    double abs_F0;
    double abs_F1;
    VectorXd dx(dim_J);
    int num_itr=0;
    double u;
    
    //cout<<"entering iteration"<<endl;
    for(int i=0;i<dim;i++)
    {
		y(i)=1;
		
	}
	//give a linear guess
	for(int i=0;i<dim;i++)x0(i)=V_left+(V_right-V_left)*double(i+1)/double(dim+1);
	
	generate_Fx(F_x0,x0_solve,n_ele,max_step,dim,h,V_left,V_right,volume,epsilon,energy_DOS,DOS,ENCUT,E_CBM,n_defect,charge, delta_n, E_defected,num_site,E_MP,temp,mu,E_perfect,E_VBM,delta_F_vib,DOS_deno,E_f,group_site);
	cout<<"Fx generated!!"<<endl;
	//cout<<"F_x0="<<F_x0<<endl;
	abs_F0=F_x0.transpose()*F_x0;
	while(((abs(y.maxCoeff())>1e-10||abs(y.minCoeff())>1e-10)&&num_itr<max_step)||num_itr==0)
	{
		num_itr+=1;
		//compute f_x0(x0)

		 //cout<<"Fx generated!!"<<endl;
		 Jacobi_sparse.setZero();
		 //compute jacobi matrix of x0
		 generate_Jacobi(Jacobi_sparse,x0_solve,n_ele,max_step,dim,h,V_left,V_right,volume,epsilon,energy_DOS,DOS,ENCUT,E_CBM,n_defect,charge, delta_n, E_defected,num_site,E_MP,temp,mu,E_perfect,E_VBM,delta_F_vib,DOS_deno,E_f,group_site);
		 cout<<"Jacobi generated!!"<<endl;
		 //cout<<Jacobi_sparse<<endl;
		 //cout<<"Jacobi genrated!!"<<endl;
	    
	 
		solver.compute(Jacobi_sparse);
		
        if(solver.info()!=Success){cout<<"solver failed!!"<<endl;return 0;}	
		dx=solver.solve(F_x0);

		//dx=Jacobi_sparse.transpose()*F_x0;
		 u=1.;
		 y=u*dx;
		 while((abs(y.maxCoeff())>1e-2)||abs(y.minCoeff())>1e-2)
		 {
			 u=u/2.;
			 y=u*dx;	 
		 }
		 x1=x0_solve-u*dx;
		generate_Fx(F_x1,x1,n_ele,max_step,dim,h,V_left,V_right,volume,epsilon,energy_DOS,DOS,ENCUT,E_CBM,n_defect,charge, delta_n, E_defected,num_site,E_MP,temp,mu,E_perfect,E_VBM,delta_F_vib,DOS_deno,E_f,group_site);
		abs_F1=F_x1.transpose()*F_x1;
		// while((abs(F_x1.maxCoeff())>abs(F_x0.maxCoeff()))&&(abs(F_x1.minCoeff())>abs(F_x0.minCoeff())))
		
		while((abs_F1>abs_F0))
		 {
			 
			 u=u/2.;
			 x1=x0_solve-u*dx;
			 generate_Fx(F_x1,x1,n_ele,max_step,dim,h,V_left,V_right,volume,epsilon,energy_DOS,DOS,ENCUT,E_CBM,n_defect,charge, delta_n, E_defected,num_site,E_MP,temp,mu,E_perfect,E_VBM,delta_F_vib,DOS_deno,E_f,group_site);
			 abs_F1=F_x1.transpose()*F_x1;
			 if(u<1e-8)continue;
		 }
		
		 //compute y
		 y=x1-x0_solve;
   
			cout<<"num_itr="<<num_itr<<endl;
		 cout<<"u="<<u<<endl;
		 cout<<"y_max="<<y.maxCoeff()<<endl;
		 cout<<"y_min="<<y.minCoeff()<<endl;
		 //cout<<endl;
		 
		cout<<"abs_F="<<abs_F1<<endl;
        cout<<endl;

        for(int i=0;i<dim_J;i++)
        {
        x0_solve(i)=x1(i);
        F_x0(i)=F_x1(i);
		}
		abs_F0=abs_F1;
		
		
	}
	//x0 gives the final solution after iteration
	for(int i=0;i<dim;i++)x0(i)=x0_solve(i);
	solve_system.open("solve_system_x0");
	solve_system<<x0;
	solve_system.close();

    return 1;
}

//Fx has dimension dim*(3+n_defect)
//input Ef is Ef at its bulk value
//pay attention to the face: the order is always from bulk to interface!! not from right to left!!
//VBC_bulk=0 is the reference for this specific side
int generate_Fx(VectorXd& F_x,VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f,int* group_site)
{
    int el_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int h_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int* group_site);
	
	double c_defect[n_defect];
	double G_defect[n_defect];
	double c_e,c_h;
	
	double epsilon0=8.8541878176e-12/(1.602176565e-19);   // F/m
	double k_T=temp*8.6173325e-5;
	
		VectorXd f_x(dim);
		for(int i=0;i<dim;i++)
		{
			    el_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+x0(i),temp,&c_e);
				h_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+x0(i),temp,&c_h);
				defect_concentration(n_ele,n_defect,charge,delta_n,E_defected,E_perfect,num_site,E_MP,E_f+x0(i),temp,mu,E_VBM,delta_F_vib,c_defect,G_defect,group_site);
			f_x(i)=-c_e+c_h;
			for(int j=0;j<n_defect;j++){f_x(i)+=charge[j]*c_defect[j];}
			f_x(i)=f_x(i)/volume;	
		}
		F_x(0)=-epsilon0*epsilon*(V_left-2*x0(0)+x0(1))/pow(h,2)-f_x(0);
		for(int i=1;i<(dim-1);i++)
		{
			F_x(i)=-epsilon0*epsilon*(x0(i-1)-2*x0(i)+x0(i+1))/pow(h,2)-f_x(i);
		}
		F_x(dim-1)=-epsilon0*epsilon*(x0(dim-2)-2*x0(dim-1)+V_right)/pow(h,2)-f_x(dim-1);

 
return 0;
}



//generate jacobi matrix analytically

int generate_Jacobi(Mat& jacobi,VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f,int* group_site)
{
    int el_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int h_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int* group_site);
	
	double c_defect[n_defect];
	double G_defect[n_defect];
	double c_e,c_h;
	
	double epsilon0=8.8541878176e-12/(1.602176565e-19);   // F/m
	
	int total_trip=3*dim-2;
	double k_T=temp*8.6173325e-5;
	jacobi.setZero();
	vector<Trip> tripletList;
    tripletList.reserve(total_trip);
	
	double T=temp;

	
	VectorXd f_x(dim);
	VectorXd df_x(dim);
	
    		for(int i=0;i<dim;i++)
		{
			el_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+x0(i),T,&c_e);
			h_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+x0(i),T,&c_h);
			defect_concentration(n_ele,n_defect,charge,delta_n,E_defected,E_perfect,num_site,E_MP,E_f+x0(i),T,mu,E_VBM,delta_F_vib,c_defect,G_defect,group_site);

			f_x(i)=-c_e+c_h;df_x(i)=-c_e/k_T-c_h/k_T;
			for(int j=0;j<n_defect;j++){f_x(i)+=charge[j]*c_defect[j];df_x(i)-=abs(charge[j])*c_defect[j]/k_T;}
			f_x(i)=f_x(i)/volume;	df_x(i)=df_x(i)/volume;
		}
		tripletList.push_back(Trip(0,0,2*epsilon*epsilon0/pow(h,2)-df_x(0)));
		tripletList.push_back(Trip(0,1,-1*epsilon*epsilon0/pow(h,2)));
		for(int i=1;i<dim-1;i++)
		{

		tripletList.push_back(Trip(i,i-1,-1*epsilon*epsilon0/pow(h,2)));
		tripletList.push_back(Trip(i,i,2*epsilon*epsilon0/pow(h,2)-df_x(i)));
		tripletList.push_back(Trip(i,i+1,-1*epsilon*epsilon0/pow(h,2)));

		}
		tripletList.push_back(Trip(dim-1,dim-1,2*epsilon*epsilon0/pow(h,2)-df_x(dim-1)));
		tripletList.push_back(Trip(dim-1,dim-2,-1*epsilon*epsilon0/pow(h,2)));
		
		jacobi.setFromTriplets(tripletList.begin(), tripletList.end());
	return 0;
}

//concentration has dimension dim*(n_defect+2)
//x0 have dimension dim
//E_field have dimension dim+1
//rho have dimemsion dim+2
int generate_rho(double* Q,VectorXd& concentration,VectorXd& rho,VectorXd& E_field,VectorXd& x0,int n_ele, int max_step,int dim,double h,double V_left, double V_right,double volume,double epsilon,double* energy_DOS, double* DOS, int ENCUT, double E_CBM, int n_defect, int* charge, int* delta_n, double* E_defected, int* num_site,double* E_MP, double temp,double* mu,double* E_perfect,double E_VBM,double* delta_F_vib,double DOS_deno,double E_f, int* group_site)
{
    int el_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int h_concentration(double* , double* , int,int , double , double , double ,double,double*);
    int defect_concentration(int n_ele, int n_defect, int* charge, int* delta_n, double* E_DFT,double* E_perfect, int* num_site,double* E_MP,double E_f, double T,double* mu,double E_VBM,double* delta_F_vib,double* c,double* G_defect, int*);
	
	double c_defect[n_defect];
	double G_defect[n_defect];
	double c_e,c_h;
	
	double epsilon0=8.8541878176e-12/(1.602176565e-19);   // F/m
	double k_T=temp*8.6173325e-5;
	double T=temp;
	double Q1,Q2;
	//first point
	el_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+V_left,T,&c_e);
    h_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+V_left,T,&c_h);
    defect_concentration(n_ele,n_defect,charge,delta_n,E_defected,E_perfect,num_site,E_MP,E_f+V_left,T,mu,E_VBM,delta_F_vib,c_defect,G_defect,group_site);
	Q1=-c_e+c_h;
	for(int j=0;j<n_defect;j++){Q1+=charge[j]*c_defect[j];}
		for(int i=0;i<dim;i++)
		{
			    el_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+x0(i),T,&c_e);
				h_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+x0(i),T,&c_h);
				defect_concentration(n_ele,n_defect,charge,delta_n,E_defected,E_perfect,num_site,E_MP,E_f+x0(i),T,mu,E_VBM,delta_F_vib,c_defect,G_defect,group_site);
			concentration(i)=c_e;
			concentration(dim+i)=c_h;
			for(int j=0;j<n_defect;j++)concentration((j+2)*dim+i)=c_defect[j];
			rho(i)=-c_e+c_h;
			for(int j=0;j<n_defect;j++){rho(i)+=charge[j]*c_defect[j];}
			rho(i)=rho(i);	
		}
	//last point
	el_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+V_right,T,&c_e);
    h_concentration(energy_DOS,DOS,ENCUT,DOS_deno,E_VBM,E_CBM,E_f+V_right,T,&c_h);
    defect_concentration(n_ele,n_defect,charge,delta_n,E_defected,E_perfect,num_site,E_MP,E_f+V_right,T,mu,E_VBM,delta_F_vib,c_defect,G_defect,group_site);
	Q2=-c_e+c_h;
	for(int j=0;j<n_defect;j++){Q2+=charge[j]*c_defect[j];}
		
		E_field(0)=-(x0(0)-V_left)/h;
		for(int i=1;i<dim;i++)E_field(i)=-(x0(i)-x0(i-1))/h;
		E_field(dim)=-(V_right-x0(dim-1))/h;
		
	double Q_tot;
	Q_tot=(Q1+Q2)*h/2.;
	for(int i=0;i<dim;i++)Q_tot+=rho(i)*h;
	for(int i=0;i<dim;i++)rho(i)=rho(i)/volume;
	
		
 *Q=Q_tot/volume;
return 0;
}


