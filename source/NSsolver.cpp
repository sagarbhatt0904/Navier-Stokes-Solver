#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "init.h"
#include "gridgen.h"
#include "metric.h"
#include "RHS.h"
#include "BC.h"
#include "contour.h"

using namespace std;

int main()
{
	
	int N= 61;                          /*Grid size`*/
	double st=1; 			  /*Stretching ratio, change this according to need*/
	double Re=100; 		/*Reynolds number*/
	double eps=0.01;   /*dissipation factor*/
	double ep=0.01;		/*residual smoothing factor*/    
	double L=30; /*Length of the domain*/
	double t=0;
	double t1=0.1;		/*Time duration, in seconds, for the study. Put time according to need or put a very large number for steady state problems*/
	double dt=0.005;	/*Time step*/
	double xmax=2*M_PI;	/* range of the grid in x-direction, change this according to the problem*/
	double ymax=2*M_PI; /* range of the grid in y-direction, change this according to the problem*/
	vector<vector<double> > x (N,vector<double>(N, 0));
	vector<vector<double> > y (N,vector<double>(N, 0));
	vector<vector<double> > xvel (N,vector<double>(N, 0));	
	vector<vector<double> > yvel (N,vector<double>(N, 0));
	vector<vector<double> > xvel1 (N,vector<double>(N, 0));
	vector<vector<double> > yvel1 (N,vector<double>(N, 0));
	vector<vector<double> > Press (N,vector<double>(N, 0));
	
	
	gridgen(N, st, xmax, ymax, x, y);  // Grid generation
	
	
	init( N, Re, t, t1, x, y, xvel, xvel1, yvel, yvel1, Press); // Initial conditions
	
	vector<vector<double> > u_new (N,vector<double>(N, 0));
	vector<vector<double> > u_new1 (N,vector<double>(N, 0));
	vector<vector<double> > u_new2 (N,vector<double>(N, 0));
	vector<vector<double> > u_new3 (N,vector<double>(N, 0));
	vector<vector<double> > u_k (N,vector<double>(N, 0));
	vector<vector<double> > u_k1 (N,vector<double>(N, 0));
	vector<vector<double> > u_k2 (N,vector<double>(N, 0));
	vector<vector<double> > u_k3 (N,vector<double>(N, 0));
	vector<vector<double> > u_n (N,vector<double>(N, 0));
	vector<vector<double> > u_n1 (N,vector<double>(N, 0));
	vector<vector<double> > u_n2 (N,vector<double>(N, 0));
	vector<vector<double> > u_n3 (N,vector<double>(N, 0));
	vector<vector<double> > u_old (N,vector<double>(N, 0));
	vector<vector<double> > u_old1 (N,vector<double>(N, 0));
	vector<vector<double> > u_old2 (N,vector<double>(N, 0));
	vector<vector<double> > u_old3 (N,vector<double>(N, 0));

	vector<vector<double> > v_new (N,vector<double>(N, 0));
	vector<vector<double> > v_new1 (N,vector<double>(N, 0));
	vector<vector<double> > v_new2 (N,vector<double>(N, 0));
	vector<vector<double> > v_new3 (N,vector<double>(N, 0));
	vector<vector<double> > v_k (N,vector<double>(N, 0));
	vector<vector<double> > v_k1 (N,vector<double>(N, 0));
	vector<vector<double> > v_k2 (N,vector<double>(N, 0));
	vector<vector<double> > v_k3 (N,vector<double>(N, 0));
	vector<vector<double> > v_n (N,vector<double>(N, 0));
	vector<vector<double> > v_n1 (N,vector<double>(N, 0));
	vector<vector<double> > v_n2 (N,vector<double>(N, 0));
	vector<vector<double> > v_n3 (N,vector<double>(N, 0));
	vector<vector<double> > v_old (N,vector<double>(N, 0));
	vector<vector<double> > v_old1 (N,vector<double>(N, 0));
	vector<vector<double> > v_old2 (N,vector<double>(N, 0));
	vector<vector<double> > v_old3 (N,vector<double>(N, 0));

	vector<vector<double> > p_new (N,vector<double>(N, 0));
	vector<vector<double> > p_new1 (N,vector<double>(N, 0));
	vector<vector<double> > p_new2 (N,vector<double>(N, 0));
	vector<vector<double> > p_new3 (N,vector<double>(N, 0));
	vector<vector<double> > p_k (N,vector<double>(N, 0));
	vector<vector<double> > p_k1 (N,vector<double>(N, 0));
	vector<vector<double> > p_k2 (N,vector<double>(N, 0));
	vector<vector<double> > p_k3 (N,vector<double>(N, 0));
	vector<vector<double> > p_n (N,vector<double>(N, 0));
	vector<vector<double> > p_n1 (N,vector<double>(N, 0));
	vector<vector<double> > p_n2 (N,vector<double>(N, 0));
	vector<vector<double> > p_n3 (N,vector<double>(N, 0));
	vector<vector<double> > p_old (N,vector<double>(N, 0));
	vector<vector<double> > p_old1 (N,vector<double>(N, 0));
	vector<vector<double> > p_old2 (N,vector<double>(N, 0));
	vector<vector<double> > p_old3 (N,vector<double>(N, 0));
	vector<vector<double> > JC (N,vector<double>(N, 0));
	vector<vector<double> > ex (N,vector<double>(N, 0));
	vector<vector<double> > ey (N,vector<double>(N, 0));
	vector<vector<double> > zx (N,vector<double>(N, 0));
	vector<vector<double> > zy (N,vector<double>(N, 0));
	vector<vector<double> > rus (N,vector<double>(N, 0));
	vector<vector<double> > rvs (N,vector<double>(N, 0));
	vector<vector<double> > rcs (N,vector<double>(N, 0));
	vector<vector<double> > rho1 (N,vector<double>(N, 0));
	vector<vector<double> > rho2(N,vector<double>(N, 0));

	vector<vector<double> > RHSu11(N,vector<double>(N, 0));
	vector<vector<double> > RHSu22(N,vector<double>(N, 0));
	vector<vector<double> > RHSu33(N,vector<double>(N, 0));
	vector<vector<double> > RHSv11(N,vector<double>(N, 0));
	vector<vector<double> > RHSv22(N,vector<double>(N, 0));
	vector<vector<double> > RHSv33(N,vector<double>(N, 0));
	vector<vector<double> > RHSu(N,vector<double>(N, 0));
	vector<vector<double> > RHSv(N,vector<double>(N, 0));


	metric( N, x, y,zx,zy,ex,ey,JC);   //Metric Calculation

	// Initializing all the velocity variables and dtau
	u_new1=xvel;
	u_new2=xvel;
	u_new3=xvel;
	u_new=xvel;
	u_k1=xvel;
	u_k2=xvel;
	u_k3=xvel;
	u_k=xvel;
	u_n1=xvel;
	u_n2=xvel;
	u_n3=xvel;
	u_n=xvel;
	u_old1=xvel;
	u_old2=xvel;
	u_old3=xvel;
	u_old=xvel;
	v_new1=yvel;
	v_new2=yvel;
	v_new3=yvel;
	v_new=yvel;
	v_k1=yvel;
	v_k2=yvel;
	v_k3=yvel;
	v_k=yvel;
	v_n1=yvel;
	v_n2=yvel;
	v_n3=yvel;
	v_n=yvel;
	v_old1=yvel;
	v_old2=yvel;
	v_old3=yvel;
	v_old=yvel;
	p_new1=Press;
	p_new2=Press;
	p_new3=Press;
	p_new=Press;
	p_n=Press;
	double dtau=0.00005;

	double dj=1;

    for(double ti=0; ti<t1;ti+=dt)
    {
		
		//while( dj>=0.00001)
		for (int lop = 0; lop < 100; ++lop)
		{
		 	 

		 	/*Fourth order Runge-Kutta*/

		    RHS(N,JC,u_k,v_k,Press,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// First step of RK
		    {    for (int j =1; j<N-1; j++)
		        {
		            RHSu[i][j]=(((-3*u_k[i][j]+4*u_n[i][j]-u_old[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv[i][j]=(((-3*v_k[i][j]+4*v_n[i][j]-v_old[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new1[i][j]=Press[i][j]+0.25*(dtau*rcs[i][j]);
		            u_new1[i][j]=u_k[i][j]+0.25*(dtau*RHSu[i][j]);
		            v_new1[i][j]=v_k[i][j]+0.25*(dtau*RHSv[i][j]);
		        }
		    }
		    BC(N,u_new1,v_new1,p_new1);			// BC after first step RK
		    
		    RHS(N,JC,u_new1,v_new1,p_new1,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2);  
		    
		    for (int i =1; i<N-1; i++)			// Second step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {
		            RHSu11[i][j]=(((-3*u_k1[i][j]+4*u_n1[i][j]-u_old1[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv11[i][j]=(((-3*v_k1[i][j]+4*v_n1[i][j]-v_old1[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new2[i][j]=Press[i][j]+0.33*(dtau*rcs[i][j]);
		            u_new2[i][j]=u_k[i][j]+0.33*(dtau*RHSu11[i][j]);
		            v_new2[i][j]=v_k[i][j]+0.33*(dtau*RHSv11[i][j]);
		        }
		    }
		    BC(N,u_new2,v_new2,p_new2);			// BC after second step RK
		    
		    RHS(N,JC,u_new2,v_new2,p_new2,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// Third step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {   
		            RHSu22[i][j]=(((-3*u_k2[i][j]+4*u_n2[i][j]-u_old2[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv22[i][j]=(((-3*v_k2[i][j]+4*v_n2[i][j]-v_old2[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new3[i][j]=Press[i][j]+0.5*(dtau*rcs[i][j]);
		            u_new3[i][j]=u_k[i][j]+0.5*(dtau*RHSu22[i][j]);
		            v_new3[i][j]=v_k[i][j]+0.5*(dtau*RHSv22[i][j]);
		        }
		    }
		    BC(N,u_new3,v_new3,p_new3);			// BC after third step RK
		    
		    RHS(N,JC,u_new3,v_new3,p_new3,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// Fourth step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {   
		            RHSu33[i][j]=(((-3*u_k3[i][j]+4*u_n3[i][j]-u_old3[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv33[i][j]=(((-3*v_k3[i][j]+4*v_n3[i][j]-v_old3[i][j])/(2*dt))+rvs[i][j]);
		            
		            p_new[i][j]=Press[i][j]+(dtau*rcs[i][j]);
		            u_new[i][j]=u_k[i][j]+(dtau*RHSu33[i][j]);
		            v_new[i][j]=v_k[i][j]+(dtau*RHSv33[i][j]);
		        }
		    }
		    BC(N,u_new,v_new,p_new);			// BC after fourth step RK`

		   
		    // Updating old values
     
		    u_old1=u_n1;
		    u_n1=u_k1;
		    u_k1=u_new1;
		    u_old2=u_n2;
		    u_n2=u_k2;
		    u_k2=u_new2;
		    u_old3=u_n3;
		    u_n3=u_k3;
		    u_k3=u_new3;
		    u_old=u_n;
		    u_n=u_k;
		    u_k=u_new;
		    v_old1=v_n1;
		    v_n1=v_k1;
		    v_k1=v_new1;
		    v_old2=v_n2;
		    v_n2=v_k2;
		    v_k2=v_new2;
		    v_old3=v_n3;
		    v_n3=v_k3;
		    v_k3=v_new3;
		    v_old=v_n;
		    v_n=v_k;
		    v_k=v_new;		    
		    Press=p_new;
			 		   
		    
		    
		}  	
    }

   	/* Writing Data to file */
    contour(u_new,N,"xvel.vtk","U");
    contour(v_new,N,"yvel.vtk","V");
    contour(p_new,N,"Press.vtk","P");
	return 0;
}
