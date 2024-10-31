#include <iostream>


//using namespace std;
#include <stdio.h>


#include "interpolation/interpolation.hpp"

#include "Vector/vector_dist.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "Solvers/petsc_solver.hpp"
#include "interpolation/mp4_kernel.hpp"
#include "Solvers/petsc_solver_AMG_report.hpp"
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "VCluster/VCluster.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"

/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/

                                      //rho     U           P    chi      dU      drho   rho_old  U_old        

                                       //drho   dU      dP   
//typedef aggregate< double,VectorS<3,double>,double,double,VectorS<3, double>,double, double, VectorS<3,double>> property_type;
typedef aggregate< double,VectorS<3,double>,double,double,VectorS<3, double>,double, double, VectorS<3,double>> property_type;



typedef vector_dist_ws<3,double,property_type> particles_type;

typedef vector_dist_subset<3,double,property_type> particles_subset_type;

void *vectorGlobal,*vectorGlobal_bulk,*vectorGlobal_boundary,*vectorGlobal_inlet,*vectorGlobal_outlet;


constexpr unsigned int rho = 0;
constexpr unsigned int U = 1;
constexpr unsigned int P = 2;
constexpr unsigned int chi = 3;
constexpr unsigned int dU = 4;
constexpr unsigned int drho=5;
constexpr unsigned int rho_pennum=6;
constexpr unsigned int U_pennum=7;
//constexpr unsigned int tauA=8;


constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;
int res;
double t = 0,Dt, endTime,Machnum;

double spacing[3];

/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/

double rho_o    =993.3;// kg/m3
double P_o      =100000;//397320;
double P_ref    =100000;
double U_o      =0.001;//0.3;//0.0001;//0.3       MA always at 0.01;
double U_s      =0.0;
double U_out    = 0.0001;
double mu       = 0.0073;  //0.0021;//0.000691;//1.7;
double cs       =0.1;//120;//0.0066;//30.0;
double Tini     =1;
double Gamma    =7;


double  total_rho=0, rhoMax=0, pMax=0;
double Umax=0, CFLmax=0;  
/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************

struct state_type_4d_ofp
{   state_type_4d_ofp()
    {
        }
        typedef size_t size_type;
        typedef int is_state_vector;
        aggregate<texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>> data;
        size_t size() const
        {return data.get<0>().size(); }
        void resize(size_t n)
        {
            data.get<0>().resize(n);
            data.get<1>().resize(n);
            data.get<2>().resize(n);
            data.get<3>().resize(n);
        }
};

namespace boost {
    namespace numeric {
        namespace odeint {
            template<>
            struct is_resizeable<state_type_4d_ofp> {
                typedef boost::true_type type;
                static const bool value = type::value;
            };
            template<>
            struct vector_space_norm_inf<state_type_4d_ofp>
            {
                typedef double result_type;
            };
        }
    }
}

**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/
void init_porous(particles_type & pr, const Box<3,double> & domain, const int & res)
{
    double x0, y0, x1, y1,z0,z1,z2;
        x0 = domain.getLow(0)+1*spacing[0];//4
        y0 = domain.getLow(1)+1*spacing[1];//    5
	z0 = domain.getLow(2)+16*spacing[2];//4  8 16 

        x1 = domain.getHigh(0)-1*spacing[0];//14  7
        y1 =domain.getHigh(1)-1*spacing[1];//6    75z1 =domain.getHigh(2)-12*spacing[2];//14
	z1 =domain.getHigh(2)-1*spacing[2];// 6   12 24 48
	//z2 =domain.getHigh(2)-4*spacing[2];//8  25 45
	

    auto it = pr.getDomainIterator();
    while (it.isNext())
        {   	

		
        	auto key = it.get();
                double xp = pr.getPos(key)[x];
                double yp = pr.getPos(key)[y];
		double zp = pr.getPos(key)[z];
		


		pr.template getProp<rho>(key)  = rho_o;
               // pr.template getProp<P>(key)    =(rho_o*cs*cs/Gamma)*(pow((pr.template getProp<rho>(key)/rho_o),Gamma));
                //pr.template getProp<P>(key)    = P_o;//rho_o*(cs*cs);//(rho_o*cs*cs/Gamma)*(pow((pr.template getProp<rho>(key)/rho_o),Gamma));
		P_o=rho_o*(cs*cs);// (rho_o*cs*cs/Gamma)*(pow((pr.template getProp<rho>(key)/rho_o),Gamma));
		//P_o= P_o+(rho_o*cs*cs/Gamma)*(pow((pr.template getProp<rho>(key)/rho_o),Gamma)-1);
		
		pr.template getProp<P>(key)= P_o;

		pr.template getProp<U>(key)[x] = U_s;
                pr.template getProp<U>(key)[y] = U_s;
                pr.template getProp<U>(key)[z] = U_s;
		pr.template getProp<U_pennum>(key)[x]=U_s;
		pr.template getProp<U_pennum>(key)[y]=U_s;
		pr.template getProp<U_pennum>(key)[z]=U_s;
		pr.template getProp<rho_pennum>(key)=rho_o;

                        pr.setSubset(key,0);

/*		
//front
               if(xp >=x1 )//&& zp>=z0 && zp<= z1 && yp>=y0 && yp<=y1 )//else if(xp >4.64)
                {
                        pr.setSubset(key,3);
		//	pr.template getProp<U>(key)[z] = U_s;
			pr.template getProp<chi>(key) =1.0;

                }
//back
                if(xp <=x0)//(xp ==0 && zp>=z0 && zp<=z1 && yp>=y0 && yp<=y1 )
                {
                        pr.setSubset(key,0);
		//	pr.template getProp<U>(key)[z] = U_s;
			pr.template getProp<chi>(key) =1.0;
                }

//top
                 if(yp >=y1 )//  && zp>=z0 && zp<=z1 && xp>=x0 && xp<=x1 )
                {
                        pr.setSubset(key,0);
		//	pr.template getProp<U>(key)[z] = U_s;
			pr.template getProp<chi>(key) =1.0;
                }
//bottom
                 if(yp<=y0)//yp ==0 && zp>=z0 && zp<= z1 && xp>=x0 && xp<=x1 )
                {
                        pr.setSubset(key,0);
		//	pr.template getProp<U>(key)[z] = U_s;
			pr.template getProp<chi>(key) =1.0;
                }
*/	

//SOURCE
                 if(zp <=z0 && 	pr.template getProp<chi>(key) ==0.0 )//&& xp>x0 && xp<x1 && yp>y0 && yp<y1)
                {
                        pr.setSubset(key,1);
			//pr.template getProp<U>(key)[z] =U_o;
			pr.template getProp<P>(key) =P_o+P_o*1.9;
			pr.template getProp<chi>(key) =0.0;
                }
//outlet
                 if (zp>=z1 && 	pr.template getProp<chi>(key) ==0.0 )// && xp>x0 && xp<x1  && yp>y0 && yp<y1  )
                {
                        pr.setSubset(key,2);
			pr.template getProp<chi>(key) =0.0;
//			pr.template getProp<U>(key)[z] =0.001*U_o;
	       	}



//BULK
                //else
                //{
                 //       pr.setSubset(key,0);
                //}

       		++it;
    	}  
    
 
}
/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/
void flip_chi(particles_type & fluid_pr) 
{
double sphere_radius=0.5;

double radius1 = -0.5*sphere_radius;
double radius2 =  0.5*sphere_radius;
double centerx = 3.0;
double centery = 3.0;
double centerz = 1.5;
double rho;



        auto it = fluid_pr.getDomainIterator();
        while (it.isNext())
        {
                auto key = it.get();
	 if( fluid_pr.template getProp<chi>(key) >=0.52)//0.82 )//&& fluid_pr.template getProp<chi>(key)<0.9 )
		{
            		fluid_pr.template getProp<chi>(key) =1.0;
		}
       	     else
             	{
            		fluid_pr.template getProp<chi>(key) =0.0;
        
		}

             ++it;
        }

}
/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/
void diagnostics(particles_type & pr)
{
    Vcluster<> & ranks = create_vcluster();
    double Umag =0, lCFL=0; //total_rho=0, rhoMax=0, pMax=0;
//  double Umax=0;
CFLmax=0;  
Umax=0;
    auto it = pr.getDomainIterator();
        while (it.isNext())
    {
        auto key = it.get();
        Umag = sqrt(pow(pr.template getProp<U>(key)[x],2) +pow (pr.template getProp<U>(key)[y],2) + pow(pr.template getProp<U>(key)[z],2));
	lCFL= pr.template getProp<U>(key)[x]*Dt/spacing[0] + pr.template getProp<U>(key)[y]*Dt/spacing[1] + pr.template getProp<U>(key)[z]*Dt/spacing[3]; 
        Umax =std::max(Umag,Umax);
	CFLmax =std::max(lCFL,CFLmax);
        total_rho = total_rho + pr.template getProp<rho>(key);
        rhoMax = std::max (pr.template getProp<rho>(key) , rhoMax); 
        ++it;
    }
   // Umax =sqrt(Umax);
    ranks.max(Umax);
    ranks.max(CFLmax);
    ranks.sum(total_rho);
    ranks.max(rhoMax);
    ranks.execute();
    //Umax =sqrt(Umax);
	Machnum =Umax/cs;
//    cs = Umax / 0.01;//Machnum;
    
}

/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/
template<typename DX,typename DY,typename DZ,typename DXX,typename DXY,typename DXZ,typename DYY,typename DYZ,typename DZZ>
//template<typename DX,typename DY,typename DZ>
struct RHSFunctor
{
    	DX &Dx;
        DY &Dy;
        DZ &Dz;

        DXX &Dxx;
        DXY &Dxy;
        DXZ &Dxz;
        DYY &Dyy;
        DYZ &Dyz;
        DZZ &Dzz;

        //Constructor
        RHSFunctor(DX &Dx,DY &Dy,DZ &Dz,DXX &Dxx,DXY &Dxy,DXZ &Dxz,DYY &Dyy,DYZ &Dyz,DZZ &Dzz):Dx(Dx),Dy(Dy),Dz(Dz),Dxx(Dxx),Dxy(Dxy),Dxz(Dxz),Dyy(Dyy),Dyz(Dyz),Dzz(Dzz)
        //RHSFunctor(DX &Dx,DY &Dy,DZ &Dz):Dx(Dx),Dy(Dy),Dz(Dz)
        {}
//This is a computation of a stage
void operator()( const state_type_4d_ofp &X , state_type_4d_ofp &dxdt , const double t ) const
        {      
                timer tm_solve;
                tm_solve.start();
                particles_type &pr= *(particles_type *) vectorGlobal;
                particles_subset_type &pr_bulk= *(particles_subset_type *) vectorGlobal_bulk;
                particles_subset_type &pr_inlet= *(particles_subset_type *) vectorGlobal_inlet;
                particles_subset_type &pr_outlet= *(particles_subset_type *) vectorGlobal_outlet;
        particles_subset_type &pr_boundary= *(particles_subset_type *) vectorGlobal_boundary;

//**********************************************************************                 
                auto P_fld = getV<P>(pr);
                auto U_fld = getV<U>(pr);
                auto chi_fld = getV<chi>(pr);
                auto rho_fld = getV<rho>(pr);
                auto dU_fld = getV<dU>(pr);
//                auto drho_fld = getV<drho>(pr);
                auto Upenn = getV<U_pennum>(pr);
                //auto tau = getV<tauA>(pr);
                auto rhopenn = getV<rho_pennum>(pr);

		auto U_bulk = getV<U>(pr_bulk);
                auto dU_bulk = getV<dU>(pr_bulk);
                auto rho_bulk = getV<rho>(pr_bulk);
                auto drho_bulk = getV<drho>(pr_bulk);
                auto P_bulk = getV<P>(pr_bulk);

                auto dP_bulk = getV<drho>(pr_bulk);
		auto dP_fld = getV<drho>(pr);
                auto dP_outlet= getV<drho>(pr_outlet);
                auto dP_inlet = getV<drho>(pr_inlet);
/*
		auto drho_outlet=getV<drho>(pr_outlet);
		auto drho_inlet= getV<drho>(pr_inlet);
*/
		auto dU_outlet = getV<dU>(pr_outlet);
		auto dU_inlet = getV<dU>(pr_inlet);
                pr.ghost_get<P,U,rho>(SKIP_LABELLING);
//**********************************************************************
      
		texp_v<double> dP_x=Dx(P_fld),dP_y=Dy(P_fld),dP_z=Dz(P_fld),
                	DivM = Dx(rho_fld*U_fld[x])+Dy(rho_fld*U_fld[y])+Dz(rho_fld*U_fld[z]),
			DivM1 =rho_fld*(Dx(U_fld[x])+Dy(U_fld[y])+Dz(U_fld[z])),
			DivM2 =U_fld[x]*Dx(rho_fld)+U_fld[y]*Dy(rho_fld)+U_fld[z]*Dz(rho_fld),
			dzUz=Dz(U_fld[z]),
		        dxxUx=Dxx(U_fld[x]),dxxUy=Dxx(U_fld[y]),dxxUz=Dxx(U_fld[z]),
                	dxyUx=Dxy(U_fld[x]),dxyUy=Dxy(U_fld[y]),dxyUz=Dxy(U_fld[z]),
                	dxzUx=Dxz(U_fld[x]),dxzUy=Dxz(U_fld[y]),dxzUz=Dxz(U_fld[z]),
                	dyyUx=Dyy(U_fld[x]),dyyUy=Dyy(U_fld[y]),dyyUz=Dyy(U_fld[z]),
                	dyzUx=Dyz(U_fld[x]),dyzUy=Dyz(U_fld[y]),dyzUz=Dyz(U_fld[z]),
                	dzzUx=Dzz(U_fld[x]),dzzUy=Dzz(U_fld[y]),dzzUz=Dzz(U_fld[z]),
                	divU= Dx(U_fld[x])+ Dy(U_fld[y])+ Dz(U_fld[z]),
		        divUx = Dx(U_fld[x])*U_fld[x]+ Dy(U_fld[x])*U_fld[y]+ Dz(U_fld[x])*U_fld[z],
		        divUy = Dx(U_fld[y])*U_fld[x]+ Dy(U_fld[y])*U_fld[y]+ Dz(U_fld[y])*U_fld[z],
		        divUz = Dx(U_fld[z])*U_fld[x]+ Dy(U_fld[z])*U_fld[y]+ Dz(U_fld[z])*U_fld[z],
			DivP  = U_fld[x]*Dx(P_fld)+ U_fld[y]*Dy(P_fld)+ U_fld[z]*Dz(P_fld),
			dP_xx=Dxx(P_fld),dP_yy=Dyy(P_fld),dP_zz=Dzz(P_fld);

  		auto tau_xx = mu*(2*dxxUx -2/3*(dxxUx + dxyUy + dxzUz));
                auto tau_xy = mu*(dyyUx+ dxyUy);
                auto tau_xz = mu*(dzzUx+dxzUz);

                auto tau_yx = mu*(dxyUx+ dxxUy);
                auto tau_yy = mu*(2*dyyUy -2/3*(dxyUx + dyyUy + dyzUz));
                auto tau_yz = mu*(dzzUy+dyzUz);

                auto tau_zx = mu*(dxzUx+dxxUz);
                auto tau_zy = mu*(dyzUy+dyyUz);
                auto tau_zz = mu*(2*dzzUz -2/3*(dxzUx + dyzUy + dzzUz));
//**********************************************************************
		 

		dU_bulk[x]= -divUx + ( -dP_x + tau_xx + tau_xy + tau_xz)/rho_fld;//-chi_fld/Dt*(-Upenn[x]+U_fld[x]);  
               	dU_bulk[y]= -divUy + ( -dP_y + tau_yx + tau_yy + tau_yz)/rho_fld;//-chi_fld/Dt*(-Upenn[y]+U_fld[y]);
		dU_bulk[z]= -divUz + ( -dP_z + tau_zx + tau_zy + tau_zz)/rho_fld;//-chi_fld/Dt*(-Upenn[z]+U_fld[z]);  

                 dP_bulk =  -DivP-rho_fld*cs*cs*divU + 0.5*spacing[0]*cs/8* (dP_xx+dP_yy+dP_zz);
                //dP_inlet = -rho_fld*cs*cs*divU + 0.5*spacing[0]*cs/8* (dP_xx+dP_yy);//+dP_zz);
                //dP_outlet= -rho_fld*cs*cs*divU + 0.5*spacing[0]*cs/8* (dP_xx+dP_yy);//+dP_zz);
		
//**********************************************************************
//********************************************************************** 
/*	texp_v<double>
			DivM1_o =rho_fld*(Dz(U_fld[z])),
			DivM2_o =U_fld[z]*Dz(rho_fld),
 
		        dxxUx_o =Dxx(U_fld[x]),dxxUy_o =Dxx(U_fld[y]),dxxUz_o =Dxx(U_fld[z]),
                	dxyUx_o =Dxy(U_fld[x]),dxyUy_o =Dxy(U_fld[y]),dxyUz_o =Dxy(U_fld[z]),
                	dxzUx_o =Dxz(U_fld[x]),dxzUy_o =Dxz(U_fld[y]),dxzUz_o =Dxz(U_fld[z]),
                	dyyUx_o =Dyy(U_fld[x]),dyyUy_o =Dyy(U_fld[y]),dyyUz_o =Dyy(U_fld[z]),
                	dyzUx_o =Dyz(U_fld[x]),dyzUy_o =Dyz(U_fld[y]),dyzUz_o =Dyz(U_fld[z]),
                	dzzUx_o =Dzz(U_fld[x]),dzzUy_o =Dzz(U_fld[y]),dzzUz_o =Dzz(U_fld[z]),
		        divUz_o = Dx(U_fld[z])*U_fld[x]+ Dy(U_fld[z])*U_fld[y]+ Dz(U_fld[z])*U_fld[z];
		
		
		auto tau_zx_o = mu*(dxzUx_o+dxxUz_o);
//		auto tau_zx_o = mu*(dxxUz_o);
	        auto tau_zy_o = mu*(dyzUy_o+dyyUz_o);
//		auto tau_zy_o = mu*(dyyUz);
                auto tau_zz_o = mu*(2*dzzUz_o -2/3*(dxzUx_o + dyzUy_o + dzzUz_o));
*/
//********************************************************************** 
//		dU_outlet[z] = 	-(U_fld[z] *Dz(U_fld[z])+U_fld[x] *Dz(U_fld[z])+U_fld[y] *Dz(U_fld[z]));
		//dU_outlet[x] =  0;//-U_fld[z] *Dz(U_fld[z]);
		//dU_outlet[y] =  0;-U_fld[y] *Dz(U_fld[z]);
		//dU_outlet[z] = 	-divUz_o + ( -dP_z + mu*(Dx(tau[z][x]) + Dy(tau[z][y]) + Dz(tau[z][z])))/rho_o;
		//dU_outlet[x] = -divUx+( -dP_x + tau_xy + tau_xz )/rho_fld-chi_fld/dt*(-Upenn[x]+U_fld[x]);  
		//dU_outlet[y] = -divUy+( -dP_y + tau_yx + tau_yz )/rho_fld-chi_fld/dt*(-Upenn[y]+U_fld[y]);  
//		dU_outlet[z] = -(U_fld[z] * Dz(U_fld[z]));//
		
		//dU_outlet[z] = -divUz+(   tau_zx + tau_zy)/rho_fld -chi_fld/dt*(-Upenn[z]+U_fld[z]);		
		//dU_inlet[z]  = -divUz+(   tau_zx + tau_zy)/rho_fld -chi_fld/dt*(-Upenn[z]+U_fld[z]);// -U_o*dzUz;	

                //drho_outlet = - (DivM);//+DivM2_o);
	
//********************************************************************** 
           //RHS for ODEint.
                dxdt.data.get<0>()=dU_fld[x];
                dxdt.data.get<1>()=dU_fld[y];
                dxdt.data.get<2>()=dU_fld[z];
                dxdt.data.get<3>()=dP_fld;
                tm_solve.stop();
                //std::cout << "Time:to BUILD the SYSTEM  " << tm_solve.getwct()  << std::endl;

    }

};

//    dU_bulk[x]= -divUx-U_fld[x]*divU+ (-dP_x + tau_xx + tau_xy + tau_xz)/rho_fld;// - chi_fld/0.00001*(-Upenn[x]+U_fld[x]))/rho_fld;
//	dU_bulk[y]=-divUy-U_fld[y]*divU+ (-dP_y + tau_yx + tau_yy + tau_yz)/rho_fld;// - chi_fld/0.00001*(-Upenn[y]+U_fld[y]))/rho_fld;
//	dU_bulk[z]=-divUz-U_fld[z]*divU+ (-dP_z + tau_zx + tau_zy + tau_zz)/rho_fld;// - chi_fld/0.00001*(-Upenn[z]+U_fld[z]))/rho_fld;
//	drho_bulk = - (DivM+ ( 1/0.001 -1) *chi_fld*DivM);//   -1*(( 1.0 +  (penParam *chi_fld))*DivM);
//               drho_bulk=-rho_fld*DivM-(rho_o)*chi_fld*DivM;

/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/
template<typename DX,typename DY,typename DZ,typename DXX,typename DXY,typename DXZ,typename DYY,typename DYZ,typename DZZ>
//template<typename DX,typename DY,typename DZ>
struct ObserverFunctor
{
    	DX &Dx;
        DY &Dy;
        DZ &Dz;

        DXX &Dxx;
        DXY &Dxy;
        DXZ &Dxz;
        DYY &Dyy;
        DYZ &Dyz;
        DZZ &Dzz;

    //Constructor
        int ctr;
        double t_old;
    //Constructor
        ObserverFunctor(DX &Dx,DY &Dy,DZ &Dz,DXX &Dxx,DXY &Dxy,DXZ &Dxz,DYY &Dyy,DYZ &Dyz,DZZ &Dzz):Dx(Dx),Dy(Dy),Dz(Dz),Dxx(Dxx),Dxy(Dxy),Dxz(Dxz),Dyy(Dyy),Dyz(Dyz),Dzz(Dzz)
       // ObserverFunctor(DX &Dx,DY &Dy,DZ &Dz):Dx(Dx),Dy(Dy),Dz(Dz)
        { 
            ctr = 0;
            t_old = -Dt;
        }
    void operator() (state_type_4d_ofp &X, double t)
        { 

        timer tm_solve;
        tm_solve.start();

        particles_type &pr= *(particles_type *) vectorGlobal;
        particles_subset_type &pr_bulk= *(particles_subset_type *) vectorGlobal_bulk;
        particles_subset_type &pr_inlet= *(particles_subset_type *) vectorGlobal_inlet;
        particles_subset_type &pr_outlet= *(particles_subset_type *) vectorGlobal_outlet;
        particles_subset_type &pr_boundary= *(particles_subset_type *) vectorGlobal_boundary;


//	penalise_interpolate(pr);
	 
//Aliasing the position and properties.
        auto Pos = getV<PROP_POS>(pr);
        auto U_fld = getV<U>(pr);
        auto rho_fld = getV<rho>(pr);
        auto chi_fld = getV<chi>(pr);
       	auto Upenn = getV<U_pennum>(pr);
	auto rhopenn = getV<rho_pennum>(pr); 
        auto Pos_bulk = getV<PROP_POS>(pr_bulk);
        auto U_bulk = getV<U>(pr_bulk);
        auto rho_bulk = getV<rho>(pr_bulk);
        auto P_fld = getV<P>(pr);
        auto P_bulk = getV<P>(pr_bulk);

        auto U_inlet = getV<U>(pr_inlet);
        auto P_inlet= getV<P>(pr_inlet);
        auto rho_inlet = getV<rho>(pr_inlet);
 	auto U_outlet = getV<U>(pr_outlet);
        auto P_outlet= getV<P>(pr_outlet);
        auto rho_outlet = getV<rho>(pr_outlet);

	

	U_fld[x] = X.data.get<0>()+( Upenn[x] - X.data.get<0>()) * chi_fld;
        U_fld[y] = X.data.get<1>()+( Upenn[y] - X.data.get<1>()) * chi_fld;
	U_fld[z] = X.data.get<2>()+( Upenn[z] - X.data.get<2>()) * chi_fld;
        //rho_fld  = X.data.get<3>()+( rhopenn - X.data.get<3>()) * chi_fld;
	
       P_fld = X.data.get<3>();//+( P_o - X.data.get<3>()) * chi_fld; 
	
//********************************************************************************
//We would like to move the particles after t=0. (Remember, odeint calls the observer before taking the step.)
                if (t != 0) {
//copy back the state after the time step into data structure. This is required as we are going to move the particles and the distributed state can be resized correctly (by copying back after map). Also this expression preserves the boundary condition on concentration.
 
                   U_fld[x] = X.data.get<0>()+( Upenn[x] - X.data.get<0>()) * chi_fld;
                   U_fld[y] = X.data.get<1>()+( Upenn[y] - X.data.get<1>()) * chi_fld;
                   U_fld[z] = X.data.get<2>()+( Upenn[z] - X.data.get<2>()) * chi_fld;
 	          // rho_fld  = X.data.get<3>()+( rhopenn - X.data.get<3>()) * chi_fld;
       		P_fld = X.data.get<3>();//+( P_o - X.data.get<3>()) * chi_fld;
	       P_outlet=P_o;	
	//	rho_fld = P_fld/(cs*cs);
	//	rho_fld= rho_fld +(rho_o - rho_fld)*chi_fld;
		//rho_inlet=rho_o;
//******************************************************************************** 

        diagnostics(pr);
//******************************************************************************** 

//UPDATE the PRESSURE 


//******************************************************************************** 

//Euler step for moving particles
//
//                    Pos_bulk = Pos + dt * U_fld;
//Map and ghost_et is required after moving particles.
//                    pr.map();
//pr.ghost_get<0>();

//Updating the subset and operators based on new positions
/*                   pr_bulk.update();
                    pr_inlet.update();
                    Dx.update(pr);
                    Dy.update(pr);
                    Dz.update(pr);
                    Dxy.update(pr);
                    Dxz.update(pr);
                    Dyz.update(pr);
                    Dxx.update(pr);
                    Dyy.update(pr);
                    Dzz.update(pr);
*/
//Sinince we did Map, we assign the Odeint state again.
                    X.data.get<0>() = U_fld[x];
                    X.data.get<1>() = U_fld[y];
                    X.data.get<2>() = U_fld[z];
                    X.data.get<3>() = P_fld;//rho_fld;
                    ctr++;

//counting the step number    
            }
            tm_solve.stop();
            Vcluster<> & ranks  =create_vcluster();                   
            if(ctr % 100  == 0)
            {   
		if(ranks.getProcessUnitID()==0) 
		{
			std::cout <<"dt="<<Dt<<"cs= "<<cs<<" Umax= "<<Umax<<" CFL ="<< CFLmax<<" Ma = "<<Machnum<< "\n";        
              		std::cout << "Time:to calc observerand update " << tm_solve.getwct()<<" @step"  <<ctr<< std::endl;

		}
		if(ctr%5000 ==0)
		{
			pr.deleteGhost();
               		pr.write_frame("PDE_sol", ctr);
			pr.save("particles_restart.hdf5");
		}

            }
                pr.ghost_get<0>();




    }
};


/*
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*/
int main(int argc, char* argv[])
{
    openfpm_init(&argc,&argv);
    	Vcluster<> & ranks = create_vcluster();
        timer tt2;
        tt2.start();
	res = int(std::atof(argv[1])) ;
        
        endTime=std::atof(argv[2]);
        Dt=std::atof(argv[3]);
	double timeTOL=std::atof(argv[4]);

	std::string hdf5_file   = "exparticles_",buf(argv[1]);


    	buf.append(".hdf5");
    	hdf5_file=hdf5_file+buf;


       	long int szu[3] ={res,res, 265};//91
	size_t sz[] = {(size_t)szu[0],(size_t)szu[1],(size_t)szu[2]};

	Box<3,double> domain({0.0,0.0,0.00},{0.002572862599889584 , 0.002572862599889584 , 0.0077055000932919815});


	size_t bc[3]={PERIODIC,PERIODIC,NON_PERIODIC};

	spacing[0]=domain.getHigh(0)/ (res );
    	spacing[1]=domain.getHigh(1)/ (res);
        spacing[2]=domain.getHigh(2)/ (265-1);


	double rCut =(3.1 * spacing[2]);
    	Ghost<3,double> ghost(rCut);

//**********************************************************************
    	particles_type particles(0,domain,bc,ghost);

    	particles.setPropNames({"00-Rho","01-U","02-P","03-Chi","04-dU","05-dRho","06-rho_pen","07-U_pen"});

//**********************************************************************
        if (ranks.getProcessUnitID() == 0)
        std::cout << "Loading the file \n"; 
        particles.load(hdf5_file);
        particles.map();
        particles.ghost_get<0>();  

//**********************************************************************

    if (ranks.getProcessUnitID() == 0)
    std::cout << "Initialising the IC and BC \n";

    flip_chi(particles);
    init_porous(particles,domain,res);

//**********************************************************************

//Creating Subset with id 0 (0s are Initialized inside init_porous)
    	particles_subset_type particles_bulk(particles,0);
   	particles_subset_type particles_inlet(particles,1);
        particles_subset_type particles_outlet(particles,2);
        particles_subset_type particles_boundary(particles,3);

//**********************************************************************
//We create aliases for referring to the BULK FLOW  properties.
        auto P_bulk   = getV<P>(particles_bulk);
        auto U_bulk   = getV<U>(particles_bulk);
        auto chi_bulk = getV<chi>(particles_bulk);
        auto rho_bulk = getV<rho>(particles_bulk);
    
        auto P_fld   = getV<P>(particles);
        auto U_fld   = getV<U>(particles);
        auto chi_fld = getV<chi>(particles);
        auto rho_fld = getV<rho>(particles);
        auto Upenn = getV<U_pennum>(particles);
        auto rhopenn = getV<rho_pennum>(particles);
        auto dU_fld = getV<dU>(particles);
        auto drho_fld = getV<drho>(particles);
        Upenn=U_s;
        rhopenn=rho_o;
        dU_fld=0;
        drho_fld=0;
//**********************************************************************
//        particles.deleteGhost();

//**********************************************************************
        //particles.write("Init",BINARY);
 //       particles.ghost_get();

//**********************************************************************
//Creating the operators
        timer tm_solve;
        tm_solve.start();
      if (ranks.getProcessUnitID() == 0)
        std::cout<<"Getting the Derivtives\n";
        Derivative_x  Dx (particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "1 \n"
        Derivative_y  Dy (particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "2 \n";
        Derivative_z  Dz (particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "3 \n";

        Derivative_xx Dxx(particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "4 \n";
        Derivative_yy Dyy(particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "5 \n";
        Derivative_zz Dzz(particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "6 \n"
        Derivative_xy Dxy(particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "7 \n"
        Derivative_xz Dxz(particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "8 \n";
        Derivative_yz Dyz(particles, 3, rCut,3.1,support_options::RADIUS);
       // std::cout << "9 \n";

//	Gradient  Grad (particles, 2, rCut);//,sampling_factor2, support_options::RADIUS);
  //      Laplacian Lap(particles, 2, rCut);//,sampling_factor, support_options::RADIUS);
//	Divergence   Div (particles, 3, rCut);//,sampling_factor, support_options::RADIUS);
//	Advection    Adv (particles, 2, rCut);//,sampling_factor, support_options::RADIUS);

        tm_solve.stop();
      if (ranks.getProcessUnitID() == 0)
        std::cout << "Time:to calc Derivatives " << tm_solve.getwct()  << std::endl;

//**********************************************************************
//Casting the global pointers for Odeint
        vectorGlobal=(void *) &particles;
        vectorGlobal_bulk=(void *) &particles_bulk;
        vectorGlobal_boundary=(void *) &particles_boundary;
        vectorGlobal_inlet=(void *) &particles_inlet;
        vectorGlobal_outlet=(void *) &particles_outlet;


//**********************************************************************
//Creating Odeint RK4, See ODEINT for more steppers.
//boost::numeric::odeint::adams_bashforth_moulton<2, state_type_4d_ofp,double,state_type_4d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp > abm;

//boost::numeric::odeint::euler< state_type_4d_ofp,double,state_type_4d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp > abm;

boost::numeric::odeint::adaptive_adams_bashforth_moulton<2, state_type_4d_ofp,double,state_type_4d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp > abmA;

//The method Odeint_rk4 from Odeint, requires system (a function which computes RHS of the PDE), an instance of the Compute RHS functor. We create the System with the correct types and parameteres for the operators as declared before.
        RHSFunctor<Derivative_x,Derivative_y,Derivative_z,Derivative_xx,Derivative_xy,Derivative_xz,Derivative_yy,Derivative_yz,Derivative_zz> System(Dx,Dy,Dz,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
        //RHSFunctor<Derivative_x,Derivative_y,Derivative_z> System(Dx,Dy,Dz);

//Since we are using Odeint to control the time steps, we also create a observer instance. Which also updates the position via an euler step for moving thr particles.
        ObserverFunctor<Derivative_x,Derivative_y,Derivative_z,Derivative_xx,Derivative_xy,Derivative_xz,Derivative_yy,Derivative_yz,Derivative_zz> ObserveAndUpdate(Dx,Dy,Dz,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
       // ObserverFunctor<Derivative_x,Derivative_y,Derivative_z> ObserveAndUpdate(Dx,Dy,Dz);

//Furhter, odeint needs data in a 4 dimesional state type "state_type_4d_ofp", we create one and fill in the initial condition.
        state_type_4d_ofp X;
    
//Since we created a 4d state_type we initialize the two fields in the object data using the method get.
        X.data.get<0>() = U_fld[0];
        X.data.get<1>() = U_fld[1];
        X.data.get<2>() = U_fld[2];
        X.data.get<3>() = P_fld;//rho_fld;

//**********************************************************************
    tm_solve.start();
        if (ranks.getProcessUnitID() == 0)
            std::cout << "HERE WE GO For Time integration!!!! \n";
    

        std::vector<double> inter_times; // vector to store intermediate time steps taken by odeint.

//	size_t steps = boost::numeric::odeint::integrate_const(abm,System,X,0.0,endTime,Dt,ObserveAndUpdate);
	size_t steps =boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(timeTOL,timeTOL,abmA),System,X,0.0,endTime,Dt,ObserveAndUpdate);
        tm_solve.stop();
        std::cout << "Time: " << tm_solve.getwct()  << std::endl;
 
        std::cout << "Time steps: " << steps << std::endl;

        U_fld[x]=X.data.get<0>();
        U_fld[y]=X.data.get<1>();
        U_fld[z]=X.data.get<2>();
        P_fld=X.data.get<3>();

//**********************************************************************
        Dx.deallocate(particles);
        Dy.deallocate(particles);
        Dz.deallocate(particles);

	Dxy.deallocate(particles);
        Dxz.deallocate(particles);
        Dyz.deallocate(particles);
        Dxx.deallocate(particles);
        Dyy.deallocate(particles);
        Dzz.deallocate(particles);


//**********************************************************************

//particles.write("KILL", VTK_WRITER | FORMAT_BINARY );
//particles.getDecomposition().write("dec");

    openfpm_finalize();
exit(0);
}

