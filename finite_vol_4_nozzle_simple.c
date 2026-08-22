#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <err.h>
#include <time.h>

// ---------------------------------------------------------------------
// Milestone 4, revisited: nozzle / contraction geometry, with a real
// pressure-correction (projection / SIMPLE-family) step.
//
// This file starts from the validated milestone-2/3 code
// (finite_vol_2_boundary_conditions.c + finite_vol_3_carreau_small.c:
// ghost-cell BCs, Carreau shear-thinning viscosity via the full stress
// divergence) and adds two things:
//
//   1. GEOMETRY: a sudden 2:1 contraction, implemented as a solid-cell
//      mask inside the same rectangular grid (no new data structures --
//      solid cells are just regular grid cells that are excluded from
//      the momentum/pressure update and held at vx=vy=0, no-slip).
//
//   2. MASS CONSERVATION: the old update_pressure()/b_src formula
//      (-dxvx^2 - dyvy^2 + 2*dxvy*dyvx) is REPLACED. That formula is the
//      classical pressure-Poisson source derived by taking the
//      divergence of the momentum equation *assuming* div(v)=0 already
//      holds -- it never actually measures whether the velocity field
//      it's given satisfies continuity, so at a geometry change nothing
//      forces the flow rate to match between the wide section and the
//      throat.
//
//      The fix is a predictor/corrector (projection-method) step, which
//      is the direct explicit-timestepping analogue of SIMPLE:
//        - predict_velocity(): advance vx,vy with advection + viscous
//          stress ONLY (no pressure) -> vxstar, vystar. This field does
//          NOT satisfy continuity in general.
//        - the actual, local mass imbalance of (vxstar,vystar) becomes
//          the source term of a standard pressure-Poisson equation:
//              lap(p) = div(v*)/dt
//        - correct_velocity(): vx = vxstar - dt*dp/dx (same for vy).
//          Because p was solved specifically to cancel the imbalance in
//          vxstar,vystar, the corrected field satisfies discrete
//          continuity in every cell -- including at the throat, where it
//          forces vx to double to match the halved cross-section.
//
//      See simple_method_nozzle_tutorial.md for the full derivation.
// ---------------------------------------------------------------------

#define Nx 100
#define Ny 20
#define NXG (Nx+2)
#define NYG (Ny+2)
#define Ntot (NXG*NYG)

// Contraction geometry: cells x in [1,CONTRACTION_X] are full height
// (wide section); cells x in [CONTRACTION_X+1,Nx] are only open for
// y in [THROAT_Y_LO,THROAT_Y_HI] (throat), everything else in that
// x-range is a solid wall cell. THROAT height = Ny/2 -> 2:1 area ratio.
#define CONTRACTION_X 50
#define THROAT_Y_LO (Ny/4+1)     // 6
#define THROAT_Y_HI (3*Ny/4)     // 15

double dt=0.006;
double U_in=0.1;     // inlet velocity (Dirichlet)
double P_out=0.0;    // outlet reference pressure (Dirichlet)

// Carreau model parameters (same convention/placeholders as milestone 5/6)
double CARREAU_MU0   = 1.0;
double CARREAU_MUINF = 0.1;
double CARREAU_GAMMA = 300.0;
double CARREAU_Q     = 0.5;

double p[Ntot], vx[Ntot], vy[Ntot];
double vxstar[Ntot], vystar[Ntot];
int E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];
int solid[Ntot];

double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];

double CDxVx[Ntot], CDxVy[Ntot], CDyVx[Ntot], CDyVy[Ntot];
double Dmag[Ntot];
double eta[Ntot];

double b_src[Ntot];
double SOR_OMEGA=1.7;

int site2index(int x,int y)
{
    return x + y*NXG;
}

void geometry()
{
    for (int i=0;i<Ntot;i++) dx[i]=dy[i]=1.0;
}

void init_control_points()
{
    for (int x=0;x<NXG;x++){
        for (int y=0;y<NYG;y++){
            int i=site2index(x,y);
            E[i]=site2index(x+1<NXG?x+1:x, y);
            W[i]=site2index(x-1>=0?x-1:x, y);
            N[i]=site2index(x, y+1<NYG?y+1:y);
            S[i]=site2index(x, y-1>=0?y-1:y);
        }
    }
}

void init_solid_mask()
{
    for (int i=0;i<Ntot;i++) solid[i]=0;
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (x>CONTRACTION_X && (y<THROAT_Y_LO || y>THROAT_Y_HI))
                solid[i]=1;
        }
    }
}

void init_fields()
{
    for (int i=0;i<Ntot;i++){
        p[i]=0.0; vx[i]=0.0; vy[i]=0.0;
        vxstar[i]=0.0; vystar[i]=0.0;
        eta[i]=CARREAU_MU0;
    }
}

// Extend a (ux,uy) velocity pair onto the ghost ring + solid cells using
// the same physical BCs each time: Dirichlet inlet, zero-gradient
// outlet, no-slip top/bottom walls and no-slip solid (contraction)
// walls. Used for both the real field (vx,vy) and the predicted field
// (vxstar,vystar) -- the BCs themselves don't involve pressure, so they
// apply identically to either.
void apply_velocity_bc(double *ux,double *uy)
{
    for (int y=1;y<=Ny;y++){
        int g=site2index(0,y);
        ux[g]=U_in; uy[g]=0.0;
    }
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y);
        int in1=site2index(Nx,y);
        ux[g]=ux[in1]; uy[g]=uy[in1];
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0);
        int in1=site2index(x,1);
        ux[g]=-ux[in1]; uy[g]=-uy[in1];
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1);
        int in1=site2index(x,Ny);
        ux[g]=-ux[in1]; uy[g]=-uy[in1];
    }
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]){ ux[i]=0.0; uy[i]=0.0; }
        }
    }
}

// Pressure BCs: zero-gradient at inlet, Dirichlet at outlet, zero-gradient
// at top/bottom walls and at every solid (contraction) wall cell.
void apply_pressure_bc()
{
    for (int y=1;y<=Ny;y++){
        int g=site2index(0,y);
        int in1=site2index(1,y);
        p[g]=p[in1];
    }
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y);
        p[g]=P_out;
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0);
        int in1=site2index(x,1);
        p[g]=p[in1];
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1);
        int in1=site2index(x,Ny);
        p[g]=p[in1];
    }
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (!solid[i]) continue;
            double psum=0.0; int pcnt=0;
            int nb[4]={E[i],W[i],N[i],S[i]};
            for (int k=0;k<4;k++){
                if (!solid[nb[k]]){ psum+=p[nb[k]]; pcnt++; }
            }
            if (pcnt>0) p[i]=psum/pcnt;
        }
    }
}

void apply_boundary_conditions()
{
    apply_velocity_bc(vx,vy);
    apply_pressure_bc();
}

void setup()
{
    geometry();
    init_control_points();
    init_solid_mask();
    init_fields();
    apply_boundary_conditions();
}

// Face-interpolated velocities, needed for the advection terms in the
// predictor step. Skips solid cells (never used as update targets).
void interpolate_velocity_faces()
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]) continue;
            vxe[i]=(vx[E[i]]+vx[i])/2;
            vxw[i]=(vx[W[i]]+vx[i])/2;
            vxn[i]=(vx[N[i]]+vx[i])/2;
            vxs[i]=(vx[S[i]]+vx[i])/2;
            vye[i]=(vy[E[i]]+vy[i])/2;
            vyw[i]=(vy[W[i]]+vy[i])/2;
            vyn[i]=(vy[N[i]]+vy[i])/2;
            vys[i]=(vy[S[i]]+vy[i])/2;
        }
    }
}

// Face-interpolated pressures, needed for the corrector step. Computed
// AFTER the pressure Poisson solve has converged.
void interpolate_pressure_faces()
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]) continue;
            pe[i]=(p[E[i]]+p[i])/2;
            pw[i]=(p[W[i]]+p[i])/2;
            pn[i]=(p[N[i]]+p[i])/2;
            ps[i]=(p[S[i]]+p[i])/2;
        }
    }
}

void compute_shear_and_viscosity()
{
    for (int i=0;i<Ntot;i++){
        CDxVx[i]=(vx[E[i]]-vx[W[i]])/(2.0*dx[i]);
        CDxVy[i]=(vy[E[i]]-vy[W[i]])/(2.0*dx[i]);
        CDyVx[i]=(vx[N[i]]-vx[S[i]])/(2.0*dy[i]);
        CDyVy[i]=(vy[N[i]]-vy[S[i]])/(2.0*dy[i]);
    }
    for (int i=0;i<Ntot;i++){
        double Dxx=CDxVx[i];
        double Dyy=CDyVy[i];
        double Dxy=0.5*(CDyVx[i]+CDxVy[i]);
        Dmag[i]=sqrt(Dxx*Dxx+Dyy*Dyy+2.0*Dxy*Dxy);
        double gd=CARREAU_GAMMA*Dmag[i];
        eta[i]=CARREAU_MUINF+(CARREAU_MU0-CARREAU_MUINF)*pow(1.0+gd*gd,(CARREAU_Q-1.0)/2.0);
    }
}

// Predictor: advection + Carreau viscous stress divergence, NO pressure.
// Result (vxstar,vystar) generally does NOT satisfy continuity -- that's
// expected and is exactly the field the pressure-correction step below
// is built to fix.
void predict_velocity()
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]){ vxstar[i]=0.0; vystar[i]=0.0; continue; }

            double dxvx=(vxe[i]-vxw[i])/dx[i];
            double dxvy=(vye[i]-vyw[i])/dx[i];
            double dyvx=(vxn[i]-vxs[i])/dy[i];
            double dyvy=(vyn[i]-vys[i])/dy[i];

            double dxvxe=2*(vx[E[i]]-vx[i])/(dx[i]+dx[E[i]]);
            double dxvye=2*(vy[E[i]]-vy[i])/(dx[i]+dx[E[i]]);
            double dxvxw=2*(-vx[W[i]]+vx[i])/(dx[i]+dx[W[i]]);
            double dxvyw=2*(-vy[W[i]]+vy[i])/(dx[i]+dx[W[i]]);
            double dyvxn=2*(vx[N[i]]-vx[i])/(dy[i]+dy[N[i]]);
            double dyvyn=2*(vy[N[i]]-vy[i])/(dy[i]+dy[N[i]]);
            double dyvxs=2*(-vx[S[i]]+vx[i])/(dy[i]+dy[S[i]]);
            double dyvys=2*(-vy[S[i]]+vy[i])/(dy[i]+dy[S[i]]);

            double eta_e=0.5*(eta[i]+eta[E[i]]);
            double eta_w=0.5*(eta[i]+eta[W[i]]);
            double eta_n=0.5*(eta[i]+eta[N[i]]);
            double eta_s=0.5*(eta[i]+eta[S[i]]);

            double dyvx_e=0.5*(CDyVx[i]+CDyVx[E[i]]);
            double dyvx_w=0.5*(CDyVx[i]+CDyVx[W[i]]);
            double dxvy_n=0.5*(CDxVy[i]+CDxVy[N[i]]);
            double dxvy_s=0.5*(CDxVy[i]+CDxVy[S[i]]);

            double tau_xx_e=2.0*eta_e*dxvxe;
            double tau_xx_w=2.0*eta_w*dxvxw;
            double tau_xy_n=eta_n*(dyvxn+dxvy_n);
            double tau_xy_s=eta_s*(dyvxs+dxvy_s);
            double visc_x=(tau_xx_e-tau_xx_w)/dx[i]+(tau_xy_n-tau_xy_s)/dy[i];

            double tau_yy_n=2.0*eta_n*dyvyn;
            double tau_yy_s=2.0*eta_s*dyvys;
            double tau_xy_e=eta_e*(dxvye+dyvx_e);
            double tau_xy_w=eta_w*(dxvyw+dyvx_w);
            double visc_y=(tau_xy_e-tau_xy_w)/dx[i]+(tau_yy_n-tau_yy_s)/dy[i];

            vxstar[i]=vx[i]+dt*(-vx[i]*dxvx-vy[i]*dyvx+visc_x);
            vystar[i]=vy[i]+dt*(-vx[i]*dxvy-vy[i]*dyvy+visc_y);
        }
    }
}

// Mass-imbalance source term for the pressure Poisson equation:
// b_src = div(v*) / dt, computed from the PREDICTED velocity field.
// This is the quantity that replaces the old ad hoc -dxvx^2-dyvy^2+... term.
void compute_b_src()
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]){ b_src[i]=0.0; continue; }
            double vxe_s=(vxstar[E[i]]+vxstar[i])/2;
            double vxw_s=(vxstar[W[i]]+vxstar[i])/2;
            double vyn_s=(vystar[N[i]]+vystar[i])/2;
            double vys_s=(vystar[S[i]]+vystar[i])/2;
            double div=(vxe_s-vxw_s)/dx[i]+(vyn_s-vys_s)/dy[i];
            b_src[i]=div/dt;
        }
    }
}

// Standard 5-point discrete Poisson solve for p, via SOR, with source
// b_src built from the actual predicted-velocity divergence above.
int solve_pressure_sor(double inner_tol,int max_inner,double *resid_out)
{
    compute_b_src();
    int k; double resid=1e300;
    for (k=0;k<max_inner;k++){
        resid=0.0;
        for (int x=1;x<=Nx;x++){
            for (int y=1;y<=Ny;y++){
                int i=site2index(x,y);
                if (solid[i]) continue;
                double idx2=1.0/(dx[i]*dx[i]);
                double idy2=1.0/(dy[i]*dy[i]);
                double p_gs=((p[E[i]]+p[W[i]])*idx2+(p[N[i]]+p[S[i]])*idy2-b_src[i])
                            /(2.0*idx2+2.0*idy2);
                double p_old=p[i];
                p[i]=(1.0-SOR_OMEGA)*p_old+SOR_OMEGA*p_gs;
                double d=fabs(p[i]-p_old);
                if (d>resid) resid=d;
            }
        }
        apply_pressure_bc();
        if (resid<inner_tol){ k++; break; }
    }
    *resid_out=resid;
    return k;
}

// Corrector: vx = vxstar - dt*dp/dx (and same for vy), using the
// pressure field that was just solved to cancel div(v*).
void correct_velocity()
{
    interpolate_pressure_faces();
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]){ vx[i]=0.0; vy[i]=0.0; continue; }
            vx[i]=vxstar[i]-dt*(pe[i]-pw[i])/dx[i];
            vy[i]=vystar[i]-dt*(pn[i]-ps[i])/dy[i];
        }
    }
}

double measure_velocity_l2()
{
    double s=0.0;
    for (int x=1;x<=Nx;x++)
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]) continue;
            s+=vx[i]*vx[i]+vy[i]*vy[i];
        }
    return sqrt(s);
}

// Volumetric flow rate through a vertical cross-section at column x
// (sum of vx*dy over the open/fluid cells there). This is the direct
// mass-conservation check: Q(x) should be constant along the duct.
double flow_rate_at(int x)
{
    double Q=0.0;
    for (int y=1;y<=Ny;y++){
        int i=site2index(x,y);
        if (!solid[i]) Q+=vx[i]*dy[i];
    }
    return Q;
}

double max_abs_divergence()
{
    double m=0.0;
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (solid[i]) continue;
            double vxe_=(vx[E[i]]+vx[i])/2, vxw_=(vx[W[i]]+vx[i])/2;
            double vyn_=(vy[N[i]]+vy[i])/2, vys_=(vy[S[i]]+vy[i])/2;
            double div=(vxe_-vxw_)/dx[i]+(vyn_-vys_)/dy[i];
            if (fabs(div)>m) m=fabs(div);
        }
    }
    return m;
}

void write_p_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %g %d\n",x-1,y-1,p[i],solid[i]);
        }
        fprintf(F,"\n");
    }
    fclose(F);
}

void write_v_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %e %e %d\n",x-1,y-1,vx[i],vy[i],solid[i]);
        }
        fprintf(F,"\n");
    }
    fclose(F);
}

void write_flow_rate_profile(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    fprintf(F,"# x  Q(x)\n");
    for (int x=1;x<=Nx;x++)
        fprintf(F,"%d %g\n",x-1,flow_rate_at(x));
    fclose(F);
}

void solver()
{
    double tolerance=1e-7;
    int it=0, max_it=15000;
    double vnorm_old=0.0, check=1.0;

    double inner_tol=1e-6;
    int max_inner=150;

    FILE *L=fopen("convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence_log.dat");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  max_abs_div  Q_wide  Q_throat\n");

    do{
        apply_boundary_conditions();
        compute_shear_and_viscosity();
        interpolate_velocity_faces();

        predict_velocity();
        apply_velocity_bc(vxstar,vystar);

        double sor_resid;
        int sor_iters=solve_pressure_sor(inner_tol,max_inner,&sor_resid);

        correct_velocity();
        apply_boundary_conditions();

        double vnorm=measure_velocity_l2();
        if (it>0 && vnorm_old>1e-12) check=fabs(vnorm-vnorm_old)/vnorm_old;
        vnorm_old=vnorm;

        if (it%200==0 || it==max_it-1){
            double mdiv=max_abs_divergence();
            double Qw=flow_rate_at(25), Qt=flow_rate_at(75);
            fprintf(L,"%d %g %g %d %g %g %g %g\n",it,vnorm,check,sor_iters,sor_resid,mdiv,Qw,Qt);
            printf("it %d  |v| %.6e  check %.3e  (SOR %d, resid %.2e)  max|div| %.2e  Q_wide %.5f  Q_throat %.5f\n",
                   it,vnorm,check,sor_iters,sor_resid,mdiv,Qw,Qt);
        }
        it++;
    } while(check>tolerance && it<max_it);

    printf("Stopped after %d iterations, check=%.3e\n",it,check);
    fclose(L);
}

int main(int argc,char **argv)
{
    setup();
    write_p_2d("start_pressure.dat");
    write_v_2d("start_velocity.dat");
    solver();
    write_p_2d("final_pressure.dat");
    write_v_2d("final_velocity.dat");
    write_flow_rate_profile("flow_rate_profile.dat");
    return 0;
}
