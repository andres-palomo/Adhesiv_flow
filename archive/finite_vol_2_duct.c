/* =====================================================================
   MILESTONE 2 -- boundary conditions: giving the flow a physical shape
   =====================================================================

   Built forward from finite_vol_1_baseline.c. The baseline solved the
   right equations (Korzec's finite-volume momentum update and pressure
   Poisson/Jacobi treatment) on the wrong domain: a periodic torus with
   no real inlet, outlet, or walls. This file keeps the same finite-volume
   machinery and replaces the domain with an actual duct:

     - a ghost-cell ring around the Nx x Ny physical interior (standard
       technique for boundary conditions in finite-volume codes -- see
       tutorial.tex, Section 4);
     - Dirichlet inlet velocity (fluid is pushed in at U_in);
     - zero-gradient outlet velocity with a Dirichlet reference pressure
       (fluid is free to leave without the outlet forcing a profile);
     - no-slip top/bottom walls, enforced by mirroring the ghost value
       so the face-averaged velocity at the wall is exactly zero.

   The pressure Poisson equation (tutorial.tex Eq. 22-25) is now solved
   to convergence via SOR every outer step, rather than the single
   relaxation pass in the baseline -- see solve_pressure_sor() below.
   This is still the *legacy* pressure formula: it is derived by
   assuming div(v)=0 already holds, rather than measuring the actual
   divergence and driving it to zero.

   That assumption is NOT harmless, even here, with no geometry change
   at all. A converged run of this file settles at max|div(v)| ~ 4e-2
   everywhere in the interior (it does not shrink further with more
   iterations -- see the mdiv column in the convergence log), and the
   flow rate Q(x) = sum_y vx(x,y)*dy drops from ~1.67 near the inlet to
   ~0.70 near the outlet, instead of staying flat at U_in*Ny = 2.0 as
   incompressible continuity requires. So the legacy formula's failure
   to conserve mass is already large on a plain duct; a real geometry
   change (finite_vol_3_seringe.c) makes it worse and easier to see,
   but it is not the origin of the problem. finite_vol_4_seringe_SIMPLE.c
   is where this actually gets fixed, by replacing this formula with a
   real projection step.

   Viscosity is Newtonian and constant (MU_CONST); shear-thinning comes
   in at finite_vol_5_duct_SIMPLE_shearthinning.c.

   Build:  gcc -O2 -o finite_vol_2_duct finite_vol_2_duct.c -lm
   Run:    ./finite_vol_2_duct
   ===================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <err.h>

/* ---------------------------------------------------------------------
   Grid. Nx=100, Ny=20 is the resolution used across every file from
   here on (finite_vol_2 through finite_vol_6), so results are directly
   comparable across the whole story.
   --------------------------------------------------------------------- */
#define Nx 100
#define Ny 20
#define NXG (Nx+2)
#define NYG (Ny+2)
#define Ntot (NXG*NYG)

double dt = 0.001;
double U_in  = 0.1;   /* inlet velocity (Dirichlet)            */
double P_out = 0.0;   /* outlet reference pressure (Dirichlet) */
double MU_CONST = 1.0; /* Newtonian dynamic viscosity           */

double p[Ntot], vx[Ntot], vy[Ntot];
int    E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];
double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];

double b_src[Ntot];
double SOR_OMEGA = 1.7; /* 1.0 = plain Gauss-Seidel, (1,2) = SOR */

/* Cell index from (x,y). No wraparound -- x,y must lie in
   [0,NXG-1]x[0,NYG-1]; the ghost ring is what keeps E/W/N/S valid at
   the domain edges without the baseline's periodic wrap. */
int site2index(int x, int y)
{
    return x + y*NXG;
}

void geometry(void)
{
    for (int i=0;i<Ntot;i++) dx[i]=dy[i]=1.0;
}

void init_control_points(void)
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

void init_fields(void)
{
    for (int i=0;i<Ntot;i++){ p[i]=0.0; vx[i]=0.0; vy[i]=0.0; }
}

/* Boundary conditions on the ghost ring -- tutorial.tex Section 4.
   Called once before the loop starts and again after every velocity
   update, since the ghost values depend on the (changing) interior
   values next to them. */
void apply_boundary_conditions(void)
{
    /* Inlet (x=0): Dirichlet velocity, zero-gradient pressure */
    for (int y=1;y<=Ny;y++){
        int g=site2index(0,y);
        vx[g]=U_in;
        vy[g]=0.0;
        int in1=site2index(1,y);
        p[g]=p[in1];
    }

    /* Outlet (x=Nx+1): zero-gradient velocity, Dirichlet pressure */
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y);
        int in1=site2index(Nx,y);
        vx[g]=vx[in1];
        vy[g]=vy[in1];
        p[g]=P_out;
    }

    /* Bottom wall (y=0): no-slip via mirrored ghost value
       (face average between ghost and first interior row = 0) */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0);
        int in1=site2index(x,1);
        vx[g]=-vx[in1];
        vy[g]=-vy[in1];
        p[g]=p[in1];
    }

    /* Top wall (y=Ny+1): no-slip via mirrored ghost value */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1);
        int in1=site2index(x,Ny);
        vx[g]=-vx[in1];
        vy[g]=-vy[in1];
        p[g]=p[in1];
    }
}

void setup(void)
{
    geometry();
    init_control_points();
    init_fields();
    apply_boundary_conditions();
}

void interpolate_faces(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            pe[i]=(p[E[i]]+p[i])/2;
            pw[i]=(p[W[i]]+p[i])/2;
            pn[i]=(p[N[i]]+p[i])/2;
            ps[i]=(p[S[i]]+p[i])/2;
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

/* Explicit momentum update -- tutorial.tex Eq. 18, i.e. Korzec Eq.
   3.59: advection + pressure gradient + Newtonian viscous diffusion,
   all evaluated at time t. */
void update_velocity(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
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

            double visc_x = MU_CONST*((dxvxe-dxvxw)/dx[i] + (dyvxn-dyvxs)/dy[i]);
            double visc_y = MU_CONST*((dxvye-dxvyw)/dx[i] + (dyvyn-dyvys)/dy[i]);

            vx[i]+=dt*(-vx[i]*dxvx-vy[i]*dyvx-(pe[i]-pw[i])/dx[i]+visc_x);
            vy[i]+=dt*(-vx[i]*dxvy-vy[i]*dyvy-(pn[i]-ps[i])/dy[i]+visc_y);
        }
    }
}

/* Legacy pressure-Poisson source -- tutorial.tex Eq. 21, obtained by
   assuming div(v)=0 already holds. Fine here: this is a plain duct
   with no geometry change, so nothing forces a local mass imbalance
   the way a contraction would (see finite_vol_3_seringe.c). */
int solve_pressure_sor(double inner_tol, int max_inner, double *resid_out)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            double dxvx=(vxe[i]-vxw[i])/dx[i];
            double dxvy=(vye[i]-vyw[i])/dx[i];
            double dyvx=(vxn[i]-vxs[i])/dy[i];
            double dyvy=(vyn[i]-vys[i])/dy[i];
            b_src[i]=-dxvx*dxvx-dyvy*dyvy+2*dxvy*dyvx;
        }
    }

    int k; double resid=1e300;
    for (k=0;k<max_inner;k++){
        resid=0.0;
        for (int x=1;x<=Nx;x++){
            for (int y=1;y<=Ny;y++){
                int i=site2index(x,y);
                double pe_i=(p[E[i]]+p[i])/2;
                double pw_i=(p[W[i]]+p[i])/2;
                double pn_i=(p[N[i]]+p[i])/2;
                double ps_i=(p[S[i]]+p[i])/2;
                double val=-1.0*b_src[i]*dx[i]*dx[i]*dy[i]*dy[i]
                           +(pn_i+ps_i)*dx[i]*dx[i]+(pe_i+pw_i)*dy[i]*dy[i];
                double p_gs=(0.5*val)/(dx[i]*dx[i]+dy[i]*dy[i]);
                double p_old=p[i];
                p[i]=(1.0-SOR_OMEGA)*p_old+SOR_OMEGA*p_gs;
                double d=fabs(p[i]-p_old);
                if (d>resid) resid=d;
            }
        }
        apply_boundary_conditions();
        if (resid<inner_tol){ k++; break; }
    }
    *resid_out=resid;
    return k;
}

double measure_velocity_l2(void)
{
    double s=0.0;
    for (int x=1;x<=Nx;x++)
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            s+=vx[i]*vx[i]+vy[i]*vy[i];
        }
    return sqrt(s);
}

/* Volumetric flow rate through the vertical cross-section at column x
   (sum of vx*dy over the cells there). Mass conservation requires this
   to be constant along the duct, and equal to U_in*Ny at the inlet. */
double flow_rate_at(int x)
{
    double Q=0.0;
    for (int y=1;y<=Ny;y++){
        int i=site2index(x,y);
        Q+=vx[i]*dy[i];
    }
    return Q;
}

/* Maximum cell-wise |div(v)| over the interior -- the direct, local
   measure of how badly continuity is being violated. */
double max_abs_divergence(void)
{
    double m=0.0;
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            double vxe_=(vx[E[i]]+vx[i])/2, vxw_=(vx[W[i]]+vx[i])/2;
            double vyn_=(vy[N[i]]+vy[i])/2, vys_=(vy[S[i]]+vy[i])/2;
            double div=(vxe_-vxw_)/dx[i]+(vyn_-vys_)/dy[i];
            if (fabs(div)>m) m=fabs(div);
        }
    }
    return m;
}

void write_flow_rate_profile(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    fprintf(F,"# x  Q(x)\n");
    for (int x=1;x<=Nx;x++) fprintf(F,"%d %g\n",x-1,flow_rate_at(x));
    fclose(F);
}

void write_p_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %g\n",x-1,y-1,p[i]);
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
            fprintf(F,"%d %d %e %e\n",x-1,y-1,vx[i],vy[i]);
        }
        fprintf(F,"\n");
    }
    fclose(F);
}

void solver(void)
{
    double tolerance=1e-6;
    int it=0, max_it=100000;
    double vnorm_old=0.0, check=1.0;

    double inner_tol=1e-6;
    int max_inner=500;

    FILE *L=fopen("finite_vol_2_duct_convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence log");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  mdiv  Q_in  Q_mid  Q_out\n");

    do{
        apply_boundary_conditions();
        interpolate_faces();

        double sor_resid;
        int sor_iters=solve_pressure_sor(inner_tol,max_inner,&sor_resid);

        interpolate_faces();      /* refresh pe/pw/pn/ps from the converged pressure */
        update_velocity();
        apply_boundary_conditions();

        double vnorm=measure_velocity_l2();
        if (it>0 && vnorm_old>1e-12) check=fabs(vnorm-vnorm_old)/vnorm_old;
        vnorm_old=vnorm;

        if (it%2000==0 || it==max_it-1){
            double mdiv=max_abs_divergence();
            double q_in=flow_rate_at(1), q_mid=flow_rate_at(Nx/2), q_out=flow_rate_at(Nx);
            fprintf(L,"%d %g %g %d %g %g %g %g %g\n",
                    it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_mid,q_out);
            printf("it %d  |v| %.6e  check %.3e  (SOR %d, resid %.2e)  mdiv %.3e  Q(in/mid/out) %.4f/%.4f/%.4f\n",
                   it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_mid,q_out);
        }
        it++;
    } while (check>tolerance && it<max_it);

    printf("Stopped after %d iterations, check=%.3e (converged=%s)\n",
           it,check,(check<=tolerance)?"yes":"NO -- hit iteration cap");
    fclose(L);
}

int main(void)
{
    setup();
    write_p_2d("finite_vol_2_duct_start_pressure.dat");
    write_v_2d("finite_vol_2_duct_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_2_duct_final_pressure.dat");
    write_v_2d("finite_vol_2_duct_final_velocity.dat");
    write_flow_rate_profile("finite_vol_2_duct_flow_rate_profile.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}
