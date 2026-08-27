/* =====================================================================
   PLAIN DUCT + SIMPLE + NEWTONIAN -- the missing fair baseline

   This file fills a gap the project's original six milestones didn't
   need but a fair report comparison does: finite_vol_2_duct.c (plain
   duct, Newtonian) still uses the *legacy* pressure formula, which
   Table tab:case-comparison in tutorial.tex shows loses 58% of the
   flow rate along a plain duct with no geometry change at all
   (Q_in=1.673, Q_out=0.704, target 2.0). That makes any velocity-
   profile comparison against it confounded: differences downstream
   are dominated by the solver's mass-conservation bug, not by
   rheology. finite_vol_4_seringe_SIMPLE.c already fixed this for the
   syringe geometry (Newtonian + SIMPLE); this file is the same fix,
   Newtonian + SIMPLE, but on the plain duct with no solid mask at
   all -- exactly finite_vol_5_duct_SIMPLE_shearthinning.c's geometry
   and SIMPLE machinery, with the power-law viscosity closure removed
   and replaced by the same constant MU_CONST Newtonian viscosity used
   in finite_vol_2_duct.c and finite_vol_4_seringe_SIMPLE.c.

   With this file, finite_vol_duct_SIMPLE.c (Newtonian) can be compared
   directly against finite_vol_5_duct_SIMPLE_shearthinning.c (power-law)
   on the plain duct -- both SIMPLE, both mass-conserving -- giving the
   duct-side counterpart to the syringe-side comparison already fixed
   in make_figures.py (finite_vol_4_seringe_SIMPLE vs
   finite_vol_6_seringe_SIMPLE_shearthinning, at x=90).

   Build:  gcc -O2 -o finite_vol_duct_SIMPLE finite_vol_duct_SIMPLE.c -lm
   Run:    ./finite_vol_duct_SIMPLE
   ===================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <err.h>

#define Nx 100
#define Ny 20
#define NXG (Nx+2)
#define NYG (Ny+2)
#define Ntot (NXG*NYG)

double dt = 0.0001;
double U_in  = 0.1;
double P_out = 0.0;
double MU_CONST = 1.0;   /* same value used in finite_vol_2_duct.c / finite_vol_4_seringe_SIMPLE.c */

double p[Ntot], vx[Ntot], vy[Ntot];
double vx_star[Ntot], vy_star[Ntot];
int    E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];

double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];
double vxe_s[Ntot], vxw_s[Ntot], vxn_s[Ntot], vxs_s[Ntot];
double vye_s[Ntot], vyw_s[Ntot], vyn_s[Ntot], vys_s[Ntot];

double b_src[Ntot];
double SOR_OMEGA = 1.7;

int site2index(int x, int y) { return x + y*NXG; }

void geometry(void) { for (int i=0;i<Ntot;i++) dx[i]=dy[i]=1.0; }

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
    for (int i=0;i<Ntot;i++){ p[i]=0.0; vx[i]=0.0; vy[i]=0.0; vx_star[i]=0.0; vy_star[i]=0.0; }
}

/* Same velocity BCs as every other file: Dirichlet inlet, zero-gradient
   outlet, mirrored no-slip top/bottom walls (tutorial.tex Section 4). */
void apply_velocity_bc(double *ux, double *uy)
{
    for (int y=1;y<=Ny;y++){ int g=site2index(0,y); ux[g]=U_in; uy[g]=0.0; }
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y), in1=site2index(Nx,y);
        ux[g]=ux[in1]; uy[g]=uy[in1];
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0), in1=site2index(x,1);
        ux[g]=-ux[in1]; uy[g]=-uy[in1];
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1), in1=site2index(x,Ny);
        ux[g]=-ux[in1]; uy[g]=-uy[in1];
    }
}

/* Pressure BCs isolated from velocity BCs, same as finite_vol_4/5/6:
   zero-gradient inlet/walls, Dirichlet outlet. */
void apply_pressure_bc(void)
{
    for (int y=1;y<=Ny;y++) p[site2index(0,y)]=p[site2index(1,y)];
    for (int y=1;y<=Ny;y++) p[site2index(Nx+1,y)]=P_out;
    for (int x=1;x<=Nx;x++) p[site2index(x,0)]=p[site2index(x,1)];
    for (int x=1;x<=Nx;x++) p[site2index(x,Ny+1)]=p[site2index(x,Ny)];
}

void setup(void)
{
    geometry();
    init_control_points();
    init_fields();
    apply_velocity_bc(vx,vy);
    apply_pressure_bc();
}

void interp_v_faces(double *ux, double *uy,
                     double *uxe, double *uxw, double *uxn, double *uxs,
                     double *uye, double *uyw, double *uyn, double *uys)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            uxe[i]=(ux[E[i]]+ux[i])/2; uxw[i]=(ux[W[i]]+ux[i])/2;
            uxn[i]=(ux[N[i]]+ux[i])/2; uxs[i]=(ux[S[i]]+ux[i])/2;
            uye[i]=(uy[E[i]]+uy[i])/2; uyw[i]=(uy[W[i]]+uy[i])/2;
            uyn[i]=(uy[N[i]]+uy[i])/2; uys[i]=(uy[S[i]]+uy[i])/2;
        }
    }
}

void interp_p_faces(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            pe[i]=(p[E[i]]+p[i])/2; pw[i]=(p[W[i]]+p[i])/2;
            pn[i]=(p[N[i]]+p[i])/2; ps[i]=(p[S[i]]+p[i])/2;
        }
    }
}

/* Step 1: predict v* with advection + Newtonian viscous diffusion only,
   no pressure -- identical formula to finite_vol_4_seringe_SIMPLE.c's
   predict_velocity(), just without the solid-cell skip. */
void predict_velocity(void)
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

            vx_star[i]=vx[i]+dt*(-vx[i]*dxvx-vy[i]*dyvx+visc_x);
            vy_star[i]=vy[i]+dt*(-vx[i]*dxvy-vy[i]*dyvy+visc_y);
        }
    }
}

/* Step 2: real Poisson source, div(v*)/dt -- measured, not assumed. */
void compute_b_src(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            double div=(vxe_s[i]-vxw_s[i])/dx[i]+(vyn_s[i]-vys_s[i])/dy[i];
            b_src[i]=div/dt;
        }
    }
}

/* Step 3: SOR solve of Lap(phi)=b_src, same stencil as every other file. */
int solve_pressure_sor(double inner_tol, int max_inner, double *resid_out)
{
    int k; double resid=1e300;
    for (k=0;k<max_inner;k++){
        resid=0.0;
        for (int x=1;x<=Nx;x++){
            for (int y=1;y<=Ny;y++){
                int i=site2index(x,y);
                double pe_i=(p[E[i]]+p[i])/2, pw_i=(p[W[i]]+p[i])/2;
                double pn_i=(p[N[i]]+p[i])/2, ps_i=(p[S[i]]+p[i])/2;
                double val=-1.0*b_src[i]*dx[i]*dx[i]*dy[i]*dy[i]
                           +(pn_i+ps_i)*dx[i]*dx[i]+(pe_i+pw_i)*dy[i]*dy[i];
                double p_gs=(0.5*val)/(dx[i]*dx[i]+dy[i]*dy[i]);
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

/* Step 4: project v* back toward div(v)=0. */
void correct_velocity(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            vx[i]=vx_star[i]-dt*(pe[i]-pw[i])/dx[i];
            vy[i]=vy_star[i]-dt*(pn[i]-ps[i])/dy[i];
        }
    }
}

double measure_velocity_l2(void)
{
    double s=0.0;
    for (int x=1;x<=Nx;x++) for (int y=1;y<=Ny;y++){ int i=site2index(x,y); s+=vx[i]*vx[i]+vy[i]*vy[i]; }
    return sqrt(s);
}

double flow_rate_at(int x)
{
    double Q=0.0;
    for (int y=1;y<=Ny;y++){ int i=site2index(x,y); Q+=vx[i]*dy[i]; }
    return Q;
}

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
        for (int y=1;y<=Ny;y++){ int i=site2index(x,y); fprintf(F,"%d %d %g\n",x-1,y-1,p[i]); }
        fprintf(F,"\n");
    }
    fclose(F);
}

void write_v_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){ int i=site2index(x,y); fprintf(F,"%d %d %e %e\n",x-1,y-1,vx[i],vy[i]); }
        fprintf(F,"\n");
    }
    fclose(F);
}

void solver(void)
{
    double tolerance=1e-6;
    int it=0, max_it=2000000;
    double vnorm_old=0.0, check=1.0;
    /* Require at least MIN_TIME of physical (simulation) time before the check above is
       trusted. Reason: check is a step-to-step comparison, so at a very small dt it can
       look "converged" simply because one step is a tiny slice of physical time -- long
       before the flow has actually developed -- not because it has reached steady state. */
    double MIN_TIME=10.0;

    double inner_tol=1e-6;
    int max_inner=1000;

    FILE *L=fopen("finite_vol_duct_SIMPLE_convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence log");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  mdiv  Q_in  Q_mid  Q_out\n");

    do{
        apply_velocity_bc(vx,vy);
        interp_v_faces(vx,vy,vxe,vxw,vxn,vxs,vye,vyw,vyn,vys);

        predict_velocity();
        apply_velocity_bc(vx_star,vy_star);
        interp_v_faces(vx_star,vy_star,vxe_s,vxw_s,vxn_s,vxs_s,vye_s,vyw_s,vyn_s,vys_s);

        compute_b_src();
        double sor_resid;
        int sor_iters=solve_pressure_sor(inner_tol,max_inner,&sor_resid);
        interp_p_faces();

        correct_velocity();
        apply_velocity_bc(vx,vy);

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
    } while ((check>tolerance || (double)it*dt<MIN_TIME) && it<max_it);

    printf("Stopped after %d iterations, check=%.3e (converged=%s)\n",
           it,check,(check<=tolerance)?"yes":"NO -- hit iteration cap");
    fclose(L);
}

int main(void)
{
    setup();
    write_p_2d("finite_vol_duct_SIMPLE_start_pressure.dat");
    write_v_2d("finite_vol_duct_SIMPLE_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_duct_SIMPLE_final_pressure.dat");
    write_v_2d("finite_vol_duct_SIMPLE_final_velocity.dat");
    write_flow_rate_profile("finite_vol_duct_SIMPLE_flow_rate_profile.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}