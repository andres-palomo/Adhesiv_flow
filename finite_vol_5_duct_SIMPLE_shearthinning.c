/* =====================================================================
   MILESTONE 5 -- shear-thinning viscosity (power law) on the plain duct

   Built forward from finite_vol_4_seringe_SIMPLE.c, but on the plain
   duct geometry (no solid mask, no constriction -- that combination is
   finite_vol_6_seringe_SIMPLE_shearthinning.c). Two things change from
   the Newtonian SIMPLE file:

     1. Viscosity is no longer a constant MU_CONST. It is the power-law
        (Ostwald-de Waele) apparent viscosity, tutorial.tex Eq. 40:
              eta_a = K * D_eff^(n-1),   D_eff = sqrt(|D|^2 + eps^2)
        with |D| the rate-of-strain invariant built from the local
        velocity gradients (tutorial.tex Section 7.1, sourced from
        Yadav (2026), Section 2.3.5 "Power Law Fluids"), and D_eff the
        zero-shear-regularized version that keeps eta_a finite where
        the flow is locally at rest (tutorial.tex Eq. 43-44).

     2. Because eta_a now varies from cell to cell, the viscous term in
        the momentum update can no longer collapse to a scalar
        Laplacian; it has to stay in the full stress-divergence form
        div(2*eta_a*D) (tutorial.tex Eq. 46-47), built here from
        face-averaged viscosities and, for the tangential derivatives
        that need a value at a face that isn't in the E/W/N/S stencil,
        the average of the two cell-centred tangential gradients that
        share that face -- exactly as tutorial.tex Section 7.5 states,
        since this grid has no diagonal-neighbour stencil.

   The nonlinear coupling eta_a(v)-v is resolved by Picard (lagged)
   iteration, tutorial.tex Eq. 48-49: eta_a is evaluated once per outer
   step from the *current* velocity field, held fixed while that step's
   predict/project/correct solves for the *next* velocity, which then
   produces the eta_a used on the following step. This is exactly the
   outer time-marching loop already in the code, with a viscosity
   update folded into it -- no separate inner loop is needed.

   The SIMPLE / projection pressure treatment itself (predict, measure
   div(v*), solve, correct) is unchanged from finite_vol_4 -- that part
   does not need to know anything about how v* was computed, so nothing
   about it has to be re-derived for the non-Newtonian case.

   Build:  gcc -O2 -o finite_vol_5_duct_SIMPLE_shearthinning finite_vol_5_duct_SIMPLE_shearthinning.c -lm
   Run:    ./finite_vol_5_duct_SIMPLE_shearthinning
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

double dt = 0.001;
double U_in  = 0.1;
double P_out = 0.0;

/* Power-law parameters -- same (K,n) as used for the Case 2/Case 3
   viscosity-profile comparison already in tutorial.tex Figure 2. */
double POWERLAW_K = 0.5;
double POWERLAW_N = 0.6;
double POWERLAW_EPSILON = 1e-6;

double p[Ntot], vx[Ntot], vy[Ntot];
double vx_star[Ntot], vy_star[Ntot];
int    E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];

double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];
double vxe_s[Ntot], vxw_s[Ntot], vxn_s[Ntot], vxs_s[Ntot];
double vye_s[Ntot], vyw_s[Ntot], vyn_s[Ntot], vys_s[Ntot];

/* Cell-centred velocity gradients (needed both for D and, for the
   tangential terms, by neighbouring cells when building face stresses)
   and the resulting apparent viscosity. */
double dxvx_c[Ntot], dyvy_c[Ntot], dxvy_c[Ntot], dyvx_c[Ntot];
double eta[Ntot], shear_mag[Ntot];

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
    for (int i=0;i<Ntot;i++){
        p[i]=0.0; vx[i]=0.0; vy[i]=0.0; vx_star[i]=0.0; vy_star[i]=0.0;
        eta[i]=POWERLAW_K*pow(POWERLAW_EPSILON,POWERLAW_N-1.0);
    }
}

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

void apply_pressure_bc(void)
{
    for (int y=1;y<=Ny;y++) p[site2index(0,y)]=p[site2index(1,y)];
    for (int y=1;y<=Ny;y++) p[site2index(Nx+1,y)]=P_out;
    for (int x=1;x<=Nx;x++) p[site2index(x,0)]=p[site2index(x,1)];
    for (int x=1;x<=Nx;x++) p[site2index(x,Ny+1)]=p[site2index(x,Ny)];
}

/* Ghost-ring values for eta itself, so a face average involving a
   ghost cell (E/W/N/S of a boundary interior cell) is well defined --
   zero-gradient is the natural choice, matching the pressure BC. */
void apply_eta_bc(void)
{
    for (int y=1;y<=Ny;y++){ eta[site2index(0,y)]=eta[site2index(1,y)]; eta[site2index(Nx+1,y)]=eta[site2index(Nx,y)]; }
    for (int x=1;x<=Nx;x++){ eta[site2index(x,0)]=eta[site2index(x,1)]; eta[site2index(x,Ny+1)]=eta[site2index(x,Ny)]; }
}

void setup(void)
{
    geometry();
    init_control_points();
    init_fields();
    apply_velocity_bc(vx,vy);
    apply_pressure_bc();
    apply_eta_bc();
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

/* tutorial.tex Section 7.1: rate-of-strain tensor at cell centres from
   the face-averaged velocity field, and Eq. 40/43-44: regularized
   power-law apparent viscosity. |D| here uses the "paper convention"
   sqrt(D:D) (no factor of 2) -- see tutorial.tex Eq. 45's discussion. */
void compute_viscosity(double *ux, double *uy,
                        double *uxe, double *uxw, double *uxn, double *uxs,
                        double *uye, double *uyw, double *uyn, double *uys)
{
    (void)ux; (void)uy;
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            double Dxx=(uxe[i]-uxw[i])/dx[i];
            double Dyy=(uyn[i]-uys[i])/dy[i];
            double dxvy=(uye[i]-uyw[i])/dx[i];
            double dyvx=(uxn[i]-uxs[i])/dy[i];
            double Dxy=0.5*(dyvx+dxvy);

            dxvx_c[i]=Dxx; dyvy_c[i]=Dyy; dxvy_c[i]=dxvy; dyvx_c[i]=dyvx;

            double Dmag=sqrt(Dxx*Dxx+Dyy*Dyy+2.0*Dxy*Dxy);
            double Deff=sqrt(Dmag*Dmag+POWERLAW_EPSILON*POWERLAW_EPSILON);
            shear_mag[i]=Dmag;
            eta[i]=POWERLAW_K*pow(Deff,POWERLAW_N-1.0);
        }
    }
    apply_eta_bc();
}

/* Step 1 (Picard-lagged): advance velocity with advection + the full
   non-Newtonian stress-divergence term div(2*eta_a*D), tutorial.tex
   Eq. 46-47 -- no pressure gradient, same as the Newtonian predictor.
   eta[] and the tangential centre-gradients (dxvy_c, dyvx_c) must
   already be current, from compute_viscosity(). */
void predict_velocity_nonnewtonian(void)
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

            double eta_e=0.5*(eta[i]+eta[E[i]]);
            double eta_w=0.5*(eta[i]+eta[W[i]]);
            double eta_n=0.5*(eta[i]+eta[N[i]]);
            double eta_s=0.5*(eta[i]+eta[S[i]]);

            /* tangential derivatives at faces = average of the two
               cell-centred tangential gradients sharing that face */
            double dyvx_e=0.5*(dyvx_c[i]+dyvx_c[E[i]]);
            double dyvx_w=0.5*(dyvx_c[i]+dyvx_c[W[i]]);
            double dxvy_n=0.5*(dxvy_c[i]+dxvy_c[N[i]]);
            double dxvy_s=0.5*(dxvy_c[i]+dxvy_c[S[i]]);

            double tau_xx_e=2*eta_e*dxvxe;
            double tau_xx_w=2*eta_w*dxvxw;
            double tau_xy_n=eta_n*(dyvxn+dxvy_n);
            double tau_xy_s=eta_s*(dyvxs+dxvy_s);

            double tau_yy_n=2*eta_n*dyvyn;
            double tau_yy_s=2*eta_s*dyvys;
            double tau_xy_e=eta_e*(dxvye+dyvx_e);
            double tau_xy_w=eta_w*(dxvyw+dyvx_w);

            double visc_x=(tau_xx_e-tau_xx_w)/dx[i]+(tau_xy_n-tau_xy_s)/dy[i];
            double visc_y=(tau_xy_e-tau_xy_w)/dx[i]+(tau_yy_n-tau_yy_s)/dy[i];

            vx_star[i]=vx[i]+dt*(-vx[i]*dxvx-vy[i]*dyvx+visc_x);
            vy_star[i]=vy[i]+dt*(-vx[i]*dxvy-vy[i]*dyvy+visc_y);
        }
    }
}

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

/* Shear-rate and viscosity profile writer, in the same column layout
   as the old reference files' *_shear_eta.dat (x y |D| eta) so
   make_figures.py's Figure 2 keeps working unchanged. */
void write_shear_eta(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    fprintf(F,"# x  y  |D|  eta_a\n");
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %g %g\n",x-1,y-1,shear_mag[i],eta[i]);
        }
    }
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
    int it=0, max_it=100000;
    double vnorm_old=0.0, check=1.0;
    double inner_tol=1e-6;
    int max_inner=500;

    FILE *L=fopen("finite_vol_5_duct_SIMPLE_shearthinning_convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence log");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  mdiv  Q_in  Q_mid  Q_out  eta_min  eta_max\n");

    do{
        apply_velocity_bc(vx,vy);
        interp_v_faces(vx,vy,vxe,vxw,vxn,vxs,vye,vyw,vyn,vys);

        /* Picard: eta_a^(k) from the *current* (already-converged-so-far)
           velocity field, then lagged/held fixed through this step's solve. */
        compute_viscosity(vx,vy,vxe,vxw,vxn,vxs,vye,vyw,vyn,vys);

        predict_velocity_nonnewtonian();
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
            double eta_min=1e300, eta_max=-1e300;
            for (int x=1;x<=Nx;x++) for (int y=1;y<=Ny;y++){
                double e=eta[site2index(x,y)];
                if (e<eta_min) eta_min=e; if (e>eta_max) eta_max=e;
            }
            fprintf(L,"%d %g %g %d %g %g %g %g %g %g %g\n",
                    it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_mid,q_out,eta_min,eta_max);
            printf("it %d  |v| %.6e  check %.3e  (SOR %d, resid %.2e)  mdiv %.3e  Q(in/mid/out) %.4f/%.4f/%.4f  eta[%.3f,%.3f]\n",
                   it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_mid,q_out,eta_min,eta_max);
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
    write_p_2d("finite_vol_5_duct_SIMPLE_shearthinning_start_pressure.dat");
    write_v_2d("finite_vol_5_duct_SIMPLE_shearthinning_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_5_duct_SIMPLE_shearthinning_final_pressure.dat");
    write_v_2d("finite_vol_5_duct_SIMPLE_shearthinning_final_velocity.dat");
    write_flow_rate_profile("finite_vol_5_duct_SIMPLE_shearthinning_flow_rate_profile.dat");
    write_shear_eta("finite_vol_5_duct_SIMPLE_shearthinning_shear_eta.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}
