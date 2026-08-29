/* fetched from project for figure regeneration - see project doc for full comments */

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

double POWERLAW_K = 1.0;
double POWERLAW_N = 0.5;
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

/* -----------------------------------------------------------------------
   Explicit-diffusion stability check -- the "diffusion/Fourier number",
   not to be confused with the advective Courant number U*dt/dx (that one
   checks out fine here; it stays under 1e-3 for every case tried). The
   viscous term in predict_velocity_nonnewtonian() is a standard 2D
   central-difference discretization of div(2*eta*D), and an explicit
   scheme on that stencil is only stable if

       Fo = eta_max * dt / dx^2  <=  0.25

   eta_max has a hard upper bound that's knowable before the solver ever
   runs: eta = K*D_eff^(n-1) is maximized exactly where D_eff is smallest,
   i.e. D_eff=POWERLAW_EPSILON, giving eta_zero_shear = K*epsilon^(n-1) --
   the same value init_fields() already uses as the starting eta
   everywhere. So Fo can (and should) be checked once, up front, from
   (K, n, epsilon, dt, dx) alone -- no need to run anything first.

   Found empirically while tracking down a wrong-looking result: K=0.5,
   n=0.6 gives Fo=0.126 and the run is clean; K=1.0, n=0.6 gives Fo=0.251
   -- just over the line -- and shows exactly the signature of marginal
   instability (max|div(v)| plateaus around 1e-3-2e-3 instead of
   continuing down toward 1e-4/1e-5, and the resulting velocity profile
   comes out qualitatively wrong: power-law reads as MORE parabolic than
   Newtonian in the syringe throat, the opposite of the textbook
   shear-thinning blunting). K=1.0, n=0.4 gives Fo=3.98, ~16x over, and
   degrades much further. This check exists so that's caught before
   spending a run and a plot on it, not after. */
void check_stability(void)
{
    double eta_zero_shear = POWERLAW_K * pow(POWERLAW_EPSILON, POWERLAW_N - 1.0);
    double h2 = dx[site2index(1,1)] * dy[site2index(1,1)];   /* dx==dy here, kept general */
    double Fo = eta_zero_shear * dt / h2;
    double dt_max_stable = 0.25 * h2 / eta_zero_shear;

    printf("--- explicit-diffusion stability check ---\n");
    printf("  eta_zero_shear (K*eps^(n-1), eta's hard upper bound) = %.4f\n", eta_zero_shear);
    printf("  Fo = eta_zero_shear * dt / dx^2                      = %.4f  (need <= 0.25)\n", Fo);
    if (Fo > 0.25) {
        printf("  *** UNSTABLE: dt=%.6g is %.2fx above the safe limit.\n", dt, Fo/0.25);
        printf("  *** Use dt <= %.6g instead (or raise POWERLAW_EPSILON / lower K) before trusting this run.\n",
               dt_max_stable);
    } else {
        printf("  OK: dt=%.6g is within the safe limit (dt_max_stable=%.6g, margin %.1fx).\n",
               dt, dt_max_stable, dt_max_stable/dt);
    }
    printf("-------------------------------------------\n");
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

            double DdotD=Dxx*Dxx+Dyy*Dyy+2.0*Dxy*Dxy;
            double gamma_dot=sqrt(2.0*DdotD);
            double Deff=sqrt(gamma_dot*gamma_dot+POWERLAW_EPSILON*POWERLAW_EPSILON);
            shear_mag[i]=gamma_dot;
            eta[i]=POWERLAW_K*pow(Deff,POWERLAW_N-1.0);
        }
    }
    apply_eta_bc();
}

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

void write_shear_eta(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    fprintf(F,"# x  y  gamma_dot  eta_a\n");
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
    int it=0, max_it=2000000;
    double vnorm_old=0.0, check=1.0;
    /* Require at least MIN_TIME of physical (simulation) time before the check above is
       trusted. Reason: check is a step-to-step comparison, so at a very small dt it can
       look "converged" simply because one step is a tiny slice of physical time -- long
       before the flow has actually developed -- not because it has reached steady state. */
    double MIN_TIME=10.0;
    double inner_tol=1e-6;
    int max_inner=1000;

    FILE *L=fopen("finite_vol_5_duct_SIMPLE_shearthinning_convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence log");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  mdiv  Q_in  Q_mid  Q_out  eta_min  eta_max\n");

    do{
        apply_velocity_bc(vx,vy);
        interp_v_faces(vx,vy,vxe,vxw,vxn,vxs,vye,vyw,vyn,vys);

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
    } while ((check>tolerance || (double)it*dt<MIN_TIME) && it<max_it);

    printf("Stopped after %d iterations, check=%.3e (converged=%s)\n",
           it,check,(check<=tolerance)?"yes":"NO -- hit iteration cap");
    fclose(L);
}

int main(void)
{
    setup();
    check_stability();
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