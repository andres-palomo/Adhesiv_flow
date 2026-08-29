/* =====================================================================
   Syringe geometry, power-law shear-thinning viscosity, SIMPLE method.

   Build:  gcc -O2 -o finite_vol_6_syringe_SIMPLE_shearthinning finite_vol_6_syringe_SIMPLE_shearthinning.c -lm
   Run:    ./finite_vol_6_syringe_SIMPLE_shearthinning
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

/* Syringe geometry: one contraction */
#define CONTRACTION_X (Nx/2)
#define THROAT_Y_LO (Ny/4+1)
#define THROAT_Y_HI (3*Ny/4)

double dt = 0.0001;
double U_in  = 0.1;
double P_out = 0.0;

double POWERLAW_K = 1.0;        /* power-law consistency index K */
double POWERLAW_N = 0.5;        /* power-law flow-behavior index n (n<1 => shear-thinning) */
double POWERLAW_EPSILON = 1e-6; /* regularization floor on the strain rate, avoids eta -> infinity */

double p[Ntot], vx[Ntot], vy[Ntot];          /* current (corrected, divergence-free) pressure and velocity */
double vx_star[Ntot], vy_star[Ntot];         /* predicted velocity v* (based on momentum-only / before pressure correction) */
int    E[Ntot], W[Ntot], N[Ntot], S[Ntot];   /* neighbor index of each cell */
double dx[Ntot], dy[Ntot];                   /* cell size (uniform, = 1.0 everywhere) */
int    is_solid[Ntot];                       /* 1 = cell is a wall,
                                                0 = cell is a fluid */

double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];              /* p interpolated to the E/W/N/S faces */
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];          /* vx interpolated to the E/W/N/S faces */
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];          /* vy interpolated to the E/W/N/S faces */
double vxe_s[Ntot], vxw_s[Ntot], vxn_s[Ntot], vxs_s[Ntot];  /* vx_star interpolated to the E/W/N/S faces */
double vye_s[Ntot], vyw_s[Ntot], vyn_s[Ntot], vys_s[Ntot];  /* vy_star interpolated to the E/W/N/S faces */

double dxvx_c[Ntot], dyvy_c[Ntot], dxvy_c[Ntot], dyvx_c[Ntot];  /* cell-centered strain-rate components */

double eta[Ntot];        /* apparent viscosity from the power law */
double shear_mag[Ntot];  /* scalar shear rate gamma_dot at each cell*/

double b_src[Ntot];      /* pressure-Poisson right-hand side: div(v*)/dt, actually measured */
double SOR_OMEGA = 1.7;  /* omega of SOR method */

/* grid coordinate to latice (1D array index). */
int site2index(int x, int y) { return x + y*NXG; }

/* cell size. */
void geometry(void) { for (int i=0;i<Ntot;i++) dx[i]=dy[i]=1.0; }

/* Precompute each cell's E/W/N/S neighbor index */
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

/* Mark the solid cells. */
void init_solid_mask(void)
{
    for (int i=0;i<Ntot;i++) is_solid[i]=0;
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            if (x>=CONTRACTION_X && (y<THROAT_Y_LO || y>THROAT_Y_HI))
                is_solid[site2index(x,y)]=1;
        }
    }
}

/* Start from rest, zero pressure, and initialize eta everywhere
   (including inside the solid mask) to its zero-shear-rate value
   K*epsilon^(n-1). */
void init_fields(void)
{
    for (int i=0;i<Ntot;i++){
        p[i]=0.0; vx[i]=0.0; vy[i]=0.0; vx_star[i]=0.0; vy_star[i]=0.0;
        eta[i]=POWERLAW_K*pow(POWERLAW_EPSILON,POWERLAW_N-1.0);
    }
}

/* Velocity at ghost-cell BC, by using mirrored (no-slip) top/bottom walls. */
void apply_velocity_bc(double *ux, double *uy)
{
    /* inlet velocity */
    for (int y=1;y<=Ny;y++){ 
        int g=site2index(0,y); 
        ux[g]=U_in;
        uy[g]=0.0; 
        }
    /* inlet gradient 0 */
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y), in1=site2index(Nx,y);
        ux[g]=ux[in1]; 
        uy[g]=uy[in1];
    }

    /* bottom wall no slip */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0), in1=site2index(x,1);
        ux[g]=-ux[in1]; 
        uy[g]=-uy[in1];
    }

    /* top wall no slip */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1), in1=site2index(x,Ny);
        ux[g]=-ux[in1];
        uy[g]=-uy[in1];
    }
}

/* No-slip at the fluid/solid interface */
void apply_solid_bc_v(double *ux, double *uy)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (!is_solid[i]) continue; /* only act on solid cells */
            int nb[4]={E[i],W[i],N[i],S[i]};
            double uxs=0.0, uys=0.0; int n=0;
            /* average out velocities of neighboring fluid cell(s) and create the appropriate BC */
            for (int k=0;k<4;k++){
                int j=nb[k];
                if (!is_solid[j]){ uxs+=-ux[j]; uys+=-uy[j]; n++; }
            }
            if (n>0){ ux[i]=uxs/n; uy[i]=uys/n; }
            else     { ux[i]=0.0;  uy[i]=0.0; }
        }
    }
}

/* Pressure ghost-cell BC: zero-gradient at inlet/walls, Dirichlet at the outlet. */
void apply_pressure_bc(void)
{
    /* zero gradient at inlet*/
    for (int y=1;y<=Ny;y++) p[site2index(0,y)]=p[site2index(1,y)];
    /* P out*/
    for (int y=1;y<=Ny;y++) p[site2index(Nx+1,y)]=P_out;
    /* zero gradient at bottom wall*/
    for (int x=1;x<=Nx;x++) p[site2index(x,0)]=p[site2index(x,1)];
    /* zero gradient at top wall*/
    for (int x=1;x<=Nx;x++) p[site2index(x,Ny+1)]=p[site2index(x,Ny)];
}

/* Zero-gradient pressure at the fluid/solid interface */
void apply_solid_bc_p(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (!is_solid[i]) continue;
            int nb[4]={E[i],W[i],N[i],S[i]};
            double ps_=0.0; int n=0;
            for (int k=0;k<4;k++){
                int j=nb[k];
                if (!is_solid[j]){ ps_+=p[j]; n++; }
            }
            if (n>0) p[i]=ps_/n;
        }
    }
}

/* Zero-gradient viscosity ghost-cell BC */
void apply_eta_bc(void)
{
    for (int y=1;y<=Ny;y++){ eta[site2index(0,y)]=eta[site2index(1,y)]; eta[site2index(Nx+1,y)]=eta[site2index(Nx,y)]; }
    for (int x=1;x<=Nx;x++){ eta[site2index(x,0)]=eta[site2index(x,1)]; eta[site2index(x,Ny+1)]=eta[site2index(x,Ny)]; }
}

/* Zero-gradient viscosity at the fluid/solid interface. */
void apply_solid_bc_eta(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (!is_solid[i]) continue;
            int nb[4]={E[i],W[i],N[i],S[i]};
            double es=0.0; int n=0;
            for (int k=0;k<4;k++){
                int j=nb[k];
                if (!is_solid[j]){ es+=eta[j]; n++; }
            }
            if (n>0) eta[i]=es/n;
        }
    }
}

/* Build the grid */
void setup(void)
{
    geometry();
    init_control_points();
    init_solid_mask();
    init_fields();
    apply_velocity_bc(vx,vy);
    apply_solid_bc_v(vx,vy);
    apply_pressure_bc();
    apply_solid_bc_p();
    apply_eta_bc();
    apply_solid_bc_eta();
}

/* Explicit-diffusion stability check for the chosen (K, n, epsilon, dt, dx)
   the viscous term is only stable explicitly if Fo = eta_max*dt/dx^2 stays below 0.25,
   and eta_max = K*epsilon^(n-1) is known up front. */
void check_stability(void)
{
    double eta_zero_shear = POWERLAW_K * pow(POWERLAW_EPSILON, POWERLAW_N - 1.0);
    double h2 = dx[site2index(1,1)] * dy[site2index(1,1)];  
    double Fo = eta_zero_shear * dt / h2;
    double dt_max_stable = 0.25 * h2 / eta_zero_shear;

    printf("-------------------------------------------\n");
    printf("Explicit-diffusion stability check \n");
    printf("  eta_zero_shear (K*eps^(n-1), eta's hard upper bound) = %.4f\n", eta_zero_shear);
    printf("  Fo = eta_zero_shear * dt / dx^2                      = %.4f  (need <= 0.25)\n", Fo);
    if (Fo > 0.25) {
        printf("  *** UNSTABLE: dt=%.6g is %.2fx above the safe limit.\n", dt, Fo/0.25);
        printf("  *** Use dt <= %.6g instead (or raise POWERLAW_EPSILON / lower K) before trusting this run.\n",
               dt_max_stable);
    } 
    else {
        printf("  OK: dt=%.6g is within the safe limit (dt_max_stable=%.6g, margin %.1fx).\n",
               dt, dt_max_stable, dt_max_stable/dt);
    }
    printf("-------------------------------------------\n");
}

/* Average a velocity field (ux,uy) with its E/W/N/S neighbors. */
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

/* Average pressure with its E/W/N/S neighbors.*/
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

/* Power-law closure (skipped for solid cells): 
    build the strain-rate tensor D from the current face-interpolated velocity
    Use its magnitude gamma_dot = sqrt(2 D:D)
    Regularizations with POWERLAW_EPSILON */
void compute_viscosity(double *uxe, double *uxw, double *uxn, double *uxs,
                        double *uye, double *uyw, double *uyn, double *uys)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (is_solid[i]) continue; /* only act on fluid cells*/
            double Dxx=(uxe[i]-uxw[i])/dx[i];
            double Dyy=(uyn[i]-uys[i])/dy[i];
            double dxvy=(uye[i]-uyw[i])/dx[i];
            double dyvx=(uxn[i]-uxs[i])/dy[i];
            double Dxy=0.5*(dyvx+dxvy);

            /* store it - used on SIMPLE method */
            dxvx_c[i]=Dxx; 
            dyvy_c[i]=Dyy; 
            dxvy_c[i]=dxvy; 
            dyvx_c[i]=dyvx;

            double DdotD=Dxx*Dxx+Dyy*Dyy+2.0*Dxy*Dxy;
            double gamma_dot=sqrt(2.0*DdotD);
            double Deff=sqrt(gamma_dot*gamma_dot+POWERLAW_EPSILON*POWERLAW_EPSILON);
            shear_mag[i]=gamma_dot;
            eta[i]=POWERLAW_K*pow(Deff,POWERLAW_N-1.0);
        }
    }

    /* ensure BC are kept */
    apply_eta_bc();
    apply_solid_bc_eta();
}

/* SIMPLE method */
void predict_velocity_nonnewtonian(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            /* set solid cell velocity to 0 */
            if (is_solid[i]){ vx_star[i]=0.0; vy_star[i]=0.0; continue; }

            /* central points - central derivaties */
            double dxvx=(vxe[i]-vxw[i])/dx[i];
            double dxvy=(vye[i]-vyw[i])/dx[i];
            double dyvx=(vxn[i]-vxs[i])/dy[i];
            double dyvy=(vyn[i]-vys[i])/dy[i];

            /* face points - central derivaties */
            double dxvxe=2*(vx[E[i]]-vx[i])/(dx[i]+dx[E[i]]);
            double dxvye=2*(vy[E[i]]-vy[i])/(dx[i]+dx[E[i]]);
            double dxvxw=2*(-vx[W[i]]+vx[i])/(dx[i]+dx[W[i]]);
            double dxvyw=2*(-vy[W[i]]+vy[i])/(dx[i]+dx[W[i]]);
            double dyvxn=2*(vx[N[i]]-vx[i])/(dy[i]+dy[N[i]]);
            double dyvyn=2*(vy[N[i]]-vy[i])/(dy[i]+dy[N[i]]);
            double dyvxs=2*(-vx[S[i]]+vx[i])/(dy[i]+dy[S[i]]);
            double dyvys=2*(-vy[S[i]]+vy[i])/(dy[i]+dy[S[i]]);

            /* eta averaged onto each face from the two cells sharing it */
            double eta_e=0.5*(eta[i]+eta[E[i]]);
            double eta_w=0.5*(eta[i]+eta[W[i]]);
            double eta_n=0.5*(eta[i]+eta[N[i]]);
            double eta_s=0.5*(eta[i]+eta[S[i]]);

            /* off-diagonal strain-rate components averaged onto the same faces */
            double dyvx_e=0.5*(dyvx_c[i]+dyvx_c[E[i]]);
            double dyvx_w=0.5*(dyvx_c[i]+dyvx_c[W[i]]);
            double dxvy_n=0.5*(dxvy_c[i]+dxvy_c[N[i]]);
            double dxvy_s=0.5*(dxvy_c[i]+dxvy_c[S[i]]);

            /* viscous stress tau = 2*eta*D at each face */
            double tau_xx_e=2*eta_e*dxvxe;
            double tau_xx_w=2*eta_w*dxvxw;
            double tau_xy_n=eta_n*(dyvxn+dxvy_n);
            double tau_xy_s=eta_s*(dyvxs+dxvy_s);

            double tau_yy_n=2*eta_n*dyvyn;
            double tau_yy_s=2*eta_s*dyvys;
            double tau_xy_e=eta_e*(dxvye+dyvx_e);
            double tau_xy_w=eta_w*(dxvyw+dyvx_w);

            /* net viscous force = divergence of the stress tensor */
            double visc_x=(tau_xx_e-tau_xx_w)/dx[i]+(tau_xy_n-tau_xy_s)/dy[i];
            double visc_y=(tau_xy_e-tau_xy_w)/dx[i]+(tau_yy_n-tau_yy_s)/dy[i];

            vx_star[i]=vx[i]+dt*(-vx[i]*dxvx-vy[i]*dyvx+visc_x);
            vy_star[i]=vy[i]+dt*(-vx[i]*dxvy-vy[i]*dyvy+visc_y);
        }
    }
}

/* SIMPLE step 2 -- pressure-correction
    b_src = div(v*)/dt, the actually measured divergence of the predicted field. */
void compute_b_src(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (is_solid[i]){ b_src[i]=0.0; continue; }
            double div=(vxe_s[i]-vxw_s[i])/dx[i]+(vyn_s[i]-vys_s[i])/dy[i];
            b_src[i]=div/dt;
        }
    }
}

/* SIMPLE step 3 -- solve Lap(p) = b_src by SOR, 
skipped for solid cells. */
int solve_pressure_sor(double inner_tol, int max_inner, double *resid_out)
{
    int k; 
    double resid=1e300;
    for (k=0;k<max_inner;k++){
        resid=0.0;
        for (int x=1;x<=Nx;x++){
            for (int y=1;y<=Ny;y++){
                int i=site2index(x,y);
                if (is_solid[i]) continue; /* skip solid */
                
                /* pressure at the faces*/
                double pe_i=(p[E[i]]+p[i])/2;
                double pw_i=(p[W[i]]+p[i])/2;
                double pn_i=(p[N[i]]+p[i])/2;
                double ps_i=(p[S[i]]+p[i])/2;

                double val=-1.0*b_src[i]*dx[i]*dx[i]*dy[i]*dy[i]
                           +(pn_i+ps_i)*dx[i]*dx[i]+(pe_i+pw_i)*dy[i]*dy[i];
                
                /* p - gauss seidel - 5-point Laplacian using the face-averaged pressures*/
                double p_gs=(0.5*val)/(dx[i]*dx[i]+dy[i]*dy[i]);
                double p_old=p[i];
                p[i]=(1.0-SOR_OMEGA)*p_old+SOR_OMEGA*p_gs;
                double d=fabs(p[i]-p_old);
                if (d>resid) resid=d;
            }
        }
        /* ensure bc conditions*/
        apply_pressure_bc();
        apply_solid_bc_p();
        if (resid<inner_tol){ k++; break; }
    }
    *resid_out=resid;
    return k;
}

/* SIMPLE step 4 - velocity correction based on pressure
v = v* - dt*grad(p); 
solid cells are pinned to zero. */
void correct_velocity(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (is_solid[i]){ vx[i]=0.0; vy[i]=0.0; continue; }
            vx[i]=vx_star[i]-dt*(pe[i]-pw[i])/dx[i];
            vy[i]=vy_star[i]-dt*(pn[i]-ps[i])/dy[i];
        }
    }
}

/* L2 norm of the velocity field. */
double measure_velocity_l2(void)
{
    double s=0.0;
    for (int x=1;x<=Nx;x++) {
        for (int y=1;y<=Ny;y++){ 
            int i=site2index(x,y); s+=vx[i]*vx[i]+vy[i]*vy[i]; 
        }
    }
    return sqrt(s);
}

/* Volumetric flow rate through x (fluid cells only). 
Since this study is just a flow in x direction, this is enoguh to check if it is preserved.*/

double flow_rate_at(int x)
{
    double Q=0.0;
    for (int y=1;y<=Ny;y++){ 
        int i=site2index(x,y); 
        if (!is_solid[i]) Q+=vx[i]*dy[i]; }
    return Q;
}

/* Maximum cell-wise |div(v)| over the fluid interior. */
double max_abs_divergence(void)
{
    double m=0.0;
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (is_solid[i]) continue;
            double vxe_=(vx[E[i]]+vx[i])/2, vxw_=(vx[W[i]]+vx[i])/2;
            double vyn_=(vy[N[i]]+vy[i])/2, vys_=(vy[S[i]]+vy[i])/2;
            double div=(vxe_-vxw_)/dx[i]+(vyn_-vys_)/dy[i];
            if (fabs(div)>m) m=fabs(div);
        }
    }
    return m;
}

/* Write Q(x) (flow rate per column) to a text file for plotting. */
void write_flow_rate_profile(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    fprintf(F,"# x  Q(x)\n");
    for (int x=1;x<=Nx;x++) fprintf(F,"%d %g\n",x-1,flow_rate_at(x));
    fclose(F);
}

/* Write the local shear rate, apparent viscosity and solid mask at
   every cell, for plotting the shear-thinning profile across the
   constriction. */
void write_shear_eta(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    fprintf(F,"# x  y  gamma_dot  eta_a  is_solid\n");
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %g %g %d\n",x-1,y-1,shear_mag[i],eta[i],is_solid[i]);
        }
    }
    fclose(F);
}

/* Write the 2D pressure field, plus the solid mask, for plotting. */
void write_p_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){ int i=site2index(x,y); fprintf(F,"%d %d %g %d\n",x-1,y-1,p[i],is_solid[i]); }
        fprintf(F,"\n");
    }
    fclose(F);
}

/* Write the 2D velocity field, plus the solid mask, for plotting. */
void write_v_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){ int i=site2index(x,y); fprintf(F,"%d %d %e %e %d\n",x-1,y-1,vx[i],vy[i],is_solid[i]); }
        fprintf(F,"\n");
    }
    fclose(F);
}

/* Main SIMPLE (projection) time-stepping loop, non-Newtonian + solid
   mask. This is the union of finite_vol_4_syringe_SIMPLE.c's solver
   (SIMPLE + solid mask) and
   finite_vol_5_duct_SIMPLE_shearthinning.c's solver (SIMPLE +
   power-law viscosity):

     1) refresh velocity BC + solid BC on (vx,vy), face-interpolate;
     2) compute_viscosity(): rebuild eta(x,y) from the current velocity
        field via the power law (fluid cells only);
     3) predict_velocity_nonnewtonian(): momentum-only step using that
        eta -> vx_star, vy_star (solid cells pinned to 0);
     4) BC + solid BC + face-interpolate the predictor;
     5) compute_b_src(): measured div(v*)/dt (fluid cells only);
     6) solve_pressure_sor(): solve for the pressure field, with the
        solid-cell pressure BC re-applied every sweep;
     7) interp_p_faces() on the converged pressure;
     8) correct_velocity(): project v* back onto div(v)=0 (solid cells
        pinned to 0);
     9) refresh velocity BC + solid BC on the corrected (vx,vy);
     10) convergence check (only trusted after MIN_TIME of simulated
         time); log |v|, SOR residual, max|div(v)|, the flow rate at
         four stations (inlet, pre-constriction, throat, outlet) and
         the eta range every 2000 steps.
*/
void solver(void)
{
    double tolerance=1e-6;
    int it=0, max_it=2000000;
    double vnorm_old=0.0, check=1.0;
    /* Require at least MIN_TIME of physical (simulation) time before the check above is
       trusted. 
       This is due to validation of results where not accurate.
       I noticed that it needs more time to reach steady-flow, even though convergence could be considered reached...
       A very small dt it can look "converged" simply because one step is a tiny slice of physical 
       at the same time very small dt is needed to allow precise carrying of information. */
    double MIN_TIME=10.0;
    double inner_tol=1e-6;
    int max_inner=1000;

    FILE *L=fopen("finite_vol_6_syringe_SIMPLE_shearthinning_convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence log");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  mdiv  Q_in  Q_pre  Q_throat  Q_out  eta_min  eta_max\n");

    do{
        apply_velocity_bc(vx,vy);
        apply_solid_bc_v(vx,vy);
        interp_v_faces(vx,vy,vxe,vxw,vxn,vxs,vye,vyw,vyn,vys);

        compute_viscosity(vxe,vxw,vxn,vxs,vye,vyw,vyn,vys);

        predict_velocity_nonnewtonian();
        apply_velocity_bc(vx_star,vy_star);
        apply_solid_bc_v(vx_star,vy_star);
        interp_v_faces(vx_star,vy_star,vxe_s,vxw_s,vxn_s,vxs_s,vye_s,vyw_s,vyn_s,vys_s);

        compute_b_src();
        double sor_resid;
        int sor_iters=solve_pressure_sor(inner_tol,max_inner,&sor_resid);
        interp_p_faces();

        correct_velocity();
        apply_velocity_bc(vx,vy);
        apply_solid_bc_v(vx,vy);

        double vnorm=measure_velocity_l2();
        if (it>0 && vnorm_old>1e-12) check=fabs(vnorm-vnorm_old)/vnorm_old;
        vnorm_old=vnorm;

        if (it%2000==0 || it==max_it-1){
            double mdiv=max_abs_divergence();
            double q_in=flow_rate_at(1), q_pre=flow_rate_at(CONTRACTION_X-1),
                   q_throat=flow_rate_at(CONTRACTION_X+1), q_out=flow_rate_at(Nx);
            double eta_min=1e300, eta_max=-1e300;
            for (int x=1;x<=Nx;x++) for (int y=1;y<=Ny;y++){
                int i=site2index(x,y);
                if (is_solid[i]) continue;
                double e=eta[i];
                if (e<eta_min) eta_min=e;
                if (e>eta_max) eta_max=e;
            }
            fprintf(L,"%d %g %g %d %g %g %g %g %g %g %g %g\n",
                    it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_pre,q_throat,q_out,eta_min,eta_max);
            printf("it %d  |v| %.6e  check %.3e  (SOR %d, resid %.2e)  mdiv %.3e  Q(in/pre/throat/out) %.4f/%.4f/%.4f/%.4f  eta[%.3f,%.3f]\n",
                   it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_pre,q_throat,q_out,eta_min,eta_max);
        }
        it++;
    } while ((check>tolerance || (double)it*dt<MIN_TIME) && it<max_it);

    printf("Stopped after %d iterations, check=%.3e (converged=%s)\n",
           it,check,(check<=tolerance)?"yes":"NO -- hit iteration cap");
    fclose(L);
}

/* Set up the grid + constriction, print the stability check, dump the
   starting fields, run to steady state, then dump the final fields,
   the flow rate profile and the shear-rate/viscosity map. */
int main(void)
{
    setup();
    check_stability();
    write_p_2d("finite_vol_6_syringe_SIMPLE_shearthinning_start_pressure.dat");
    write_v_2d("finite_vol_6_syringe_SIMPLE_shearthinning_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_6_syringe_SIMPLE_shearthinning_final_pressure.dat");
    write_v_2d("finite_vol_6_syringe_SIMPLE_shearthinning_final_velocity.dat");
    write_flow_rate_profile("finite_vol_6_syringe_SIMPLE_shearthinning_flow_rate_profile.dat");
    write_shear_eta("finite_vol_6_syringe_SIMPLE_shearthinning_shear_eta.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}