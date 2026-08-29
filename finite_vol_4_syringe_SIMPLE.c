/* =====================================================================
   Syringe geometry, Newtonian viscosity, SIMPLE (projection) method.

   Build:  gcc -O2 -o finite_vol_4_syringe_SIMPLE finite_vol_4_syringe_SIMPLE.c -lm
   Run:    ./finite_vol_4_syringe_SIMPLE
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

double dt = 0.001;
double U_in  = 0.1;
double P_out = 0.0;
double MU_CONST = 1.0;  /* constant Newtonian viscosity */

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

double b_src[Ntot];      /* pressure-Poisson right-hand side: div(v*)/dt, actually measured */
double SOR_OMEGA = 1.7;  /* omega of SOR method */

/* grid coordinate to latice (1D array index). */
int site2index(int x, int y) { return x + y*NXG; }

/* cell size. */
void geometry(void)
{
    for (int i=0;i<Ntot;i++) dx[i]=dy[i]=1.0;
}

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

/* Start from rest, zero pressure (both the real field and the
   predictor v_star). */
void init_fields(void)
{
    for (int i=0;i<Ntot;i++){ p[i]=0.0; vx[i]=0.0; vy[i]=0.0; vx_star[i]=0.0; vy_star[i]=0.0; }
}

/* Velocity at ghost-cell BC, by using mirrored (no-slip) top/bottom
   walls. Called for both (vx,vy) and (vx_star,vy_star). */
void apply_velocity_bc(double *ux, double *uy)
{
    /* inlet velocity */
    for (int y=1;y<=Ny;y++){
        int g=site2index(0,y);
        ux[g]=U_in; uy[g]=0.0;
    }
    /* inlet gradient 0 */
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y);
        int in1=site2index(Nx,y);
        ux[g]=ux[in1]; uy[g]=uy[in1];
    }
    /* bottom wall no slip */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0);
        int in1=site2index(x,1);
        ux[g]=-ux[in1]; uy[g]=-uy[in1];
    }
    /* top wall no slip */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1);
        int in1=site2index(x,Ny);
        ux[g]=-ux[in1]; uy[g]=-uy[in1];
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

/* SIMPLE method */
void predict_velocity(void)
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

            double visc_x = MU_CONST*((dxvxe-dxvxw)/dx[i] + (dyvxn-dyvxs)/dy[i]);
            double visc_y = MU_CONST*((dxvye-dxvyw)/dx[i] + (dyvyn-dyvys)/dy[i]);

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
    for (int x=1;x<=Nx;x++)
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            s+=vx[i]*vx[i]+vy[i]*vy[i];
        }
    return sqrt(s);
}

/* Volumetric flow rate through x (fluid cells only). Since this study
   is just a flow in x direction, this is enoguh to check if it is
   preserved. */
double flow_rate_at(int x)
{
    double Q=0.0;
    for (int y=1;y<=Ny;y++){
        int i=site2index(x,y);
        if (is_solid[i]) continue;
        Q+=vx[i]*dy[i];
    }
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

/* Write the 2D pressure field, plus the solid mask, for plotting. */
void write_p_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %g %d\n",x-1,y-1,p[i],is_solid[i]);
        }
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
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %e %e %d\n",x-1,y-1,vx[i],vy[i],is_solid[i]);
        }
        fprintf(F,"\n");
    }
    fclose(F);
}

/* Main SIMPLE (projection) time-stepping loop:

     1) refresh velocity BC + solid BC on (vx,vy), then face-
        interpolate them;
     2) predict_velocity(): explicit momentum step without pressure ->
        vx_star, vy_star;
     3) apply the same BCs to the predictor and face-interpolate it too
        (needed to measure its divergence in the next step);
     4) compute_b_src(): b_src = div(v*)/dt, the real divergence the
        pressure correction has to cancel;
     5) solve_pressure_sor(): solve Lap(p) = b_src for the pressure
        field, to convergence;
     6) interp_p_faces(): face-interpolate the converged pressure;
     7) correct_velocity(): v = v* - dt*grad(p) -- this is the
        projection step that makes v (approximately) divergence-free;
     8) refresh velocity BC + solid BC on the corrected (vx,vy);
     9) check the L2-norm convergence criterion (only trusted once at
        least MIN_TIME of simulated time has elapsed -- see comment on
        MIN_TIME below); every 2000 steps log |v|, the SOR residual,
        max|div(v)| and the flow rate at four stations so mass
        conservation can be tracked across the constriction.
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

    FILE *L=fopen("finite_vol_4_syringe_SIMPLE_convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence log");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  mdiv  Q_in  Q_pre  Q_throat  Q_out\n");

    do{
        apply_velocity_bc(vx,vy);
        apply_solid_bc_v(vx,vy);
        interp_v_faces(vx,vy,vxe,vxw,vxn,vxs,vye,vyw,vyn,vys);

        predict_velocity();
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
            fprintf(L,"%d %g %g %d %g %g %g %g %g %g\n",
                    it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_pre,q_throat,q_out);
            printf("it %d  |v| %.6e  check %.3e  (SOR %d, resid %.2e)  mdiv %.3e  Q(in/pre/throat/out) %.4f/%.4f/%.4f/%.4f\n",
                   it,vnorm,check,sor_iters,sor_resid,mdiv,q_in,q_pre,q_throat,q_out);
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
    write_p_2d("finite_vol_4_syringe_SIMPLE_start_pressure.dat");
    write_v_2d("finite_vol_4_syringe_SIMPLE_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_4_syringe_SIMPLE_final_pressure.dat");
    write_v_2d("finite_vol_4_syringe_SIMPLE_final_velocity.dat");
    write_flow_rate_profile("finite_vol_4_syringe_SIMPLE_flow_rate_profile.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}
