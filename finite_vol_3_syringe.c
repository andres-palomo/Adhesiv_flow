/* =====================================================================
   Syringe geometry (symmetric constriction), Newtonian viscosity,
   legacy (non-projection) pressure solve.

   Build:  gcc -O2 -o finite_vol_3_syringe finite_vol_3_syringe.c -lm
   Run:    ./finite_vol_3_syringe
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

double p[Ntot], vx[Ntot], vy[Ntot];          /* pressure and velocity in each cell */
int    E[Ntot], W[Ntot], N[Ntot], S[Ntot];   /* neighbor index of each cell */
double dx[Ntot], dy[Ntot];                   /* cell size (uniform, = 1.0 everywhere) */
int    is_solid[Ntot];                       /* 1 = cell is a wall,
                                                 0 = cell is a fluid */

double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];       /* p interpolated to the E/W/N/S faces */
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];   /* vx interpolated to the E/W/N/S faces */
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];   /* vy interpolated to the E/W/N/S faces */

double b_src[Ntot];      /* pressure-Poisson right-hand side, assuming div(v)=0 already holds */
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

/* Start from rest, zero pressure. */
void init_fields(void)
{
    for (int i=0;i<Ntot;i++){ p[i]=0.0; vx[i]=0.0; vy[i]=0.0; }
}

/* Ghost-cell BC on velocity AND pressure together (this file solves
   both with one pressure formula, so there's no separate velocity/
   pressure BC split like the SIMPLE files use). */
void apply_boundary_conditions(void)
{
    /* inlet velocity, zero gradient pressure */
    for (int y=1;y<=Ny;y++){
        int g=site2index(0,y);
        vx[g]=U_in; vy[g]=0.0;
        p[g]=p[site2index(1,y)];
    }
    /* outlet gradient 0, pressure fixed at P_out */
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y);
        int in1=site2index(Nx,y);
        vx[g]=vx[in1]; vy[g]=vy[in1];
        p[g]=P_out;
    }
    /* bottom wall no slip */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0);
        int in1=site2index(x,1);
        vx[g]=-vx[in1]; vy[g]=-vy[in1];
        p[g]=p[in1];
    }
    /* top wall no slip */
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1);
        int in1=site2index(x,Ny);
        vx[g]=-vx[in1]; vy[g]=-vy[in1];
        p[g]=p[in1];
    }
}

/* No-slip at the fluid/solid interface (velocity AND pressure
   together, since this file doesn't split them). */
void apply_solid_bc(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (!is_solid[i]) continue; /* only act on solid cells */
            int nb[4]={E[i],W[i],N[i],S[i]};
            double vxs=0.0, vys=0.0, ps=0.0; int n=0;
            /* average out velocities of neighboring fluid cell(s) and create the appropriate BC */
            for (int k=0;k<4;k++){
                int j=nb[k];
                if (!is_solid[j]){ vxs+=-vx[j]; vys+=-vy[j]; ps+=p[j]; n++; }
            }
            if (n>0){ vx[i]=vxs/n; vy[i]=vys/n; p[i]=ps/n; }
            else     { vx[i]=0.0;  vy[i]=0.0; }
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
    apply_boundary_conditions();
    apply_solid_bc();
}

/* Average each cell with its E/W/N/S neighbors. */
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

/* Explicit velocity update: advection + pressure gradient + Newtonian
   viscous diffusion, all in one step; skipped for solid cells. */
void update_velocity(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (is_solid[i]) continue;
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

/* Pressure Poisson solve, legacy formula (assumes div(v)=0), skipped
   for solid cells. */
int solve_pressure_sor(double inner_tol, int max_inner, double *resid_out)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (is_solid[i]) { b_src[i]=0.0; continue; }
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
        apply_boundary_conditions();
        apply_solid_bc();
        if (resid<inner_tol){ k++; break; }
    }
    *resid_out=resid;
    return k;
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

/* Main time-stepping loop -- same legacy (non-projection) scheme as
   finite_vol_2_duct.c, with the solid mask applied at every step:

     1) refresh ghost cells + solid-cell BC, then face-interpolated
        values, for the current (p, vx, vy);
     2) solve the pressure Poisson equation to convergence (SOR);
     3) re-interpolate p to the faces;
     4) explicit update of vx, vy (already includes the pressure
        gradient from step 3);
     5) refresh ghost + solid BC again;
     6) check convergence; every 2000 steps log |v|, the SOR residual,
        max|div(v)| and the flow rate at four stations (inlet, just
        before the constriction, just into the throat, outlet) to
        track how mass conservation behaves across the area change.
*/
void solver(void)
{
    double tolerance=1e-6;
    int it=0, max_it=100000;
    double vnorm_old=0.0, check=1.0;

    double inner_tol=1e-6;
    int max_inner=500;

    FILE *L=fopen("finite_vol_3_syringe_convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence log");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual  mdiv  Q_in  Q_pre  Q_throat  Q_out\n");

    do{
        apply_boundary_conditions();
        apply_solid_bc();
        interpolate_faces();

        double sor_resid;
        int sor_iters=solve_pressure_sor(inner_tol,max_inner,&sor_resid);

        interpolate_faces();
        update_velocity();
        apply_boundary_conditions();
        apply_solid_bc();

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
    } while (check>tolerance && it<max_it);

    printf("Stopped after %d iterations, check=%.3e (converged=%s)\n",
           it,check,(check<=tolerance)?"yes":"NO -- hit iteration cap");
    fclose(L);
}

int main(void)
{
    setup();
    write_p_2d("finite_vol_3_syringe_start_pressure.dat");
    write_v_2d("finite_vol_3_syringe_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_3_syringe_final_pressure.dat");
    write_v_2d("finite_vol_3_syringe_final_velocity.dat");
    write_flow_rate_profile("finite_vol_3_syringe_flow_rate_profile.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}
