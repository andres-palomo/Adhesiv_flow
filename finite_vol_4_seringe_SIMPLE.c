/* =====================================================================
   MILESTONE 4 -- the SIMPLE method: actually conserving mass

   Built forward from finite_vol_3_seringe.c. Same syringe geometry
   (solid-cell mask, half-height constriction from CONTRACTION_X on),
   same Newtonian constant viscosity. What changes is the pressure
   treatment: the legacy formula (assume div(v)=0, solve a Poisson
   equation that never actually looks at the divergence) is replaced
   by a real projection / SIMPLE-family step:

     1. predict_velocity(): advance velocity with advection + viscous
        diffusion only, no pressure gradient at all, producing an
        intermediate field v* that is generally NOT divergence-free.
     2. compute_b_src(): b_src = div(v*)/dt -- the actual, measured
        divergence of v*, not an assumed one.
     3. solve_pressure_sor(): solve the standard pressure-Poisson
        equation Lap(p) = b_src for the correction pressure.
     4. correct_velocity(): v = v* - dt*grad(p), which drives the
        velocity field back toward div(v) = 0.

   This is tutorial.tex Section 6 (Eq. 26-31) implemented directly.
   finite_vol_3_seringe.c's own numbers are the baseline to beat: there,
   Q(x) fell from ~1.67 near the inlet to ~0.24 near the outlet, and
   max|div(v)| sat at ~3.9e-2 no matter how long the run went. The
   printouts at the end of this file are the same measurements, for
   the same geometry, with the fix applied.

   Build:  gcc -O2 -o finite_vol_4_seringe_SIMPLE finite_vol_4_seringe_SIMPLE.c -lm
   Run:    ./finite_vol_4_seringe_SIMPLE
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

#define CONTRACTION_X (Nx/2)
#define THROAT_Y_LO (Ny/4+1)
#define THROAT_Y_HI (3*Ny/4)

double dt = 0.001;
double U_in  = 0.1;
double P_out = 0.0;
double MU_CONST = 1.0;

double p[Ntot], vx[Ntot], vy[Ntot];
double vx_star[Ntot], vy_star[Ntot];
int    E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];
int    is_solid[Ntot];

/* face-interpolated real velocity and pressure */
double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];
/* face-interpolated predicted (star) velocity, used only to build b_src */
double vxe_s[Ntot], vxw_s[Ntot], vxn_s[Ntot], vxs_s[Ntot];
double vye_s[Ntot], vyw_s[Ntot], vyn_s[Ntot], vys_s[Ntot];

double b_src[Ntot];
double SOR_OMEGA = 1.7;

int site2index(int x, int y) { return x + y*NXG; }

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

void init_fields(void)
{
    for (int i=0;i<Ntot;i++){ p[i]=0.0; vx[i]=0.0; vy[i]=0.0; vx_star[i]=0.0; vy_star[i]=0.0; }
}

/* Velocity-only boundary conditions, parameterized so the exact same
   physical rules (Dirichlet inlet, zero-gradient outlet, no-slip
   walls) apply to the real field and to the intermediate field v*. */
void apply_velocity_bc(double *ux, double *uy)
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
}

/* No-slip mirroring at the solid mask, for whichever velocity field
   (real or star) is passed in -- same idea as finite_vol_3_seringe.c's
   apply_solid_bc(), split so it can be reused for both fields. */
void apply_solid_bc_v(double *ux, double *uy)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (!is_solid[i]) continue;
            int nb[4]={E[i],W[i],N[i],S[i]};
            double uxs=0.0, uys=0.0; int n=0;
            for (int k=0;k<4;k++){
                int j=nb[k];
                if (!is_solid[j]){ uxs+=-ux[j]; uys+=-uy[j]; n++; }
            }
            if (n>0){ ux[i]=uxs/n; uy[i]=uys/n; }
            else     { ux[i]=0.0;  uy[i]=0.0; }
        }
    }
}

/* Pressure-only boundary conditions: zero-gradient at inlet and walls,
   Dirichlet at the outlet -- same rules as before, just isolated from
   the velocity BCs since the projection step handles them separately. */
void apply_pressure_bc(void)
{
    for (int y=1;y<=Ny;y++) p[site2index(0,y)]=p[site2index(1,y)];
    for (int y=1;y<=Ny;y++) p[site2index(Nx+1,y)]=P_out;
    for (int x=1;x<=Nx;x++) p[site2index(x,0)]=p[site2index(x,1)];
    for (int x=1;x<=Nx;x++) p[site2index(x,Ny+1)]=p[site2index(x,Ny)];
}

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

/* Step 1: advance velocity with advection + Newtonian viscous
   diffusion only -- no pressure gradient. Uses the real field's face
   values (vxe..vys), which must already be current. Result goes into
   vx_star/vy_star, leaving vx/vy untouched until correct_velocity(). */
void predict_velocity(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (is_solid[i]){ vx_star[i]=0.0; vy_star[i]=0.0; continue; }
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

/* Step 2: the real pressure-Poisson source, measured directly from
   the (generally non-solenoidal) predicted field v*. */
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

/* Step 3: standard pressure-Poisson solve, same SOR sweep and stencil
   as the earlier files -- only b_src's meaning has changed. */
int solve_pressure_sor(double inner_tol, int max_inner, double *resid_out)
{
    int k; double resid=1e300;
    for (k=0;k<max_inner;k++){
        resid=0.0;
        for (int x=1;x<=Nx;x++){
            for (int y=1;y<=Ny;y++){
                int i=site2index(x,y);
                if (is_solid[i]) continue;
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
        apply_pressure_bc();
        apply_solid_bc_p();
        if (resid<inner_tol){ k++; break; }
    }
    *resid_out=resid;
    return k;
}

/* Step 4: project v* back onto (approximately) the divergence-free
   manifold using the pressure correction just solved for. */
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
            fprintf(F,"%d %d %g %d\n",x-1,y-1,p[i],is_solid[i]);
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
            fprintf(F,"%d %d %e %e %d\n",x-1,y-1,vx[i],vy[i],is_solid[i]);
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

    FILE *L=fopen("finite_vol_4_seringe_SIMPLE_convergence_log.dat","w");
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
    } while (check>tolerance && it<max_it);

    printf("Stopped after %d iterations, check=%.3e (converged=%s)\n",
           it,check,(check<=tolerance)?"yes":"NO -- hit iteration cap");
    fclose(L);
}

int main(void)
{
    setup();
    write_p_2d("finite_vol_4_seringe_SIMPLE_start_pressure.dat");
    write_v_2d("finite_vol_4_seringe_SIMPLE_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_4_seringe_SIMPLE_final_pressure.dat");
    write_v_2d("finite_vol_4_seringe_SIMPLE_final_velocity.dat");
    write_flow_rate_profile("finite_vol_4_seringe_SIMPLE_flow_rate_profile.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}
