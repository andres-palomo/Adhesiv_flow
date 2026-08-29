/* =====================================================================
   MILESTONE 3 -- the syringe geometry: a symmetric strangulation

   Built forward from finite_vol_2_duct.c. Same ghost-cell duct, same
   Newtonian constant viscosity, same *legacy* pressure-Poisson formula
   (still assuming div(v)=0 rather than measuring it). The only change
   is geometric: for x >= CONTRACTION_X the duct is symmetrically
   strangled down to half its height, like the barrel of a syringe
   necking into its needle. That is done with a solid-cell mask: cells
   in the blocked band are excluded from the momentum and pressure
   updates and held at v=0, with a mirrored no-slip treatment (the same
   trick as the outer wall ghost cells, extended to whichever interior
   faces border the fluid) at the interface -- see apply_solid_bc().

   This file is where the legacy pressure formula's mass-conservation
   problem, already present on the plain duct (finite_vol_2_duct.c), is
   supposed to become visible. It does: with a real area change forcing
   the issue, Q(x) is not just off by a roughly constant amount, it
   drops sharply right at the constriction and never recovers, which is
   the numbers this file prints and logs. Nothing here fixes it -- that
   is finite_vol_4_seringe_SIMPLE.c's job.

   Build:  gcc -O2 -o finite_vol_3_seringe finite_vol_3_seringe.c -lm
   Run:    ./finite_vol_3_seringe
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

/* Constriction geometry: for x >= CONTRACTION_X, only the band
   [THROAT_Y_LO, THROAT_Y_HI] (half the channel height, centered) stays
   open; everything else in that band is solid. */
#define CONTRACTION_X (Nx/2)
#define THROAT_Y_LO (Ny/4+1)
#define THROAT_Y_HI (3*Ny/4)

double dt = 0.001;
double U_in  = 0.1;
double P_out = 0.0;
double MU_CONST = 1.0;

double p[Ntot], vx[Ntot], vy[Ntot];
int    E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];
double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];
int    is_solid[Ntot];

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

/* Solid-cell mask for the syringe constriction. A physical interior
   cell (x,y) in [1,Nx]x[1,Ny] is solid if x is past the contraction
   start and y falls outside the throat band. Ghost cells are never
   marked solid -- the outer-wall/inlet/outlet BCs already own those. */
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
    for (int i=0;i<Ntot;i++){ p[i]=0.0; vx[i]=0.0; vy[i]=0.0; }
}

void apply_boundary_conditions(void)
{
    for (int y=1;y<=Ny;y++){
        int g=site2index(0,y);
        vx[g]=U_in; vy[g]=0.0;
        p[g]=p[site2index(1,y)];
    }
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y);
        int in1=site2index(Nx,y);
        vx[g]=vx[in1]; vy[g]=vy[in1];
        p[g]=P_out;
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0);
        int in1=site2index(x,1);
        vx[g]=-vx[in1]; vy[g]=-vy[in1];
        p[g]=p[in1];
    }
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1);
        int in1=site2index(x,Ny);
        vx[g]=-vx[in1]; vy[g]=-vy[in1];
        p[g]=p[in1];
    }
}

/* No-slip at the interior faces where fluid meets the solid mask,
   extending the same ghost-mirroring idea used at the outer walls:
   a solid cell's velocity is set to minus the (average of the) mirror
   of its fluid neighbor(s), so the face-averaged velocity between a
   fluid cell and this solid cell is zero. Pressure is carried across
   with a zero-gradient average instead. Solid cells with no fluid
   neighbor (deep inside the blocked band) are never read by anything
   that matters, so they are simply pinned to zero. */
void apply_solid_bc(void)
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            if (!is_solid[i]) continue;
            int nb[4]={E[i],W[i],N[i],S[i]};
            double vxs=0.0, vys=0.0, ps=0.0; int n=0;
            for (int k=0;k<4;k++){
                int j=nb[k];
                if (!is_solid[j]){ vxs+=-vx[j]; vys+=-vy[j]; ps+=p[j]; n++; }
            }
            if (n>0){ vx[i]=vxs/n; vy[i]=vys/n; p[i]=ps/n; }
            else     { vx[i]=0.0;  vy[i]=0.0; }
        }
    }
}

void setup(void)
{
    geometry();
    init_control_points();
    init_solid_mask();
    init_fields();
    apply_boundary_conditions();
    apply_solid_bc();
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
        apply_boundary_conditions();
        apply_solid_bc();
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

    FILE *L=fopen("finite_vol_3_seringe_convergence_log.dat","w");
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
    write_p_2d("finite_vol_3_seringe_start_pressure.dat");
    write_v_2d("finite_vol_3_seringe_start_velocity.dat");
    solver();
    write_p_2d("finite_vol_3_seringe_final_pressure.dat");
    write_v_2d("finite_vol_3_seringe_final_velocity.dat");
    write_flow_rate_profile("finite_vol_3_seringe_flow_rate_profile.dat");

    double q_in=flow_rate_at(1), q_out=flow_rate_at(Nx);
    printf("Final flow rate: Q(inlet)=%.6f  Q(outlet)=%.6f  (expected %.6f = U_in*Ny)\n",
           q_in,q_out,U_in*Ny);
    printf("Final max|div(v)| = %.6e\n",max_abs_divergence());
    return 0;
}
