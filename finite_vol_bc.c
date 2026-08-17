#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <err.h>
#include <time.h>

// ---------------------------------------------------------------------
// Milestone 2: inlet / outlet / no-slip wall boundary conditions
//
// Change from the original: the domain now has a ghost-cell ring around
// the Nx x Ny physical interior. Total grid is (Nx+2) x (Ny+2).
//   x = 0        -> inlet ghost column
//   x = Nx+1      -> outlet ghost column
//   y = 0        -> bottom wall ghost row
//   y = Ny+1      -> top wall ghost row
//   x in [1,Nx], y in [1,Ny] -> physical interior cells (the ones we update)
//
// site2index no longer wraps -- indices are computed directly, and the
// ghost ring is what keeps E/W/N/S lookups valid at the domain edges.
// ---------------------------------------------------------------------

#define Nx 10
#define Ny 10
#define NXG (Nx+2)
#define NYG (Ny+2)
#define Ntot (NXG*NYG)

double dt=0.001;
double U_in=0.1;     // inlet velocity (Dirichlet)
double P_out=0.0;    // outlet reference pressure (Dirichlet)

double p[Ntot], vx[Ntot], vy[Ntot];
int E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];
double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];

int site2index(int x,int y)
{
    return x + y*NXG;   // no wraparound -- x,y must be in [0,NXG-1]x[0,NYG-1]
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

void init_fields()
{
    for (int i=0;i<Ntot;i++){ p[i]=0.0; vx[i]=0.0; vy[i]=0.0; }
}

// Enforce boundary conditions on the ghost ring. Called once before the
// loop starts and again after every velocity update, since the ghost
// values depend on the (changing) interior values next to them.
void apply_boundary_conditions()
{
    // Inlet (x=0): Dirichlet on velocity, vx = U_in, vy = 0
    for (int y=1;y<=Ny;y++){
        int g=site2index(0,y);
        vx[g]=U_in;
        vy[g]=0.0;
        // zero-gradient pressure at inlet
        int in1=site2index(1,y);
        p[g]=p[in1];
    }

    // Outlet (x=Nx+1): zero-gradient velocity, Dirichlet pressure
    for (int y=1;y<=Ny;y++){
        int g=site2index(Nx+1,y);
        int in1=site2index(Nx,y);
        vx[g]=vx[in1];
        vy[g]=vy[in1];
        p[g]=P_out;
    }

    // Bottom wall (y=0): no-slip via mirrored ghost value
    // (face average between ghost and first interior row = 0)
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,0);
        int in1=site2index(x,1);
        vx[g]=-vx[in1];
        vy[g]=-vy[in1];
        p[g]=p[in1];
    }

    // Top wall (y=Ny+1): no-slip via mirrored ghost value
    for (int x=1;x<=Nx;x++){
        int g=site2index(x,Ny+1);
        int in1=site2index(x,Ny);
        vx[g]=-vx[in1];
        vy[g]=-vy[in1];
        p[g]=p[in1];
    }
}

void setup()
{
    geometry();
    init_control_points();
    init_fields();
    apply_boundary_conditions();
}

void interpolate_faces()
{
    // only interior cells need their faces interpolated
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

void update_velocity()
{
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            double dxvx=(vxe[i]-vxw[i])/dx[i];
            double dxvy=(vye[i]-vyw[i])/dx[i];   // fixed: was vxw[i] in the original
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
            vx[i]+=dt*(-vx[i]*dxvx-vy[i]*dyvx-(pe[i]-pw[i])/dx[i]+(dxvxe-dxvxw)/dx[i]+(dyvxn-dyvxs)/dx[i]);
            vy[i]+=dt*(-vx[i]*dxvy-vy[i]*dyvy-(pn[i]-ps[i])/dy[i]+(dxvye-dxvyw)/dx[i]+(dyvyn-dyvys)/dx[i]);
        }
    }
}

double b_src[Ntot];   // pressure Poisson source term, frozen for the duration of the inner solve
double SOR_OMEGA=1.7; // relaxation factor: 1.0 = plain Gauss-Seidel, (1,2) = SOR

// Solve the pressure Poisson equation to convergence via SOR/Gauss-Seidel,
// holding the velocity field (and hence b_src) fixed.
// Returns the number of inner sweeps used; writes final residual to *resid_out.
int solve_pressure_sor(double inner_tol,int max_inner,double *resid_out)
{
    // b_src computed once from the current (frozen) velocity field
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            double dxvx=(vxe[i]-vxw[i])/dx[i];
            double dxvy=(vye[i]-vyw[i])/dx[i];   // fixed: was vxw[i] in the original
            double dyvx=(vxn[i]-vxs[i])/dy[i];
            double dyvy=(vyn[i]-vys[i])/dy[i];
            b_src[i]=-dxvx*dxvx-dyvy*dyvy+2*dxvy*dyvx;
        }
    }

    int k;
    double resid=1e300;
    for (k=0;k<max_inner;k++){
        resid=0.0;
        for (int x=1;x<=Nx;x++){
            for (int y=1;y<=Ny;y++){
                int i=site2index(x,y);
                // face pressures use the latest p values available -- this is what
                // makes it Gauss-Seidel rather than Jacobi (sweeping x then y,
                // west/south neighbours in this sweep are already updated)
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
        // ghost pressure values feed into pe/pw/pn/ps for boundary-adjacent
        // interior cells, so refresh them every sweep
        apply_boundary_conditions();
        if (resid<inner_tol) { k++; break; }
    }
    *resid_out=resid;
    return k;
}

double measure_velocity_l2()
{
    double s=0.0;
    for (int x=1;x<=Nx;x++)
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            s+=vx[i]*vx[i]+vy[i]*vy[i];
        }
    return sqrt(s);
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

void solver()
{
    double tolerance=1e-6;
    int it=0, max_it=20000;
    double vnorm_old=0.0, check=1.0;

    double inner_tol=1e-8;
    int max_inner=2000;

    FILE *O=fopen("v_t.dat","w");
    if(!O) err(1,"Could not create v_t.dat");
    // report-ready log: outer iteration, |v|, check, SOR sweeps used, SOR residual
    FILE *L=fopen("convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence_log.dat");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual\n");

    do{
        apply_boundary_conditions();
        interpolate_faces();               // freezes velocity faces for this step's b_src

        double sor_resid;
        int sor_iters=solve_pressure_sor(inner_tol,max_inner,&sor_resid);

        interpolate_faces();               // refresh pe/pw/pn/ps from the converged pressure
        update_velocity();
        apply_boundary_conditions();

        double vnorm=measure_velocity_l2();
        fprintf(O,"%d %g\n",it,vnorm);
        if (it>0 && vnorm_old>1e-12) check=fabs(vnorm-vnorm_old)/vnorm_old;
        fprintf(L,"%d %g %g %d %g\n",it,vnorm,check,sor_iters,sor_resid);
        vnorm_old=vnorm;
        it++;
        if (it%500==0)
            printf("it %d  |v| %.6e  check %.3e  (SOR: %d sweeps, resid %.2e)\n",
                   it,vnorm,check,sor_iters,sor_resid);
    } while(check>tolerance && it<max_it);

    printf("Stopped after %d iterations, check=%.3e\n",it,check);
    fclose(O);
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
    return 0;
}
