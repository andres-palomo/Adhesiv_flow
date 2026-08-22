#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <err.h>
#include <time.h>

// ---------------------------------------------------------------------
// Milestones 5 & 6: local shear rates + Carreau shear-thinning viscosity
//
// Same plain-channel geometry and BCs as before (no nozzle). What's new:
//   - cell-centered strain-rate tensor D_ij = 0.5*(dvi/dxj + dvj/dxi),
//     and the scalar invariant |D| = sqrt(D:D)  (paper convention, not
//     the rheology-standard sqrt(2D:D) -- see memory/derivation notes)
//   - eta(|D|) via the Carreau model:
//         eta = mu_inf + (mu_0-mu_inf) * (1 + (gamma*|D|)^2)^((q-1)/2)
//   - the viscous term in the momentum equation is now the full stress
//     divergence div(2*eta*D), not eta*laplacian(v) -- required because
//     eta varies in space, so the two forms are no longer equivalent.
//
// Discretization notes:
//   - face-NORMAL derivatives (e.g. dvx/dx at the east face) are exact,
//     using the existing direct neighbor-difference stencils.
//   - face-TANGENTIAL derivatives (e.g. dvy/dy at the east face, needed
//     for |D| there) aren't directly available at that face, so they're
//     approximated by averaging the cell-centered tangential gradient
//     between the two cells sharing the face. This only needs the
//     existing E/W/N/S connectivity -- no diagonal-neighbor stencil.
//   - the pressure Poisson source term (b_src) is left as-is, i.e. still
//     derived assuming constant viscosity. Re-deriving it for spatially
//     varying eta would add extra grad(eta) terms; out of scope for now
//     and noted as a documented simplification.
// ---------------------------------------------------------------------

#define Nx 100
#define Ny 10
#define NXG (Nx+2)
#define NYG (Ny+2)
#define Ntot (NXG*NYG)

double dt=0.02;      // raised from 0.001 -- Courant/diffusive stability checks
                     // showed >10x margin at this value (see dimensionless_checks.py);
                     // needed since a 10x longer duct needs proportionally longer
                     // to reach steady state, and the old dt would be impractically slow
double U_in=0.1;     // inlet velocity (Dirichlet)
double P_out=0.0;    // outlet reference pressure (Dirichlet)

// Carreau model parameters (paper convention: mu_0, mu_inf, gamma, q, |D|)
// PLACEHOLDER VALUES -- replace with the paper's fitted parameters before
// drawing quantitative conclusions. Setting mu_0=mu_inf recovers the
// Newtonian (constant-viscosity) case, useful as a sanity check that this
// reduces back to the milestone-3-validated behaviour.
double CARREAU_MU0   = 1.0;   // zero-shear-rate viscosity
double CARREAU_MUINF = 0.1;   // infinite-shear-rate viscosity
double CARREAU_GAMMA = 300.0; // time-constant parameter (paper's "gamma")
double CARREAU_Q     = 0.5;   // power-law index (q<1 : shear-thinning)

double p[Ntot], vx[Ntot], vy[Ntot];
int E[Ntot], W[Ntot], N[Ntot], S[Ntot];
double dx[Ntot], dy[Ntot];
double pe[Ntot], pw[Ntot], pn[Ntot], ps[Ntot];
double vxe[Ntot], vxw[Ntot], vxn[Ntot], vxs[Ntot];
double vye[Ntot], vyw[Ntot], vyn[Ntot], vys[Ntot];

// cell-centered velocity gradients, strain-rate invariant, and viscosity
double CDxVx[Ntot], CDxVy[Ntot], CDyVx[Ntot], CDyVy[Ntot];
double Dmag[Ntot];
double eta[Ntot];

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
    for (int i=0;i<Ntot;i++){
        p[i]=0.0; vx[i]=0.0; vy[i]=0.0;
        eta[i]=CARREAU_MU0;  // start at zero-shear viscosity everywhere
    }
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

// Milestone 5: local shear rate.  Milestone 6: Carreau viscosity from it.
// Computed once per outer iteration, on the velocity field as it stands
// at the start of that iteration (lagged, like b_src) -- eta is then held
// fixed while update_velocity() uses it.
void compute_shear_and_viscosity()
{
    // cell-centered gradients via direct central differences. Loop over
    // the WHOLE grid (including ghosts) since E/W/N/S are always valid
    // (clamped at the array edges in init_control_points), and interior
    // cells next to a ghost need that ghost's gradient for the face
    // averaging step below.
    for (int i=0;i<Ntot;i++){
        CDxVx[i]=(vx[E[i]]-vx[W[i]])/(2.0*dx[i]);
        CDxVy[i]=(vy[E[i]]-vy[W[i]])/(2.0*dx[i]);
        CDyVx[i]=(vx[N[i]]-vx[S[i]])/(2.0*dy[i]);
        CDyVy[i]=(vy[N[i]]-vy[S[i]])/(2.0*dy[i]);
    }

    for (int i=0;i<Ntot;i++){
        double Dxx=CDxVx[i];
        double Dyy=CDyVy[i];
        double Dxy=0.5*(CDyVx[i]+CDxVy[i]);
        Dmag[i]=sqrt(Dxx*Dxx+Dyy*Dyy+2.0*Dxy*Dxy);   // |D| = sqrt(D:D), paper convention
        double gd=CARREAU_GAMMA*Dmag[i];
        eta[i]=CARREAU_MUINF+(CARREAU_MU0-CARREAU_MUINF)*pow(1.0+gd*gd,(CARREAU_Q-1.0)/2.0);
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
            // face-normal derivatives (exact)
            double dxvxe=2*(vx[E[i]]-vx[i])/(dx[i]+dx[E[i]]);
            double dxvye=2*(vy[E[i]]-vy[i])/(dx[i]+dx[E[i]]);
            double dxvxw=2*(-vx[W[i]]+vx[i])/(dx[i]+dx[W[i]]);
            double dxvyw=2*(-vy[W[i]]+vy[i])/(dx[i]+dx[W[i]]);
            double dyvxn=2*(vx[N[i]]-vx[i])/(dy[i]+dy[N[i]]);
            double dyvyn=2*(vy[N[i]]-vy[i])/(dy[i]+dy[N[i]]);
            double dyvxs=2*(-vx[S[i]]+vx[i])/(dy[i]+dy[S[i]]);
            double dyvys=2*(-vy[S[i]]+vy[i])/(dy[i]+dy[S[i]]);

            // face viscosities (averaged from the cell-centered Carreau eta)
            double eta_e=0.5*(eta[i]+eta[E[i]]);
            double eta_w=0.5*(eta[i]+eta[W[i]]);
            double eta_n=0.5*(eta[i]+eta[N[i]]);
            double eta_s=0.5*(eta[i]+eta[S[i]]);

            // face-tangential derivatives (approximated: average of the two
            // cell-centered gradients sharing that face -- see header note)
            double dyvx_e=0.5*(CDyVx[i]+CDyVx[E[i]]);
            double dyvx_w=0.5*(CDyVx[i]+CDyVx[W[i]]);
            double dxvy_n=0.5*(CDxVy[i]+CDxVy[N[i]]);
            double dxvy_s=0.5*(CDxVy[i]+CDxVy[S[i]]);

            // full stress divergence div(2*eta*D), replacing the old
            // constant-viscosity Laplacian form
            double tau_xx_e=2.0*eta_e*dxvxe;
            double tau_xx_w=2.0*eta_w*dxvxw;
            double tau_xy_n=eta_n*(dyvxn+dxvy_n);
            double tau_xy_s=eta_s*(dyvxs+dxvy_s);
            double visc_x=(tau_xx_e-tau_xx_w)/dx[i]+(tau_xy_n-tau_xy_s)/dy[i];

            double tau_yy_n=2.0*eta_n*dyvyn;
            double tau_yy_s=2.0*eta_s*dyvys;
            double tau_xy_e=eta_e*(dxvye+dyvx_e);
            double tau_xy_w=eta_w*(dxvyw+dyvx_w);
            double visc_y=(tau_xy_e-tau_xy_w)/dx[i]+(tau_yy_n-tau_yy_s)/dy[i];

            vx[i]+=dt*(-vx[i]*dxvx-vy[i]*dyvx-(pe[i]-pw[i])/dx[i]+visc_x);
            vy[i]+=dt*(-vx[i]*dxvy-vy[i]*dyvy-(pn[i]-ps[i])/dy[i]+visc_y);
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

void write_shear_eta_2d(const char *name)
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    fprintf(F,"# x y |D| eta\n");
    for (int x=1;x<=Nx;x++){
        for (int y=1;y<=Ny;y++){
            int i=site2index(x,y);
            fprintf(F,"%d %d %e %e\n",x-1,y-1,Dmag[i],eta[i]);
        }
        fprintf(F,"\n");
    }
    fclose(F);
}

void solver()
{
    double tolerance=1e-6;
    int it=0, max_it=60000;
    double vnorm_old=0.0, check=1.0;

    double inner_tol=1e-8;
    int max_inner=6000;  // bumped for the 10x larger pressure solve (100x10 grid)

    FILE *O=fopen("v_t.dat","w");
    if(!O) err(1,"Could not create v_t.dat");
    // report-ready log: outer iteration, |v|, check, SOR sweeps used, SOR residual
    FILE *L=fopen("convergence_log.dat","w");
    if(!L) err(1,"Could not create convergence_log.dat");
    fprintf(L,"# it  |v|  check  sor_sweeps  sor_residual\n");

    do{
        apply_boundary_conditions();
        compute_shear_and_viscosity();     // milestones 5+6: |D| and eta from the current field
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
        if (it%100==0){
            printf("it %d  |v| %.6e  check %.3e  (SOR: %d sweeps, resid %.2e)\n",
                   it,vnorm,check,sor_iters,sor_resid);
            fflush(stdout);
        }
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
    write_shear_eta_2d("final_shear_eta.dat");
    return 0;
}
