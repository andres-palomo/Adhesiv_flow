#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <err.h>
#include <time.h>

#define Nx 10
#define Ny 10
int Ntot=Nx*Ny;
double dt=0.001;
//take all the physics parameter = 1

double p[Nx*Ny]; //pressure in control points
double vx[Nx*Ny]; //x-component of the velocity in control points
double vy[Nx*Ny]; //y-component of the velocity in control points
int E[Nx*Ny]; //index of the neighbour control points in the E direction
int W[Nx*Ny]; //index of the neighbour control points in the W direction
int N[Nx*Ny]; //index of the neighbour control points in the N direction
int S[Nx*Ny]; //index of the neighbour control points in the S direction
int dx[Nx*Ny]; //deltax for each cell (it is just a constant in the easy case)
int dy[Nx*Ny]; //deltay for each cell
double pe[Nx*Ny]; //p in the neighbour face in the E direction
double pw[Nx*Ny]; //p in the neighbour face in the W direction
double pn[Nx*Ny]; //p in the neighbour face in the N direction
double ps[Nx*Ny]; //p in the neighbour face in the S direction
double vxe[Nx*Ny]; //vx in the neighbour face in the E direction
double vxw[Nx*Ny]; //vx in the neighbour face in the W direction
double vxn[Nx*Ny]; //vx in the neighbour face in the N direction
double vxs[Nx*Ny]; //vx in the neighbour face in the S direction
double vye[Nx*Ny]; //vy in the neighbour face in the E direction
double vyw[Nx*Ny]; //vy in the neighbour face in the W direction
double vyn[Nx*Ny]; //vy in the neighbour face in the N direction
double vys[Nx*Ny]; //vy in the neighbour face in the S direction

// geometry is setting up the grid spacing (dx, dy) for every cell.
// In this baseline version the mesh is uniform, so it just fills dx/dy with 1,
// but it is kept as its own function so a non-uniform grid can be plugged in later.
void geometry()
{
	//setup the geometry of the system
	//here we write an easy case where all the cells have the same dimension, but from here one can generalise
	for (int i=0;i<Ntot;i++){
		dx[i]=dy[i]=1;
	}
}

// site2index is converting 2D grid coordinates (x,y) into the 1D array index used
// to store every field. It applies periodic wrap-around (modulo nx/ny) so that
// out-of-range coordinates (e.g. neighbour lookups past the last cell) fold back
// onto the opposite side of the domain.
int site2index(int x,int y,int nx,int ny)
{
    //write the coordinates inside the box, considering the boundary conditions
    int xx=(x+nx)%nx; //if x < nx, xx=x. if for example x=nx+2, xx=2
    int yy=(y+ny)%ny;
    return xx+yy*nx; //transform the 2d coordinates (xx,yy) in a 1d index
}

// init_control_points is doing the bookkeeping of the mesh connectivity: for every
// cell it computes and stores the array index of its East/West/North/South
// neighbours (E, W, N, S), using site2index and the local cell sizes.
void init_control_points()
{
	for (int x=0;x<Nx;x++){
		for (int y=0;y<Ny;y++){
			int i=site2index(x,y,Nx,Ny);
			//neighbour control points
			//each control point is in the centre of its own cell
			E[i]=site2index(x+(int)(dx[i]+dx[i+1])/2,y,Nx,Ny);
			W[i]=site2index(x-(int)(dx[i]+dx[i-1])/2,y,Nx,Ny);
			N[i]=site2index(x,y+(int)(dy[i]+dy[i+1])/2,Nx,Ny);
			S[i]=site2index(x,y-(int)(dy[i]+dy[i+1])/2,Nx,Ny);
		}
	}
}

// measure_pressure is doing a simple diagnostic: it sums the pressure field p[]
// over every cell, giving a single scalar used elsewhere to track how much the
// total pressure changes between iterations (convergence check).
double measure_pressure()
{
	double pval=0.;
	for (int i=0;i<Ntot;i++){
		pval+=p[i];
	}
	return pval;
}

// init_p is setting the initial pressure field, imposing a zero-gradient
// (Neumann) condition at the left (x=0) and right (x=Nx) boundaries by giving
// matching random values to the E/W and N/S neighbour pressures there.
void init_p()
{
	//dx p(c) = dy p(c) = 0 at the boundary
	// so p[E[i]]=p[W[i]], p[N[i]]=p[S[i]]
	srand(time(NULL));
	int x=0;
	for (int y=0;y<Ny;y++){
		int i=site2index(x,y,Nx,Ny);
		double a=0.01*(rand()/2147483647);
		p[E[i]]=a; // remember E[i] is still an index
		p[W[i]]=a;
		a=rand()/2147483647;
		p[N[i]]=a;
		p[S[i]]=a;
	}
	x=Nx;
	for (int y=0;y<Ny;y++){
		int i=site2index(x,y,Nx,Ny);
		double a=0.01*(rand()/2147483647);
		p[E[i]]=a;
		p[W[i]]=a;
		a=rand()/2147483647;
		p[N[i]]=a;
		p[S[i]]=a;
	}
}

// init_v is setting the initial velocity field: zero velocity everywhere in the
// interior of the domain, and an imposed horizontal inflow/outflow velocity
// (vx=0.1, vy=0) on the left (x=0) and right (x=Nx) boundary columns.
void init_v()
{
	for (int x=1;x<Nx-1;x++){
		for (int y=1;y<Ny-1;y++){
			int i=site2index(x,y,Nx,Ny);
			vx[i]=0.;
			vy[i]=0.;
		}
	}
	int x=0;
	for (int y=0;y<Ny;y++){
		int i=site2index(x,y,Nx,Ny);
		vx[i]=0.1;
		vy[i]=0.;
	}
	x=Nx;
	for (int y=0;y<Ny;y++){
		int i=site2index(x,y,Nx,Ny);
		vx[i]=0.1;
		vy[i]=0.;
	}

}

// setup is the top-level initialization routine: it builds the geometry, computes
// the neighbour connectivity, and initializes the pressure and velocity fields,
// in that order, before the solver runs.
void setup()
{
	geometry();
	init_control_points();
	init_p();
	//double p_init=measure_pressure(); //check
	init_v();
}

// interpolate_faces is computing face values of pressure and velocity by
// averaging each cell's value with its E/W/N/S neighbour's value. These
// interpolated face quantities (pe/pw/pn/ps, vxe/vxw/..., vye/vyw/...) are what
// the finite-volume update equations actually use.
void interpolate_faces()
{
	for (int i=0;i<Ntot;i++){
		pe[i]=(p[E[i]]+p[i])/2; //interpolate between control points to find the value on the face
		pw[i]=(p[W[i]]+p[i])/2;
		pn[i]=(p[N[i]]+p[i])/2;
		ps[i]=(p[S[i]]+p[i])/2;
		vxe[i]=(vx[E[i]]+vx[i])/2; //interpolate between control points to find the value on the face
		vxw[i]=(vx[W[i]]+vx[i])/2;
		vxn[i]=(vx[N[i]]+vx[i])/2;
		vxs[i]=(vx[S[i]]+vx[i])/2;
		vye[i]=(vy[E[i]]+vy[i])/2; //interpolate between control points to find the value on the face
		vyw[i]=(vy[W[i]]+vy[i])/2;
		vyn[i]=(vy[N[i]]+vy[i])/2;
		vys[i]=(vy[S[i]]+vy[i])/2;
	}
}

// update_velocity is advancing vx and vy by one explicit time step dt, evaluating
// the discretized momentum equation at every cell: advection terms (dxvx, dyvx,
// ...), the pressure gradient term ((pe-pw)/dx, (pn-ps)/dy), and viscous diffusion
// terms built from the face-normal derivatives (dxvxe, dxvxw, dyvxn, dyvxs, ...).
// The equation numbers in the comments refer to Korzec's CFD notes.
void update_velocity()
{
	for (int i=0;i<Ntot;i++){
	double dxvx=(vxe[i]-vxw[i])/dx[i]; //3.53
	double dxvy=(vye[i]-vxw[i])/dx[i];
	double dyvx=(vxn[i]-vxs[i])/dy[i]; //3.54
	double dyvy=(vyn[i]-vys[i])/dy[i];
	double dxvxe=2*(vx[E[i]]-vx[i])/(dx[i]+dx[E[i]]); //3.55
	// || E - c || = dx[i]/2 + dx[E[i]]/2
	double dxvye=2*(vy[E[i]]-vy[i])/(dx[i]+dx[E[i]]);
	double dxvxw=2*(-vx[W[i]]+vx[i])/(dx[i]+dx[W[i]]); //3.56
	double dxvyw=2*(-vy[W[i]]+vy[i])/(dx[i]+dx[W[i]]);
	double dyvxn=2*(vx[N[i]]-vx[i])/(dy[i]+dy[N[i]]); //3.57
	double dyvyn=2*(vy[N[i]]-vy[i])/(dy[i]+dy[N[i]]);
	double dyvxs=2*(-vx[S[i]]+vx[i])/(dy[i]+dy[S[i]]); //3.58
	double dyvys=2*(-vy[S[i]]+vy[i])/(dy[i]+dy[S[i]]);
	//see Eq.3.59 in Korzec's notes
	vx[i]+=dt*(-vx[i]*dxvx-vy[i]*dyvx-(pe[i]-pw[i])/dx[i]+(dxvxe-dxvxw)/dx[i]+(dyvxn-dyvxs)/dx[i]);
	vy[i]+=dt*(-vx[i]*dxvy-vy[i]*dyvy-(pn[i]-ps[i])/dy[i]+(dxvye-dxvyw)/dx[i]+(dyvyn-dyvys)/dx[i]);
	//printf("update v vx[i] %e pe[i] %e pw[i] %e \n",vx[i],pe[i],pw[i]);
	}
}

// update_pressure is doing the pressure-correction step of the SIMPLE-like
// algorithm: it recomputes p[i] at every cell from the local velocity-divergence
// term (b) and the interpolated neighbour face pressures, following the
// un-numbered pressure-update equation in section 3.5.5 of Korzec's notes.
void update_pressure()
{
	//see un-numbered Eq for p evolution in section 3.5.5 of Korzec's notes
	for (int i=0;i<Ntot;i++){
		double dxvx=(vxe[i]-vxw[i])/dx[i]; //3.53
		double dxvy=(vye[i]-vxw[i])/dx[i];
		double dyvx=(vxn[i]-vxs[i])/dy[i]; //3.54
		double dyvy=(vyn[i]-vys[i])/dy[i];
		double b=-dxvx*dxvx-dyvy*dyvy+2*dxvy*dyvx;
		//printf("b %e \n",b);
		p[i]=-0.125*b*dx[i]*dx[i]*dy[i]*dy[i]+(pn[i]+ps[i])*dx[i]*dx[i]+(pe[i]+pw[i])*dy[i]*dy[i];
		p[i]=(0.5*p[i])/(dx[i]*dx[i]+dy[i]*dy[i]);
	}
}

// set_guess is filling the pressure field with random values as the starting
// guess for the iterative solver in solver().
void set_guess()
{
	//tmp guess
	srand(time(NULL));
        for (int i=0;i<Ntot;i++){
		p[i]=rand()/2147483647.;
	}
}

// write_p_2d is writing the current pressure field to a text file, one line per
// cell as "x y p(x,y)", with a blank line between rows of x, so the file can be
// plotted directly (e.g. with gnuplot's pm3d/splot).
void write_p_2d(const char *name)/*{{{*/
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=0;x<Nx;x++){
	    for (int y=0;y<Ny;y++){
		    int i=site2index(x,y,Nx,Ny);
		    fprintf(F,"%d %d %g\n",x,y,p[i]);
	    }
	    fprintf(F,"\n");
    }
    fclose(F);
}

// write_v_2d is writing the current velocity field to a text file, one line per
// cell as "x y vx(x,y) vy(x,y)", with a blank line between rows of x, in the same
// plot-friendly layout as write_p_2d.
void write_v_2d(const char *name)/*{{{*/
{
    FILE *F=fopen(name,"w");
    if(!F) err(1,"Could not create %s",name);
    for (int x=0;x<Nx;x++){
	    for (int y=0;y<Ny;y++){
		    int i=site2index(x,y,Nx,Ny);
		    fprintf(F,"%d %d %e %e \n",x,y,vx[i],vy[i]);
	    }
	    fprintf(F,"\n");
    }
    fclose(F);
}

// write_p_tot is appending one line per cell index with the same scalar value p
// to an already-open file F. It is used by solver() to log the total pressure
// (measure_pressure's result) at every solver iteration into p_t.dat.
void write_p_tot(FILE *F,double p)/*{{{*/
{
    for (int i=0;i<Ntot;i++){
	    fprintf(F,"%d %g\n",i,p);
    }
}

// solver is the main iterative loop of the simulation: starting from a random
// pressure guess, it repeatedly interpolates face values, updates velocity and
// pressure, measures the total pressure before/after the step, and logs it to
// p_t.dat, stopping once the relative change in total pressure between two
// consecutive iterations drops below the tolerance (1e-6).
void solver()
{
    //See pseudocode in 3.5.5 of Korzec's notes
    double tolerance=1e-6;
    int it=0;
    set_guess();
    double pval=0.;
    double pstart=0.;
    double check=0.;
    FILE *O=fopen("p_t.dat","w");
    if(!O) err(1,"Could not create p_2d.dat");
    do{
	interpolate_faces();
	double pstart=measure_pressure(); //measure total p
        update_velocity();
	update_pressure();
        pval=measure_pressure();
	write_p_tot(O,pval);
        it++;
	check=fabs((pval-pstart)/pstart);
        printf("pstart %.10e pval %.10e pval-pstart %e (pval-pstart)/pstart %e \n",pstart,pval,fabs(pval-pstart),check);
    } while(check>tolerance);

    fclose(O);
}


// main is the program entry point: it builds and initializes the mesh/fields
// (setup), dumps the initial pressure and velocity fields to file, runs the
// iterative solver, and finally dumps the converged pressure and velocity
// fields to file.
int main(int argc,char **argv)
{
    setup();
    write_p_2d("start_pressure.dat");
    write_v_2d("start_velocity.dat");
    solver();
    write_p_2d("final_pressure.dat");
    write_v_2d("final_velocity.dat");
}