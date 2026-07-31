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

void geometry()
{
	//setup the geometry of the system
	//here we write an easy case where all the cells have the same dimension, but from here one can generalise
	for (int i=0;i<Ntot;i++){
		dx[i]=dy[i]=1;
	}
}

int site2index(int x,int y,int nx,int ny)
{
    //write the coordinates inside the box, considering the boundary conditions
    int xx=(x+nx)%nx; //if x < nx, xx=x. if for example x=nx+2, xx=2
    int yy=(y+ny)%ny;
    return xx+yy*nx; //transform the 2d coordinates (xx,yy) in a 1d index 
}

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

double measure_pressure()
{
	double pval=0.;
	for (int i=0;i<Ntot;i++){
		pval+=p[i];
	}
	return pval;
}

void init_p()
{
	//dx p(c) = dy p(c) = 0 at the boundary
        // so p[E[i]]=p[W[i]], p[N[i]]=p[S[i]]
        srand(time(NULL));
        int x=0;
        for (int y=0;y<Ny;y++){
                int i=site2index(x,y,Nx,Ny);
                double a=0.01*(rand()/2147483647);
                p[E[i]]=a;
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

void setup()
{
	geometry();
	init_control_points();
	init_p();
	//double p_init=measure_pressure(); //check
	init_v();
}

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

void update_velocity()
{
}

void update_pressure()
{
}

void set_guess()
{
}

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

void write_p_tot(FILE *F,double p)/*{{{*/
{
    for (int i=0;i<Ntot;i++){
	    fprintf(F,"%d %g\n",i,p);
    }
}

void solver()
{
}


int main(int argc,char **argv)
{
    setup();
    write_p_2d("start_pressure.dat");
    write_v_2d("start_velocity.dat");
    solver();
    write_p_2d("final_pressure.dat");
    write_v_2d("final_velocity.dat");
}

