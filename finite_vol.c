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

void set_guess()
{
	//tmp guess
	srand(time(NULL));
        for (int i=0;i<Ntot;i++){
		p[i]=rand()/2147483647.;
	}
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


int main(int argc,char **argv)
{
    setup();
    write_p_2d("start_pressure.dat");
    write_v_2d("start_velocity.dat");
    solver();
    write_p_2d("final_pressure.dat");
    write_v_2d("final_velocity.dat");
}

