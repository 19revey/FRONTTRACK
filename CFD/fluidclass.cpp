//
//  fluidclass.cpp
//  CFD
//
//  Created by YDuan on 9/30/16.
//  Copyright © 2016 Duan. All rights reserved.
//
#include <string.h>
#include "fluidclass.h"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std ;



void fluidclass::timeloop(){
    for (int is=0;is<nstep;is++)
    {
        time=is;
        //cout<<ny<<endl;
        copy2d(un, u, nx+1, ny+2);
        copy2d(vn, v, nx+2, ny+1);
        copy2d(rn, r, nx+2, ny+2);
        copy2d(mn, m, nx+2, ny+2);
        copy1d(xfn, xf, 10*Nf+2);
        copy1d(yfn, yf, 10*Nf+2);

        for (int substep=1; substep<=2; substep++) {
            
            find_surface_tension();
            
            calculate_boundary_velocity();
            
            solve_velocity_adv_diff();
            //solve_velocity_without_pressure();          //temporary u velocity
            
            solve_pressure();                           //compute source term and the coefficient for p(i,j)
            
            correct_velocity_after_pressure();          //correct velocity
            
            
            //advect_density();
            
            
            advect_front();

            distribute_gradient();
            update_viscosity();
            
            
            
            
        }
        
        
        /////
        for (int i=0; i<=nx; i++) {
            for (int j=0; j<=ny+1; j++) {
                u[i][j]=(u[i][j]+un[i][j])/2;
            }
        }
        for (int i=0; i<=nx+1; i++) {
            for (int j=0; j<=ny; j++) {
                v[i][j]=(v[i][j]+vn[i][j])/2;
            }
        }
        for (int i=0; i<=nx+1; i++) {
            for (int j=0; j<=ny+1; j++) {
                r[i][j]=(r[i][j]+rn[i][j])/2;
                m[i][j]=(m[i][j]+mn[i][j])/2;
            }
        }
        for (int i=0; i<10*Nfv+2; i++) {
            xf[i]=(xf[i]+xfn[i])/2;
            yf[i]=(yf[i]+yfn[i])/2;
        }//////
        add_points_to_front();
        
        if((int)time % print_interval==0) {
            
            print_vtk();
            print_dem();
        }
        
    }
    
    //cout<<u[25][10]<<endl;
}

void fluidclass::initialize_fluid(){
    
    rho1=1;rho2=2;m0=0.01;rro=rho1;
    rad=0.15;xc=0.5;yc=0.7;
    
    u=make2dmem(nx+1, ny+2,0);
    v=make2dmem(nx+2, ny+1,0);
    p=make2dmem(nx+2, ny+2,0);
    ut=make2dmem(nx+1,ny+2,0);
    vt=make2dmem(nx+2, ny+1,0);
    tmp1=make2dmem(nx+2, ny+2,0);
    uu=make2dmem(nx+1, ny+1,0);
    vv=make2dmem(nx+1, ny+1,0);
    tmp2=make2dmem(nx+2, ny+2,0);
    r=make2dmem(nx+2, ny+2, rho1);
    
    
    //pressure from last iteration
    p_old=make2dmem(nx+2, ny+2, 0);
    rt=make2dmem(nx+2, ny+2, 0);
    ro=make2dmem(nx+2, ny+2,0);
    
    //set the grid
    dx=lx/nx;dy=ly/ny;
    x=new double[nx+2];
    y=new double[ny+2];
    for (int i=0; i<=nx+1; i++) {
        x[i]=dx*(i-0.5);
    }
    for (int j=0; j<=ny+1; j++) {
        y[j]=dy*(j-0.5);
    }
    
    

    //set front//////////////////////////////////////////////
    Nf=100;Nfv=Nf;                                                //
    xf=make1dmem(10*Nf+2, 0); yf=make1dmem(10*Nf+2, 0);    //
    uf=make1dmem(10*Nf+2, 0); vf=make1dmem(10*Nf+2, 0);    //
    xfold=make1dmem(10*Nf+2, 0);                           //
    yfold=make1dmem(10*Nf+2, 0);                           //
    for (int l=0; l<=Nf+1; l++) {                          //
        xf[l]=xc-rad*sin(2.0*M_PI*(l)/Nf);                 //
        yf[l]=yc+rad*cos(2.0*M_PI*(l)/Nf);                 //
    }                                                      //
    fx=make2dmem(nx+2, ny+2, 0);                           //
    fy=make2dmem(nx+2, ny+2, 0);                           //
    r_old=make2dmem(nx+2, ny+2, 0);                        //
    /////////////////////////////////////////////////////////
    
    
    //tensor part/////////////////////////////////////////////////
    m1=0.01;m2=0.05;sigma=10;                                   //
    un=make2dmem(nx+1, ny+2, 0);                                //
    vn=make2dmem(nx+2, ny+1, 0);                                //
    m=make2dmem(nx+2, ny+2, m1);                                //
    rn=make2dmem(nx+2, ny+2, 0); mn=make2dmem(nx+2, ny+2, 0);   //
    tx=make1dmem(10*Nf+2, 0); ty=make1dmem(10*Nf+2, 0);         // tangent vectors
    xfn=make1dmem(10*Nf+2, 0); yfn=make1dmem(10*Nf+2, 0);
    //set density
    for (int i=1; i<=nx; i++)
        for(int j=1;j<=ny;j++)
            if( square(x[i]-xc) + square(y[j]-yc) <rad*rad)
            {
                r[i][j]=rho2;m[i][j]=m2;
            }
    
    //////////////////////////////////////////////////////////////

};

void fluidclass::calculate_boundary_velocity(){
    //tangential velocity at boundaries
    for (int index=0;index<=nx;index++)
    {
        u[index][0]=2*usouth-u[index][1];
        u[index][ny+1]=2*unorth-u[index][ny];
    }
    for (int index=0;index<=ny;index++)
    {
        v[0][index]=2*vwest-v[1][index];
        v[nx+1][index]=2*veast-v[nx][index];
    }
    
    
}

void fluidclass::solve_velocity_without_pressure(){
    for (int i=1; i<=nx-1; i++)
        for (int j=1;j<=ny;j++)
        {
            double Advx=0.25 * (square(u[i+1][j]+u[i][j])-square(u[i][j]+u[i-1][j]))/dx + 0.25*( (u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j]) - (u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]) )/dy;
            double Difx=( ( u[i+1][j]-2*u[i][j]+u[i-1][j] )/square(dx)+ ( u[i][j+1]-2*u[i][j]+u[i][j-1] )/square(dy))*m0;
            ut[i][j]=u[i][j]+dt*(-Advx+Difx/(0.5*(r[i+1][j]+r[i][j])) +(1.0-rro/(0.5*(r[i+1][j]+r[i][j])))*gx );
        }
    for (int i=1;i<=nx;i++)
        for(int j=1;j<=ny-1;j++)
        {
            double Advy=0.25 * ((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i-1][j+1]+u[i-1][j])*(v[i][j]+v[i-1][j]))/dx+0.25*(square(v[i][j+1]+v[i][j])-square(v[i][j]+v[i][j-1]))/dy;
            double Dify=((v[i+1][j]-2*v[i][j]+v[i-1][j])/square(dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/square(dy))*m0;
            vt[i][j]=v[i][j]+dt*(-Advy+Dify/(0.5*(r[i][j+1]+r[i][j])) +(1.0-rro/(0.5*(r[i][j+1]+r[i][j])))*gy);
        }
}

int fluidclass::check_p_convg(){

    for (int i=0; i<nx+2; i++) {
        for (int j=0; j<ny+2; j++) {
            if (square(p[i][j]-p_old[i][j])>square(maxError)) {
                return 0;
            }
        }
    }
    return 1;
}

void fluidclass::solve_pressure(){
    double lrg=1000;

    copy2d(rt, r, (nx+2),(ny+2));

    for(int i=0;i<=nx+1;i++) {rt[i][0]=lrg; rt[i][ny+1]=lrg;} //set up boundary p
    for(int j=0;j<=ny+1;j++) {rt[0][j]=lrg; rt[nx+1][j]=lrg;} //set up boundary p
    for (int i=1; i<=nx; i++) {
        for (int j=1; j<=ny; j++) {
            tmp1[i][j]=(0.5/dt)* ( (ut[i][j]-ut[i-1][j])/dx+(vt[i][j]-vt[i][j-1])/dy  );
            tmp2[i][j]=1.0/(   1.0/square(dx)* (1.0/(rt[i+1][j]+rt[i][j])+1.0/(rt[i-1][j]+rt[i][j])) + 1.0/square(dy)*(   1.0/(rt[i][j+1]+rt[i][j])+1.0/(rt[i][j-1]+rt[i][j]) )   );
        }
    }
    for (int it=0; it<=maxit; it++) {

        copy2d(p_old, p, (nx+2),(ny+2));
        //solve pressure from iteration
        for (int i=1; i<=nx; i++) {
            for (int j=1; j<=ny; j++) {
                p[i][j]=beta*tmp2[i][j]*( 1.0/square(dx)*( p[i+1][j]/(rt[i+1][j]+rt[i][j])+p[i-1][j]/(rt[i][j]+rt[i-1][j]) )+ 1.0/square(dy)*( p[i][j+1]/(rt[i][j+1]+rt[i][j])+p[i][j-1]/(rt[i][j-1]+rt[i][j]) ) -tmp1[i][j])+(1-beta)*p[i][j];
            }
        }
        //check error to break
        if(check_p_convg())
        {   //cout<<it<<endl;
            break;}
        
    }

}

void fluidclass::correct_velocity_after_pressure(){
    for (int i=1; i<=nx-1; i++) {
        for (int j=1; j<=ny; j++) {
            u[i][j]=ut[i][j]-dt*(2.0/dx)*(p[i+1][j]-p[i][j])/(r[i+1][j]+r[i][j]);
        }
    }
    for (int i=1; i<=nx; i++) {
        for (int j=1; j<=ny-1; j++) {
            v[i][j]=vt[i][j]-dt*(2.0/dy)*(p[i][j+1]-p[i][j])/(r[i][j+1]+r[i][j]);
        }
    }
    
}

void fluidclass::advect_density(){

    copy2d(ro,r,(nx+2),(ny+2));
    for (int i=1; i<=nx; i++) {
        for (int j=1; j<=ny; j++) {
            r[i][j]=ro[i][j]-(0.5*dt/dx)*( u[i][j]*(ro[i+1][j]+ro[i][j])-u[i-1][j]*(ro[i-1][j]+ro[i][j]) ) - (0.5*dt/dy)*( v[i][j]*(ro[i][j+1]+ro[i][j])-v[i][j-1]*(ro[i][j-1]+ro[i][j]) ) + m0*dt/dx/dx*(ro[i+1][j]-2*ro[i][j]+ro[i-1][j]) + m0*dt/dy/dy*(ro[i][j+1]-2*ro[i][j]+ro[i][j-1]);
        }
    }

}

void fluidclass::print_vtk(){
    int i,j;
    
    for (i=0; i<=nx; i++) {
        for (j=0; j<=ny; j++) {
            uu[i][j]=0.5*( u[i][j+1]+u[i][j] );
            vv[i][j]=0.5*( v[i+1][j]+v[i][j] );
        }
    }
    
    FILE* fp;
    char fname[40];
    sprintf(fname, "fluid.%d.vtk",(int)time);
    fp=fopen(fname,"w");
    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"example\nASCII\n");
    fprintf(fp,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp,"DIMENSIONS  %d %d 1\n",(nx+1), (ny+1));
    fprintf(fp,"ORIGIN   0 0 0\n");
    fprintf(fp,"SPACING   %e %e 0\n",1.0/(nx),1.0/(ny));
    fprintf(fp,"POINT_DATA %d\n",(nx+1)*(ny+1)*1);
    
    // to print scalar
    fprintf(fp,"SCALARS pressure float 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for(j=0;j<=ny;j++){
        for(i=0;i<=nx;i++){
            fprintf(fp,"%.8e\n", p[i][j]);
        }
    }
    // add more scalars if needed
    fprintf(fp,"SCALARS density float 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for(j=0;j<=ny;j++){
        for(i=0;i<=nx;i++){
            fprintf(fp,"%.8e\n", r[i][j]);
        }
    }
    
    // to print vector
    fprintf(fp,"VECTORS velocity float\n");
    for(j=0;j<=ny;j++){
        for(i=0;i<=nx;i++){
            fprintf(fp,"%.8e %.8e  0.\n",uu[i][j],vv[i][j]);
        }
    }
    
    fclose(fp);
    
}
void fluidclass::print_dem(){
    FILE *fvtk;
    char fname1[40];
    sprintf(fname1, "point.%d.vtk",(int)time);
    fvtk=fopen(fname1,"w");
    //..
    
    
    fprintf(fvtk,"# vtk DataFile Version 2.0\n");
    fprintf(fvtk,"Particle Tracking: ...\n");
    fprintf(fvtk,"ASCII\n\n");
    fprintf(fvtk,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fvtk,"POINTS %d double\n",Nfv);
    
    for(int j=1;j<=Nfv;j++){
        fprintf(fvtk,"%.4e %.4e 0\n",xf[j]+0.5*dx,yf[j]+0.5*dy);
    }
    
    fprintf(fvtk,"\nPOINT_DATA %d\n", Nfv);
    fprintf(fvtk,"SCALARS diameter double 1\n");
    fprintf(fvtk,"LOOKUP_TABLE DEFAULT\n");
    
    for(int j=1;j<=Nfv;j++){
        fprintf(fvtk,"%.5e\n",2*rad);
    }
    
    fclose(fvtk);
}



/////********************************************/////
void fluidclass::advect_front(){

    int ip=0;int jp=0;
    double ax,ay;
    //*****//calculate front velocity
    for (int lfn=1; lfn<=Nfv; lfn++) {
        //swap u
        ip=floor(xf[lfn]/dx);                         //u velocity id in x dir
        jp=floor((yf[lfn]+0.5*dy)/dy);                //u velocity id in y dir
        ax=xf[lfn]/dx-ip;
        ay=(yf[lfn]+0.5*dy)/dy-jp;
        uf[lfn]=(1-ax)*(1-ay)*u[ip][jp]+ax*(1-ay)*u[ip+1][jp]+(1-ax)*ay*u[ip][jp+1]+ax*ay*u[ip+1][jp+1];
        //swap v
        ip=floor((xf[lfn]+0.5*dx)/dx);
        jp=floor(yf[lfn]/dy);
        ax=(xf[lfn]+0.5*dx)/dx-ip;
        ay=yf[lfn]/dy-jp;
        vf[lfn]=(1-ax)*(1-ay)*v[ip][jp]+ax*(1-ay)*v[ip+1][jp]+(1-ax)*ay*v[ip][jp+1]+ax*ay*v[ip+1][jp+1];
    }
    //*****//move front
    for (int i=1; i<=Nfv; i++) {
        xf[i]=xf[i]+dt*uf[i];
        yf[i]=yf[i]+dt*vf[i];
    }
    xf[0]=xf[Nfv];
    yf[0]=yf[Nfv];
    xf[Nfv+1]=xf[1];
    yf[Nfv+1]=yf[1];
    
    



}

/////********************************************/////
void fluidclass::add_points_to_front(){
    copy1d(xfold, xf, 10*Nf+2);
    copy1d(yfold, yf, 10*Nf+2);
    int j=0; double ds;
    for (int lfn=1; lfn<=Nfv; lfn++) {
        ds=sqrt( square((xfold[lfn]-xf[j])/dx)+square((yfold[lfn]-yf[j])/dy) );
        if (ds>0.5) {
            j=j+1;
            xf[j]=0.5*(xfold[lfn]+xf[j-1]);
            yf[j]=0.5*(yfold[lfn]+yf[j-1]);
            j=j+1;
            xf[j]=xfold[lfn]; yf[j]=yfold[lfn];
            /////////////******
        

            /////////////******
        }
        else if (ds<0.25){j=j;}
        else {
            j=j+1;
            xf[j]=xfold[lfn]; yf[j]=yfold[lfn];
        }
    }
    Nfv=j;
    xf[0]=xf[Nfv];
    yf[0]=yf[Nfv];
    xf[Nfv+1]=xf[1];
    yf[Nfv+1]=yf[1];
    for (int i=Nfv+2; i<10*Nf+2; i++) {
        xf[i]=0;yf[i]=0;
    }
    

}

void fluidclass::distribute_gradient(){
    //**normal vector
    double nfx,nfy;
    int ip,jp;
    //ip=0;jp=0;
    double ax,ay;
    //**set fx & fy to zero
    for (int i=0; i<=nx+1; i++) {
        for (int j=0; j<=ny+1; j++) {
            fx[i][j]=0;
            fy[i][j]=0;
        }
    }
    for (int lfn=1; lfn<=Nfv; lfn++) {
        nfx=-0.5*(yf[lfn+1]-yf[lfn-1])*(rho2-rho1);
        nfy=0.5*(xf[lfn+1]-xf[lfn-1])*(rho2-rho1);
        ip=floor(xf[lfn]/dx);
        jp=floor((yf[lfn]+0.5*dy)/dy);
        ax=xf[lfn]/dx-ip;
        ay=(yf[lfn]+0.5*dy)/dy-jp;
        fx[ip][jp]=fx[ip][jp]+(1-ax)*(1-ay)*nfx/dx/dy;
        fx[ip+1][jp]=fx[ip+1][jp]+ax*(1-ay)*nfx/dx/dy;
        fx[ip][jp+1]=fx[ip][jp+1]+(1-ax)*ay*nfx/dx/dy;
        fx[ip+1][jp+1]=fx[ip+1][jp+1]+ax*ay*nfx/dx/dy;
        
        ip=floor((xf[lfn]+0.5*dx)/dx);
        jp=floor(yf[lfn]/dy);
        ax=(xf[lfn]+0.5*dx)/dx-ip;
        ay=yf[lfn]/dy-jp;
        fy[ip][jp]=fy[ip][jp]+(1-ax)*(1-ay)*nfy/dx/dy;
        fy[ip+1][jp]=fy[ip+1][jp]+ax*(1-ay)*nfy/dx/dy;
        fy[ip][jp+1]=fy[ip][jp+1]+(1-ax)*ay*nfy/dx/dy;
        fy[ip+1][jp+1]=fy[ip+1][jp+1]+ax*ay*nfy/dx/dy;
        
    }
    
    //set rho=rho1
    for (int i=0; i<nx+2; i++) {
        for (int j=0; j<ny+2; j++) {
            r[i][j]=0+rho1;
        }
    }
    for (int iter=1; iter<=maxit; iter++) {
        copy2d(r_old, r, nx+2, ny+2);
        for (int i=1; i<=nx; i++) {
            for (int j=1; j<=ny; j++) {
                r[i][j]=0.25*(r[i+1][j]+r[i-1][j]+r[i][j+1]+r[i][j-1]+dx*fx[i-1][j]-dx*fx[i][j]+dy*fy[i][j-1]-dy*fy[i][j]);
            }
        }
        if (check_r_convg()) {
            break;
        }
    }
    
    
}

int fluidclass::check_r_convg(){
    
    for (int i=0; i<nx+2; i++) {
        for (int j=0; j<ny+2; j++) {
            if (square(r[i][j]-r_old[i][j])>square(maxError)) {
                return 0;
            }
        }
    }
    return 1;
}

void fluidclass::find_surface_tension(){
    double ds;
    double nfx,nfy;
    int ip,jp;
    double ax,ay;
    for (int lfn=0; lfn<=Nfv; lfn++) {
        ds=sqrt( square(xf[lfn+1]-xf[lfn])+square(yf[lfn+1]-yf[lfn]) );
        tx[lfn]=( xf[lfn+1]-xf[lfn] )/ds;
        ty[lfn]=( yf[lfn+1]-yf[lfn] )/ds;
    }
    tx[Nfv+1]=tx[1];ty[Nfv+1]=ty[1];
    for (int lfn=1; lfn<=Nfv; lfn++) {
        nfx=sigma*( tx[lfn]-tx[lfn-1] );
        nfy=sigma*( ty[lfn]-ty[lfn-1] );
        ip=floor(xf[lfn]/dx);
        jp=floor((yf[lfn]+0.5*dy)/dy);
        ax=xf[lfn]/dx-ip;
        ay=(yf[lfn]+0.5*dy)/dy-jp;
        fx[ip][jp]=fx[ip][jp]+(1-ax)*(1-ay)*nfx/dx/dy;
        fx[ip+1][jp]=fx[ip+1][jp]+ax*(1-ay)*nfx/dx/dy;
        fx[ip][jp+1]=fx[ip][jp+1]+(1-ax)*ay*nfx/dx/dy;
        fx[ip+1][jp+1]=fx[ip+1][jp+1]+ax*ay*nfx/dx/dy;
        
        ip=floor((xf[lfn]+0.5*dx)/dx);
        jp=floor(yf[lfn]/dy);
        ax=(xf[lfn]+0.5*dx)/dx-ip;
        ay=yf[lfn]/dy-jp;
        fy[ip][jp]=fy[ip][jp]+(1-ax)*(1-ay)*nfy/dx/dy;
        fy[ip+1][jp]=fy[ip+1][jp]+ax*(1-ay)*nfy/dx/dy;
        fy[ip][jp+1]=fy[ip][jp+1]+(1-ax)*ay*nfy/dx/dy;
        fy[ip+1][jp+1]=fy[ip+1][jp+1]+ax*ay*nfy/dx/dy;
    }
    //tangent velocity at boundaries
    
}

void fluidclass::solve_velocity_adv_diff(){
    for(int i=1; i<=nx-1;i++)
        for (int j=1;j<=ny;j++)
        {
            double Advx=0.25 * (square(u[i+1][j]+u[i][j])-square(u[i][j]+u[i-1][j]))/dx + 0.25*( (u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j]) - (u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]) )/dy ;

            //double Difx=( ( u[i+1][j]-2*u[i][j]+u[i-1][j] )/square(dx)+ ( u[i][j+1]-2*u[i][j]+u[i][j-1] )/square(dy))*m0;
            ut[i][j]=u[i][j]+dt*(-Advx + fx[i][j]/(0.5*(r[i+1][j]+r[i][j]))+(1.0-rro/(0.5*(r[i+1][j]+r[i][j])))*gx );
        }

    for (int i=1;i<=nx;i++)
        for(int j=1;j<=ny-1;j++)
        {
            double Advy=0.25 * ((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i-1][j+1]+u[i-1][j])*(v[i][j]+v[i-1][j]))/dx+0.25*(square(v[i][j+1]+v[i][j])-square(v[i][j]+v[i][j-1]) )/dy ;
            //double Dify=((v[i+1][j]-2*v[i][j]+v[i-1][j])/square(dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/square(dy))*m0;
            vt[i][j]=v[i][j]+dt*(-Advy + fy[i][j]/(0.5*(r[i][j+1]+r[i][j]))+(1.0-rro/(0.5*(r[i][j+1]+r[i][j])))*gy);
        }
    
    for (int i=1; i<=nx-1; i++)
        for (int j=1; j<=ny; j++) {
            double Difx=1.0/dx*2*(m[i+1][j]*(u[i+1][j]-u[i][j])/dx-m[i][j]*(u[i][j]-u[i-1][j])/dy) +1/dy*(0.25*(m[i][j]+m[i+1][j]+m[i+1][j+1]+m[i][j+1])*( (u[i][j+1]-u[i][j])/dy+(v[i+1][j]-v[i][j])/dx) - 0.25*(m[i][j]+m[i+1][j]+m[i+1][j-1]+m[i][j-1])*( (u[i][j]-u[i][j-1])/dy+ (v[i+1][j-1]-v[i][j-1])/dx) );
            ut[i][j]=ut[i][j]+dt*Difx/(0.5*(r[i+1][j]+r[i][j]));
        }
    for (int i=1;i<=nx;i++)
        for(int j=1;j<=ny-1;j++)
        {
            double Dify=1.0/dy*2*(m[i][j+1]*(v[i][j+1]-v[i][j])/dy-m[i][j]*(v[i][j]-v[i][j-1])/dy) +1/dx*(0.25*(m[i][j]+m[i+1][j]+m[i+1][j+1]+m[i][j+1])*((u[i][j+1]-u[i][j])/dy+(v[i+1][j]-v[i][j])/dx) - 0.25*(m[i][j]+m[i][j+1]+m[i-1][j+1]+m[i-1][j])*(1/dy*(u[i-1][j+1]-u[i-1][j])+1/dx*(v[i][j]-v[i-1][j])) );
            vt[i][j]=vt[i][j]+dt*Dify/(0.5*(r[i][j+1]+r[i][j]));
        }
}


void fluidclass::update_viscosity(){
    for (int i=1; i<=nx; i++) {
        for (int j=1; j<=ny; j++) {
            m[i][j]=m1+(m2-m1)*(r[i][j]-rho1)/(rho2-rho1);
        }
    }
}

    
        
    










