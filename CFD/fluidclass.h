//
//  fluidclass.h
//  CFD
//
//  Created by YDuan on 9/30/16.
//  Copyright © 2016 Duan. All rights reserved.
//
#include "para.h"
#ifndef fluidclass_h
#define fluidclass_h
class fluidclass: public param{
protected:
    double rho1,rho2,m0,rro;
    double rad,xc,yc;
    double** r;
    double** u;double** v;double** p;
    double** ut;double** vt;double** tmp1;
    double** uu;double** vv;double** tmp2;
    double* x;double* y;
    double** rt;
    double** p_old;
    double** ro;

    //front tracking information//
    int Nf,Nfv;                     //
    double* xf;double* yf;      //front position
    double* uf;double* vf;      //front velocity
    double* tx;double* ty;      //tangent vectors
    double** fx; double** fy;   //gradient near front
    double** r_old;             //
    double* xfold;double* yfold;//
    //////////////////////////////
    
    
    //surface tension information//
    double m1,m2;                //viscosity
    double sigma;                //
    double** un;double** vn;     // second order???
    double** m;                  // viscosity field
    double** rn;double** mn;     //
    double* xfn;double* yfn;
    
public:
    void initialize_fluid();
    void timeloop();
    void calculate_boundary_velocity();
    void solve_velocity_without_pressure();
    void solve_pressure();
    int check_p_convg();
    void correct_velocity_after_pressure();
    void advect_density();
    void print_vtk();
    void print_dem();
    
    
    //front tracking call functions (called in time loop)//
    void advect_front();                                 //
    void add_points_to_front();                          //
    void distribute_gradient();                          //
    int check_r_convg();                                 //
    //void construct_density();                          //
    ///////////////////////////////////////////////////////
    
    //surface tension/////////////////////////////////////////
    void find_surface_tension();                            //
    void solve_velocity_adv_diff();
    void update_viscosity();

    
    ~fluidclass(){};
    
};

#endif /* fluidclass_h */
