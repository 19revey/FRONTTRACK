//
//  para.h
//  CFD
//
//  Created by YDuan on 9/29/16.
//  Copyright © 2016 Duan. All rights reserved.
//

#ifndef para_h
#define para_h


class param{
protected:
    //domain size
    double lx;
    double ly;
    double gx,gy;
    double dx,dy;
    double unorth,usouth,veast,vwest;
    double time;
    //numerical variables
    int nx;
    int ny;
    double dt;
    int nstep;
    int maxit;
    double maxError;
    double beta;
    int print_interval;
    //arrays

    
public:
    double* make1dmem(int arraySize, double initial_value);
    double** make2dmem(int arraySizeX, int arraySizeY,double initial_value);
    param();
    double square(double a){return a*a;}
    void copy1d(double* a,double* b,int n);
    void copy2d(double** a,double** b,int nx,int ny);

    ~param(){};
};




#endif /* para_h */
