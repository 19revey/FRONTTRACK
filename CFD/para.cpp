//
//  para.cpp
//  CFD
//
//  Created by YDuan on 9/30/16.
//  Copyright © 2016 Duan. All rights reserved.
//

#include <stdio.h>
#include "para.h"


param::param(){
    lx=1.0;ly=1.0;gx=0;gy=-100;
    unorth=0;usouth=0;veast=0;vwest=0;time=0;
    nx=64;ny=64;dt=0.000625;
    nstep=800;
    maxit=200;maxError=0.001;beta=1.2;
    print_interval=20;
    //zero various arrays
    
    
}

double* param::make1dmem(int arraySize, double initial_value){
    int i;
    double* theArray;
    theArray=new double[arraySize];
    for(i=0;i<arraySize;i++){
        theArray[i]=initial_value;
    }
    return theArray;
}

double** param::make2dmem(int arraySizeX, int arraySizeY, double initial_value){
    int i,j;
    double** theArray;
    theArray = new double*[arraySizeX];
    for(i=0;i<arraySizeX;i++){
        theArray[i]=new double[arraySizeY];
        for(j=0;j<arraySizeY;j++)
            theArray[i][j]=initial_value;
    }
    
    return theArray;
};

void param::copy1d(double* a, double* b, int n){
    for (int i=0; i<n; i++) {
        a[i]=b[i];
    }
}

void param::copy2d(double** a, double** b, int nx, int ny){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            a[i][j]=b[i][j];
        }
    }
}



