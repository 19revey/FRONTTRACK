//
//  main.cpp
//  CFD
//
//  Created by YDuan on 9/29/16.
//  Copyright © 2016 Duan. All rights reserved.
//

#include <iostream>
#include "para.h"
#include "fluidclass.h"



int main(int argc, const char * argv[]) {
    // insert code here...
    fluidclass A;
    A.initialize_fluid();
    A.timeloop();
    //A.print_vtk();
    std::cout << "Hello, World!\n";
    
    return 0;
}
