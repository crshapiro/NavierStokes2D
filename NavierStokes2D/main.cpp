#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include "NavierStokes2D.h"
#include "Multigrid.h"
#include "Shape2D.h"
using namespace std;

void strouhal(NavierStokes2D& NS, double L, double H, double A, double d, double time)
{
    ofstream f;
    f.open("output/strouhal.bin", std::ios::out | std::ios::binary | std::ios::app);
    f.write(reinterpret_cast<const char*>(&time), sizeof(double));
    double v = NS.get_v(d+A/2+2*A, H/2);
    f.write(reinterpret_cast<const char*>(&v), sizeof(double));
    f.close();
}

int main(int argc, const char * argv[])
{
    // Make output directory if needed
    system("mkdir output");
    
    
    // parse input arguments
    if (argc < 8)
    {
        cout << "Usage is NavierStokes2D Nx Ny L H Re nStep nDisp nPrint nPrintStart (shape) (d) (A)." << endl;
        cin.get();
        return -1;
    }
    double Nx = stod(argv[1]);
    double Ny = stod(argv[2]);
    double L = stod(argv[3]);
    double H = stod(argv[4]);
    double Re = stod(argv[5]);
    size_t nStep = stoi(argv[6]);
    size_t nDisp = stoi(argv[7]);
    size_t nPrint = stoi(argv[8]);
    size_t nPrintStart = stoi(argv[9]);

    string shape("none");
    double d(L/2), A(1);
    if (argc > 10)
    {
        if (argc < 13)
        {
            cout << "Usage is NavierStokes2D Nx Ny L H Re nStep nDisp nPrint nPrintStart (shape) (d) (A)." << endl;
            cout << "Shape, d, and A must be included together." << endl;
            cin.get();
            return -1;
        }
        shape = argv[10];
        d = stod(argv[11]);
        A = stod(argv[12]);
    }
    
    // Output input arguments
    cout << "Nx = " << Nx << endl;
    cout << "Ny = " << Ny << endl;
    cout << "L = " << L << endl;
    cout << "H = " << H << endl;
    cout << "Re = " << Re << endl;
    cout << "nStep = " << nStep << endl;
    cout << "nDisp = " << nDisp << endl;
    cout << "nPrint = " << nPrint << endl;
    cout << "nPrintStart = " << nPrintStart << endl;
    cout << "shape = " << shape << endl;
    cout << "d = " << d << endl;
    cout << "A = " << A << endl;
    
    // Allocate boundary conditions
    vector<double> ul(Ny,0.0), ut(Nx,0.0), ur(Ny,0.0), ub(Nx,0.0);
    vector<double> vl(Ny,0.0), vt(Nx,0.0), vr(Ny,0.0), vb(Nx,0.0);
    for (size_t j = 0; j < ul.size(); j++)
    {
        double dy = H/Ny;
        double y = j*dy + dy/2;
        ul[j] = 4*y/(H*H)*(H-y);
    }
    
    // Types of boundary conditions
    bool leftIsDirichlet(true), topIsDirichlet(true), rightIsDirichlet(false), bottomIsDirichlet(true);
    
    // Construct solver
    NavierStokes2D NS = NavierStokes2D(ul, ut, ur, ub, vl, vt, vr, vb, leftIsDirichlet, topIsDirichlet, rightIsDirichlet, bottomIsDirichlet, H, L, Re, true);
    
    // Apply internal boundary conditions
    if (shape == "square")
    {
        NS.setInternalBoundary(Square(d, H/2, A));
        cout << "Adding square to flow." << endl;
    }
    else if (shape == "circle")
    {
        NS.setInternalBoundary(Circle(d, H/2, A/2));
        cout << "Adding circle to flow." << endl;
    }
    else
    {
        cout << "No shape added to flow." << endl;
    }
    
    // Print initial grid information and initial conditions
    NS.printGrid("output/");
    NS.print("output/");
    
    cout << "Staring simulation...." << endl;
    // Simulate until number of steps is reached
    for (size_t n = 1; n <= nStep; n++)
    {
        strouhal(NS, L, H, A, d, NS.step());
        if (!(n % nDisp))                           NS.displayInfo();
        if (!(n % nPrint) && n > nPrintStart)       NS.print("output/");

    }

    cout << "Simulation complete." << endl;
    
    // Multigrid Poisson test with Neumann boundary conditions
    /*MultigridController *MGC;
    MGC = new Fcycle(Nx, Ny);
    
    vector<vector<double>> R(Nx+2, vector<double>(Ny+2, 0.0));
    double vol(0.0);
    for (size_t j = 1; j<= Ny; j++)
    {
        for (size_t i = 1; i <= Nx; i++)
        {
            double dx = L/Nx;
            double x = i*dx - dx/2;
            double dy = H/Ny;
            double y = j*dy - dy/2;
            R[i][j] = -2*cos(x)*cos(y);
            vol += R[i][j]*dx*dy;
        }
        cout << vol << endl;
    }
    
    vector<double> ul(Ny,0.0), ut(Nx,0.0), ur(Ny,0.0), ub(Nx,0.0);
    
    Multigrid MG(L, H, ul, ut, ur, ub, false, false, false, false, MGC);
    MG.setRHS(R);
    MG.reset();
    MG.solve();
    //.reset();
    MG.print("MultigridResults.out");*/
    
    return 0;
}
