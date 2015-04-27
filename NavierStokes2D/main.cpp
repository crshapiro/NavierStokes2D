#include <iostream>
#include <vector>
#include <math.h>
#include "NavierStokes2D.h"
#include "Multigrid.h"
#include "Shape2D.h"
using namespace std;

int main(int argc, const char * argv[]) {
    
    // parse input arguments
    if (argc < 6)
    {
        cout << "Usage is NavierStokes2D Nx Ny Re nsteps shape";
        cin.get();
        return -1;
    }
    double Nx = stod(argv[1]);
    double Ny = stod(argv[2]);
    double Re = stod(argv[3]);
    size_t nsteps = stoi(argv[4]);
    string shape = argv[5];
    
    // domain values
    double H = 8.0;
    double L = 16.0;
    
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
    if (shape == "square") NS.setInternalBoundary(Square(4, 4, 1));
    if (shape == "circle") NS.setInternalBoundary(Circle(4, 4, 0.5));
    NS.printGrid();
    NS.print();
    for (size_t n = 1; n <= nsteps; n++)
    {
        NS.step();
        if ( !(n % 100) )
        {
            NS.displayInfo();
            NS.print();
        }
    }
    NS.step();
    NS.step();

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
