#ifndef __NavierStokes2D__Multigrid__
#define __NavierStokes2D__Multigrid__

#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include "MultigridController.h"
#include "Shape2D.h"
using std::vector;
using std::invalid_argument;
using std::string;
using std::to_string;
using std::ofstream;
using std::endl;
using std::cout;

class Multigrid
{
public:
    Multigrid(double Lx, double Ly, vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub, bool leftIsDirichlet, bool topIsDirichlet, bool rightIsDirichlet, bool bottomIsDirichlet, MultigridController *MGC);
    Multigrid(double Lx, double Ly, vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub, bool leftIsDirichlet, bool topIsDirichlet, bool rightIsDirichlet, bool bottomIsDirichlet, size_t level_number, Multigrid *lp, MultigridController *MGC);
    ~Multigrid();
    void solve();
    void iterate();
    void smooth();
    void resetZero();
    void setBoundaries(vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub);
    void reset();
    void setRHS(vector<vector<double>> RHS_val);
    void print(string fileName);
    
    void maxIter(size_t val);
    unsigned int maxIter();
    void tol(double val);
    double tol();
    
    void setInternalBoundary(const Shape2D& obj);
    
    double elem(size_t i, size_t j);
    
protected:
    // Grid
    size_t Nx, Ny;
    double Lx, Ly;
    double dx, dy;
    vector<vector<double>> u;
    vector<vector<double>> R;
    
    // Coordinates
    vector<vector<double>> x, y;
    
    // Internal boundary
    vector<vector<double>> isSolid;
    
    // Error
    vector<vector<double>> e;
    double e_sum;
    double error();
    
    // Boundary condititions
    double leftIsDirichlet, topIsDirichlet, rightIsDirichlet, bottomIsDirichlet;
    vector<vector<double>> atLeftBoundary, atTopBoundary, atRightBoundary, atBottomBoundary;
    bool allNeumann;
    void subtractMean();
    
    // Points to controller and stack information
    MultigridController *controller;
    Multigrid *sublevel;
    Multigrid *superlevel;
    
    // Level number (bottom to top)
    size_t nlev;
    
    // Description
    void computeResidual();
    void restriction();
    void prolongation();
    
    //Iteration values
    unsigned int n = 0;                  // current iteration number
    
    // Retrieval of adjacent values
    inline double ddx_east_p(size_t i, size_t j);
    inline double ddy_north_p(size_t i, size_t j);
    inline double ddx_west_p(size_t i, size_t j);
    inline double ddy_south_p(size_t i, size_t j);
    inline double ddx_east_r(size_t i, size_t j);
    inline double ddy_north_r(size_t i, size_t j);
    inline double ddx_west_r(size_t i, size_t j);
    inline double ddy_south_r(size_t i, size_t j);
    inline double ddx_east(size_t i, size_t j);
    inline double ddy_north(size_t i, size_t j);
    inline double ddx_west(size_t i, size_t j);
    inline double ddy_south(size_t i, size_t j);

private:
    unsigned int _maxIter = 1E6;        // maximum number of iterations
    double _tol = 1E-6;                 // solution tolerance on error
    double _large_value = 1E100;        // large number to use for first error on reset
    
};

#endif /* defined(__NavierStokes2D__Multigrid__) */
