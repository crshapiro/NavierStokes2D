#ifndef __NavierStokes2D__NavierStokes2D__
#define __NavierStokes2D__NavierStokes2D__

#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include "TridiagonalMatrix.h"
#include "Multigrid.h"
#include "MultigridController.h"
#include "Shape2D.h"
using std::vector;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::ofstream;
using std::to_string;

class NavierStokes2D
{
public:
    NavierStokes2D(vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub, vector<double> vl, vector<double> vt, vector<double> vr, vector<double> vb, bool leftIsDirichlet, bool topIsDirichlet, bool rightIsDirichlet, bool bottomIsDirichlet, double H, double L, double Re, bool applyOutletBoundary);
    ~NavierStokes2D();
    double step();

    void setInternalBoundary(const Shape2D& obj);

    void print(string fileName);
    void displayInfo();
    void printGrid(string fileName);
    
    double get_u(const double& xi, const double& yi);
    double get_v(const double& xi, const double& yi);
    double get_p(const double& xi, const double& yi);
    
private:
    // Velocities
    vector<vector<double>> u, v;
    vector<vector<double>> ustar, vstar;
    
    // Grid
    vector<vector<double>> x, y;
    
    // Find nearest index
    size_t getNearestIndex_i(const double& xi);
    size_t getNearestIndex_j(const double &yi);

    // CFL condition
    inline double uMagnitude(size_t i, size_t j);
    void CFL();
    double dt;
    
    // Time
    double time;
    size_t n;
    
    // Internal boundary
    vector<vector<double>> isSolid;
    
    // Pressure
    Multigrid MG;
    inline double divergence(size_t i, size_t j);
    inline double ustar_e(size_t i, size_t j);
    inline double ustar_w(size_t i, size_t j);
    inline double vstar_n(size_t i, size_t j);
    inline double vstar_s(size_t i, size_t j);
    
    // Divergence
    double divergence();
    inline double ue(size_t i, size_t j);
    inline double uw(size_t i, size_t j);
    inline double vn(size_t i, size_t j);
    inline double vs(size_t i, size_t j);
    
    // Streamlines
    Multigrid psi;
    inline double ve(size_t i, size_t j);
    inline double vw(size_t i, size_t j);
    inline double un(size_t i, size_t j);
    inline double us(size_t i, size_t j);
    

    // Domain values
    size_t Nx, Ny;
    double H, L;
    double dx, dy;
    double Re;
    
    // Boundary condititions
    double leftIsDirichlet, topIsDirichlet, rightIsDirichlet, bottomIsDirichlet;
    vector<vector<double>> atLeftBoundary, atTopBoundary, atRightBoundary, atBottomBoundary;
    bool applyOutletBoundary;

    // Convective terms
    inline double uu_e(size_t i, size_t j);
    inline double uu_w(size_t i, size_t j);
    inline double uv_n(size_t i, size_t j);
    inline double uv_s(size_t i, size_t j);
    inline double uv_e(size_t i, size_t j);
    inline double uv_w(size_t i, size_t j);
    inline double vv_n(size_t i, size_t j);
    inline double vv_s(size_t i, size_t j);
    
    // Indexing
    inline size_t xVectorIndex(size_t i, size_t j);
    inline size_t xVectorIndex_i(size_t k);
    inline size_t xVectorIndex_j(size_t k);
    inline size_t yVectorIndex(size_t i, size_t j);
    inline size_t yVectorIndex_i(size_t k);
    inline size_t yVectorIndex_j(size_t k);
    
    // Approximate Factorization solutions by TDMA
    void yTDMA();
    void xTDMA();
    
    // Viscous term
    void Ly(size_t i, size_t j, double& as, double& ap, double& an, double &R, const double dt, const double& un, const double& us);
    void Lx(size_t i, size_t j, double& aw, double& ap, double& ae, double &R, const double dt, const double& ue, const double& uw);
    
};





#endif /* defined(__NavierStokes2D__NavierStokes2D__) */
