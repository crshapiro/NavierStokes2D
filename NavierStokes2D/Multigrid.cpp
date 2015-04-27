#include "Multigrid.h"
#include <iostream>

Multigrid::Multigrid(double Lx, double Ly, vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub, bool leftIsDirichlet, bool topIsDirichlet, bool rightIsDirichlet, bool bottomIsDirichlet, MultigridController *MGC) :  Multigrid::Multigrid(Lx, Ly, ul, ut, ur, ub, leftIsDirichlet, topIsDirichlet, rightIsDirichlet, bottomIsDirichlet, MGC->numberOfLevels(), nullptr, MGC)
{ return; }

Multigrid::Multigrid(double Lx, double Ly, vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub, bool leftIsDirichlet, bool topIsDirichlet, bool rightIsDirichlet, bool bottomIsDirichlet,  size_t level_number, Multigrid *lp, MultigridController *MGC)
{
    
    // Check that N is a power of 2
    if (Nx % 2 != 0 || Nx % 2 != 0)
        throw invalid_argument("Nx and Ny must be a power of 2");
    
    // Assign inputs
    nlev = level_number - 1;
    controller = MGC;
    superlevel = lp;
    
    // Assign domain values
    Nx = ub.size();
    Ny = ul.size();
    this->Lx = Lx;
    this->Ly = Ly;
    dx = Lx/Nx;
    dy = Ly/Ny;
    
    // Boundary condition booleans
    this->leftIsDirichlet = (double)leftIsDirichlet;
    this->topIsDirichlet = (double)topIsDirichlet;
    this->rightIsDirichlet = (double)rightIsDirichlet;
    this->bottomIsDirichlet = (double)bottomIsDirichlet;
    
    if (!leftIsDirichlet && !topIsDirichlet && !rightIsDirichlet && !bottomIsDirichlet)
        allNeumann = true;
    else
        allNeumann = false;
    
    // Boundary flags
    atLeftBoundary = vector<vector<double>>(Nx+2,vector<double>(Ny+2, 0.0));
    atTopBoundary = vector<vector<double>>(Nx+2,vector<double>(Ny+2, 0.0));
    atRightBoundary = vector<vector<double>>(Nx+2,vector<double>(Ny+2, 0.0));
    atBottomBoundary = vector<vector<double>>(Nx+2,vector<double>(Ny+2, 0.0));
    for (size_t i = 1; i < Nx+1; i++)
    {
        atTopBoundary[i][Ny] = 1.0;
        atBottomBoundary[i][1] = 1.0;
    }
    for (size_t j = 1; j < Ny+1; j++)
    {
        atLeftBoundary[1][j] = 1.0;
        atRightBoundary[Nx][j] = 1.0;
    }
    
    // Construct grid with buffer lines
    u = vector<vector<double>>(Nx+2);
    for (size_t i = 0; i < Nx+2; i++)
        u[i] = vector<double>(Ny+2,0.0);
    
    setBoundaries(ul, ut, ur, ub);
    
    // Allocate right hand side and error term with same same size as u
    R = vector<vector<double>>(Nx+2);
    e = vector<vector<double>>(Nx+2);
    for (size_t i = 0; i < Nx+2; i++)
    {
        R[i] = vector<double>(Ny+2,0.0);
        e[i] = vector<double>(Ny+2,0.0);
    }
    
    // Generate coordinate values
    x = vector<vector<double>>(Nx+2, vector<double>(Ny+2));
    y = vector<vector<double>>(Nx+2, vector<double>(Ny+2));
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j<= Ny; j++)
        {
            x[i][j] = i*dx - dx/2;
            y[i][j] = j*dy - dy/2;
        }
    
    // Internal boundary
    isSolid = vector<vector<double>>(Nx+2, vector<double>(Ny+2, 0.0));
    
    // Create sublevel
    if (nlev == 0)
        sublevel = nullptr;
    else
        sublevel = new Multigrid(Lx, Ly, vector<double>(Ny/2, 0.0), vector<double>(Nx/2, 0.0), vector<double>(Ny/2, 0.0), vector<double>(Nx/2, 0.0), leftIsDirichlet, topIsDirichlet, rightIsDirichlet, bottomIsDirichlet, nlev, this, MGC);

}

Multigrid::~Multigrid()
{ delete sublevel; }

double Multigrid::elem(size_t i, size_t j)
{
    return u[i][j];
}

void Multigrid::setInternalBoundary(const Shape2D& obj)
{
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
            if (obj.inVolume(x[i][j], y[i][j]))
                isSolid[i][j] = 1.0;
    //if (sublevel != nullptr)
        //sublevel->setInternalBoundary(obj);
}

void Multigrid::setBoundaries(vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub)
{
    
    // Check that sizes are consistent
    if ( ul.size() != ur.size() )
        throw invalid_argument("boundary condition vectors must be the same size");
    if ( ut.size() != ub.size() )
        throw invalid_argument("boundary condition vectors must be the same size");
    
    // Store boundary conditions in buffer lines
    for (size_t i = 1; i < Nx+1; i++)
    {
        u[i][0] = ub[i-1];
        u[i][Ny+1] = ut[i-1];
    }
    for (size_t j = 1; j < Ny+1; j++)
    {
        u[0][j] = ul[j-1];
        u[Nx+1][j] = ur[j-1];
    }
}

void Multigrid::print(string fileName)
{
    ofstream f;
    f.open(fileName);
    for (size_t j = 1; j <= Ny; j++)
    {
        for (size_t i = 1; i < Nx; i++)
            f << u[i][j] << "\t";
        f << u[Nx][j] << endl;
    }
    f.close();
}

void Multigrid::setRHS(vector<vector<double>> RHS_val)
{
    // Really need to check that boundaries are 0.
    if (RHS_val.size() != Nx + 2 )
        throw invalid_argument("RHS must have Nx + 2 by Ny + 2 elements");
    if (RHS_val[1].size() != Ny + 2 )
        throw invalid_argument("RHS must have Nx + 2 by Ny + 2 elements");
    
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
            R[i][j] = RHS_val[i][j];
    
}

void Multigrid::reset()
{
    n = 0;
    e_sum = _large_value;
    if (sublevel != nullptr)
        sublevel->reset();
}

void Multigrid::solve()
{
    smooth();
    computeResidual();
    while (error() > tol() & n < maxIter())
        iterate();
    if (n == maxIter())
        cout << "Multigrid reached maximum number of iterations" << endl;
}

double Multigrid::error()
{ return e_sum; }

void Multigrid::iterate()
{
    int short o = controller->getCommand();
    if (o < 0) {
        if (superlevel != nullptr) smooth();
        sublevel->resetZero();
        restriction();
        sublevel->iterate();
    } else if (o == 0) {
        smooth();
    } else{
        if (superlevel != nullptr) {
            smooth();
            prolongation();
            superlevel->iterate();
        } else {
            smooth();
            if (allNeumann) subtractMean();
            controller->reset();
        }
    }
    //cout << "Hello" << endl;
}

void Multigrid::subtractMean()
{
    double um = 0;
    double N = 0;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            um += u[i][j] * (1.0 - isSolid[i][j]);
            N += 1.0 - isSolid[i][j];
        }
    
    um/= N;
    
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
            u[i][j] -= um * (1.0 - isSolid[i][j]);
}

void Multigrid::resetZero()
{
    for (size_t i = 1; i <= Nx; i++)
        for(size_t j = 1; j <= Ny; j++)
            u[i][j] = 0;
}

void Multigrid::smooth()
{
    double a1, a2, a3, a4;
    double b1, b2, b3, b4;
    double a, b;
    for (int i = 1; i <= Nx; i++)
        for (int j = 1; j <= Ny; j++)
        {
            a1 = ddx_east_p(i, j);
            a2 = ddx_west_p(i, j);
            a3 = ddy_north_p(i, j);
            a4 = ddy_south_p(i, j);
            b1 = ddx_east_r(i, j);
            b2 = ddx_west_r(i, j);
            b3 = ddy_north_r(i, j);
            b4 = ddy_south_r(i, j);
            a = (a1 - a2)/dx + (a3 - a4)/dy + isSolid[i][j];
            b = (b1 - b2)/dx + (b3 - b4)/dy;
            u[i][j] = (R[i][j]/a - b/a)*(1.0 - isSolid[i][j]);
        }
    
    
    computeResidual();
    n++;
    return;
}

void Multigrid::maxIter(size_t val)
{ _maxIter = (int)val; }

unsigned int Multigrid::maxIter()
{ return _maxIter; }

void Multigrid::tol(double val)
{
    if (val > 0)
        _tol = val;
    else
        throw invalid_argument("tolerance must be positive");
}

double Multigrid::tol()
{ return _tol; }

void Multigrid::restriction()
{
    for (size_t i = 1; i <= Nx/2; i++)
        for (size_t j = 1; j <= Ny/2; j++)
            sublevel->R[i][j] = -0.25*(e[2*i-1][2*j-1] + e[2*i-1][2*j] + e[2*i][2*j-1] + e[2*i][2*j])*(1.0 - sublevel->isSolid[i][j]);
}

void Multigrid::prolongation()
{
    // Correction term for superlevel
    vector<vector<double>> uc(superlevel->Nx + 2, vector<double>(superlevel->Ny + 2, 0.0));
    
    // Bi-linearly interpolate lower level to higher level.
    // This is simplified because prolongation is only used for lower levels with zero boundary conditions.
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            double ut, ub;
            
            // Lower left
            ut = ((0.75*u[i][j] + 0.25*u[i-1][j])*(1.0 - atLeftBoundary[i][j])
            + 0.5*u[i][j]*atLeftBoundary[i][j]*leftIsDirichlet
            + u[i][j]*atLeftBoundary[i][j]*(1.0-leftIsDirichlet))
            * (1.0 - isSolid[i-1][j]) + u[i][j]*isSolid[i-1][j];
            
            ub = ((0.75*u[i][j-1] + 0.25*u[i-1][j-1])*(1.0 - atLeftBoundary[i][j])
            + 0.5*u[i][j-1]*atLeftBoundary[i][j]*leftIsDirichlet
            + u[i][j-1]*atLeftBoundary[i][j]*(1.0-leftIsDirichlet))
            * (1.0 - isSolid[i-1][j]) + u[i][j-1]*isSolid[i-1][j];
            
            uc[2*i-1][2*j-1] = ((0.75*ut + 0.25*ub)*(1.0 - atBottomBoundary[i][j])
            + 0.5*ut*atBottomBoundary[i][j]*bottomIsDirichlet
            + ut*atBottomBoundary[i][j]*(1.0 - bottomIsDirichlet))
            * (1.0 - isSolid[i][j-1]) + ut*isSolid[i][j-1];
            
            // Upper left
            ut = ((0.75*u[i][j+1] + 0.25*u[i-1][j+1])*(1.0 - atLeftBoundary[i][j])
            + 0.5*u[i][j+1]*atLeftBoundary[i][j]*leftIsDirichlet
            + u[i][j+1]*atLeftBoundary[i][j]*(1.0 - leftIsDirichlet))
            * (1.0 - isSolid[i-1][j]) + u[i][j+1]*isSolid[i-1][j];
            
            ub = ((0.75*u[i][j] + 0.25*u[i-1][j])*(1.0 - atLeftBoundary[i][j])
            + 0.5*u[i][j]*atLeftBoundary[i][j]*leftIsDirichlet
            + u[i][j]*atLeftBoundary[i][j]*(1.0 - leftIsDirichlet))
            * (1.0 - isSolid[i-1][j]) + u[i][j]*isSolid[i-1][j];
            
            uc[2*i-1][2*j] = ((0.75*ub + 0.25*ut)*(1.0 - atTopBoundary[i][j])
            + 0.5*ub*atTopBoundary[i][j]*topIsDirichlet
            + ub*atTopBoundary[i][j]*(1.0 - topIsDirichlet))
            * (1.0 - isSolid[i][j+1]) + ub*isSolid[i][j+1];
            
            // Lower right
            ut = ((0.75*u[i][j] + 0.25*u[i+1][j])*(1.0 - atRightBoundary[i][j])
            + 0.5*u[i][j]*atRightBoundary[i][j]*rightIsDirichlet
            + u[i][j]*atRightBoundary[i][j]*(1.0 - rightIsDirichlet))
            * (1.0 - isSolid[i+1][j]) + u[i][j]*isSolid[i-1][j];
            
            ub = ((0.75*u[i][j-1] + 0.25*u[i+1][j-1])*(1.0 - atRightBoundary[i][j])
            + 0.5*u[i][j-1]*atRightBoundary[i][j]*rightIsDirichlet
            + u[i][j-1]*atRightBoundary[i][j]*(1.0 - rightIsDirichlet))
            * (1.0 - isSolid[i+1][j]) + u[i][j-1]*isSolid[i-1][j];
            
            uc[2*i][2*j-1] = ((0.75*ut + 0.25*ub)*(1.0 - atBottomBoundary[i][j])
            + 0.5*ut*atBottomBoundary[i][j]*bottomIsDirichlet
            + ut*atBottomBoundary[i][j]*(1.0 - bottomIsDirichlet))
            * (1.0 - isSolid[i][j-1]) + ut*isSolid[i][j-1];
            
            // Upper right
            ut = ((0.75*u[i][j+1] + 0.25*u[i+1][j+1])*(1.0 - atRightBoundary[i][j])
            + 0.5*u[i][j+1]*atRightBoundary[i][j]*rightIsDirichlet
            + u[i][j+1]*atRightBoundary[i][j]*(1.0 - rightIsDirichlet))
            * (1.0 - isSolid[i+1][j]) + u[i][j+1]*isSolid[i-1][j];
            
            ub = ((0.75*u[i][j] + 0.25*u[i+1][j])*(1.0 - atRightBoundary[i][j])
            + 0.5*u[i][j]*atRightBoundary[i][j]*rightIsDirichlet
            + u[i][j]*atRightBoundary[i][j]*(1.0 - rightIsDirichlet))
            * (1.0 - isSolid[i+1][j]) + u[i][j]*isSolid[i-1][j];
            
            uc[2*i][2*j] = ((0.75*ub + 0.25*ut)*(1.0 - atTopBoundary[i][j])
            + 0.5*ub*atTopBoundary[i][j]*topIsDirichlet
            + ub*atTopBoundary[i][j]*(1.0 - topIsDirichlet))
            * (1.0 - isSolid[i][j+1]) + ub*isSolid[i][j+1];
        }
    
    // Inject to higher level
    for (int i = 1; i <= 2*Nx; i++)
        for (int j = 1; j <= 2*Ny; j++)
            superlevel->u[i][j] += uc[i][j]*(1.0 - superlevel->isSolid[i][j]);
}

void Multigrid::computeResidual()
{
    e_sum = 0.0;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            e[i][j] = (ddx_east(i, j) - ddx_west(i, j))/dx + (ddy_north(i, j) - ddy_south(i, j))/dy - R[i][j];
            e[i][j] *= 1.0 - isSolid[i][j];
            e_sum += e[i][j]*e[i][j];
        }
    
    e_sum /= Nx;
    e_sum /= Ny;
    e_sum = sqrt(e_sum);
    return;
}

inline double Multigrid::ddx_east_p(size_t i, size_t j)
{
    return (-1.0/dx*(1.0-atRightBoundary[i][j])
    + -0.5/dx*rightIsDirichlet*atRightBoundary[i][j])*(1.0 - isSolid[i+1][j]);
    //+ u[i+1][j]*(1.0-rightIsDirichlet)*atRightBoundary[i][j];
}

inline double Multigrid::ddy_north_p(size_t i, size_t j)
{
    return (-1.0/dy*(1.0-atTopBoundary[i][j])
    + -0.5/dy*topIsDirichlet*atTopBoundary[i][j])*(1.0 - isSolid[i][j+1]);
    //+ u[i][j+1]*(1.0-topIsDirichlet)*atTopBoundary[i][j];
}

inline double Multigrid::ddx_west_p(size_t i, size_t j)
{
    return (1.0/dx*(1.0-atLeftBoundary[i][j])
    + 0.5/dx*leftIsDirichlet*atLeftBoundary[i][j])*(1.0 - isSolid[i-1][j]);
    //+ u[i-1][j]*(1.0-leftIsDirichlet)*atLeftBoundary[i][j];
}

inline double Multigrid::ddy_south_p(size_t i, size_t j)
{
    return (1.0/dy*(1.0-atBottomBoundary[i][j])
    + 0.5/dy*bottomIsDirichlet*atBottomBoundary[i][j])*(1.0 - isSolid[i][j-1]);
    //+ u[i][j-1]*(1.0-bottomIsDirichlet)*atBottomBoundary[i][j];
}

inline double Multigrid::ddx_east_r(size_t i, size_t j)
{
    return (u[i+1][j]/dx*(1.0-atRightBoundary[i][j])
    + 0.5*u[i+1][j]/dx*rightIsDirichlet*atRightBoundary[i][j]
    + u[i+1][j]*(1.0-rightIsDirichlet)*atRightBoundary[i][j])*(1.0 - isSolid[i+1][j]);
}

inline double Multigrid::ddy_north_r(size_t i, size_t j)
{
    return (u[i][j+1]/dy*(1.0-atTopBoundary[i][j])
    + 0.5*u[i][j+1]/dy*topIsDirichlet*atTopBoundary[i][j]
    + u[i][j+1]*(1.0-topIsDirichlet)*atTopBoundary[i][j])*(1.0 - isSolid[i][j+1]);
}

inline double Multigrid::ddx_west_r(size_t i, size_t j)
{
    return (-u[i-1][j]/dx*(1.0-atLeftBoundary[i][j])
    + -0.5*u[i-1][j]/dx*leftIsDirichlet*atLeftBoundary[i][j]
    + u[i-1][j]*(1.0-leftIsDirichlet)*atLeftBoundary[i][j])*(1.0 - isSolid[i-1][j]);
}

inline double Multigrid::ddy_south_r(size_t i, size_t j)
{
    return (-u[i][j-1]/dy*(1.0-atBottomBoundary[i][j])
    + -0.5*u[i][j-1]/dy*bottomIsDirichlet*atBottomBoundary[i][j]
    + u[i][j-1]*(1.0-bottomIsDirichlet)*atBottomBoundary[i][j])*(1.0 - isSolid[i][j-1]);
}

inline double Multigrid::ddx_east(size_t i, size_t j)
{
    return (1.0/dx*(u[i+1][j] - u[i][j])*(1.0-atRightBoundary[i][j])
    + 0.5/dx*(u[i+1][j] - u[i][j])*atRightBoundary[i][j]*rightIsDirichlet
    + u[i+1][j]*atRightBoundary[i][j]*(1.0-rightIsDirichlet)) * (1.0 - isSolid[i+1][j]);
}

inline double Multigrid::ddy_north(size_t i, size_t j)
{
    return (1.0/dy*(u[i][j+1] - u[i][j])*(1.0-atTopBoundary[i][j])
    + 0.5/dy*(u[i][j+1] - u[i][j])*atTopBoundary[i][j]*topIsDirichlet
    + u[i][j+1]*atTopBoundary[i][j]*(1.0-topIsDirichlet)) * (1.0 - isSolid[i][j+1]);
}

inline double Multigrid::ddx_west(size_t i, size_t j)
{
    return (1.0/dx*(u[i][j] - u[i-1][j])*(1.0-atLeftBoundary[i][j])
    + 0.5/dx*(u[i][j] - u[i-1][j])*atLeftBoundary[i][j]*leftIsDirichlet
    + u[i-1][j]*atLeftBoundary[i][j]*(1.0-leftIsDirichlet)) * (1.0 - isSolid[i-1][j]);
}

inline double Multigrid::ddy_south(size_t i, size_t j)
{
    return (1.0/dy*(u[i][j] - u[i][j-1])*(1.0-atBottomBoundary[i][j])
    + 0.5/dy*(u[i][j] - u[i][j-1])*atBottomBoundary[i][j]*bottomIsDirichlet
    + u[i][j-1]*atBottomBoundary[i][j]*(1.0-bottomIsDirichlet)) * (1.0 - isSolid[i][j-1]);
}

