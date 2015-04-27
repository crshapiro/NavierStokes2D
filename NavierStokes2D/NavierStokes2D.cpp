#include "NavierStokes2D.h"

NavierStokes2D::NavierStokes2D(vector<double> ul, vector<double> ut, vector<double> ur, vector<double> ub, vector<double> vl, vector<double> vt, vector<double> vr, vector<double> vb, bool leftIsDirichlet, bool topIsDirichlet, bool rightIsDirichlet, bool bottomIsDirichlet, double H, double L, double Re, bool applyOutletBoundary) : MG(L, H, ul, ut, ur, ub, false, false, false, false, new Fcycle(ub.size(), ul.size()))
{
    // Check that sizes are consistent
    if ( ul.size() != ur.size() || vl.size() != vr.size() || ul.size() != vl.size() )
        throw invalid_argument("boundary condition vectors must be the same size");
    if ( ut.size() != ub.size() || vt.size() != vb.size() || ut.size() != vb.size() )
        throw invalid_argument("boundary condition vectors must be the same size");
    
    // Assign domain values
    Nx = ub.size();
    Ny = ul.size();
    this->H = H;
    this->L = L;
    dx = L/Nx;
    dy = H/Ny;
    this->Re = Re;
    
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

    // Boundary condition booleans
    this->leftIsDirichlet = (double)leftIsDirichlet;
    this->topIsDirichlet = (double)topIsDirichlet;
    this->rightIsDirichlet = (double)rightIsDirichlet;
    this->bottomIsDirichlet = (double)bottomIsDirichlet;
    this->applyOutletBoundary = applyOutletBoundary;
    
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
    v = vector<vector<double>>(Nx+2);
    for (size_t i = 0; i < Nx+2; i++)
    {
        u[i] = vector<double>(Ny+2,0.0);
        v[i] = vector<double>(Ny+2,0.0);
    }

    // Store boundary conditions in buffer lines
    for (size_t i = 1; i < Nx+1; i++)
    {
        u[i][0] = ub[i-1];
        v[i][0] = vb[i-1];
        u[i][Ny+1] = ut[i-1];
        v[i][Ny+1] = vt[i-1];
    }
    for (size_t j = 1; j < Ny+1; j++)
    {
        u[0][j] = ul[j-1];
        v[0][j] = vl[j-1];
        if (applyOutletBoundary)
            u[Nx+1][j] = ul[j-1];
        else
            u[Nx+1][j] = ur[j-1];
        v[Nx+1][j] = vr[j-1];
    }
    
    // Make right boundary dirichlet if necessary
    if (applyOutletBoundary)
        this->rightIsDirichlet = 1.0;
    
    //Copy to ustar and vstar
    ustar = u;
    vstar = v;
    
    // Loosen convergence for pressure
    MG.tol(1E-6);
    MG.maxIter(100*Nx*Ny);
    
    // Set time and step number
    n = 0;
    time = 0;
    
    return;
}

NavierStokes2D::~NavierStokes2D()
{ return; }

void NavierStokes2D::displayInfo()
{
    cout << "Step " << n << "\t" << "t = " << time << endl;
    cout << "Total divergence = " << divergence() << endl;
}

double NavierStokes2D::divergence()
{
    double div = 0;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
            div += (((ue(i,j) - uw(i, j))/dx + (vn(i, j) - vs(i, j))/dy))*dx*dy;
    return div;
}

double NavierStokes2D::dt_CFL()
{
    double dt = 0.4*max(dy, dx);
    double U = 0;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
            U = max(U, uMagnitude(i, j));
    

    for (size_t j = 1; j <= Ny; j++)
    {
        if (leftIsDirichlet)
            U = max(U, uMagnitude(0, j));
        if (rightIsDirichlet)
            U = max(U, uMagnitude(Nx, j));
    }
    
    for (size_t i = 1; i <= Nx; i++)
    {
        if (bottomIsDirichlet)
            U = max(U, uMagnitude(i, 0));
        if (topIsDirichlet)
            U = max(U, uMagnitude(i, Ny));
    }
    
    return dt/U;
}

inline double NavierStokes2D::uMagnitude(size_t i, size_t j)
{ return sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]); }

double NavierStokes2D::step()
{
    // Some Parameters
    double dt = dt_CFL();
    double dx_inv = 1.0/dx;
    double dy_inv = 1.0/dy;
    
    // Convective right hand side stored in ustar
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            ustar[i][j] = u[i][j] - dt*(dx_inv*(uu_e(i, j) - uu_w(i, j)) + dy_inv*(uv_n(i, j) - uv_s(i, j)) );
            vstar[i][j] = v[i][j] - dt*(dx_inv*(uv_e(i, j) - uv_w(i, j)) + dy_inv*(vv_n(i, j) - vv_s(i, j)) );
        }
    
    // Do approximate factorization TDMA
    yTDMA(dt);
    xTDMA(dt);
    
    // Pressure Poisson equation RHS
    vector<vector<double>> R(Nx+2);
    double div;
    for (size_t i = 0; i < Nx+2; i++)
        R[i] = vector<double>(Ny+2, 0.0);
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            R[i][j] = ((ustar_e(i,j) - ustar_w(i, j))/dx + (vstar_n(i, j) - vstar_s(i, j))/dy);
            div += R[i][j]*dx*dy;
        }
    
    // Solve for pressure
    MG.setBoundaries(vector<double>(Ny, 0.0), vector<double>(Nx, 0.0), vector<double>(Ny, 0.0), vector<double>(Nx, 0.0));
    MG.setRHS(R);
    MG.reset();
    MG.solve();
    
    vector<vector<double>> dpx(Nx+2), dpy(Nx+2);
    for (size_t i = 0; i < Nx + 2; i++)
    {
        dpx[i] = vector<double>(Ny+2);
        dpy[i] = vector<double>(Ny+2);
    }
    
    double pw, pe, pn, ps;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            pe = 0.5*(MG.elem(i+1, j) + MG.elem(i, j))*(1.0 - atRightBoundary[i][j]) + MG.elem(i, j)*atRightBoundary[i][j];
            pw = 0.5*(MG.elem(i, j) + MG.elem(i-1, j))*(1.0 - atLeftBoundary[i][j]) + MG.elem(i, j)*atLeftBoundary[i][j];
            pn = 0.5*(MG.elem(i, j+1) + MG.elem(i, j))*(1.0 - atTopBoundary[i][j]) + MG.elem(i, j)*atTopBoundary[i][j];
            ps = 0.5*(MG.elem(i, j) + MG.elem(i, j-1))*(1.0 - atBottomBoundary[i][j]) + MG.elem(i, j)*atBottomBoundary[i][j];
            pe *= (1.0 - isSolid[i+1][j]);
            pw *= (1.0 - isSolid[i-1][j]);
            pn *= (1.0 - isSolid[i][j+1]);
            ps *= (1.0 - isSolid[i][j-1]);
            pe += isSolid[i+1][j]*MG.elem(i, j);
            pw += isSolid[i-1][j]*MG.elem(i, j);;
            pn += isSolid[i][j+1]*MG.elem(i, j);;
            ps += isSolid[i][j-1]*MG.elem(i, j);;
            dpx[i][j] = (pe - pw)/dx;
            dpy[i][j] = (pn - ps)/dy;
            u[i][j] = ustar[i][j] - dpx[i][j]*(1.0 - isSolid[i][j]);
            v[i][j] = vstar[i][j] - dpy[i][j]*(1.0 - isSolid[i][j]);
        }
    
    // Fix boundaries for u and ustar
    if (applyOutletBoundary)
    {
        double massOut = 0.0;
        double massIn = 0.0;
        for (size_t j = 1; j <= Ny; j++)
        {
            massIn += u[0][j];
            massOut += u[Nx][j];
        }
        double massDiff = (massOut - massIn)/Ny;
        for (size_t j = 1; j <= Ny; j++)
        {
            u[Nx][j] -= massDiff;
            u[Nx+1][j] = u[Nx][j];
            v[Nx+1][j] = v[Nx][j];
            ustar[Nx+1][j] = u[Nx][j];
            vstar[Nx+1][j] = v[Nx][j];
        }
    }
    
    // Advance time and return
    n++;
    time += dt;
    return time;
}

void NavierStokes2D::setInternalBoundary(const Shape2D& obj)
{
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
            if (obj.inVolume(x[i][j], y[i][j]))
                isSolid[i][j] = 1.0;
    MG.setInternalBoundary(obj);
}

inline double NavierStokes2D::ustar_e(size_t i, size_t j)
{
    return (0.5*(ustar[i+1][j] + ustar[i][j])*(1.0 - atRightBoundary[i][j])
    + ustar[i+1][j]*rightIsDirichlet*atRightBoundary[i][j]
    + (0.5*dx*ustar[i+1][j] + ustar[i][j])*(1.0-rightIsDirichlet)*atRightBoundary[i][j])
    * (1.0 - isSolid[i+1][j]) * (1.0 - isSolid[i][j]);
}

inline double NavierStokes2D::ustar_w(size_t i, size_t j)
{
    return (0.5*(ustar[i-1][j] + ustar[i][j])*(1.0 - atLeftBoundary[i][j])
    + ustar[i-1][j]*atLeftBoundary[i][j]*leftIsDirichlet
    + (-0.5*dx*ustar[i-1][j] + ustar[i][j])*(1.0-leftIsDirichlet)*atLeftBoundary[i][j])
    * (1.0 - isSolid[i-1][j]) * (1.0 - isSolid[i][j]);
}

inline double NavierStokes2D::vstar_n(size_t i, size_t j)
{
    return (0.5*(vstar[i][j+1] + vstar[i][j])*(1.0 - atTopBoundary[i][j])
    + vstar[i][j+1]*atTopBoundary[i][j]*topIsDirichlet
    + (0.5*dy*vstar[i][j+1] + vstar[i][j])*(1.0-topIsDirichlet)*atTopBoundary[i][j])
    * (1.0 - isSolid[i][j+1]) * (1.0 - isSolid[i][j]);
}

inline double NavierStokes2D::vstar_s(size_t i, size_t j)
{
    return (0.5*(vstar[i][j-1] +  vstar[i][j])*(1.0 - atBottomBoundary[i][j])
    + vstar[i][j-1]*atBottomBoundary[i][j]*bottomIsDirichlet
    + (-0.5*dy*vstar[i][j-1] + vstar[i][j])*(1.0-bottomIsDirichlet)*atBottomBoundary[i][j])
    * (1.0 - isSolid[i][j-1]) * (1.0 - isSolid[i][j]);
}

inline double NavierStokes2D::ue(size_t i, size_t j)
{
    return (0.5*(u[i+1][j] + u[i][j])*(1.0 - atRightBoundary[i][j])
    + u[i+1][j]*rightIsDirichlet*atRightBoundary[i][j]
    + (0.5*dx*u[i+1][j] + u[i][j])*(1.0-rightIsDirichlet)*atRightBoundary[i][j])
    * (1.0 - isSolid[i+1][j]) * (1.0 - isSolid[i][j]);
}

inline double NavierStokes2D::uw(size_t i, size_t j)
{
    return (0.5*(u[i-1][j] + u[i][j])*(1.0 - atLeftBoundary[i][j])
    + u[i-1][j]*atLeftBoundary[i][j]*leftIsDirichlet
    + (-0.5*dx*u[i-1][j] + u[i][j])*(1.0-leftIsDirichlet)*atLeftBoundary[i][j])
    * (1.0 - isSolid[i-1][j]) * (1.0 - isSolid[i][j]);
}

inline double NavierStokes2D::vn(size_t i, size_t j)
{
    return (0.5*(v[i][j+1] + v[i][j])*(1.0 - atTopBoundary[i][j])
    + v[i][j+1]*atTopBoundary[i][j]*topIsDirichlet
    + (0.5*dy*v[i][j+1] + v[i][j])*(1.0-topIsDirichlet)*atTopBoundary[i][j])
    * (1.0 - isSolid[i][j+1]) * (1.0 - isSolid[i][j]);
}

inline double NavierStokes2D::vs(size_t i, size_t j)
{
    return (0.5*(v[i][j-1] +  v[i][j])*(1.0 - atBottomBoundary[i][j])
    + v[i][j-1]*atBottomBoundary[i][j]*bottomIsDirichlet
    + (-0.5*dy*v[i][j-1] + v[i][j])*(1.0-bottomIsDirichlet)*atBottomBoundary[i][j])
    * (1.0 - isSolid[i][j-1]) * (1.0 - isSolid[i][j]);
}

void NavierStokes2D::print()
{
    ofstream fu, fv;
    fu.open("u_" + to_string(n) + ".bin", std::ios::out | std::ios::binary);
    fv.open("v_" + to_string(n) + ".bin", std::ios::out | std::ios::binary);
    for (size_t j = 1; j <= Ny; j++)
    {
        for (size_t i = 1; i < Nx; i++)
        {
            fu.write(reinterpret_cast<const char*>(&u[i][j]), sizeof(double));
            fv.write(reinterpret_cast<const char*>(&v[i][j]), sizeof(double));
        }
        fu.write(reinterpret_cast<const char*>(&u[Nx][j]), sizeof(double));
        fv.write(reinterpret_cast<const char*>(&v[Nx][j]), sizeof(double));
    }
    fu.close();
    fv.close();
}

void NavierStokes2D::printGrid()
{
    ofstream fx, fy, fb;
    fx.open("x.bin", std::ios::out | std::ios::binary);
    fy.open("y.bin", std::ios::out | std::ios::binary);
    fb.open("obj.bin", std::ios::out | std::ios::binary);
    for (size_t j = 1; j <= Ny; j++)
    {
        for (size_t i = 1; i < Nx; i++)
        {
            fx.write(reinterpret_cast<const char*>(&x[i][j]), sizeof(double));
            fy.write(reinterpret_cast<const char*>(&y[i][j]), sizeof(double));
            fb.write(reinterpret_cast<const char*>(&isSolid[i][j]), sizeof(double));
        }
        fx.write(reinterpret_cast<const char*>(&x[Nx][j]), sizeof(double));
        fy.write(reinterpret_cast<const char*>(&y[Nx][j]), sizeof(double));
        fb.write(reinterpret_cast<const char*>(&isSolid[Nx][j]), sizeof(double));
    }
    fx.close();
    fy.close();
    fb.close();
}

void NavierStokes2D::yTDMA(const double dt)
{
    // Flip ustar and vstar to right hand side with the right indexing
    vector<double> Rx(Nx*Ny), Ry(Nx*Ny);
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = yVectorIndex(i, j);
            Rx[k] = ustar[i][j];
            Ry[k] = vstar[i][j];
        }
    
    // Initialize TDMA vectors
    vector<double> a(Nx*Ny, 0.0);
    vector<double> b(Nx*Ny, 1.0);
    vector<double> c(Nx*Ny, 0.0);
    
    // Find ustar
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = yVectorIndex(i, j);
            Ly(i, j, a[k], b[k], c[k], Rx[k], dt, u[i][j+1], u[i][j-1]);
        }
    // for solid cells, set R, a, c = 0 and b = 1;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = yVectorIndex(i, j);
            a[k] *= (1.0 - isSolid[i][j]);
            c[k] *= (1.0 - isSolid[i][j]);
            Rx[k] *= (1.0 - isSolid[i][j]);
            b[k] *= (1.0 - isSolid[i][j]);
            b[k] += isSolid[i][j];
        }
    TridiagonalMatrix Tx = TridiagonalMatrix(a, b, c);
    Tx.solve(Rx);
    
    // Find vstar
    fill(a.begin(), a.end(), 0.0);
    fill(b.begin(), b.end(), 1.0);
    fill(c.begin(), c.end(), 0.0);
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = yVectorIndex(i, j);
            Ly(i, j, a[k], b[k], c[k], Ry[k], dt, v[i][j+1], v[i][j-1]);
        }
    // for solid cells, set R, a, c = 0 and b = 1;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = yVectorIndex(i, j);
            a[k] *= (1.0 - isSolid[i][j]);
            c[k] *= (1.0 - isSolid[i][j]);
            Ry[k] *= (1.0 - isSolid[i][j]);
            b[k] *= (1.0 - isSolid[i][j]);
            b[k] += isSolid[i][j];
        }
    TridiagonalMatrix Ty = TridiagonalMatrix(a, b, c);
    Ty.solve(Ry);

    // Store these in ustar and vstar
    for (size_t k = 0; k < Rx.size(); k++)
    {
        size_t i = yVectorIndex_i(k);
        size_t j = yVectorIndex_j(k);
        ustar[i][j] = Rx[k];
        vstar[yVectorIndex_i(k)][yVectorIndex_j(k)] = Ry[k];
    }
    return;
}

void NavierStokes2D::xTDMA(const double dt)
{
    // Flip ustar and vstar to right hand side with the right indexing
    vector<double> Rx(Nx*Ny), Ry(Nx*Ny);
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = xVectorIndex(i, j);
            Rx[k] = ustar[i][j];
            Ry[k] = vstar[i][j];
        }
    
    // Initialize TDMA vectors
    vector<double> a(Nx*Ny, 0.0);
    vector<double> b(Nx*Ny, 1.0);
    vector<double> c(Nx*Ny, 0.0);
    
    // Find ustar
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = xVectorIndex(i, j);
            Lx(i, j, a[k], b[k], c[k], Rx[k], dt, u[i+1][j], u[i-1][j]);
        }
    // for solid cells, set R, a, c = 0 and b = 1;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = xVectorIndex(i, j);
            a[k] *= (1.0 - isSolid[i][j]);
            c[k] *= (1.0 - isSolid[i][j]);
            Rx[k] *= (1.0 - isSolid[i][j]);
            b[k] *= (1.0 - isSolid[i][j]);
            b[k] += isSolid[i][j];
        }
    TridiagonalMatrix Tx = TridiagonalMatrix(a, b, c);
    Tx.solve(Rx);
    
    // Reset TDMA vectors
    fill(a.begin(), a.end(), 0.0);
    fill(b.begin(), b.end(), 1.0);
    fill(c.begin(), c.end(), 0.0);

    // Find vstar
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = xVectorIndex(i, j);
            Lx(i, j, a[k], b[k], c[k], Ry[k], dt, v[i+1][j], v[i-1][j]);
        }
    // for solid cells, set R, a, c = 0 and b = 1;
    for (size_t i = 1; i <= Nx; i++)
        for (size_t j = 1; j <= Ny; j++)
        {
            size_t k = xVectorIndex(i, j);
            a[k] *= (1.0 - isSolid[i][j]);
            c[k] *= (1.0 - isSolid[i][j]);
            Rx[k] *= (1.0 - isSolid[i][j]);
            b[k] *= (1.0 - isSolid[i][j]);
            b[k] += isSolid[i][j];
        }
    TridiagonalMatrix Ty = TridiagonalMatrix(a, b, c);
    Ty.solve(Ry);
    
    // Store these in ustar and vstar
    for (size_t k = 0; k < Rx.size(); k++)
    {
        size_t i = xVectorIndex_i(k);
        size_t j = xVectorIndex_j(k);
        ustar[i][j] = Rx[k];
        vstar[xVectorIndex_i(k)][xVectorIndex_j(k)] = Ry[k];
    }
    return;
}

inline double NavierStokes2D::uu_e(size_t i, size_t j)
{
    return (0.5*(u[i+1][j]*u[i+1][j] + u[i][j]*u[i][j])*(1.0 - atRightBoundary[i][j])
    + u[i+1][j]*u[i+1][j]*rightIsDirichlet*atRightBoundary[i][j]
    + (0.5*dx*u[i+1][j] + u[i][j])*(0.5*dx*u[i+1][j] + u[i][j])*(1.0-rightIsDirichlet)*atRightBoundary[i][j])
    * (1.0 - isSolid[i+1][j]);
}

inline double NavierStokes2D::uu_w(size_t i, size_t j)
{
    return (0.5*(u[i-1][j]*u[i-1][j] + u[i][j]*u[i][j])*(1.0 - atLeftBoundary[i][j])
    + u[i-1][j]*u[i-1][j]*atLeftBoundary[i][j]*leftIsDirichlet
    + (-0.5*dx*u[i-1][j] + u[i][j])*(-0.5*dx*u[i-1][j] + u[i][j])*(1.0-leftIsDirichlet)*atLeftBoundary[i][j])
    * (1.0 - isSolid[i-1][j]);
}

inline double NavierStokes2D::uv_n(size_t i, size_t j)
{
    return (0.5*(u[i][j+1]*v[i][j+1] + u[i][j]*v[i][j])*(1.0 - atTopBoundary[i][j])
    + u[i][j+1]*v[i][j+1]*atTopBoundary[i][j]*topIsDirichlet
    + (0.5*dy*u[i][j+1] + u[i][j])*(0.5*dy*v[i][j+1] + v[i][j])*(1.0-topIsDirichlet)*atTopBoundary[i][j])
    * (1.0 - isSolid[i][j+1]);
}

inline double NavierStokes2D::uv_s(size_t i, size_t j)
{
    return (0.5*(u[i][j-1]*v[i][j-1] + u[i][j]*v[i][j])*(1.0 - atBottomBoundary[i][j])
    + u[i][j-1]*v[i][j-1]*atBottomBoundary[i][j]*bottomIsDirichlet
    + (-0.5*dy*u[i][j-1] + u[i][j])*(-0.5*dy*v[i][j-1] + v[i][j])*(1.0-bottomIsDirichlet)*atBottomBoundary[i][j])
    * (1.0 - isSolid[i][j-1]);
}

inline double NavierStokes2D::uv_e(size_t i, size_t j)
{
    return (0.5*(u[i+1][j]*v[i+1][j] + u[i][j]*v[i][j])*(1.0 - atRightBoundary[i][j])
    + u[i+1][j]*v[i+1][j]*rightIsDirichlet*atRightBoundary[i][j]
    + (0.5*dx*u[i+1][j] + u[i][j])*(0.5*dx*v[i+1][j] + v[i][j])*(1.0-rightIsDirichlet)*atRightBoundary[i][j])
    * (1.0 - isSolid[i+1][j]);
}

inline double NavierStokes2D::uv_w(size_t i, size_t j)
{
    return (0.5*(u[i-1][j]*v[i-1][j] + u[i][j]*v[i][j])*(1.0 - atLeftBoundary[i][j])
    + u[i-1][j]*v[i-1][j]*atLeftBoundary[i][j]*leftIsDirichlet
    + (-0.5*dx*u[i-1][j] + u[i][j])*(-0.5*dx*v[i-1][j] + v[i][j])*(1.0-leftIsDirichlet)*atLeftBoundary[i][j])
    * (1.0 - isSolid[i-1][j]);
}

inline double NavierStokes2D::vv_n(size_t i, size_t j)
{
    return (0.5*(v[i][j+1]*v[i][j+1] + v[i][j]*v[i][j])*(1.0 - atTopBoundary[i][j])
    + v[i][j+1]*v[i][j+1]*atTopBoundary[i][j]*topIsDirichlet
    + (0.5*dy*v[i][j+1] + v[i][j])*(0.5*dy* v[i][j+1] + v[i][j])*(1.0-topIsDirichlet)*atTopBoundary[i][j])
    * (1.0 - isSolid[i][j+1]);
}

inline double NavierStokes2D::vv_s(size_t i, size_t j)
{
    return (0.5*(v[i][j-1]*v[i][j-1] +  v[i][j]*v[i][j])*(1.0 - atBottomBoundary[i][j])
    + v[i][j-1]*v[i][j-1]*atBottomBoundary[i][j]*bottomIsDirichlet
    + (-0.5*dy* v[i][j-1] + v[i][j])*(-0.5*dy* v[i][j-1] + v[i][j])*(1.0-bottomIsDirichlet)*atBottomBoundary[i][j])
    * (1.0 - isSolid[i][j-1]);
}

inline size_t NavierStokes2D::yVectorIndex(size_t i, size_t j)
{ return (i-1)*Ny + j-1; }

inline size_t NavierStokes2D::yVectorIndex_i(size_t k)
{ return k / Ny + 1; }

inline size_t NavierStokes2D::yVectorIndex_j(size_t k)
{ return k % Ny + 1; }

inline size_t NavierStokes2D::xVectorIndex(size_t i, size_t j)
{ return i-1 + (j-1)*Nx; }

inline size_t NavierStokes2D::xVectorIndex_i(size_t k)
{ return k % Nx + 1; }

inline size_t NavierStokes2D::xVectorIndex_j(size_t k)
{ return k / Nx + 1; }

void NavierStokes2D::Ly(size_t i, size_t j, double& as, double& ap, double& an, double &R, const double dt, const double& un,const double& us)
{
    double pre = -dt/Re;
    double dy_inv = 1/dy;
    double dy_inv2 = dy_inv*dy_inv;
    
    /******** coefficients  for (dudy)_n ********/
    // If the point is not on the boudary, there is a contribution to an and ap
    an += pre*dy_inv2*(1.0 - atTopBoundary[i][j])*(1.0 - isSolid[i][j+1]);
    ap -= pre*dy_inv2*(1.0 - atTopBoundary[i][j])*(1.0 - isSolid[i][j+1]);

    // If it is on a Neumann boundary, move the boundary condition to the right hand side
    R -= pre*dy_inv*atTopBoundary[i][j]*(1.0 - topIsDirichlet)*un;

    // If it is on a dirichlet boundary, there is a contribution to the right hand side and ap
    ap -= 2.0*pre*dy_inv2*atTopBoundary[i][j]*topIsDirichlet;
    R -= 2.0*pre*dy_inv2*atTopBoundary[i][j]*topIsDirichlet*un;
    
    // If it is a solid boundary, there is a contribution to ap
    ap -= 2.0*pre*dy_inv2*isSolid[i][j+1];
    
    /******** coefficients  for (dudy)_s ********/
    // If the point is not on the boudary, there is a contribution to an and ap
    as += pre*dy_inv2*(1.0 - atBottomBoundary[i][j])*(1.0 - isSolid[i][j-1]);
    ap -= pre*dy_inv2*(1.0 - atBottomBoundary[i][j])*(1.0 - isSolid[i][j-1]);

    // If it is on a Neumann boundary, move the boundary condition to the right hand side
    R += pre*dy_inv*atBottomBoundary[i][j]*(1.0 - bottomIsDirichlet)*us;
    
    // If it is on a dirichlet boundary, there is a contribution to the right hand side and ap
    ap -= 2.0*pre*dy_inv2*atBottomBoundary[i][j]*bottomIsDirichlet;
    R -= 2.0*pre*dy_inv2*atBottomBoundary[i][j]*bottomIsDirichlet*us;
    
    // If it is a solid boundary, there is a contribution to ap
    ap -= 2.0*pre*dy_inv2*isSolid[i][j-1];
}

void NavierStokes2D::Lx(size_t i, size_t j, double& aw, double& ap, double& ae, double &R, const double dt, const double& ue, const double& uw)
{
    double pre = -dt/Re;
    double dx_inv = 1.0/dx;
    double dx_inv2 = dx_inv*dx_inv;
    
    /******** coefficients  for (dudx)_e ********/
    // If the point is not on the boudary, there is a contribution to ae and ap
    ae += pre*dx_inv2*(1.0 - atRightBoundary[i][j])*(1.0 - isSolid[i+1][j]);
    ap -= pre*dx_inv2*(1.0 - atRightBoundary[i][j])*(1.0 - isSolid[i+1][j]);
    
    // If it is on a Neumann boundary, move the boundary condition to the right hand side
    R -= pre*dx_inv*atRightBoundary[i][j]*(1.0 - rightIsDirichlet)*ue;
    
    // If it is on a dirichlet boundary, there is a contribution to the right hand side and ap
    ap -= 2.0*pre*dx_inv2*atRightBoundary[i][j]*rightIsDirichlet;
    R -= 2.0*pre*dx_inv2*atRightBoundary[i][j]*rightIsDirichlet*ue;
    
    // If it is a solid boundary, there is a contribution to ap
    ap -= 2.0*pre*dx_inv2*isSolid[i+1][j];
    
    /******** coefficients  for (dudx)_w ********/
    // If the point is not on the boudary, there is a contribution to ae and ap
    aw += pre*dx_inv2*(1.0 - atLeftBoundary[i][j])*(1.0 - isSolid[i-1][j]);
    ap -= pre*dx_inv2*(1.0 - atLeftBoundary[i][j])*(1.0 - isSolid[i-1][j]);
    
    // If it is on a Neumann boundary, move the boundary condition to the right hand side
    R += pre*dx_inv*atLeftBoundary[i][j]*(1.0 - leftIsDirichlet)*uw;
    
    // If it is on a dirichlet boundary, there is a contribution to the right hand side and ap
    ap -= 2.0*pre*dx_inv2*atLeftBoundary[i][j]*leftIsDirichlet;
    R -= 2.0*pre*dx_inv2*atLeftBoundary[i][j]*leftIsDirichlet*uw;
    
    // If it is a solid boundary, there is a contribution to ap
    ap -= 2.0*pre*dx_inv2*isSolid[i-1][j];
}



