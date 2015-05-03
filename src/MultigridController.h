#ifndef __NavierStokes2D__MultigridController__
#define __NavierStokes2D__MultigridController__

#include <stdio.h>
#include <vector>
#include "math.h"
using std::vector;

class MultigridController
{
public:
    MultigridController(size_t Nx, size_t Ny);
    virtual ~MultigridController() = 0;
    short int getCommand();
    void reset();
    size_t numberOfLevels();
    
protected:
    vector<short int> control_order;
    vector<short int>::iterator it;
    size_t nlev;
    virtual void initialize() = 0;
};

class Vcycle : public MultigridController
{
public:
    Vcycle(size_t Nx, size_t Ny);
    ~Vcycle();
    
protected:
    void initialize();
};

class Wcycle : public MultigridController
{
public:
    Wcycle(size_t Nx, size_t Ny);
    ~Wcycle();
    
protected:
    void initialize();
};

class Fcycle : public MultigridController
{
public:
    Fcycle(size_t Nx, size_t Ny);
    ~Fcycle();
    
protected:
    void initialize();
};

#endif /* defined(__NavierStokes2D__MultigridController__) */
