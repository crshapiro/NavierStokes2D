#include "MultigridController.h"

MultigridController::MultigridController(size_t Nx, size_t Ny)
{
    nlev = floor(log2(std::min(Nx, Ny)));
}

MultigridController::~MultigridController()
{ return; }

void MultigridController::reset()
{ it = control_order.begin(); }

short int MultigridController::getCommand()
{ return *it++; }

size_t MultigridController::numberOfLevels()
{ return nlev; }

Vcycle::Vcycle(size_t Nx, size_t Ny) : MultigridController::MultigridController(Nx, Ny)
{ initialize(); }

Vcycle::~Vcycle()
{ return; }

void Vcycle::initialize()
{
    // Create call order of levels
    for (int i = 0; i < nlev - 1; i++)
        control_order.push_back(-1);
    for (int i = 0; i < nlev; i++)
        control_order.push_back(1);
    
    it = control_order.begin();
}

Fcycle::Fcycle(size_t Nx, size_t Ny) : MultigridController::MultigridController(Nx, Ny)
{ initialize(); }

Fcycle::~Fcycle()
{ return; }

void Fcycle::initialize()
{
    // Create call order of levels
    for (int i = 0; i < nlev - 1; i++)
        control_order.push_back(-1);
    for (int i = 0; i < nlev - 2; i++)
    {
        for (int j = 0; j < i + 1; j++)
            control_order.push_back(1);
        for (int j = 0; j < i + 1; j++)
            control_order.push_back(-1);
    }
    for (int i = 0; i < nlev; i++)
        control_order.push_back(1);
    
    it = control_order.begin();
}

Wcycle::Wcycle(size_t Nx, size_t Ny) : MultigridController::MultigridController(Nx, Ny)
{ initialize(); }

Wcycle::~Wcycle()
{ return; }

void Wcycle::initialize()
{
    // Create call order of levels
    for (int i = 0; i < nlev - 1; i++)
        control_order.push_back(-1);
    for (int i = 0; i < nlev - 2; i++)
    {
        for (int j = 0; j < i + 1; j++)
            control_order.push_back(1);
        for (int j = 0; j < i + 1; j++)
            control_order.push_back(-1);
    }
    for (int i = (int)nlev - 4; i >= 0; i--)
    {
        for (int j = 0; j < i + 1; j++)
            control_order.push_back(1);
        for (int j = 0; j < i + 1; j++)
            control_order.push_back(-1);
    }
    for (int i = 0; i < nlev; i++)
        control_order.push_back(1);
    
    it = control_order.begin();
}