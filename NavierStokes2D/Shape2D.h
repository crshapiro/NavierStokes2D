#ifndef __NavierStokes2D__Shape2D__
#define __NavierStokes2D__Shape2D__

#include <stdio.h>
#include <vector>
#include <math.h>
#include <stdexcept>
using std::invalid_argument;

class Shape2D
{
public:
    Shape2D(double xc, double yc);
    ~Shape2D();
    virtual bool inVolume(double x, double y) const = 0;

protected:
    double xc, yc;
};

class Square : public Shape2D
{
public:
    Square(double xc, double yc, double s);
    ~Square();
    bool inVolume(double x, double y) const;
    
protected:
    double s;
};

class Circle : public Shape2D
{
public:
    Circle(double xc, double yc, double R);
    ~Circle();
    bool inVolume(double x, double y) const;
    
protected:
    double R;
};

#endif /* defined(__NavierStokes2D__Shape2D__) */
