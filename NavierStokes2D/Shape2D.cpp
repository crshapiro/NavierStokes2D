#include "Shape2D.h"

Shape2D::Shape2D(double xc, double yc)
{
    this->xc = xc;
    this->yc = yc;
}

Shape2D::~Shape2D()
{ return; }

Square::Square(double xc, double yc, double s) : Shape2D(xc, yc)
{
    if (s < 0)
        throw invalid_argument("Length of side must be positive");
    this->s = s;
}

Square::~Square()
{ return; }

bool Square::inVolume(double x, double y) const
{ return (x <= xc + s/2 && x >= xc - s/2 && y <= yc + s/2 && y >= yc - s/2); }

Circle::Circle(double xc, double yc, double R) : Shape2D(xc, yc)
{
    if (R < 0)
        throw invalid_argument("Radius must be positive");
    this->R = R;
}

Circle::~Circle()
{ return; }

bool Circle::inVolume(double x, double y) const
{ return (sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) ) <= R); }