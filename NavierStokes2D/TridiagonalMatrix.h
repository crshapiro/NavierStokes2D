#ifndef __Poisson2D__TridiagonalMatrix__
#define __Poisson2D__TridiagonalMatrix__

#include <stdio.h>
#include <vector>
using namespace std;

class TridiagonalMatrix
{
public:
    TridiagonalMatrix();
    TridiagonalMatrix(const vector<double>& a, const vector<double>& b, const vector<double>& c);
    ~TridiagonalMatrix();
    void solve(vector<double> &d);
    vector<double> superdiagonal();
    vector<double> diagonal();
    vector<double> subdiagonal();
    
protected:
    unsigned int N;
    vector<double> a, b, c;
    vector<double> bf, m;
};

#endif /* defined(__Poisson2D__TridiagonalMatrix__) */
