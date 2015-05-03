#include "TridiagonalMatrix.h"

TridiagonalMatrix::TridiagonalMatrix(){ N = 0; }

TridiagonalMatrix::TridiagonalMatrix(const vector<double>& a, const vector<double>& b, const vector<double>& c)
{
    N = (unsigned int)a.size();
    if (b.size() == N && c.size() == N)
    {
        this->a = a;
        this->b = b;
        this->c = c;
    }
    else
        throw invalid_argument("a, b, and c must be the same size");
    
    // Resize vectors
    bf.resize(N);
    m.resize(N);
    
    // Perform the forward matrix evaluations on construction
    bf[0] = b[0];
    m[0] = 0;
    for (int i = 1; i < N; i++)
    {
        m[i] = a[i]/bf[i-1];
        bf[i] = b[i] - m[i]*c[i-1];
    }
}

TridiagonalMatrix::~TridiagonalMatrix() { return; }

void TridiagonalMatrix::solve(vector<double> &d)
{
    if (d.size() != N)
        throw invalid_argument("u must be the same size as the matrix");
    
    // Forward sweep
    for (unsigned int i = 1; i < N; i++)
        d[i] -= m[i]*d[i-1];
    
    // backward sweep
    d[N-1] = d[N-1]/bf[N-1];
    for (int i = N-2; i >= 0; i--)
        d[i] = (d[i] - c[i]*d[i+1])/bf[i];
}

vector<double> TridiagonalMatrix::superdiagonal() { return c; };

vector<double> TridiagonalMatrix::diagonal() { return b; };

vector<double> TridiagonalMatrix::subdiagonal() { return a; };