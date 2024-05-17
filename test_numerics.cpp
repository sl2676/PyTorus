#include <iostream>
#include <cmath>
#include <algorithm>

typedef double* Pdbl;
typedef double** PPdbl;
typedef int* Pint;

// Function prototypes
int LUDecomposition(PPdbl a, const int n, Pint indx, int& d);
void LUSolution(PPdbl a, const int n, const Pint indx, Pdbl b);

int main() {
    const int n = 3;
    double matrix[n][n] = {
        {4, 12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    };
    double* a[n];
    for(int i = 0; i < n; ++i)
        a[i] = matrix[i];

    double b[n] = {1, 2, 3};

    int indx[n];
    int d;

    int res = LUDecomposition(a, n, indx, d);

    if(res == 0) {
        LUSolution(a, n, indx, b);

        std::cout << "Solution vector x:\n";
        for(int i = 0; i < n; ++i)
            std::cout << b[i] << " ";
        std::cout << "\n";
    } else {
        std::cout << "LU decomposition failed.\n";
    }

	if(res == 0) {
        std::cout << "LU decomposition succeeded.\n";
        std::cout << "LU matrix:\n";
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j)
                std::cout << a[i][j] << " ";
            std::cout << "\n";
        }
        std::cout << "Permutation vector:\n";
        for(int i = 0; i < n; ++i)
            std::cout << indx[i] << " ";
        std::cout << "\n";
    } else {
        std::cout << "LU decomposition failed.\n";
    }

    return 0;
}

int LUDecomposition(PPdbl a, const int n, Pint indx, int& d) {
    const double tiny = 1.e-20;
    int i, imax, j, k;
    double big, dum, sum;
    Pdbl vv = new double[n];

    d = 1;
    for(i = 0; i < n; i++) {
        big = 0.;
        for(j = 0; j < n; j++) big = std::max(big, std::abs(a[i][j]));
        if(big == 0.) return -1; // Singular matrix
        vv[i] = 1. / big;
    }
    for(j = 0; j < n; j++) {
        imax = j;
        for(i = 0; i < j; i++) {
            sum = a[i][j];
            for(k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.;
        for(i = j; i < n; i++) {
            sum = a[i][j];
            for(k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if(j != imax) {
            for(k = 0; k < n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if(a[j][j] == 0.) a[j][j] = tiny;
        if(j < n - 1) {
            dum = 1. / a[j][j];
            for(i = j + 1; i < n; i++) a[i][j] *= dum;
        }
    }
    delete[] vv;
    return 0;
}

void LUSolution(PPdbl a, const int n, const Pint indx, Pdbl b) {
    int i, ii = -1, ip, j;
    double sum;
    for(i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if(ii >= 0) {
            for(j = ii; j < i; j++) sum -= a[i][j] * b[j];
        } else if(sum) {
            ii = i;
        }
        b[i] = sum;
    }
    for(i = n - 1; i >= 0; i--) {
        sum = b[i];
        for(j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

