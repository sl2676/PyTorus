#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

typedef double** PPdbl;
typedef double* Pdbl;

void tred2(PPdbl a, const int n, Pdbl d, Pdbl e, const char EV) {
    int l, k, j, i;
    double scale, hh, h, g, f;
    for (i = n - 1; i > 0; i--) {
        l = i - 1;
        h = scale = 0.;
        if (l > 0) {
            for (k = 0; k <= l; k++)
                scale += fabs(a[i][k]);
            if (scale == 0.)
                e[i] = a[i][l];
            else {
                for (k = 0; k <= l; k++) {
                    a[i][k] /= scale;
                    h += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = (f >= 0.) ? -sqrt(h) : sqrt(h);
                e[i] = scale * g;
                h -= f * g;
                a[i][l] = f - g;
                f = 0.;
                for (j = 0; j <= l; j++) {
                    if (EV) a[j][i] = a[i][j] / h;
                    for (k = 0, g = 0.; k <= j; k++)
                        g += a[j][k] * a[i][k];
                    for (k = j + 1; k <= l; k++)
                        g += a[k][j] * a[i][k];
                    e[j] = g / h;
                    f += e[j] * a[i][j];
                }
                hh = f / (h + h);
                for (j = 0; j <= l; j++) {
                    f = a[i][j];
                    e[j] = g = e[j] - hh * f;
                    for (k = 0; k <= j; k++)
                        a[j][k] -= f * e[k] + g * a[i][k];
                }
            }
        } else
            e[i] = a[i][l];
        d[i] = h;
    }
    d[0] = e[0] = 0.;
    if (EV) {
        for (i = 0; i < n; i++) {
            l = i - 1;
            if (d[i]) {
                for (j = 0; j <= l; j++) {
                    for (k = 0, g = 0.; k <= l; k++)
                        g += a[i][k] * a[k][j];
                    for (k = 0; k <= l; k++)
                        a[k][j] -= g * a[k][i];
                }
            }
            d[i] = a[i][i];
            a[i][i] = 1.0;
            for (j = 0; j <= l; j++)
                a[j][i] = a[i][j] = 0.;
        }
    } else
        for (i = 0; i < n; i++)
            d[i] = a[i][i];
}

void printMatrix(PPdbl a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << a[i][j] << " ";
        }
        cout << endl;
    }
}

void printArray(Pdbl arr, int n) {
    for (int i = 0; i < n; i++) {
        cout << setw(10) << arr[i] << " ";
    }
    cout << endl;
}

int main() {
    const int n = 4;
    double a_data[n][n] = {
        {4, 1, -2, 2},
        {1, 2, 0, 1},
        {-2, 0, 3, -2},
        {2, 1, -2, -1}
    };

    PPdbl a = new Pdbl[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) {
            a[i][j] = a_data[i][j];
        }
    }

    double d[n], e[n];
    char EV = 1; // Compute eigenvectors

    cout << "Original matrix:" << endl;
    printMatrix(a, n);

    tred2(a, n, d, e, EV);

    cout << "\nTridiagonal matrix (stored in 'a'):" << endl;
    printMatrix(a, n);

    cout << "\nDiagonal elements (d):" << endl;
    printArray(d, n);

    cout << "\nOff-diagonal elements (e):" << endl;
    printArray(e, n);

    for (int i = 0; i < n; i++) {
        delete[] a[i];
    }
    delete[] a;
	cout << "EV " << EV << " Value" << endl;
    return 0;
}

