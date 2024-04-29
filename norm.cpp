#include <iostream>
#include <cmath>

template<class T, int N>
class Vector {
private:
    T elements[N];

public:
    Vector() {
        for (int i = 0; i < N; ++i)
            elements[i] = 0;
    }

    Vector(const T (&arr)[N]) {
        for (int i = 0; i < N; ++i)
            elements[i] = arr[i];
    }

    T norm() const {
        T x = elements[0] * elements[0];
        for (int i = 1; i < N; i++)
            x += elements[i] * elements[i];
        return x;
    }
};

template<class T, int N>
inline T norm(const Vector<T,N>& V) {
    return V.norm();
}

int main() {
    // Example usage
    double arr[] = {1.0, 1.0, 1.0};
    Vector<double, 3> vec(arr);
    std::cout << "Norm of vector: " << norm(vec) << std::endl;

    return 0;
}

