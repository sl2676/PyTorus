#include "../../Torus/src/utils/Compress.h"

std::string put_double_wrapper(const std::vector<float>& x) {
    std::ostringstream to;
    int C = 0;
    for (size_t k = 0; k < x.size(); ++k, ++C) {
        static_cast<void (*)(const float, std::ostream&)>(put)(x[k], to);
        if (C >= 15) {
            C = -1;
            to << '\n';
        }
    }
    return to.str();
}
std::string put_three_wrapper(const std::vector<std::vector<std::vector<float>>>& x) {
    std::ostringstream to;
    int C = 0;
    for (const auto& layer : x) {
        for (const auto& row : layer) {
            for (const auto& val : row) {
                put(val, to);
                if (C >= 15) {
                    C = -1;
                    to << '\n';
                }
                ++C;
            }
        }
    }
    return to.str();
}
std::string put_single_wrapper(const float x) {
    std::ostringstream to;
    put(x, to); 
    return to.str();
}

std::string getFORTRAN_two_Wrapper(const std::vector<std::vector<float>>& x) {
    std::ostringstream to;
    int C = 0;
    int n[2] = {static_cast<int>(x.size()), x.empty() ? 0 : static_cast<int>(x[0].size())};

    for (int j = 0; j < n[1]; ++j) {
        for (int i = 0; i < n[0]; ++i, ++C) {
            if (x[i].size() != n[1]) {
                throw std::invalid_argument("Inconsistent number of columns in 2D array input.");
            }
            put(x[i][j], to);
            if (C >= 14 && (i < n[0] - 1 || j < n[1] - 1)) {
                C = -1;
                to << '\n';
            }
        }
    }
    return to.str();
}
std::string getFORTRAN_three_Wrapper(const std::vector<std::vector<std::vector<float>>>& x) {
    std::ostringstream to;
    int n[3] = {static_cast<int>(x.size()), 
                x.empty() ? 0 : static_cast<int>(x[0].size()), 
                x.empty() || x[0].empty() ? 0 : static_cast<int>(x[0][0].size())};
    int C = 0;
    for (int k = 0; k < n[2]; ++k) {
        for (int j = 0; j < n[1]; ++j) {
            for (int i = 0; i < n[0]; ++i, ++C) {
                put(x[i][j][k], to); 
                if (C >= 14 && (i < n[0] - 1 || j < n[1] - 1 || k < n[2] - 1)) {
                    C = -1;
                    to << '\n';
                }
            }
        }
    }
    return to.str();
}
std::vector<std::vector<float>> getFORTRAN_two_routine_Wrapper(const std::string& inputData, const std::vector<int>& n) {
    if (n.size() != 2) {
        throw std::invalid_argument("Dimension array must have exactly two elements.");
    }

    std::istringstream in(inputData);
    int nArray[2] = {n[0], n[1]};
    float** x;
    
    Alloc2D(x, nArray);

    get(x, nArray, in);

    std::vector<std::vector<float>> result(n[0], std::vector<float>(n[1]));
    for (int i = 0; i < n[0]; ++i) {
        for (int j = 0; j < n[1]; ++j) {
            result[i][j] = x[i][j];
        }
    }

    Free2D(x);
    return result;
}
std::vector<std::vector<std::vector<float>>> getFORTRAN_three_routine_Wrapper(const std::string& inputData, const std::vector<int>& dimensions) {
    if (dimensions.size() != 3) {
        throw std::invalid_argument("Dimensions array must have exactly three elements.");
    }
    
    int n[3] = {dimensions[0], dimensions[1], dimensions[2]};
    float*** x;
    Alloc3D(x, n);

    std::istringstream inStream(inputData);
    get(x, n, inStream);

    std::vector<std::vector<std::vector<float>>> result(n[0], std::vector<std::vector<float>>(n[1], std::vector<float>(n[2])));
    for (int i = 0; i < n[0]; ++i) {
        for (int j = 0; j < n[1]; ++j) {
            for (int k = 0; k < n[2]; ++k) {
                result[i][j][k] = x[i][j][k];
            }
        }
    }

    Free3D(x);
    return result;
}

std::string compress_wrapper(const float x) {
    char s[5];
    compress(x, s);
    return std::string(s, 5);
}

float uncompress_wrapper(const std::string& input) {
    if (input.length() != 5) {
        throw std::invalid_argument("Input string must be exactly 5 characters long.");
    }
    const char s[5] = {input[0], input[1], input[2], input[3], input[4]};
    return uncompress(s);
}

std::string Put_one_Wrapper(const std::vector<float>& x) {
    std::ostringstream to;
    int C = 0; 
    for (unsigned long k = 0; k < x.size(); ++k) {
        float xk = x[k];
        if (xk == 0.f) {
            if (C > 79) {
                to << '\n';
                C = 0;
            }
            to << char(95); 
            C++;
        } else {
            if (C > 75) {
                to << '\n';
                C = 0;
            }
            put(xk, to);
            C += 5; 
        }
    }
    return to.str();
}

std::string Put_two_Wrapper(const std::vector<std::vector<float>>& x) {
    std::ostringstream to;
    int C = 0; 
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < x[i].size(); ++j) {
            float xk = x[i][j];
            if (xk == 0.f) {
                if (C > 79) {
                    to << '\n';
                    C = 0;
                }
                to << char(95); 
                C++;
            } else {
                if (C > 75) {
                    to << '\n';
                    C = 0;
                }
                put(xk, to); 
                C += 5; 
            }
        }
    }
    return to.str();
}

std::string Put_three_Wrapper(const std::vector<std::vector<std::vector<float>>>& x) {
    std::ostringstream to;
    int C = 0; 
    for (const auto& layer : x) {
        for (const auto& row : layer) {
            for (const auto& val : row) {
                if (val == 0.f) {
                    if (C > 79) {
                        to << '\n';
                        C = 0;
                    }
                    to << char(95); 
                    C++;
                } else {
                    if (C > 75) {
                        to << '\n';
                        C = 0;
                    }
                    put(val, to); 
                    C += 5; 
                }
            }
        }
    }
    return to.str();
}


namespace py = pybind11;
void init_compress(py::module_ &torus) {
	torus.def("compress", &compress_wrapper);
    torus.def("compress", py::overload_cast<const float*, char*, const int>(&compress));
	torus.def("compress", py::overload_cast<const double*, char*, const int>(&compress));
	torus.def("uncompress", &uncompress_wrapper);
	torus.def("uncompress", py::overload_cast<char*, float*, const int>(&uncompress));
	torus.def("uncompress", py::overload_cast<char*, double*, const int>(&uncompress));
	torus.def("put", &put_single_wrapper);
	torus.def("put", &put_double_wrapper);	
	torus.def("put", &put_three_wrapper);
	torus.def("Put", &Put_one_Wrapper);
	torus.def("Put", &Put_two_Wrapper);
	torus.def("Put", &Put_three_Wrapper);

	torus.def("getFORTRAN", &getFORTRAN_two_Wrapper);
	torus.def("getFORTRAN", &getFORTRAN_three_Wrapper);
	torus.def("getFORTRAN", &getFORTRAN_two_routine_Wrapper);
	torus.def("getFORTRAN", &getFORTRAN_three_routine_Wrapper);

}


