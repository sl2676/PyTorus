#pragma once
#include <iomanip> 
#include <iostream>
#include <sstream> 
#include <type_traits>
#include <variant>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <stdio.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include <numeric>
#include <cmath>
#include <complex>
#include <algorithm>
#include <sstream>
#include "../../Torus/src/utils/Vector.h"
#include "../../Torus/src/utils/Compress.h"
namespace py = pybind11;

class PyMatrix;

template<typename T> struct Promote { using Type = T; };
template<> struct Promote<int> { using Type = double; };
template<> struct Promote<double> { using Type = std::complex<double>; };
template<> struct Promote<std::pair<int, std::complex<double>>> { using Type = std::complex<double>; };
template<> struct Promote<std::pair<std::complex<double>, int>> { using Type = std::complex<double>; };

template<typename T, typename U>
struct PromotedType {
    using type = decltype(T() * U() - T() * U());
};

template<typename Scalar>
using Promote_t = typename Promote<Scalar>::Type;

template<typename T, typename U>
struct CommonTypeDeduction {
    using type = std::complex<double>; 
};

template<>
struct CommonTypeDeduction<int, int> {
	using type = double;
};

template<>
struct CommonTypeDeduction<double, double> {
    using type = double;
};
template<>
struct CommonTypeDeduction<double, int> {
	using type = double;
};

template<>
struct CommonTypeDeduction<int, double> {
	using type = double;
};

template<>
struct CommonTypeDeduction<std::complex<double>, int> {
	using type = std::complex<double>;
};

template<>
struct CommonTypeDeduction<std::complex<double>, double> {
	using type = std::complex<double>;
};

template<>
struct CommonTypeDeduction<std::complex<double>, std::complex<double>>  {
	using type = std::complex<double>;
};

template<>
struct CommonTypeDeduction<int, std::complex<double>> {
	using type = std::complex<double>;
};

template<>
struct CommonTypeDeduction<double, std::complex<double>> {
	using type = std::complex<double>;
};

template<typename T, typename U>
struct DotProductResultType {
    using type = typename CommonTypeDeduction<T, U>::type;
};

template<typename T, typename U>
using CommonType = typename CommonTypeDeduction<T, U>::type;


template<typename T, typename U>
using DotProductResultType_t = typename DotProductResultType<T, U>::type;


struct BaseVector;

template<typename T, typename PromotedType>
auto convertVector(const BaseVector* vector);


// Base Vector struct parent struct
struct BaseVector {
    virtual ~BaseVector() = default;
    
	virtual void setItem(size_t index, py::object value) = 0;
	virtual void serialize(std::ostream& os) const = 0;
	virtual void serialize2(std::ostream& os) const = 0;
	
	virtual std::unique_ptr<BaseVector> add(const BaseVector* other) const = 0;
	virtual std::unique_ptr<BaseVector> subtract(const BaseVector* other) const = 0;
	virtual std::unique_ptr<BaseVector> divide(const BaseVector* other) const = 0;
	virtual std::unique_ptr<BaseVector> cross(const BaseVector* other) const = 0;
	virtual std::unique_ptr<BaseVector> conj() const = 0;
	virtual std::unique_ptr<BaseVector> arg() const = 0;
	virtual std::unique_ptr<BaseVector> imag() const = 0;
	virtual std::unique_ptr<BaseVector> real() const = 0;
	virtual std::unique_ptr<BaseVector> multiply_elements(const BaseVector* other) const = 0;
	
	virtual py::object dot_product(const BaseVector* other) const = 0;
	virtual py::object dot(const BaseVector* other) const = 0;
	virtual py::object getItem(size_t index) const = 0;
	virtual py::object toPython() const = 0;

	virtual size_t size() const = 0;
};

template<typename T, typename U>
auto dotProduct(const std::vector<T>& vec1, const std::vector<U>& vec2) {
    using CommonResultType = typename CommonTypeDeduction<T, U>::type;

    if (vec1.size() != vec2.size()) {
        throw std::runtime_error("Vector sizes do not match for dot product");
    }

    CommonResultType result = CommonResultType{};

    for (size_t i = 0; i < vec1.size(); ++i) {
        result += CommonResultType(vec1[i]) * CommonResultType(vec2[i]);
    }

    return result;
}

// Typed Vector struct for operations and child of base vector
template<typename T>
struct TypedVector : BaseVector {
     std::vector<T> vec;

    TypedVector() = default;
    TypedVector(const std::initializer_list<T>& init) : vec(init) {}

    explicit TypedVector(std::vector<T> v) : vec(std::move(v)) {}


    TypedVector(const Vector<T, 3>& customVec) {
        for (int i = 0; i < 3; ++i) {
            vec.push_back(customVec[i]); 
        }
    }
	
	void serialize2(std::ostream& os) const override {
    	if constexpr (std::is_same_v<T, double>) {
        	double* A = new double[this->vec.size()];
        	std::copy(this->vec.begin(), this->vec.end(), A);
        	::put(A, this->vec.size(), os);
        	delete[] A;
    	} else {
        	serialize(os); 
    	}
	}
	
	void serialize(std::ostream& os) const override {
        for (const auto& value : vec) {
            if constexpr (std::is_same_v<T, std::complex<double>>) {
                double realPart = value.real();
                double imagPart = value.imag();
                put(realPart, os); // Serialize the real part
                put(imagPart, os); // Serialize the imaginary part
            } else {
                put(static_cast<double>(value), os);
            }
        }
    }

	size_t size() const override {
        return vec.size();
    }


	template<typename Func>
    auto apply(Func&& func) const -> std::unique_ptr<BaseVector> {
        using ResultTypeRaw = decltype(func(std::declval<T>()));
        using ResultType = typename std::decay<ResultTypeRaw>::type;
        using CommonType = typename CommonTypeDeduction<T, ResultType>::type;

        std::vector<CommonType> resultVec;
        resultVec.reserve(vec.size());

        for (const auto& element : vec) {
            resultVec.push_back(static_cast<CommonType>(func(element)));
        }

        return std::make_unique<TypedVector<CommonType>>(resultVec);
    }

	std::unique_ptr<BaseVector> multiply_elements(const BaseVector* other) const override {
    	if (auto o = dynamic_cast<const TypedVector<int>*>(other)) {
        	return multiply_elements_typed(*o);
    	} else if (auto o = dynamic_cast<const TypedVector<double>*>(other)) {
        	return multiply_elements_typed(*o);
    	} else if (auto o = dynamic_cast<const TypedVector<std::complex<double>>*>(other)) {
        	return multiply_elements_typed(*o);
    	} else {
        	throw std::runtime_error("Unsupported vector type for element-wise multiplication.");
    	}
	}

	template<typename U>
	std::unique_ptr<BaseVector> multiply_elements_typed(const TypedVector<U>& other) const {
    	using CommonType = typename CommonTypeDeduction<T, U>::type;

    	if (vec.size() != other.vec.size()) {
        	throw std::runtime_error("Vector sizes do not match for element-wise multiplication.");
    	}

    	std::vector<CommonType> result(vec.size());
    	for (size_t i = 0; i < vec.size(); ++i) {
        	result[i] = CommonType(vec[i]) * CommonType(other.vec[i]);
    	}

    	return std::make_unique<TypedVector<CommonType>>(std::move(result));
	}




	std::unique_ptr<BaseVector> conj() const override {
    	if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<double>>) {
        	std::vector<T> resultVec(vec.size());
        	std::transform(vec.begin(), vec.end(), resultVec.begin(), [](const T& value) -> T {
            	return std::conj(value);
        	});
        	return std::make_unique<TypedVector<T>>(std::move(resultVec));
    	} else {
        	std::vector<std::complex<double>> resultVec(vec.size());
        	std::transform(vec.begin(), vec.end(), resultVec.begin(), [](const T& value) -> std::complex<double> {
            	return {static_cast<double>(value), 0.0f};
        	});
        	return std::make_unique<TypedVector<std::complex<double>>>(std::move(resultVec));
    	}
	}

	std::unique_ptr<BaseVector> imag() const override {
    	if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<double>>) {
        	std::vector<double> resultVec; 
        	resultVec.reserve(vec.size());
        	std::transform(vec.begin(), vec.end(), std::back_inserter(resultVec), [](const T& value) -> double {
            	return std::imag(value); 
       		});
        	return std::make_unique<TypedVector<double>>(std::move(resultVec));
    	} else {
        	std::vector<double> resultVec(vec.size(), 0.0f); 
        	return std::make_unique<TypedVector<double>>(std::move(resultVec));
    	}
	}

	std::unique_ptr<BaseVector> real() const override {
		if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<double>>) {
			std::vector<double> resultVec;
			resultVec.reserve(vec.size());
			std::transform(vec.begin(), vec.end(), std::back_inserter(resultVec), [](const T& value) -> double {
				return std::real(value);
			});
			return std::make_unique<TypedVector<double>>(std::move(resultVec));
		} else {
			std::vector<double> resultVec(vec.size(), 0.0f);
			return std::make_unique<TypedVector<double>>(std::move(resultVec));
		}
	}

	std::unique_ptr<BaseVector> arg() const override {
		if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<double>>) {
			std::vector<double> resultVec;
			resultVec.reserve(vec.size());
			std::transform(vec.begin(), vec.end(), std::back_inserter(resultVec), [](const T& value) -> double {
				return std::arg(value);
			});
			return std::make_unique<TypedVector<double>>(std::move(resultVec));
		} else {
			std::vector<double> resultVec(vec.size(), 0.0f);
			return std::make_unique<TypedVector<double>>(std::move(resultVec));
		}
	}
	
	double norm() const {
        double sum = 0;
        for (const auto& element : vec) {
            sum += std::norm(element); 
        }
		return std::sqrt(sum);
	}	
	
	py::object dot(const BaseVector* other) const override {
        if (const auto* o = dynamic_cast<const TypedVector<T>*>(other)) {
            T result = std::inner_product(vec.begin(), vec.end(), o->vec.begin(), T(0));
            return py::cast(result);
        }
        
        if (const auto* oInt = dynamic_cast<const TypedVector<int>*>(other)) {
            if constexpr (std::is_same_v<T, double>) {
                double result = std::inner_product(vec.begin(), vec.end(), oInt->vec.begin(), 0.0f,
                    std::plus<>(), [](double a, int b) { return a * b; });
                return py::cast(result);
            }
        }
        if (const auto* oFloat = dynamic_cast<const TypedVector<double>*>(other)) {
            if constexpr (std::is_same_v<T, int>) {
                double result = std::inner_product(vec.begin(), vec.end(), oFloat->vec.begin(), 0.0f,
                	std::plus<>(), [](int a, double b) { return a * b; });
                return py::cast(result);
            }
        }
        throw std::runtime_error("Incompatible vector types for dot product.");
    }
	std::unique_ptr<BaseVector> add(const BaseVector* other) const override {
        const auto* o = dynamic_cast<const TypedVector<T>*>(other);
        if (!o || o->vec.size() != this->vec.size()) {
            throw std::runtime_error("Vector sizes or types do not match for addition.");
        }
        std::vector<T> result;
        result.reserve(this->vec.size());
        for (size_t i = 0; i < this->vec.size(); ++i) {
            result.push_back(this->vec[i] + o->vec[i]);
        }
        return std::make_unique<TypedVector<T>>(result);
    }
	// added stuff
	std::unique_ptr<BaseVector> subtract(const BaseVector* other) const override {
		const auto* o = dynamic_cast<const TypedVector<T>*>(other);
		if (!o || o->vec.size() != this->vec.size()) {
			throw std::runtime_error("Vector sizes or types do not match for subtraction.");
		}
		std::vector<T> result;
		result.reserve(this->vec.size());
		for (size_t i = 0; i < this->vec.size(); ++i) {
			result.push_back(this->vec[i] + o->vec[i]);
		}
		return std::make_unique<TypedVector<T>>(result);
	}

	std::unique_ptr<BaseVector> divide(const BaseVector* other) const override {
		const auto* o = dynamic_cast<const TypedVector<T>*>(other);
		if (!o || o->vec.size() != this->vec.size()) {
			throw std::runtime_error("Vector sizes or types do not match for division.");
		}
		std::vector<T> result;
		result.reserve(this->vec.size());
		for (size_t i = 0; i < this->vec.size(); ++i) {
			result.push_back(this->vec[i] + o->vec[i]);
		}
		return std::make_unique<TypedVector<T>>(result);
	}
	
    TypedVector operator+(const TypedVector& other) const {
        if (other.vec.size() != this->vec.size()) {
            throw std::runtime_error("Vector sizes do not match for addition.");
        }
        std::vector<T> result;
        result.reserve(this->vec.size());
        for (size_t i = 0; i < this->vec.size(); ++i) {
            result.push_back(this->vec[i] + other.vec[i]);
        }
        return TypedVector(result);
    }
	// added stuff
	TypedVector operator-(const TypedVector& other) const {
		if (other.vec.size() != this->vec.size()) {
			throw std::runtime_error("Vector sizes do not match for subtraction.");
		}
		std::vector<T> result;
		result.reserve(this->vec.size());
		for (size_t i = 0; i < this->vec.size(); ++i) {
			result.push_back(this->vec[i] + other.vec[i]);
		}
		return TypedVector(result);
	}


    TypedVector operator*(T scalar) const {
        std::vector<T> result;
        result.reserve(this->vec.size());
        for (const auto& item : this->vec) {
            result.push_back(item * scalar);
        }
        return TypedVector(result);
    }

    py::object toPython() const override {
        py::tuple tuple(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            tuple[i] = py::cast(vec[i]);
        }
        return tuple;
    }

    py::object getItem(size_t index) const override {
        if (index < vec.size()) {
            return py::cast(vec[index]);
        }
        throw py::index_error();
    }

    void setItem(size_t index, py::object value) override {
        if (index >= vec.size()) {
            throw py::index_error();
        }
        vec[index] = value.cast<T>();
    }

	std::unique_ptr<BaseVector> cross(const BaseVector* other) const override {
    	if (auto o = dynamic_cast<const TypedVector<int>*>(other)) {
        	auto result = crossProduct(this, o);
        	return std::unique_ptr<BaseVector>(std::move(result));
    	} else if (auto o = dynamic_cast<const TypedVector<double>*>(other)) {
        	auto result = crossProduct(this, o);
        	return std::unique_ptr<BaseVector>(std::move(result));
    	} else if (auto o = dynamic_cast<const TypedVector<std::complex<double>>*>(other)) {
        	auto result = crossProduct(this, o);
        	return std::unique_ptr<BaseVector>(std::move(result));
    	} else {
        	throw std::runtime_error("Unsupported vector type for cross product.");
    	}
	}

	py::object dot_product(const BaseVector* other) const override{
    	if (auto o = dynamic_cast<const TypedVector<int>*>(other)) {
        	auto result = dotProduct(this->vec, o->vec);
    		return py::cast(result);
		} else if (auto o = dynamic_cast<const TypedVector<double>*>(other)) {
        	auto result = dotProduct(this->vec, o->vec); 
    		return py::cast(result);
		} else if (auto o = dynamic_cast<const TypedVector<std::complex<double>>*>(other)) {
        	auto result = dotProduct(this->vec, o->vec); 
			return py::cast(result);
    	} else {
        	throw std::runtime_error("Unsupported vector type for dot product");
    	}
	}


	template<typename Type>
	std::unique_ptr<BaseVector> addByVec(Type type) const {
		using ResultType = CommonType<T, Type>;
		std::vector<ResultType> resultVec;
	}

	template<typename Scalar>
	std::unique_ptr<BaseVector> multiplyByScalar(Scalar scalar) const {
    	using ResultType = CommonType<T, Scalar>;
    	std::vector<ResultType> resultVec;
    	std::transform(this->vec.begin(), this->vec.end(), std::back_inserter(resultVec),
        	[scalar](const T& value) -> ResultType {
            	return static_cast<ResultType>(value) * static_cast<ResultType>(scalar);
        	});
    
    	return std::make_unique<TypedVector<ResultType>>(resultVec);
	}


	template<typename Scalar>
	std::unique_ptr<BaseVector> addByScalar(Scalar scalar) const {
		using ResultType = CommonType<T, Scalar>;
		std::vector<ResultType> resultVec;
		std::transform(this->vec.begin(), this->vec.end(), std::back_inserter(resultVec),
			[scalar](const T& value) -> ResultType {
				return static_cast<ResultType>(value) + static_cast<ResultType>(scalar);
			});
		return std::make_unique<TypedVector<ResultType>>(resultVec);
	}

	template<typename Scalar>
	std::unique_ptr<BaseVector> subByScalar(Scalar scalar) const {
		using ResultType = CommonType<T, Scalar>;
		std::vector<ResultType> resultVec;
		std::transform(this->vec.begin(), this->vec.end(), std::back_inserter(resultVec),
			[scalar](const T& value) -> ResultType {
				return static_cast<ResultType>(value) - static_cast<ResultType>(scalar);
			});
		return std::make_unique<TypedVector<ResultType>>(resultVec);
	}

	template<typename Scalar>
	std::unique_ptr<BaseVector> divByScalar(Scalar scalar) const {
		using ResultType = CommonType<T, Scalar>;
		std::vector<ResultType> resultVec;
		std::transform(this->vec.begin(), this->vec.end(), std::back_inserter(resultVec),
			[scalar](const T& value) -> ResultType {
				return static_cast<ResultType>(value) / static_cast<ResultType>(scalar);
		});
		return std::make_unique<TypedVector<ResultType>>(resultVec);
	}

};
// conversion template functions for operations 

template<typename From, typename To>
std::unique_ptr<TypedVector<To>> performConversion(const TypedVector<From>* source) {
    std::vector<To> convertedVec;
    std::transform(source->vec.begin(), source->vec.end(), std::back_inserter(convertedVec),
                   [](const From& val) -> To { return static_cast<To>(val); });
    return std::make_unique<TypedVector<To>>(convertedVec);
}

template<typename From, typename To>
std::unique_ptr<TypedVector<To>> convertIfNeeded(const TypedVector<From>* source);

template<typename T>
std::unique_ptr<TypedVector<T>> convertIfNeeded(const TypedVector<T>* source) {
    return std::make_unique<TypedVector<T>>(*source);
}

double roundToDecimal(double value, int decimalPlaces) {
    double scale = std::pow(10.0, decimalPlaces);
    return std::round(value * scale) / scale;
}

template<typename From, typename To>
std::unique_ptr<TypedVector<To>> convertIfNeeded(const TypedVector<From>* source) {
    std::vector<To> convertedVec;

    if constexpr (std::is_same<From, To>::value) {
        convertedVec.assign(source->vec.begin(), source->vec.end());
    } else {
    	std::transform(source->vec.begin(), source->vec.end(), std::back_inserter(convertedVec),
            [](const From& val) -> To {
                if constexpr (std::is_same<To, std::complex<double>>::value && 
                             (std::is_same<From, int>::value || std::is_same<From, double>::value)) {
                    return std::complex<double>(val, 0); 
                } else if constexpr (std::is_same<To, double>::value && std::is_same<From, int>::value) {
                    return static_cast<double>(val); 
                } else if constexpr (std::is_same<To, int>::value && std::is_same<From, double>::value) {
                    return static_cast<int>(val); 
                } else if constexpr (std::is_same<To, double>::value && std::is_same<From, std::complex<double>>::value) {
                    return val.real(); 
                } else if constexpr (std::is_same<To, int>::value && std::is_same<From, std::complex<double>>::value) {
                    return static_cast<int>(val.real()); 
                } else {
                    static_assert(std::is_same<From, void>::value, "Unsupported conversion type.");
                }
            });
	}

    return std::make_unique<TypedVector<To>>(TypedVector<To>{convertedVec});
}



template<typename T, typename U>
std::unique_ptr<BaseVector> crossProduct(const TypedVector<T>* a, const TypedVector<U>* b) {
//    using CommonType = decltype(T() * U());
	using CommonType = typename CommonTypeDeduction<T, U>::type;    
    auto aConverted = convertIfNeeded<T, CommonType>(a); 
    auto bConverted = convertIfNeeded<U, CommonType>(b); 

    std::vector<CommonType> result(3);
    result[0] = aConverted->vec[1] * bConverted->vec[2] - aConverted->vec[2] * bConverted->vec[1];
    result[1] = aConverted->vec[2] * bConverted->vec[0] - aConverted->vec[0] * bConverted->vec[2];
    result[2] = aConverted->vec[0] * bConverted->vec[1] - aConverted->vec[1] * bConverted->vec[0];

    return std::make_unique<TypedVector<CommonType>>(result); 
}


template<typename From, typename To>
std::unique_ptr<BaseVector> convertVector(const TypedVector<From>* source);


template<typename TargetType, typename SourceType>
std::unique_ptr<TypedVector<TargetType>> convertVectorToType(const TypedVector<SourceType>* source) {
    std::vector<TargetType> convertedVec(source->vec.begin(), source->vec.end());
    return std::make_unique<TypedVector<TargetType>>(convertedVec);
}


static std::unique_ptr<TypedVector<std::complex<double>>> convertToComplex(const BaseVector* vec) {
    if (auto* vInt = dynamic_cast<const TypedVector<int>*>(vec)) {
        std::vector<std::complex<double>> resultVec(vInt->vec.size());
        std::transform(vInt->vec.begin(), vInt->vec.end(), resultVec.begin(), 
                       [](int element) { return std::complex<double>(element); });
        return std::make_unique<TypedVector<std::complex<double>>>(resultVec);
    } else if (auto* vFloat = dynamic_cast<const TypedVector<double>*>(vec)) {
        std::vector<std::complex<double>> resultVec(vFloat->vec.size());
        std::transform(vFloat->vec.begin(), vFloat->vec.end(), resultVec.begin(), 
                       [](double element) { return std::complex<double>(element); });
        return std::make_unique<TypedVector<std::complex<double>>>(resultVec);
    }
    throw std::runtime_error("Unsupported vector type for conversion to complex.");
}

static std::unique_ptr<TypedVector<double>> convertToFloat(const BaseVector* vec) {
    if (auto* vInt = dynamic_cast<const TypedVector<int>*>(vec)) {
        std::vector<double> resultVec(vInt->vec.size());
        std::transform(vInt->vec.begin(), vInt->vec.end(), resultVec.begin(), 
                       [](int element) { return static_cast<double>(element); });
        return std::make_unique<TypedVector<double>>(resultVec);
    }
    throw std::runtime_error("Unsupported vector type for conversion to float.");
}

template<typename FromType, typename ToType>
std::unique_ptr<BaseVector> convertVector(const BaseVector* source) {
    throw std::runtime_error("Conversion not implemented for these types.");
}

template<>
std::unique_ptr<BaseVector> convertVector<int, double>(const BaseVector* source) {
    auto* srcTyped = dynamic_cast<const TypedVector<int>*>(source);
    if (!srcTyped) {
        throw std::runtime_error("Source vector type mismatch.");
    }
    std::vector<double> converted(srcTyped->vec.begin(), srcTyped->vec.end());
    return std::make_unique<TypedVector<double>>(converted);
}

template<>
std::unique_ptr<BaseVector> convertVector<double, std::complex<double>>(const BaseVector* source) {
    auto* srcTyped = dynamic_cast<const TypedVector<double>*>(source);
    if (!srcTyped) {
        throw std::runtime_error("Source vector type mismatch.");
    }
    std::vector<std::complex<double>> converted(srcTyped->vec.begin(), srcTyped->vec.end());
    return std::make_unique<TypedVector<std::complex<double>>>(converted);
}

template<typename Scalar, typename T>
TypedVector<decltype(std::declval<Scalar>() + std::declval<T>())> operator+(const Scalar& scalar, const TypedVector<T>& vec) {
    using PromotedType = decltype(scalar + T());
    std::vector<PromotedType> result;
    result.reserve(vec.vec.size());
    for (const auto& element : vec.vec) {
        result.push_back(scalar + element);
    }
    return TypedVector<PromotedType>(result);
}

// PyVector class returned to python binding for class operations returns tuple
class PyVector {
public:
    std::unique_ptr<BaseVector> baseVec;

    explicit PyVector(std::unique_ptr<BaseVector> vec) : baseVec(std::move(vec)) {}

	PyVector(py::args args) {
    if (args.empty()) {
        throw std::runtime_error("Arguments cannot be empty");
    }

    bool has_int = false, has_float = false, has_bool = false;

    for (auto item : args) {
        if (py::isinstance<py::bool_>(item)) {
            has_bool = true;
        } else if (py::isinstance<py::float_>(item)) {
            has_float = true;
        } else if (py::isinstance<py::int_>(item)) {
            has_int = true;
        } else {
            has_int = false;
            has_float = false;
            has_bool = false;
            break;
        }
    }

    if (has_bool) {
        std::vector<bool> values;
        for (auto item : args) {
            values.push_back(item.cast<bool>());
        }
        std::vector<int> int_values(values.begin(), values.end());
        baseVec = std::make_unique<TypedVector<int>>(int_values); 
    } else if (has_float) {
        std::vector<double> values;
        for (auto item : args) {
            values.push_back(item.cast<double>());
        }
        baseVec = std::make_unique<TypedVector<double>>(values);
    } else if (has_int) {
        std::vector<int> values;
        for (auto item : args) {
            values.push_back(item.cast<int>());
        }
        baseVec = std::make_unique<TypedVector<int>>(values);
    } else {
        std::vector<std::complex<double>> values;
        for (auto item : args) {
            if (py::isinstance<py::int_>(item) || py::isinstance<py::float_>(item)) {
                values.push_back(std::complex<double>(item.cast<double>(), 0));
            } else {
                values.push_back(item.cast<std::complex<double>>());
            }
        }
        baseVec = std::make_unique<TypedVector<std::complex<double>>>(values);
    }
}

	PyVector(py::handle input) {
    if (py::isinstance<py::args>(input) || py::isinstance<py::list>(input) || py::isinstance<py::tuple>(input)) {
        
        bool has_int = false, has_float = false, has_bool = false;

        for (auto item : input) {
            if (py::isinstance<py::bool_>(item)) {
                has_bool = true;
            } else if (py::isinstance<py::float_>(item)) {
                has_float = true;
            } else if (py::isinstance<py::int_>(item)) {
                has_int = true;
            } else {
                has_int = false;
                has_float = false;
                has_bool = false;
                break;
            }
        }

        if (has_bool) {
            std::vector<bool> values;
            for (auto item : input) {
                values.push_back(item.cast<bool>());
            }
            std::vector<int> int_values(values.begin(), values.end()); 
            baseVec = std::make_unique<TypedVector<int>>(int_values);
        }
        else if (has_float) {
            std::vector<double> values;
            for (auto item : input) {
                values.push_back(item.cast<double>());
            }
            baseVec = std::make_unique<TypedVector<double>>(values);
        }
        else if (has_int) {
            std::vector<int> values;
            for (auto item : input) {
                values.push_back(item.cast<int>());
            }
            baseVec = std::make_unique<TypedVector<int>>(values);
        }
        else {
            std::vector<std::complex<double>> values;
            for (auto item : input) {
                if (py::isinstance<py::int_>(item) || py::isinstance<py::float_>(item)) {
                    values.push_back(std::complex<double>(item.cast<double>(), 0));
                } else {
                    values.push_back(item.cast<std::complex<double>>());
                }
            }
            baseVec = std::make_unique<TypedVector<std::complex<double>>>(values);
        }
    } else {
        throw std::runtime_error("Unsupported input type for initialization. Expected args, list, or tuple.");
    }
}


	size_t size() const {
        if (!baseVec) throw std::runtime_error("Vector is uninitialized");
        return baseVec->size();
    }

	

	std::string getType() const {
        if (!baseVec) {
            return "Uninitialized";
        }
        if (dynamic_cast<TypedVector<int>*>(baseVec.get())) {
            return "int";
        } else if (dynamic_cast<TypedVector<double>*>(baseVec.get())) {
            return "double";
        } else if (dynamic_cast<TypedVector<std::complex<double>>*>(baseVec.get())) {
            return "std::complex<double>";
        }
        return "Unknown";
    }
	
	const std::unique_ptr<BaseVector>& getBaseVector() const {
    	return baseVec;
	}
	
	template<typename T>	
	bool isType() const {
        if (!baseVec) return false; 
        auto* castedVec = dynamic_cast<TypedVector<T>*>(baseVec.get());
        return castedVec != nullptr;
    }
	

	static PyVector deserializeCompressed(const std::string& compressedData) {
        std::vector<double> decompressedValues;
        
        for (size_t i = 0; i < compressedData.length(); i += 5) {
            if (i + 5 > compressedData.length()) throw std::runtime_error("Invalid compressed data length.");
            
            char s[5];
            std::copy(compressedData.begin() + i, compressedData.begin() + i + 5, s);
            double value = uncompress(s);
            decompressedValues.push_back(value);
        }

        return PyVector(std::make_unique<TypedVector<double>>(decompressedValues));
    }

	std::string serialize() const {
    	std::ostringstream oss;
    	if (baseVec) {
        	baseVec->serialize(oss);
    	}
    	return oss.str();
	}	
	std::string serialize2() const {
		std::ostringstream oss;
		if (baseVec) {
			baseVec->serialize2(oss);
		}
		return oss.str();
	}

	template<typename T>
    void initialize(size_t size, T fillValue) {
        std::vector<T> newVec(size, fillValue);
        baseVec = std::make_unique<TypedVector<T>>(newVec);
    }

	PyVector multiply_elements(const PyVector& other) const {
    	if (!baseVec) throw std::runtime_error("Vector is uninitialized");
    	auto resultVec = baseVec->multiply_elements(other.baseVec.get());
    	return PyVector(std::move(resultVec));
	}

	PyVector arg() const {
    	if (!baseVec) throw std::runtime_error("Vector is uninitialized");
    	auto resultVec = dynamic_cast<const TypedVector<std::complex<double>>*>(baseVec.get())->arg(); 
    	return PyVector(std::move(resultVec));
	}


	PyVector conj() const {
        if (!baseVec) throw std::runtime_error("Vector is uninitialized");
        auto resultVec = baseVec->conj(); 
        return PyVector(std::move(resultVec));
    }

	PyVector real() const {
		if (!baseVec) throw std::runtime_error("Vector is uninitalized");
		auto resultVec = baseVec->real();
		return PyVector(std::move(resultVec));
	}
	PyVector imag() const {
		if (!baseVec) throw std::runtime_error("Vector is unintialized");
		auto resultVec = baseVec->imag();
		return PyVector(std::move(resultVec));
	}
	
	double norm() const {
    if (!baseVec) throw std::runtime_error("Vector is uninitialized");
    	double sum = 0;
    	if (auto* vecInt = dynamic_cast<TypedVector<int>*>(baseVec.get())) {
        	sum = vecInt->norm();
    	} else if (auto* vecFloat = dynamic_cast<TypedVector<double>*>(baseVec.get())) {
        	sum = vecFloat->norm();
    	} else if (auto* vecComplex = dynamic_cast<TypedVector<std::complex<double>>*>(baseVec.get())) {
        	sum = vecComplex->norm(); 
    	} else {
        	throw std::runtime_error("Unsupported vector type for norm calculation.");
    	}

    	return std::round(sum * std::pow(10.0, 11)) / std::pow(10.0, 11);
	}


	template<typename T, typename U>
		static std::unique_ptr<BaseVector> dot_product(const TypedVector<T>& v1, const TypedVector<U>& v2) {
    	if (v1.vec.size() != v2.vec.size()) {
        	throw std::runtime_error("Vector sizes do not match for dot product.");
    	}

    	using ResultType = typename CommonTypeDeduction<T, U>::type;

    	ResultType result = std::inner_product(v1.vec.begin(), v1.vec.end(), v2.vec.begin(), ResultType(0),
        	std::plus<>(), [](const T& a, const U& b) { return ResultType(a) * ResultType(b); });

    	return result;
	}
	
	// add operations with two vectors
	static std::unique_ptr<BaseVector> addComplexVectors(const TypedVector<std::complex<double>>* a, const TypedVector<std::complex<double>>* b) {
    	if (a == NULL || b == NULL) {
       		throw std::runtime_error("Null vector passed to addComplexVectors");
    	}
    	if (a->vec.size() != b->vec.size()) {
        	throw std::runtime_error("Vector sizes do not match");
    	}

    	std::vector<std::complex<double>> resultVec(a->vec.size());
    	for (size_t i = 0; i < a->vec.size(); ++i) {
        	resultVec[i] = a->vec[i] + b->vec[i];
    	}
    	return std::make_unique<TypedVector<std::complex<double>>>(resultVec);
	}

	static std::unique_ptr<BaseVector> addFloatVectors(const TypedVector<double>* a, const TypedVector<double>* b) {
    	if (a == NULL|| b == NULL) {
        	throw std::runtime_error("Null vector passed to addFloatVectors");
    	}
    	if (a->vec.size() != b->vec.size()) {
        	throw std::runtime_error("Vector sizes do not match");
    	}

    	std::vector<double> resultVec(a->vec.size());
    	for (size_t i = 0; i < a->vec.size(); ++i) {
        	resultVec[i] = a->vec[i] + b->vec[i];
    	}
    	return std::make_unique<TypedVector<double>>(resultVec);
	}

	static std::unique_ptr<BaseVector> addIntVectors(const TypedVector<int>* a, const TypedVector<int>* b) {
    	if (a == NULL || b == NULL) {
        	throw std::runtime_error("Null vector passed to addIntVectors");
    	}
    	if (a->vec.size() != b->vec.size()) {
        	throw std::runtime_error("Vector sizes do not match");
    	}

    	std::vector<int> resultVec(a->vec.size());
    	for (size_t i = 0; i < a->vec.size(); ++i) {
        	resultVec[i] = a->vec[i] + b->vec[i];
    	}
    	return std::make_unique<TypedVector<int>>(resultVec);
	}

	// subtraction operations with two vectors
	static std::unique_ptr<BaseVector> subComplexVectors(const TypedVector<std::complex<double>>* a, const TypedVector<std::complex<double>>* b) {
		if (a == NULL || b == NULL) {
			throw std::runtime_error("Null vector passed to subComplexVectors");
		}
		std::vector<std::complex<double>> resultVec(a->vec.size());
		for (size_t i = 0; i < a->vec.size(); ++i) {
			resultVec[i] = a->vec[i] - b->vec[i];
		}
		return std::make_unique<TypedVector<std::complex<double>>>(resultVec);
	}

	static std::unique_ptr<BaseVector> subFloatVectors(const TypedVector<double>* a, const TypedVector<double>* b) {
		if (a == NULL || b == NULL) {
			throw std::runtime_error("Null vector passed to subComplexVectors");
		}
		if (a->vec.size() != b->vec.size()) {
			throw std::runtime_error("Vector sizes do not match");
		}
		std::vector<double> resultVec(a->vec.size());
		for (size_t i = 0; i < a->vec.size(); ++i) {
			resultVec[i] = a->vec[i] - b->vec[i];
		}
		return std::make_unique<TypedVector<double>>(resultVec);
	}
	

	static std::unique_ptr<BaseVector> subIntVectors(const TypedVector<int>* a, const TypedVector<int>* b) {
		if (a == NULL || b == NULL) {
			throw std::runtime_error("Null vector passed to subIntVectors");
		}
		if (a->vec.size() != b->vec.size()) {
			throw std::runtime_error("Vector sizes do not match");
		}
		std::vector<double> resultVec(a->vec.size());
		for (size_t i = 0; i < a->vec.size(); ++i) {
			resultVec[i] = a->vec[i] - b->vec[i];
		}
		return std::make_unique<TypedVector<double>>(resultVec);
	}


	// division operations with two vectors
	static std::unique_ptr<BaseVector> divComplexVectors(const TypedVector<std::complex<double>>* a, const TypedVector<std::complex<double>>* b) {
		std::cout << a;
		std::cout << b;
		if (a == NULL || b == NULL) {
			throw std::runtime_error("Null vector passed to divComplexVectors");
		}
		if (a->vec.size() != b->vec.size()) {
			throw std::runtime_error("Vector sizes do not match");
		}
		std::vector<std::complex<double>> resultVec(a->vec.size());
		for (size_t i = 0; i < a->vec.size(); ++i) {
			resultVec[i] = a->vec[i] / b->vec[i];
		}
		return std::make_unique<TypedVector<std::complex<double>>>(resultVec);
	}
	
	static std::unique_ptr<BaseVector> divFloatVectors(const TypedVector<double>* a, const TypedVector<double>* b) {
		if (a == NULL || b == NULL) {
			throw std::runtime_error("Null vector passed to divFloatVectors");
		}
		if (a->vec.size() != b->vec.size()) {
			throw std::runtime_error("Vector sizes do not match");
		}
		std::vector<double> resultVec(a->vec.size());
		for (size_t i = 0; i < a->vec.size(); ++i) {
			resultVec[i] = a->vec[i] / b->vec[i];
		}
		return std::make_unique<TypedVector<double>>(resultVec);
	}

	static std::unique_ptr<BaseVector> divIntVectors(const TypedVector<int>* a, const TypedVector<int>* b) {
		if (a == NULL || b == NULL) {
			throw std::runtime_error("Null vector passed to divIntVectors");
		}
		if (a->vec.size() != b->vec.size()) {
			throw std::runtime_error("Vector sizes do not match");
		}
		std::vector<double> resultVec(a->vec.size());
		for (size_t i = 0; i < a->vec.size(); ++i) {
			resultVec[i] = a->vec[i] / b->vec[i];
		}
		return std::make_unique<TypedVector<double>>(resultVec);
	}

	static std::unique_ptr<TypedVector<std::complex<double>>> convertToComplexIfNeeded(const BaseVector* vec) {
    	if (auto* v = dynamic_cast<const TypedVector<std::complex<double>>*>(vec)) {
        	return std::make_unique<TypedVector<std::complex<double>>>(*v); 
    	} else if (auto* v = dynamic_cast<const TypedVector<double>*>(vec)) {
        	return convertToComplex(v);
    	} else if (auto* v = dynamic_cast<const TypedVector<int>*>(vec)) {
        	return convertToComplex(v);
    	}
    	throw std::runtime_error("Unsupported vector type for conversion to complex.");
	}

	static std::unique_ptr<TypedVector<double>> convertToFloatIfNeeded(const BaseVector* vec) {
    	if (auto* v = dynamic_cast<const TypedVector<double>*>(vec)) {
        	return std::make_unique<TypedVector<double>>(*v); 
    	} else if (auto* v = dynamic_cast<const TypedVector<int>*>(vec)) {
        	return convertToFloat(v);
    	}
    	throw std::runtime_error("Unsupported vector type for conversion to float.");
	}

	static std::unique_ptr<BaseVector> addVectors(const BaseVector* a, const BaseVector* b) {
    	if (dynamic_cast<const TypedVector<std::complex<double>>*>(a) || dynamic_cast<const TypedVector<std::complex<double>>*>(b)) {
        	auto aComplex = convertToComplexIfNeeded(a);
        	auto bComplex = convertToComplexIfNeeded(b);
        	return addComplexVectors(aComplex.get(), bComplex.get());
    	} else if (dynamic_cast<const TypedVector<double>*>(a) || dynamic_cast<const TypedVector<double>*>(b)) {
        	auto aFloat = convertToFloatIfNeeded(a);
        	auto bFloat = convertToFloatIfNeeded(b);
        	return addFloatVectors(aFloat.get(), bFloat.get());
    	} else {
        	const auto* aInt = dynamic_cast<const TypedVector<int>*>(a);
        	const auto* bInt = dynamic_cast<const TypedVector<int>*>(b);
        	return addIntVectors(aInt, bInt);
    	}
	}

	static std::unique_ptr<BaseVector> subVectors(const BaseVector* a, const BaseVector* b) {
		if (dynamic_cast<const TypedVector<std::complex<double>>*>(a) || dynamic_cast<const TypedVector<std::complex<double>>*>(b)) {
			auto aComplex = convertToComplexIfNeeded(a);
			auto bComplex = convertToComplexIfNeeded(b);
			return subComplexVectors(aComplex.get(), bComplex.get());
		} else if (dynamic_cast<const TypedVector<double>*>(a) || dynamic_cast<const TypedVector<double>*>(b)) {
			auto aFloat = convertToFloatIfNeeded(a);
			auto bFloat = convertToFloatIfNeeded(b);
			return subFloatVectors(aFloat.get(), bFloat.get());
		} else {
			const auto* aInt = dynamic_cast<const TypedVector<int>*>(a);
			const auto* bInt = dynamic_cast<const TypedVector<int>*>(b);
			return subIntVectors(aInt, bInt);
		}
	}

	static std::unique_ptr<BaseVector> divVectors(const BaseVector* a, const BaseVector* b) {
		if (dynamic_cast<const TypedVector<std::complex<double>>*>(a) || dynamic_cast<const TypedVector<std::complex<double>>*>(b)) {
			auto aComplex = convertToComplexIfNeeded(a);
			auto bComplex = convertToComplexIfNeeded(a);
			return divComplexVectors(aComplex.get(), bComplex.get());
		} else if (dynamic_cast<const TypedVector<double>*>(a) || dynamic_cast<const TypedVector<double>*>(b)) {
			auto aFloat = convertToFloatIfNeeded(a);
			auto bFloat = convertToFloatIfNeeded(b);
			return divFloatVectors(aFloat.get(), bFloat.get());
		} else {
			const auto* aInt = dynamic_cast<const TypedVector<int>*>(a);
			const auto* bInt = dynamic_cast<const TypedVector<int>*>(b);
			return divIntVectors(aInt, bInt);
		}
	}
	py::object __getitem__(int index) const {
    	if (!baseVec) throw py::index_error("Vector is uninitialized.");
		size_t adjusted_index = static_cast<size_t>(index);
    	if (index < 0) {
    		adjusted_index = baseVec->size() + index;
		}

		if (adjusted_index >= baseVec->size()) {
        	throw py::index_error("Index out of range.");
    	}

    	if (auto vecInt = dynamic_cast<const TypedVector<int>*>(baseVec.get())) {
        	return py::cast(vecInt->vec[adjusted_index]);
    	} else if (auto vecFloat = dynamic_cast<const TypedVector<double>*>(baseVec.get())) {
        	return py::cast(vecFloat->vec[adjusted_index]);
    	} else if (auto vecComplex = dynamic_cast<const TypedVector<std::complex<double>>*>(baseVec.get())) {
        	return py::cast(vecComplex->vec[adjusted_index]);
    	} else {
    	    throw std::runtime_error("Unsupported vector type for indexing.");
 	   }
	}

	

	void __setitem__(int index, py::object value) {
    	if (!baseVec) throw py::index_error("Vector is uninitialized.");
    		if (auto vecInt = dynamic_cast<TypedVector<int>*>(baseVec.get())) {
        		if (index < 0 || index >= vecInt->vec.size()) throw py::index_error("Index out of range.");
        		vecInt->vec[index] = py::cast<int>(value);
    		} else if (auto vecFloat = dynamic_cast<TypedVector<double>*>(baseVec.get())) {
        		if (index < 0 || index >= vecFloat->vec.size()) throw py::index_error("Index out of range.");
        		vecFloat->vec[index] = py::cast<double>(value);
    		} else if (auto vecComplex = dynamic_cast<TypedVector<std::complex<double>>*>(baseVec.get())) {
        		if (index < 0 || index >= vecComplex->vec.size()) throw py::index_error("Index out of range.");
        		vecComplex->vec[index] = py::cast<std::complex<double>>(value);
    		} else {
        		throw std::runtime_error("Unsupported vector type for indexing.");
    		}
	}

	
	template<typename Scalar>
    PyVector multiply(Scalar scalar) const {
        if (!baseVec) throw std::runtime_error("Vector is uninitialized");

        if (auto* intVec = dynamic_cast<const TypedVector<int>*>(baseVec.get())) {
            auto resultVec = intVec->multiplyByScalar(scalar);
            return PyVector(std::move(resultVec));
        } else if (auto* floatVec = dynamic_cast<const TypedVector<double>*>(baseVec.get())) {
			auto resultVec = floatVec->multiplyByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else if (auto* complexVec = dynamic_cast<const TypedVector<std::complex<double>>*>(baseVec.get())) {
			auto resultVec = complexVec->multiplyByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else {
            throw std::runtime_error("Scalar multiplication logic not implemented.");
        }
    }
    
	template<typename Scalar>
	PyVector add(Scalar scalar) const {
		if (!baseVec) throw std::runtime_error("Vector is uninitialized");
	
		if (auto* intVec = dynamic_cast<const TypedVector<int>*>(baseVec.get())) {
			auto resultVec = intVec->addByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else if (auto* floatVec = dynamic_cast<const TypedVector<double>*>(baseVec.get())) {
			auto resultVec = floatVec->addByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else if (auto* complexVec = dynamic_cast<const TypedVector<std::complex<double>>*>(baseVec.get())) {
			auto resultVec = complexVec->addByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else {
			throw std::runtime_error("Scalar addition logic not implemented");
		}
	}

	template<typename Scalar>
	PyVector subtract(Scalar scalar) const {
		if (!baseVec) throw std::runtime_error("Vector is unintialized");

		if (auto* intVec = dynamic_cast<const TypedVector<int>*>(baseVec.get())) {
			auto resultVec = intVec->subByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else if (auto* floatVec = dynamic_cast<const TypedVector<double>*>(baseVec.get())) {
			auto resultVec = floatVec->subByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else if (auto* complexVec = dynamic_cast<const TypedVector<std::complex<double>>*>(baseVec.get())) {
			auto resultVec = complexVec->subByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else {
			throw std::runtime_error("Scalar subtraction logic not implemented");
		}
	}

	template<typename Scalar>
	PyVector divide(Scalar scalar) const {
		if (!baseVec) throw std::runtime_error("Vector is uninitalized");
		if (auto* intVec = dynamic_cast<const TypedVector<int>*>(baseVec.get())) {
			auto resultVec = intVec->divByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else if (auto* floatVec = dynamic_cast<const TypedVector<double>*>(baseVec.get())) {
			auto resultVec = floatVec->divByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else if (auto* complexVec = dynamic_cast<const TypedVector<std::complex<double>>*>(baseVec.get())) {
			auto resultVec = complexVec->divByScalar(scalar);
			return PyVector(std::move(resultVec));
		} else {
			throw std::runtime_error("Scalar subtraction logic not implemented");
		}
	}
	friend PyMatrix operator%(const PyVector& u, const PyVector& v);


	static PyVector createFromBaseVector(BaseVector* baseVectorPtr) {
        return PyVector(std::unique_ptr<BaseVector>(baseVectorPtr));
    }

	template<typename T>
    std::vector<T> extractDataAs() const {
        throw std::runtime_error("Unsupported type for extraction");
    }

    py::object __repr__() const {
        return baseVec->toPython();
    }
};

template<>
std::vector<int> PyVector::extractDataAs<int>() const {
    if (getType() != "int") throw std::runtime_error("Attempting to extract int data from non-int vector");
    auto* vecInt = dynamic_cast<TypedVector<int>*>(baseVec.get());
    if (!vecInt) throw std::runtime_error("Failed to extract int data");
    return vecInt->vec;
}

template<>
std::vector<double> PyVector::extractDataAs<double>() const {
    if (getType() != "double") throw std::runtime_error("Attempting to extract double data from non-double vector");
    auto* vecDouble = dynamic_cast<TypedVector<double>*>(baseVec.get());
    if (!vecDouble) throw std::runtime_error("Failed to extract double data");
    return vecDouble->vec;
}


template<>
std::vector<std::complex<double>> PyVector::extractDataAs<std::complex<double>>() const {
    if (getType() != "std::complex<double>") throw std::runtime_error("Attempting to extract std::complex<double> data from non-std::complex<double> vector");
    auto* vecComplex = dynamic_cast<TypedVector<std::complex<double>>*>(baseVec.get());
    if (!vecComplex) throw std::runtime_error("Failed to extract std::complex<double> data");
    return vecComplex->vec;
}

#include "wrap_matrix.cpp"

PyMatrix operator%(const PyVector& u, const PyVector& v) {
    size_t uSize = u.size();
    size_t vSize = v.size();
    
    PyMatrix resultMatrix(uSize, vSize);  

    for (size_t i = 0; i < uSize; ++i) {
        for (size_t j = 0; j < vSize; ++j) {
            auto uVal = u.__getitem__(i).cast<double>(); 
            auto vVal = v.__getitem__(j).cast<double>(); 

            std::visit([&](auto& matrixPtr) {
                using MatrixType = typename std::decay<decltype(*matrixPtr)>::type;
                using ValueType = typename MatrixType::ValueType;

                matrixPtr->setValueAt(i, j, uVal * vVal);
            }, resultMatrix.getMatrixVariant());
        }
    }

    return resultMatrix;
}


// pybinding init function 
void init_vector(py::module_ &m) {
    py::class_<PyVector>(m, "Vector")
        //.def(py::init<py::args>())
		.def(py::init<py::handle>())
		.def(py::init<py::args>())
		.def("__len__", &PyVector::size)
        .def("__getitem__", &PyVector::__getitem__)
        .def("__setitem__", &PyVector::__setitem__)
		.def("__mod__", [](const PyVector& lhs, const PyVector& rhs) {
            return lhs % rhs; 
        }, "Calculate the outer product of two vectors")
		.def("__mul__", [](const PyVector& self, const PyVector& other){
			auto result = self.baseVec->dot_product(other.baseVec.get());
    		return result;
		}, py::is_operator()) 
		.def("__xor__", [](const PyVector& self, const PyVector& other) {
            auto resultVec = self.baseVec->cross(other.baseVec.get());
            return PyVector(std::move(resultVec));
        }, py::is_operator())

		.def("__add__", [](const PyVector& self, const PyVector& other) {
            return PyVector::createFromBaseVector(PyVector::addVectors(self.baseVec.get(), other.baseVec.get()).release());
        })

		.def("__sub__", [](const PyVector& self, const PyVector& other) {
			return PyVector::createFromBaseVector(PyVector::subVectors(self.baseVec.get(), other.baseVec.get()).release());
		})

		.def("__truediv__", [](const PyVector& self, const PyVector& other) {
			return PyVector::createFromBaseVector(PyVector::divVectors(self.baseVec.get(), other.baseVec.get()).release());
		})

		.def("__eq__", [](const PyVector& self, const PyVector& other) -> bool {
    		if (typeid(*self.baseVec) != typeid(*other.baseVec)) {
        		return false;
    		}

    		if (self.baseVec->size() != other.baseVec->size()) {
        		return false;
    		}

    		for (size_t i = 0; i < self.baseVec->size(); ++i) {
        		py::handle selfItem = self.baseVec->getItem(i);
        		py::handle otherItem = other.baseVec->getItem(i);

       			if (!selfItem.equal(otherItem)) {
            		return false; 
        			}
    			}

    			return true;
		}, py::is_operator())		

		.def("__mul__", [](const PyVector& self, py::handle scalar) {
    		if (py::isinstance<py::int_>(scalar)) {
        		return self.multiply(scalar.cast<int>());
    		} else if (py::isinstance<py::float_>(scalar)) {
        		return self.multiply(scalar.cast<double>());
    		} else if (py::isinstance<py::object>(scalar) && py::hasattr(scalar, "real") && py::hasattr(scalar, "imag")) {
        		return self.multiply(scalar.cast<std::complex<double>>());
    		} else {
        		throw std::runtime_error("Unsupported scalar type for multiplication.");
    		}	
		}, py::is_operator())

		.def("__add__", [](const PyVector& self, py::handle scalar) {
			if (py::isinstance<py::int_>(scalar)) {
				return self.add(scalar.cast<int>());
			} else if (py::isinstance<py::float_>(scalar)) {
				return self.add(scalar.cast<double>());
			} else if (py::isinstance<py::object>(scalar) && py::hasattr(scalar, "real") && py::hasattr(scalar, "imag")) {
				return self.add(scalar.cast<std::complex<double>>());
			} else {
				throw std::runtime_error("Unsupported scalar type for addition.");
			}
		}, py::is_operator()) 

		.def("__sub__", [](const PyVector& self, py::handle scalar) {
			if (py::isinstance<py::int_>(scalar)) {
				return self.subtract(scalar.cast<int>());
			} else if (py::isinstance<py::float_>(scalar)) {
				return self.subtract(scalar.cast<double>());
			} else if (py::isinstance<py::object>(scalar) && py::hasattr(scalar, "real") && py::hasattr(scalar, "imag")) {
				return self.add(scalar.cast<std::complex<double>>());
			} else {
				throw std::runtime_error("Unsupported scalar type for addition.");
			}
		}, py::is_operator())	
		.def("fill_value", &PyVector::initialize<int>, "Initalize vector with size and int fill value")
    	.def("fill_value", &PyVector::initialize<double>, "Initalize vector with size and float fill value")
		.def("fill_value", &PyVector::initialize<std::complex<double>>, "Initialize vector with size and complex<float> fill value")
		.def("multiply_elements", &PyVector::multiply_elements, "Performs element-wise multiplication with another vector")

		.def("__truediv__", [](const PyVector& self, py::handle scalar) {
			if (py::isinstance<py::int_>(scalar)) {
				return self.divide(scalar.cast<int>());
			} else if (py::isinstance<py::float_>(scalar)) {
				return self.divide(scalar.cast<double>());
			} else if (py::isinstance<py::object>(scalar) && py::hasattr(scalar, "real") && py::hasattr(scalar, "imag")) {
				return self.divide(scalar.cast<std::complex<double>>());
			} else {
				throw std::runtime_error("Unsupported scalar type for division.");
			}
		}, py::is_operator())
		.def("apply", [](PyVector& self, py::function func) {
            if (auto pIntVec = dynamic_cast<TypedVector<int>*>(self.baseVec.get())) {
                auto resultVec = std::make_unique<TypedVector<int>>(pIntVec->vec); 
                for (auto& val : resultVec->vec) {
                    py::object result = func(val);
                    val = result.cast<int>(); 
                }
                return PyVector(std::move(resultVec));
            }
            else if (auto pFloatVec = dynamic_cast<TypedVector<double>*>(self.baseVec.get())) {
                auto resultVec = std::make_unique<TypedVector<double>>(pFloatVec->vec); 
                for (auto& val : resultVec->vec) {
                    py::object result = func(val);
                    val = result.cast<double>(); 
                }
                return PyVector(std::move(resultVec));
            }
            else if (auto pComplexVec = dynamic_cast<TypedVector<std::complex<double>>*>(self.baseVec.get())) {
                auto resultVec = std::make_unique<TypedVector<std::complex<double>>>(pComplexVec->vec); 
                for (auto& val : resultVec->vec) {
                    py::object result = func(val);
                    val = result.cast<std::complex<double>>(); 
                }
                return PyVector(std::move(resultVec));
            }
            else {
                throw std::runtime_error("Unsupported vector type for apply operation");
            }
         }, "Applies a function to each element of the vector and returns a new vector of the same type with the results.")
		.def("put", &PyVector::serialize, "Serializes the vector to a string.")
		.def("put", &PyVector::serialize2, "Serializes the vector to a string.")

		.def("get", &PyVector::deserializeCompressed, "Deserializes a vector from a compressed string format.")
		

		.def("conj", &PyVector::conj, "Computes the conjugate of the vector.")
		.def("arg", &PyVector::arg, "Computes the argument of each element in the complex vector.")
		.def("norm", &PyVector::norm, "Calculate the Euclidean norm of the vector.")
		.def("real", &PyVector::real, "Computes the real of the vector")	
		.def("imag", &PyVector::imag, "Computes the imag of the vector")	
		.def("__repr__", [](const PyVector &p) -> std::string {
    		auto pyObj = p.baseVec->toPython(); 
    		if (py::isinstance<py::tuple>(pyObj)) {
        		py::tuple tuple = pyObj.cast<py::tuple>();
        		std::stringstream ss;
        		ss << '(';
        		for (size_t i = 0; i < tuple.size(); ++i) {
            		if (i > 0) ss << ", ";
            		auto element = tuple[i];
            		if (py::isinstance<py::float_>(element)) {
                		double value = element.cast<double>();
            			ss << value;
					} else {
                		ss << py::str(element).cast<std::string>();
            		}
        		}
        		if (tuple.size() == 1) ss << ','; 
        		ss << ')';
        		return ss.str();
    		}
    		return "Invalid vector";
		});
	}
