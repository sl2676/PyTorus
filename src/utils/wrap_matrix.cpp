#pragma once
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <sstream>
#include <complex>
#include <type_traits>
#include <complex>
#include <iostream>
#include <memory>
#include <variant>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "./wrap_vector.cpp" 
template<typename T>
struct TypeIndex;

template<> struct TypeIndex<int> { static constexpr size_t value = 0; };
template<> struct TypeIndex<double> { static constexpr size_t value = 1; };
template<> struct TypeIndex<std::complex<double>> { static constexpr size_t value = 2; };
template<typename T, bool isSupported>

struct ValidateType {
    static_assert(isSupported, "Unsupported type");
};

template<typename T>
T WDabs(const T& value) {
    if constexpr (std::is_same_v<T, std::complex<double>>) {
        return std::abs(value);
    } else {
        return std::abs(value);
    }
}

template<typename T>
void WDswap(T& a, T& b) {
    T temp = a;
    a = b;
    b = temp;
}

template<typename T>
struct ValidateType<T, true> {};

template<typename T>
constexpr size_t getTypeIndex() {
    if constexpr (std::is_same_v<T, int>) {
        return 0;
    } else if constexpr (std::is_same_v<T, double>) {
        return 1;
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        return 2;
    } else {
        ValidateType<T, false>{}; 
        return std::numeric_limits<size_t>::max(); 
    }
}

template<typename TargetType, typename SourceType>
TargetType convertValue(const SourceType& value) {
    return static_cast<TargetType>(value);
}


template<>
double convertValue<double, std::complex<double>>(const std::complex<double>& value) {
    return value.real(); 
}

template<>
int convertValue<int, std::complex<double>>(const std::complex<double>& value) {
    return static_cast<int>(value.real()); 
}

template<>
std::complex<double> convertValue<std::complex<double>, int>(const int& value) {
    return {static_cast<double>(value), 0.0};
}

template<>
std::complex<double> convertValue<std::complex<double>, double>(const double& value) {
    return {value, 0.0};
}


namespace py = pybind11;




template<typename T>
class MatrixImpl;

using MatrixVariant = std::variant<
    std::shared_ptr<MatrixImpl<int>>,
    std::shared_ptr<MatrixImpl<double>>,
    std::shared_ptr<MatrixImpl<std::complex<double>>>
>;



class BaseMatrix {
public:
	virtual ~BaseMatrix() {}
    virtual MatrixVariant variant() const = 0;
	virtual std::string toString() const = 0;
	virtual size_t getRows() const = 0;
	virtual size_t getCols() const = 0;
};

template<typename T>
class MatrixImpl : public BaseMatrix {
private:
    std::vector<std::vector<T>> data;
	size_t rows, cols;
public:
    /*
		CONSTRUCTORS
	*/
	using ValueType = T;
	MatrixImpl() : data() {}

	MatrixImpl(size_t rows, size_t cols, T defaultValue = T())
		: data(rows, std::vector<T>(cols, defaultValue)), rows(rows), cols(cols) {}
    MatrixImpl(const MatrixImpl& M) : data(M.data) {}
	
	size_t getRows() const override {
        return data.size();
    }

    size_t getCols() const override {
        return data[0].size();
    }
	T& getValueAt(size_t row, size_t col) {
    	return data[row][col];
	}

	T getValueAt(size_t row, size_t col) const {
    	return data[row][col];
	}

	void setValueAt(size_t row, size_t col, T value) {
    	data[row][col] = value;
	}
	void setColumnFromStdVector(size_t colIndex, const std::vector<T>& columnVector) {
        if (colIndex >= cols) {
            throw std::out_of_range("Column index is out of bounds.");
        }
        if (columnVector.size() != rows) {
            throw std::invalid_argument("Vector size does not match the number of matrix rows.");
        }

        for (size_t row = 0; row < rows; ++row) {
            data[row][colIndex] = columnVector[row];
        }
    }

	void setRowFromStdVector(size_t rowIndex, const std::vector<T>& rowVector) {
		if (rowIndex >= rows) {
			throw std::out_of_range("Row index is out of bounds.");
		}
		if (rowVector.size() != cols) {
			throw std::invalid_argument("Vector size does not match the number of matrix rows.");
		}
		for (size_t col = 0; col < cols; ++col) {
			data[rowIndex][col] = rowVector[col];
		}
	}

	MatrixVariant variant() const override {
		return std::make_shared<MatrixImpl<T>>(*this);
	}

	std::shared_ptr<MatrixImpl<T>> transpose() const {
        auto transposedMatrix = std::make_shared<MatrixImpl<T>>(cols, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                transposedMatrix->setValueAt(j, i, this->getValueAt(i, j));
            }
        }
        return transposedMatrix;
    }

	// ERROR HANDLING

	static void range_error() {
		throw std::runtime_error("Range error");
	}
	static void error(const char* message) {
		throw std::runtime_error(message);
	}
	static void division_by_zero_error() {
		throw std::runtime_error("Divison by zero error");
	}
	
	std::string toString() const override {
        std::ostringstream oss;
        for (const auto& row : data) {
            for (const auto& elem : row) {
                oss << elem << " ";
            }
            oss << "\n";
        }
        return oss.str();
    }
};



class PyMatrix {
private:
    MatrixVariant matrixVariant;
    
    
    
    explicit PyMatrix(MatrixVariant variant) : matrixVariant(std::move(variant)) {}
public: 
	size_t getRows() const {
        return std::visit([](const auto& matrixImpl) -> size_t {
            return matrixImpl->getRows();
        }, matrixVariant);
    }

    size_t getCols() const {
        return std::visit([](const auto& matrixImpl) -> size_t {
            return matrixImpl->getCols();
        }, matrixVariant);
    }

	MatrixVariant deepCopyVariant(const MatrixVariant& sourceVariant) {
	    return std::visit([](const auto& sourceMatrix) -> MatrixVariant {
    	    using MatrixType = typename std::decay<decltype(*sourceMatrix)>::type;
        	using ValueType = typename MatrixType::ValueType;
        	auto copiedMatrix = std::make_shared<MatrixImpl<ValueType>>(sourceMatrix->getRows(), sourceMatrix->getCols());
        	for (size_t row = 0; row < sourceMatrix->getRows(); ++row) {
            	for (size_t col = 0; col < sourceMatrix->getCols(); ++col) {
                	copiedMatrix->setValueAt(row, col, sourceMatrix->getValueAt(row, col));
            	}
        	}
        	return copiedMatrix;
    	}, sourceVariant);
	}

	
    template<typename TargetType>
	void promoteMatrixVariantIfNeeded() {
    	size_t currentTypeIndex = std::visit([](const auto& arg) -> size_t {
        	using MatrixType = typename std::decay_t<decltype(*arg)>;
        	return getTypeIndex<typename MatrixType::ValueType>();
    	}, matrixVariant);

    	size_t targetTypeIndex = getTypeIndex<TargetType>();

    	if (currentTypeIndex < targetTypeIndex) {
        	promoteMatrixVariant<TargetType>();
    	}
	}
	
	py::object getValueAt(int row, int col) const {
        return std::visit([row, col](const auto& arg) -> py::object {
            return py::cast(arg->getValueAt(row, col));
        }, matrixVariant);
    }
	
	template<typename TargetType>
	void promoteMatrixVariant() {
    	MatrixVariant newMatrixVariant = std::visit([](auto& oldMatrix) -> MatrixVariant {
        	using OldValueType = typename std::decay_t<decltype(*oldMatrix)>::ValueType;

        	size_t rows = oldMatrix->getRows();
        	size_t cols = oldMatrix->getCols();

        	auto newMatrix = std::make_shared<MatrixImpl<TargetType>>(rows, cols);

        	for (size_t row = 0; row < rows; ++row) {
            	for (size_t col = 0; col < cols; ++col) {
                	OldValueType oldValue = oldMatrix->getValueAt(row, col);
                	newMatrix->setValueAt(row, col, convertValue<TargetType, OldValueType>(oldValue));
            	}
        	}

        	return newMatrix;
    	}, matrixVariant);

    	matrixVariant = std::move(newMatrixVariant);
	}
	template<typename Scalar, typename Function>
	PyMatrix& elementWiseOperation(Scalar scalar, Function func) {
    	std::visit([&](auto& matrixPtr) {
        	for (size_t row = 0; row < matrixPtr->getRows(); ++row) {
            	for (size_t col = 0; col < matrixPtr->getCols(); ++col) {
                	func(matrixPtr->getValueAt(row, col), scalar);
            	}
        	}
    	}, matrixVariant);
    	return *this;
	}

	template<typename LhsMatrixPtr, typename RhsMatrixPtr>
    static MatrixVariant matrixMultiplyZ(const LhsMatrixPtr& lhsMatrix, const RhsMatrixPtr& rhsMatrix) {
        using LhsValueType = typename std::decay_t<decltype(*lhsMatrix)>::ValueType;
        using RhsValueType = typename std::decay_t<decltype(*rhsMatrix)>::ValueType;
        using ResultValueType = typename std::common_type<LhsValueType, RhsValueType>::type;

        if (lhsMatrix->getCols() != rhsMatrix->getRows()) {
            throw std::runtime_error("Matrix dimensions mismatch for multiplication.");
        }

        auto resultMatrix = std::make_shared<MatrixImpl<ResultValueType>>(lhsMatrix->getRows(), rhsMatrix->getCols());

        for (size_t i = 0; i < lhsMatrix->getRows(); ++i) {
            for (size_t j = 0; j < rhsMatrix->getCols(); ++j) {
                ResultValueType sum = ResultValueType();
                for (size_t k = 0; k < lhsMatrix->getCols(); ++k) {
                    LhsValueType lhsValue = lhsMatrix->getValueAt(i, k);
                    RhsValueType rhsValue = rhsMatrix->getValueAt(k, j);
                    if (lhsValue != LhsValueType(0) && rhsValue != RhsValueType(0)) {
                        sum += static_cast<ResultValueType>(lhsValue) * static_cast<ResultValueType>(rhsValue);
                    }
                }
                resultMatrix->setValueAt(i, j, sum);
            }
        }

        return resultMatrix;
    }
	
	


public:
	class RowProxy {
    private:
        PyMatrix& matrix;
        size_t rowIndex;

    public:
        RowProxy(PyMatrix& matrix, size_t rowIndex) : matrix(matrix), rowIndex(rowIndex) {}

		py::object operator[](size_t colIndex) const {
            if (colIndex >= matrix.getCols()) {
                throw std::out_of_range("Column index out of bounds.");
            }
            return matrix.getValueAt(rowIndex, colIndex);
        }

		void set(size_t colIndex, const py::object& value) {
        	matrix.set_value(rowIndex, colIndex, value);
    	}

		void operator=(const py::object& value) {
			for (size_t col = 0; col < matrix.getCols(); ++col) {
				matrix.set_value(rowIndex, col, value);
			}			
		}
		void operator=(const PyVector& vec) {
			matrix.setRow(rowIndex, vec);
		}
	};
	friend class RowProxy;
	PyMatrix(size_t rows, size_t cols) : matrixVariant(std::make_shared<MatrixImpl<int>>(rows, cols)) {}

	RowProxy operator[](size_t rowIndex) {
        return RowProxy(*this, rowIndex);
    }
	const MatrixVariant& getMatrixVariant() const {
        return matrixVariant;
    }


	MatrixVariant& getMatrixVariant() {
		return matrixVariant;
	} 
	
	void set_value(size_t row, size_t col, const py::object& value) {
    	
		if (py::isinstance<py::int_>(value)) {
    	} else if (py::isinstance<py::float_>(value)) {
        	promoteMatrixVariantIfNeeded<double>(); 
   		} else if (py::isinstance(value, py::module_::import("numbers").attr("Complex"))) {
        	promoteMatrixVariantIfNeeded<std::complex<double>>(); 
    	} else {
        	throw std::runtime_error("Unsupported value type for matrix assignment.");
    	}

	    std::visit([row, col, &value](auto& arg) {
            using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
            using ValueType = typename MatrixType::ValueType;

            if constexpr (std::is_same_v<ValueType, int>) {
                arg->setValueAt(row, col, value.cast<int>());
            } else if constexpr (std::is_same_v<ValueType, double>) {
                arg->setValueAt(row, col, value.cast<double>());
            } else if constexpr (std::is_same_v<ValueType, std::complex<double>>) {
                arg->setValueAt(row, col, value.cast<std::complex<double>>());
            } else {
                throw std::runtime_error("Unsupported type for assignment.");
            }
        }, matrixVariant);
    }
	py::object __getitem__(size_t rowIndex) {
    	return py::cast(RowProxy(*this, rowIndex));
	}

    void __setitem__(const py::tuple& index, const py::object& value) {
        if (py::len(index) != 2) {
            throw std::runtime_error("Index must be a tuple of row and column.");
        }
        size_t row = index[0].cast<size_t>();
        size_t col = index[1].cast<size_t>();
        this->set_value(row, col, value);
    }


	void setColumn(size_t colIndex, const PyVector& pyVec) {
    	auto vecType = pyVec.getType();
    	if (vecType == "int") {
        	promoteMatrixVariantIfNeeded<int>();
    	} else if (vecType == "double") {
        	promoteMatrixVariantIfNeeded<double>();
    	} else if (vecType == "std::complex<double>") {
        	promoteMatrixVariantIfNeeded<std::complex<double>>();
    	} else {
        	throw std::runtime_error("Unsupported PyVector type.");
    	}

    	std::visit([colIndex, &pyVec, vecType](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        
        	std::vector<ValueType> columnData;

        	for (size_t i = 0; i < pyVec.size(); ++i) {
            	if (vecType == "int") {
                	auto item = pyVec.__getitem__(i).cast<int>();
                	columnData.push_back(convertValue<ValueType>(item));
            	} else if (vecType == "double") {
                	auto item = pyVec.__getitem__(i).cast<double>();
                	columnData.push_back(convertValue<ValueType>(item));
            	} else if (vecType == "std::complex<double>") {
                	auto item = pyVec.__getitem__(i).cast<std::complex<double>>();
                	columnData.push_back(convertValue<ValueType>(item));
            	}
        	}
        
        	arg->setColumnFromStdVector(colIndex, columnData);
    	}, matrixVariant);
	}

	void setRow(size_t rowIndex, const PyVector& pyVec) {
		auto vecType = pyVec.getType();
		if (vecType == "int") {
			promoteMatrixVariantIfNeeded<int>();
		} else if (vecType == "double") {
			promoteMatrixVariantIfNeeded<double>();
		} else if (vecType == "std::complex<double>") {
			promoteMatrixVariantIfNeeded<std::complex<double>>();
		} else {
			throw std::runtime_error("Unsupported PyVector type.");
		}

		std::visit([rowIndex, &pyVec, vecType](auto& arg) {
			using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
			using ValueType = typename MatrixType::ValueType;
				
			std::vector<ValueType> rowData;
			
			for (size_t i = 0; i < pyVec.size(); ++i) {
				if (vecType == "int") {
					auto item = pyVec.__getitem__(i).cast<int>();
					rowData.push_back(convertValue<ValueType>(item));
				} else if (vecType == "double") {
					auto item = pyVec.__getitem__(i).cast<double>();
					rowData.push_back(convertValue<ValueType>(item));
				} else if (vecType == "std::complex<double>") {
					auto item = pyVec.__getitem__(i).cast<std::complex<double>>();
					rowData.push_back(convertValue<ValueType>(item));
				}
			}
			arg->setRowFromStdVector(rowIndex, rowData);
		}, matrixVariant);
	}

		
	PyVector getColumnAsPyVector(size_t columnIndex) const {
    	return std::visit([columnIndex](const auto& arg) -> PyVector {
        	using MatrixType = typename std::decay<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        	std::vector<ValueType> columnData;

        	if(columnIndex >= arg->getCols()) {
            	throw std::out_of_range("Column index out of bounds.");
        	}

        	for(size_t rowIndex = 0; rowIndex < arg->getRows(); ++rowIndex) {
            	columnData.push_back(arg->getValueAt(rowIndex, columnIndex));
        	}

        	return PyVector(std::make_unique<TypedVector<ValueType>>(std::move(columnData)));
    	}, matrixVariant);
	}

	PyVector getRowAsPyVector(size_t rowIndex) const {
    	return std::visit([rowIndex](const auto& arg) -> PyVector {
        	using MatrixType = typename std::decay<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        	std::vector<ValueType> rowData;

        	if(rowIndex >= arg->getRows()) {
            	throw std::out_of_range("Row index out of bounds.");
        	}	

        	for(size_t colIndex = 0; colIndex < arg->getCols(); ++colIndex) {
            	rowData.push_back(arg->getValueAt(rowIndex, colIndex));
        	}

        	return PyVector(std::make_unique<TypedVector<ValueType>>(std::move(rowData)));
    	}, matrixVariant);
	}

	py::object __getitem__(const py::object& key) {
        if (py::isinstance<py::int_>(key)) {
            size_t rowIndex = key.cast<size_t>();
            return py::cast(getRowAsPyVector(rowIndex));
        } else if (py::isinstance<py::tuple>(key)) {
            py::tuple indexTuple = key.cast<py::tuple>();
            size_t row = indexTuple[0].cast<size_t>();
            size_t col = indexTuple[1].cast<size_t>();
            return getValueAt(row, col);
        } else {
            throw std::runtime_error("Invalid index type.");
        }
    }

	PyMatrix& operator=(const PyMatrix& other) {
    	if (this != &other) {
        	this->matrixVariant = deepCopyVariant(other.matrixVariant);
    	}
    	return *this;
	}

	PyMatrix& operator+=(const PyMatrix& other) {
		auto commonTypeIndex = std::max(
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, matrixVariant),
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, other.matrixVariant)
		);


    	switch (commonTypeIndex) {
        	case 0: promoteMatrixVariantIfNeeded<int>(); break;
        	case 1: promoteMatrixVariantIfNeeded<double>(); break;
        	case 2: promoteMatrixVariantIfNeeded<std::complex<double>>(); break;
        	default: throw std::runtime_error("Unsupported matrix type for operation.");
    	}

    	std::visit([this](auto& arg) {
        	using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        	promoteMatrixVariantIfNeeded<MatrixType>();

        	std::visit([&arg](auto& thisArg) {
            	using ThisMatrixType = typename std::decay_t<decltype(*thisArg)>::ValueType;

            	if (arg->getRows() != thisArg->getRows() || arg->getCols() != thisArg->getCols()) {
                	throw std::runtime_error("Matrix dimensions must match for addition.");
            	}

            	for (size_t row = 0; row < thisArg->getRows(); ++row) {
                	for (size_t col = 0; col < thisArg->getCols(); ++col) {
                    	ThisMatrixType thisValue = thisArg->getValueAt(row, col);
                    	ThisMatrixType otherValue = convertValue<ThisMatrixType, MatrixType>(arg->getValueAt(row, col));
                    	thisArg->setValueAt(row, col, thisValue + otherValue);
                	}
            	}
        	}, matrixVariant);
    	}, other.matrixVariant);

    	return *this;
	}

	PyMatrix& operator-=(const PyMatrix& other) {
		auto commonTypeIndex = std::max(
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, matrixVariant),
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, other.matrixVariant)
		);


    	switch (commonTypeIndex) {
        	case 0: promoteMatrixVariantIfNeeded<int>(); break;
        	case 1: promoteMatrixVariantIfNeeded<double>(); break;
        	case 2: promoteMatrixVariantIfNeeded<std::complex<double>>(); break;
        	default: throw std::runtime_error("Unsupported matrix type for operation.");
    	}

    	std::visit([this](auto& arg) {
        	using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        	promoteMatrixVariantIfNeeded<MatrixType>();

        	std::visit([&arg](auto& thisArg) {
            	using ThisMatrixType = typename std::decay_t<decltype(*thisArg)>::ValueType;

            	if (arg->getRows() != thisArg->getRows() || arg->getCols() != thisArg->getCols()) {
                	throw std::runtime_error("Matrix dimensions must match for addition.");
            	}

            	for (size_t row = 0; row < thisArg->getRows(); ++row) {
                	for (size_t col = 0; col < thisArg->getCols(); ++col) {
                    	ThisMatrixType thisValue = thisArg->getValueAt(row, col);
                    	ThisMatrixType otherValue = convertValue<ThisMatrixType, MatrixType>(arg->getValueAt(row, col));
                    	thisArg->setValueAt(row, col, thisValue - otherValue);
                	}
            	}
        	}, matrixVariant);
    	}, other.matrixVariant);

    	return *this;
	}
	
	PyMatrix& operator/=(const PyMatrix& other) {
		auto commonTypeIndex = std::max(
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, matrixVariant),
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, other.matrixVariant)
		);


    	switch (commonTypeIndex) {
        	case 0: promoteMatrixVariantIfNeeded<int>(); break;
        	case 1: promoteMatrixVariantIfNeeded<double>(); break;
        	case 2: promoteMatrixVariantIfNeeded<std::complex<double>>(); break;
        	default: throw std::runtime_error("Unsupported matrix type for operation.");
    	}

    	std::visit([this](auto& arg) {
        	using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        	promoteMatrixVariantIfNeeded<MatrixType>();

        	std::visit([&arg](auto& thisArg) {
            	using ThisMatrixType = typename std::decay_t<decltype(*thisArg)>::ValueType;

            	if (arg->getRows() != thisArg->getRows() || arg->getCols() != thisArg->getCols()) {
                	throw std::runtime_error("Matrix dimensions must match for addition.");
            	}

            	for (size_t row = 0; row < thisArg->getRows(); ++row) {
                	for (size_t col = 0; col < thisArg->getCols(); ++col) {
                    	ThisMatrixType thisValue = thisArg->getValueAt(row, col);
                    	ThisMatrixType otherValue = convertValue<ThisMatrixType, MatrixType>(arg->getValueAt(row, col));
                    	thisArg->setValueAt(row, col, thisValue / otherValue);
                	}
            	}
        	}, matrixVariant);
    	}, other.matrixVariant);

    	return *this;
	}
	
	PyMatrix& operator*=(const PyMatrix& other) {
		auto commonTypeIndex = std::max(
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, matrixVariant),
    		std::visit([](const auto& arg) -> size_t {
        		using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        		return getTypeIndex<MatrixType>();
    		}, other.matrixVariant)
		);


    	switch (commonTypeIndex) {
        	case 0: promoteMatrixVariantIfNeeded<int>(); break;
        	case 1: promoteMatrixVariantIfNeeded<double>(); break;
        	case 2: promoteMatrixVariantIfNeeded<std::complex<double>>(); break;
        	default: throw std::runtime_error("Unsupported matrix type for operation.");
    	}

    	std::visit([this](auto& arg) {
        	using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
        	promoteMatrixVariantIfNeeded<MatrixType>();

        	std::visit([&arg](auto& thisArg) {
            	using ThisMatrixType = typename std::decay_t<decltype(*thisArg)>::ValueType;

            	if (arg->getRows() != thisArg->getRows() || arg->getCols() != thisArg->getCols()) {
                	throw std::runtime_error("Matrix dimensions must match for addition.");
            	}

            	for (size_t row = 0; row < thisArg->getRows(); ++row) {
                	for (size_t col = 0; col < thisArg->getCols(); ++col) {
                    	ThisMatrixType thisValue = thisArg->getValueAt(row, col);
                    	ThisMatrixType otherValue = convertValue<ThisMatrixType, MatrixType>(arg->getValueAt(row, col));
                    	thisArg->setValueAt(row, col, thisValue * otherValue);
                	}
            	}
        	}, matrixVariant);
    	}, other.matrixVariant);

    	return *this;
	}

	template<typename Scalar>
	PyMatrix& operator+=(const Scalar& scalar) {
    	auto addScalar = [&](auto& matrix) {
        	using MatrixValueType = typename std::decay_t<decltype(*matrix)>::ValueType;
        	if constexpr (!std::is_same_v<Scalar, MatrixValueType>) {
            	auto promotedScalar = convertValue<MatrixValueType>(scalar);
            	for (size_t row = 0; row < matrix->getRows(); ++row) {
                	for (size_t col = 0; col < matrix->getCols(); ++col) {
                    	matrix->getValueAt(row, col) += promotedScalar;
                	}
            	}
        	} else {
            	for (size_t row = 0; row < matrix->getRows(); ++row) {
                	for (size_t col = 0; col < matrix->getCols(); ++col) {
                    	matrix->getValueAt(row, col) += scalar;
                	}
            	}
        	}
    	};

    	std::visit(addScalar, matrixVariant);

    	return *this;
	}

	template<typename Scalar>
	PyMatrix& operator-=(const Scalar& scalar) {
		auto subScalar = [&](auto& matrix) {
			using MatrixValueType = typename std::decay_t<decltype(*matrix)>::ValueType;
			if constexpr (!std::is_same_v<Scalar, MatrixValueType>) {
				auto promotedScalar = convertValue<MatrixValueType>(scalar);
				for (size_t row = 0; row < matrix->getRows(); ++row) {
					for (size_t col = 0; col < matrix->getCols(); ++col) {
						matrix->getValueAt(row, col) -= promotedScalar;
					}
				}
			} else {
				for (size_t row = 0; row < matrix->getRows(); ++row) {
					for (size_t col = 0; col < matrix->getCols(); ++col) {
						matrix->getValueAt(row, col) -= scalar;
					}
				}
			}
		};
		std::visit(subScalar, matrixVariant);
		return *this;
	}

	template<typename Scalar>
	PyMatrix& operator*=(const Scalar& scalar) {
		auto mulScalar = [&](auto& matrix) {
			using MatrixValueType = typename std::decay_t<decltype(*matrix)>::ValueType;
			if constexpr (!std::is_same_v<Scalar, MatrixValueType>) {
				auto promotedScalar = convertValue<MatrixValueType>(scalar);
				for (size_t row = 0; row < matrix->getRows(); ++row) {
					for (size_t col = 0; col < matrix->getCols(); ++col) {
						matrix->getValueAt(row, col) *= promotedScalar;
					}
				}
			} else {
				for (size_t row = 0; row < matrix->getRows(); ++row) {
					for (size_t col = 0; col < matrix->getCols(); ++col) {
						matrix->getValueAt(row, col) *= scalar;
					}
				}
			}
		};
		std::visit(mulScalar, matrixVariant);
		return *this;
	}

	template<typename Scalar>
	PyMatrix& operator/=(const Scalar& scalar) {
		auto divScalar = [&](auto& matrix) {
			using MatrixValueType = typename std::decay_t<decltype(*matrix)>::ValueType;
			if constexpr (!std::is_same_v<Scalar, MatrixValueType>) {
				auto promotedScalar = convertValue<MatrixValueType>(scalar);
				for (size_t row = 0; row < matrix->getRows(); ++row) {
					for (size_t col = 0; col < matrix->getCols(); ++col) {
						matrix->getValueAt(row, col) /= promotedScalar;
					}
				}
			} else {
				for (size_t row = 0; row < matrix->getRows(); ++row) {
					for (size_t col = 0; col < matrix->getCols(); ++col) {
						matrix->getValueAt(row, col) /= scalar;
					}
				}
			}
		};
		std::visit(divScalar, matrixVariant);
		return *this;
	}
	
	PyMatrix operator+(const PyMatrix& other) const {
    	auto commonTypeIndex = std::max(
        	std::visit([](const auto& arg) -> size_t {
            	using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
            	return getTypeIndex<MatrixType>();
        	}, this->matrixVariant),
        	std::visit([](const auto& arg) -> size_t {
            	using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;
            	return getTypeIndex<MatrixType>();
        	}, other.matrixVariant)
    	);

    	PyMatrix result = *this;
    	switch (commonTypeIndex) {
        	case 0: result.promoteMatrixVariantIfNeeded<int>(); break;
        	case 1: result.promoteMatrixVariantIfNeeded<double>(); break;
        	case 2: result.promoteMatrixVariantIfNeeded<std::complex<double>>(); break;
        	default: throw std::runtime_error("Unsupported matrix type for operation.");
    	}

		std::visit([&result, commonTypeIndex](auto& arg) {
        	using MatrixType = typename std::decay_t<decltype(*arg)>::ValueType;

        	switch (commonTypeIndex) {
            	case 0: result.promoteMatrixVariant<int>(); break;
            	case 1: result.promoteMatrixVariant<double>(); break;
            	case 2: result.promoteMatrixVariant<std::complex<double>>(); break;
        	}

            std::visit([&arg](auto& thisArg) {
				using ThisMatrixType = typename std::decay_t<decltype(*thisArg)>::ValueType;

            	if (arg->getRows() != thisArg->getRows() || arg->getCols() != thisArg->getCols()) {
                	throw std::runtime_error("Matrix dimensions must match for addition.");
            	}

            	for (size_t row = 0; row < thisArg->getRows(); ++row) {
                	for (size_t col = 0; col < thisArg->getCols(); ++col) {
                    	auto thisValue = thisArg->getValueAt(row, col);
                    	auto otherValueTemp = arg->getValueAt(row, col);
                    	ThisMatrixType otherValue = convertValue<ThisMatrixType, MatrixType>(otherValueTemp);
                    	ThisMatrixType sumValue = thisValue + otherValue;
                    	thisArg->setValueAt(row, col, sumValue);
                	}
            	}
        	}, result.matrixVariant);
    	}, other.matrixVariant);

    	return result;
	}
	
	template<typename Scalar>
    PyMatrix operator+(const Scalar& scalar) const {
        PyMatrix result(*this);
        result.promoteMatrixVariantIfNeeded<Scalar>();

        auto addScalar = [&scalar](auto& matrixPtr) {
            using MatrixValueType = typename std::decay_t<decltype(*matrixPtr)>::ValueType;
            
            MatrixValueType convertedScalar = convertValue<MatrixValueType, Scalar>(scalar);

            for (size_t row = 0; row < matrixPtr->getRows(); ++row) {
                for (size_t col = 0; col < matrixPtr->getCols(); ++col) {
                    matrixPtr->getValueAt(row, col) += convertedScalar;
                }
            }
        };

        std::visit(addScalar, result.matrixVariant);
        return result;
    }	

	bool operator==(const PyMatrix& other) const {
    return std::visit([&](const auto& arg) -> bool {
        const auto& otherVariant = other.matrixVariant;

        const auto otherArgPtr = std::get_if<std::shared_ptr<std::decay_t<decltype(*arg)>>>(&otherVariant);
        if (!otherArgPtr) {
            return false; 
        }

        const auto& otherArg = *otherArgPtr;

        if (arg->getRows() != otherArg->getRows() || arg->getCols() != otherArg->getCols()) {
            return false;
        }

        for (size_t row = 0; row < arg->getRows(); ++row) {
            for (size_t col = 0; col < arg->getCols(); ++col) {
                if (!(arg->getValueAt(row, col) == otherArg->getValueAt(row, col))) {
                    return false; 
                }
            }
        }
        return true; 
    }, matrixVariant);
}
	template<typename T, typename U>
	bool compareValues(const T& a, const U& b) {
    	if constexpr (std::is_same_v<T, U>) {
        	return a == b;
    	} else {
        	static_assert(std::is_same_v<T, U>, "Trying to compare values of different types.");
        	return false;
    	}
	}

	template<typename Scalar>
	void fillColumn(size_t colIndex, const Scalar& value) {
    	promoteMatrixVariantIfNeeded<Scalar>();

    	std::visit([colIndex, &value](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;

        	ValueType convertedValue = convertValue<ValueType, Scalar>(value);

        	if (colIndex >= arg->getCols()) {
            	throw std::out_of_range("Column index out of bounds.");
        	}
        	for(size_t rowIndex = 0; rowIndex < arg->getRows(); ++rowIndex) {
            	arg->setValueAt(rowIndex, colIndex, convertedValue);
        	}
    	}, matrixVariant);
	}
	
	template<typename Scalar>
	void fillRow(size_t rowIndex, const Scalar& value) {
		promoteMatrixVariantIfNeeded<Scalar>();
		std::visit([rowIndex, &value](auto& arg) {
			using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
			using ValueType = typename MatrixType::ValueType;
			ValueType convertedValue = convertValue<ValueType, Scalar>(value);

			if (rowIndex >= arg->getRows()) {
				throw std::out_of_range("Row index out of bounds.");
			}
			for (size_t colIndex = 0; colIndex < arg->getCols(); ++colIndex) {
				arg->setValueAt(rowIndex, colIndex, convertedValue);
			}
		}, matrixVariant);
	}	

	template<typename Scalar>
	void multiplyColumn(size_t colIndex, const Scalar& scalar) {
    	promoteMatrixVariantIfNeeded<Scalar>();
    
    	std::visit([colIndex, &scalar](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        
        	ValueType convertedScalar = convertValue<ValueType, Scalar>(scalar);

        	if (colIndex >= arg->getCols()) {
            	throw std::out_of_range("Column index out of bounds.");
        	}
        	for (size_t rowIndex = 0; rowIndex < arg->getRows(); ++rowIndex) {
            	arg->setValueAt(rowIndex, colIndex, arg->getValueAt(rowIndex, colIndex) * convertedScalar);
        	}
    	}, matrixVariant);
	}

	template<typename Scalar>
	void multiplyRow(size_t rowIndex, const Scalar& scalar) {
    	promoteMatrixVariantIfNeeded<Scalar>();
    
    	std::visit([rowIndex, &scalar](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        
        	ValueType convertedScalar = convertValue<ValueType, Scalar>(scalar);

        	if (rowIndex >= arg->getRows()) {
            	throw std::out_of_range("Row index out of bounds.");
        	}
        	for (size_t colIndex = 0; colIndex < arg->getCols(); ++colIndex) {
            	arg->setValueAt(rowIndex, colIndex, arg->getValueAt(rowIndex, colIndex) * convertedScalar);
        	}
    	}, matrixVariant);
	}

	void setMatrixColumn(size_t colIndex, const PyVector& pyVec) {
    	auto vecType = pyVec.getType(); 
    	if (vecType == "int") {
        	promoteMatrixVariantIfNeeded<int>();
    	} else if (vecType == "double") {
        	promoteMatrixVariantIfNeeded<double>();
    	} else if (vecType == "std::complex<double>") {
        	promoteMatrixVariantIfNeeded<std::complex<double>>();
    	} else {
        	throw std::runtime_error("Unsupported PyVector type.");
    	}

    	std::visit([colIndex, &pyVec, vecType](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        
        	std::vector<ValueType> columnData;

        	for (size_t i = 0; i < pyVec.size(); ++i) {
            	if (vecType == "int") {
                	auto item = pyVec.__getitem__(i).cast<int>();
                	columnData.push_back(convertValue<ValueType>(item));
            	} else if (vecType == "double") {
                	auto item = pyVec.__getitem__(i).cast<double>();
                	columnData.push_back(convertValue<ValueType>(item));
            	} else if (vecType == "std::complex<double>") {
                	auto item = pyVec.__getitem__(i).cast<std::complex<double>>();
                	columnData.push_back(convertValue<ValueType>(item));
            	}
        	}
        
        	arg->setColumnFromStdVector(colIndex, columnData);
    	}, matrixVariant);
	}

	void setMatrixRow(size_t rowIndex, const PyVector& pyVec) {
		auto vecType = pyVec.getType();
		if (vecType == "int") {
			promoteMatrixVariantIfNeeded<int>();
		} else if (vecType == "double") {
			promoteMatrixVariantIfNeeded<double>();
		} else if (vecType == "std::complex<double>") {
			promoteMatrixVariantIfNeeded<std::complex<double>>();
		} else {
			throw std::runtime_error("Unsupported PyVector type");
		}

		std::visit([rowIndex, &pyVec, vecType](auto& arg) {
			using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
			using ValueType = typename MatrixType::ValueType;
	
			std::vector<ValueType> rowData;
		
			for (size_t i = 0; i < pyVec.size(); ++i) {
				if (vecType == "int") {
					auto item = pyVec.__getitem__(i).cast<int>();
					rowData.push_back(convertValue<ValueType>(item));
				} else if (vecType == "double") {
					auto item = pyVec.__getitem__(i).cast<double>();
					rowData.push_back(convertValue<ValueType>(item));
				} else if (vecType == "std::complex<double>") {
					auto item = pyVec.__getitem__(i).cast<std::complex<double>>();
					rowData.push_back(convertValue<ValueType>(item));
				}
			}
			arg->setRowFromStdVector(rowIndex, rowData);
		}, matrixVariant);
	}

	void add_to_column(size_t colIndex, const PyVector& pyVec) {
    	auto vecType = pyVec.getType();
    	if (vecType == "int") {
        	promoteMatrixVariantIfNeeded<int>();
    	} else if (vecType == "double") {
        	promoteMatrixVariantIfNeeded<double>();
    	} else if (vecType == "std::complex<double>") {
       		promoteMatrixVariantIfNeeded<std::complex<double>>();
    	} else {
        	throw std::runtime_error("Unsupported PyVector type.");
    	}

    	std::visit([colIndex, &pyVec](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        
        	if (colIndex >= arg->getCols()) {
            	throw std::out_of_range("Column index out of bounds.");
        	}
        	if (pyVec.size() != arg->getRows()) {
            	throw std::invalid_argument("PyVector size does not match the number of matrix rows.");
        	}

        	for (size_t row = 0; row < arg->getRows(); ++row) {
            	ValueType vecValue = pyVec.__getitem__(row).cast<ValueType>();
            	arg->setValueAt(row, colIndex, arg->getValueAt(row, colIndex) + vecValue);
        	}
    	}, matrixVariant);
	}

	void add_to_row(size_t rowIndex, const PyVector& pyVec) {
    	auto vecType = pyVec.getType();
    	if (vecType == "int") {
        	promoteMatrixVariantIfNeeded<int>();
    	} else if (vecType == "double") {
        	promoteMatrixVariantIfNeeded<double>();
    	} else if (vecType == "std::complex<double>") {
        	promoteMatrixVariantIfNeeded<std::complex<double>>();
    	} else {
        	throw std::runtime_error("Unsupported PyVector type.");
    	}

    	std::visit([rowIndex, &pyVec](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        
        	if (rowIndex >= arg->getRows()) {
            	throw std::out_of_range("Row index out of bounds.");
        	}
        	if (pyVec.size() != arg->getCols()) {
            	throw std::invalid_argument("PyVector size does not match the number of matrix columns.");
        	}

        	for (size_t col = 0; col < arg->getCols(); ++col) {
            	ValueType vecValue = pyVec.__getitem__(col).cast<ValueType>();
            	arg->setValueAt(rowIndex, col, arg->getValueAt(rowIndex, col) + vecValue);
        	}
    	}, matrixVariant);
	}

	void subtract_from_column(size_t colIndex, const PyVector& pyVec) {
    	auto vecType = pyVec.getType();
    	if (vecType == "int") {
        	promoteMatrixVariantIfNeeded<int>();
    	} else if (vecType == "double") {
        	promoteMatrixVariantIfNeeded<double>();
    	} else if (vecType == "std::complex<double>") {
        	promoteMatrixVariantIfNeeded<std::complex<double>>();
    	} else {
        	throw std::runtime_error("Unsupported PyVector type.");
    	}

    	std::visit([colIndex, &pyVec](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;

        	if (colIndex >= arg->getCols()) {
            	throw std::out_of_range("Column index out of bounds.");
        	}
        	if (pyVec.size() != arg->getRows()) {
            	throw std::invalid_argument("PyVector size does not match the number of matrix rows.");
        	}
        	for (size_t rowIndex = 0; rowIndex < arg->getRows(); ++rowIndex) {
            	auto item = pyVec.__getitem__(rowIndex);
            	ValueType vectorValue = item.cast<ValueType>();
            	arg->setValueAt(rowIndex, colIndex, arg->getValueAt(rowIndex, colIndex) - vectorValue);
        	}
    	}, matrixVariant);
	}
	
	void subtract_from_row(size_t rowIndex, const PyVector& pyVec) {
    	auto vecType = pyVec.getType();
    	if (vecType == "int") {
        	promoteMatrixVariantIfNeeded<int>();
    	} else if (vecType == "double") {
        	promoteMatrixVariantIfNeeded<double>();
    	} else if (vecType == "std::complex<double>") {
        	promoteMatrixVariantIfNeeded<std::complex<double>>();
    	} else {
        	throw std::runtime_error("Unsupported PyVector type.");
    	}

    	std::visit([rowIndex, &pyVec](auto& arg) {
        	using MatrixType = typename std::remove_reference<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;

        	if (rowIndex >= arg->getRows()) {
            	throw std::out_of_range("Row index out of bounds.");
        	}
        	if (pyVec.size() != arg->getCols()) {
            	throw std::invalid_argument("PyVector size does not match the number of matrix columns.");
        	}
        	for (size_t colIndex = 0; colIndex < arg->getCols(); ++colIndex) {
            	auto item = pyVec.__getitem__(colIndex);
            	ValueType vectorValue = item.cast<ValueType>();
            	arg->setValueAt(rowIndex, colIndex, arg->getValueAt(rowIndex, colIndex) - vectorValue);
        	}
    	}, matrixVariant);
	}

	void setToUnity() {
    	std::visit([](auto& arg) {
        	using MatrixType = typename std::decay_t<decltype(*arg)>;
        	using ValueType = typename MatrixType::ValueType;

        	size_t rows = arg->getRows();
        	size_t cols = arg->getCols();
        	if (rows != cols) {
            	throw std::runtime_error("setToIdentity can only be applied to square matrices.");
        	}

        	for (size_t row = 0; row < rows; ++row) {
            	for (size_t col = 0; col < cols; ++col) {
                	if (row == col) {
                    	arg->setValueAt(row, col, convertValue<ValueType>(1));
                	} else {
                    	arg->setValueAt(row, col, convertValue<ValueType>(0));
                	}
            	}
        	}
    	}, matrixVariant);
	}

	static PyMatrix unity(size_t N, const std::string& type) {
    	PyMatrix matrix(N, N); 
    	if (type == "float") {
        	matrix.promoteMatrixVariantIfNeeded<double>();
    	} else if (type == "complex") {
        	matrix.promoteMatrixVariantIfNeeded<std::complex<double>>();
    	}
    	std::visit([&](auto& arg) {
        	using MatrixType = typename std::decay<decltype(*arg)>::type;
        	using ValueType = typename MatrixType::ValueType;
        	ValueType one = static_cast<ValueType>(1);
        	for (size_t i = 0; i < N; ++i) {
            	arg->setValueAt(i, i, one);
        	}
    	}, matrix.getMatrixVariant());
    	return matrix;
	}

	py::object trace() const {
        return std::visit([](const auto& matrixPtr) -> py::object {
            using MatrixType = typename std::decay_t<decltype(*matrixPtr)>;
            using ValueType = typename MatrixType::ValueType;

            if (matrixPtr->getRows() != matrixPtr->getCols()) {
                throw std::runtime_error("Trace can only be computed for square matrices.");
            }

            ValueType traceValue = ValueType(0);
            for (size_t i = 0; i < matrixPtr->getRows(); ++i) {
                traceValue += matrixPtr->getValueAt(i, i);
            }
            
            return py::cast(traceValue);
        }, matrixVariant);
    }

	py::object __getitem__(py::object key) {
        if (py::isinstance<py::int_>(key)) {
            size_t rowIndex = key.cast<size_t>();
            return py::cast(RowProxy(*this, rowIndex));
        } else if (py::isinstance<py::tuple>(key)) {
            py::tuple indexTuple = key.cast<py::tuple>();
            size_t row = indexTuple[0].cast<size_t>();
            size_t col = indexTuple[1].cast<size_t>();
        } else {
            throw std::runtime_error("Indexing with unsupported type.");
        }
    }

	PyMatrix operator*(const PyMatrix& other) const {
    	auto matrixProduct = [](const auto& lhsMatrix, const auto& rhsMatrix) -> MatrixVariant {
        	using LhsValueType = typename std::decay_t<decltype(*lhsMatrix)>::ValueType;
        	using RhsValueType = typename std::decay_t<decltype(*rhsMatrix)>::ValueType;
        	using ResultValueType = typename std::common_type<LhsValueType, RhsValueType>::type;

        	if (lhsMatrix->getCols() != rhsMatrix->getRows()) {
            	throw std::runtime_error("Matrix dimensions mismatch for multiplication.");
        	}

        	auto resultMatrix = std::make_shared<MatrixImpl<ResultValueType>>(lhsMatrix->getRows(), rhsMatrix->getCols());

        	for (size_t i = 0; i < lhsMatrix->getRows(); ++i) {
            	for (size_t j = 0; j < rhsMatrix->getCols(); ++j) {
                	ResultValueType sum = ResultValueType();
                	for (size_t k = 0; k < lhsMatrix->getCols(); ++k) {
                	    sum += static_cast<ResultValueType>(lhsMatrix->getValueAt(i, k)) * static_cast<ResultValueType>(rhsMatrix->getValueAt(k, j));
            	    }
        	        resultMatrix->setValueAt(i, j, sum);
    	        }
	        }
	
        	return resultMatrix;
    	};

    	MatrixVariant resultVariant = std::visit([&](const auto& lhsVariant) -> MatrixVariant {
        	return std::visit([&](const auto& rhsVariant) -> MatrixVariant {
            	return matrixProduct(lhsVariant, rhsVariant);
        	}, other.matrixVariant);
    	}, this->matrixVariant);

    	return PyMatrix{resultVariant};
	}
	
	static void apl_ml(PyMatrix& A, const py::object& scalar, PyMatrix& C) {
		std::variant<int, double, std::complex<double>> typedScalar;

		if (py::isinstance<py::int_>(scalar)) {
			typedScalar = scalar.cast<int>();
			C.promoteMatrixVariantIfNeeded<int>();
		} else if (py::isinstance<py::float_>(scalar)) {
			typedScalar = scalar.cast<double>();
			C.promoteMatrixVariantIfNeeded<double>();
		} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {
			typedScalar = scalar.cast<std::complex<double>>();
			C.promoteMatrixVariantIfNeeded<std::complex<double>>();
		} else {
			throw std::runtime_error("Unsupported scalar type for operation.");
		}
	
		std::visit([&](auto& argA, auto& argC) {
			using ValueTypeA = typename std::remove_reference<decltype(argA)>::type::element_type::ValueType;
			using ValueTypeC = typename std::remove_reference<decltype(argC)>::type::element_type::ValueType;
				if (argA->getRows() != argC->getRows() || argA->getCols() != argC->getCols()) {
					throw std::runtime_error("Matrix dimensions must match for operation.");
				}

				for (size_t row = 0; row < argA->getRows(); ++row) {
					for (size_t col = 0; col < argA->getCols(); ++col) {
						std::visit([&](auto&& scalarValue) {
							if constexpr(std::is_same_v<decltype(scalarValue), ValueTypeA> && std::is_same_v<ValueTypeA, ValueTypeC>) {
								argC->setValueAt(row, col, argC->getValueAt(row, col) + argA->getValueAt(row, col) * scalarValue);
							} else {
								auto convertedScalar = convertValue<ValueTypeC>(scalarValue);
								auto convertedAValue = convertValue<ValueTypeC>(argA->getValueAt(row, col));
								argC->setValueAt(row, col, argC->getValueAt(row, col) + convertedAValue * convertedScalar);
							}
						}, typedScalar);
					}
				}
		}, A.getMatrixVariant(), C.getMatrixVariant());
	}
	static void multiplyZ(PyMatrix& A, const py::object& scalar, PyMatrix& C) {
		std::variant<int, double, std::complex<double>> typedScalar;
		if (py::isinstance<py::int_>(scalar)) {
			typedScalar = scalar.cast<int>();
			C.promoteMatrixVariantIfNeeded<int>();
		} else if (py::isinstance<py::float_>(scalar)) {
			typedScalar = scalar.cast<double>();
			C.promoteMatrixVariantIfNeeded<double>();
		} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {
			typedScalar = scalar.cast<std::complex<double>>();
			C.promoteMatrixVariantIfNeeded<std::complex<double>>();
		} else {
			throw std::runtime_error("Unsupported scalar type for operation.");
		}

		std::visit([&](auto& argA, auto& argC) {
			using ValueTypeA = typename std::remove_reference<decltype(argA)>::type::element_type::ValueType;	
			using ValueTypeC = typename std::remove_reference<decltype(argC)>::type::element_type::ValueType;
			
				if (argA->getRows() != argC->getRows() || argA->getCols() != argC->getCols()) {
					throw std::runtime_error("Matrix dimensions must match for operation.");;
				}

				for (size_t row = 0; row < argA->getRows(); ++row) {
					for (size_t col = 0; col < argA->getCols(); ++col) {
						std::visit([&](auto&& scalarValue) {
							if constexpr(std::is_same_v<decltype(scalarValue), ValueTypeA> && std::is_same_v<ValueTypeA, ValueTypeC>) {
								argC->setValueAt(row, col, argA->getValueAt(row, col) * scalarValue);
							} else {
								auto convertedScalar = convertValue<ValueTypeC>(scalarValue);
								auto convertedAValue = convertValue<ValueTypeC>(argA->getValueAt(row, col));
								argC->setValueAt(row, col, convertedAValue * convertedScalar);
							}
						}, typedScalar);
					}
				}
		}, A.getMatrixVariant(), C.getMatrixVariant());
	}
	static void as_ml_ml(const PyMatrix& A, const PyMatrix& B, const py::object& scalar, PyMatrix& C) {
		std::variant<int, double, std::complex<double>> typedScalar;
		if (py::isinstance<py::int_>(scalar)) {
			typedScalar = scalar.cast<int>();
			C.promoteMatrixVariantIfNeeded<int>();
		} else if (py::isinstance<py::float_>(scalar)) {
			typedScalar = scalar.cast<double>();
			C.promoteMatrixVariantIfNeeded<double>();
		} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {
			typedScalar = scalar.cast<std::complex<double>>();
			C.promoteMatrixVariantIfNeeded<std::complex<double>>();
		} else {
			throw std::runtime_error("Unsupported scalar type for operation.");
		}

		std::visit([&](auto& argA, auto& argB, auto& argC) {
			using ValueTypeA = typename std::decay_t<decltype(argA)>::element_type::ValueType;
			using ValueTypeB = typename std::decay_t<decltype(argB)>::element_type::ValueType;
			using ValueTypeC = typename std::decay_t<decltype(argC)>::element_type::ValueType;
			
			if (argA->getRows() != argB->getRows() || argA->getRows() != argC->getRows() ||
				argA->getCols() != argB->getCols() || argA->getCols() != argC->getCols()) {
					throw std::runtime_error("Matrix dimensions must match for operation.");
			}
			
			for (size_t row = 0; row < argA->getRows(); ++row) {
				for (size_t col = 0; col < argA->getCols(); ++col) {
					std::visit([&](auto&& scalarValue) {
							auto convertedScalar = convertValue<ValueTypeC>(scalarValue);
							auto convertedAValue = convertValue<ValueTypeC>(argA->getValueAt(row, col));
							auto convertedBValue = convertValue<ValueTypeC>(argB->getValueAt(row, col));
							argC->setValueAt(row, col, convertedAValue * convertedBValue * convertedScalar);
					}, typedScalar);
				}
			}
		
		}, A.getMatrixVariant(), B.getMatrixVariant(), C.getMatrixVariant());
	}
	PyVector operator*(const PyVector& vec) const {
    	return std::visit([&vec](const auto& matrixPtr) -> PyVector {
        	using MatrixType = typename std::decay<decltype(*matrixPtr)>::type;
        	using ValueType = typename MatrixType::ValueType;

        	if (vec.size() != matrixPtr->getCols()) {
            	throw std::runtime_error("Dimension mismatch between matrix and vector.");
        	}

        	std::vector<ValueType> result(matrixPtr->getRows(), ValueType(0));

        	for (size_t row = 0; row < matrixPtr->getRows(); ++row) {
            	for (size_t col = 0; col < matrixPtr->getCols(); ++col) {
                	auto vecValue = vec.__getitem__(col).cast<ValueType>();
                	result[row] += matrixPtr->getValueAt(row, col) * vecValue;
            	}
        	}
		
    	    return PyVector(std::make_unique<TypedVector<ValueType>>(std::move(result)));
    	}, matrixVariant);
	}

	PyVector multiplyLeft(const PyVector& vec) const {
        PyMatrix transposedMatrix = this->transpose();
        return transposedMatrix * vec;
    }			

	friend PyVector operator*(const PyVector& vec, const PyMatrix& mat) {
        return mat.transpose() * vec;  
    }

	

	static void multiply(const PyMatrix& A, const PyMatrix& B, PyMatrix& C) {
        auto resultVariant = std::visit([&](const auto& argA, const auto& argB) -> MatrixVariant {
            using MatrixTypeA = typename std::decay_t<decltype(*argA)>::ValueType;
            using MatrixTypeB = typename std::decay_t<decltype(*argB)>::ValueType;
            using ResultValueType = typename std::common_type<MatrixTypeA, MatrixTypeB>::type;

            if (argA->getCols() != argB->getRows()) {
                throw std::runtime_error("Matrix dimensions mismatch for multiplication.");
            }

            auto resultMatrix = std::make_shared<MatrixImpl<ResultValueType>>(argA->getRows(), argB->getCols());

            for (size_t i = 0; i < argA->getRows(); ++i) {
                for (size_t j = 0; j < argB->getCols(); ++j) {
                    ResultValueType sum = ResultValueType();
                    for (size_t k = 0; k < argA->getCols(); ++k) {
                        sum += static_cast<ResultValueType>(argA->getValueAt(i, k)) * static_cast<ResultValueType>(argB->getValueAt(k, j));
                    }
                    resultMatrix->setValueAt(i, j, sum);
                }
            }

            return resultMatrix;
        }, A.getMatrixVariant(), B.getMatrixVariant());

        C = PyMatrix(resultVariant);
    }
	
	static void multiplyZ(const PyMatrix& A, const PyMatrix& B, PyMatrix& C) {
        C = std::visit([&](const auto& lhsMatrixVariant) -> PyMatrix {
            using LhsMatrixType = typename std::decay_t<decltype(*lhsMatrixVariant)>;
            using LhsValueType = typename LhsMatrixType::ValueType;

            return std::visit([&](const auto& rhsMatrixVariant) -> PyMatrix {
                using RhsMatrixType = typename std::decay_t<decltype(*rhsMatrixVariant)>;
                using RhsValueType = typename RhsMatrixType::ValueType;
                using ResultValueType = std::common_type_t<LhsValueType, RhsValueType>;

                if (lhsMatrixVariant->getCols() != rhsMatrixVariant->getRows()) {
                    throw std::runtime_error("Matrix dimensions mismatch for multiplication.");
                }

                auto resultMatrix = std::make_shared<MatrixImpl<ResultValueType>>(lhsMatrixVariant->getRows(), rhsMatrixVariant->getCols());

                for (size_t i = 0; i < lhsMatrixVariant->getRows(); ++i) {
                    for (size_t j = 0; j < rhsMatrixVariant->getCols(); ++j) {
                        ResultValueType sum = ResultValueType();
                        for (size_t k = 0; k < lhsMatrixVariant->getCols(); ++k) {
                            LhsValueType lhsValue = lhsMatrixVariant->getValueAt(i, k);
                            RhsValueType rhsValue = rhsMatrixVariant->getValueAt(k, j);
                            if (lhsValue != LhsValueType(0) && rhsValue != RhsValueType(0)) {
                                sum += static_cast<ResultValueType>(lhsValue) * static_cast<ResultValueType>(rhsValue);
                            }
                        }
                        resultMatrix->setValueAt(i, j, sum);
                    }
                }

                return PyMatrix{resultMatrix};
            }, B.matrixVariant);
        }, A.matrixVariant);
	}
		
	PyMatrix transpose() const {
        MatrixVariant transposedVariant = std::visit([](const auto& matrixImpl) -> MatrixVariant {
            return matrixImpl->transpose();
        }, matrixVariant);
        return PyMatrix(transposedVariant);
    }

	PyMatrix operator!() const {
        MatrixVariant transposedVariant = std::visit([](const auto& matrixImpl) -> MatrixVariant {
            return matrixImpl->transpose(); 
        }, matrixVariant);
        return PyMatrix(transposedVariant);
    }


	static void GaussInvert(PyMatrix& matrix) {
    	bool hasComplex = std::visit([](const auto& arg) -> bool {
        	using MatType = typename std::decay_t<decltype(*arg)>::ValueType;
        	return std::is_same<MatType, std::complex<double>>::value;
    	}, matrix.getMatrixVariant());

    	if (hasComplex) {
        	matrix.promoteMatrixVariantIfNeeded<std::complex<double>>();
    	} else {
        	matrix.promoteMatrixVariantIfNeeded<double>();
    	}

    	auto& variant = matrix.getMatrixVariant();
    	auto matrixPtr = std::get_if<std::shared_ptr<MatrixImpl<double>>>(&variant);
    	auto matrixPtrComplex = std::get_if<std::shared_ptr<MatrixImpl<std::complex<double>>>>(&variant);

    	size_t N;
    	std::vector<std::vector<std::complex<double>>> augmented;

    	if (matrixPtr) {
        	auto& A = *matrixPtr;
        	N = A->getRows();
        	augmented.resize(N, std::vector<std::complex<double>>(2 * N));
        	for (size_t i = 0; i < N; ++i) {
            	for (size_t j = 0; j < N; ++j) {
                	augmented[i][j] = A->getValueAt(i, j);
            	}
            	augmented[i][N + i] = 1.0;
        	}
    	} else if (matrixPtrComplex) {
        	auto& A = *matrixPtrComplex;
        	N = A->getRows();
        	augmented.resize(N, std::vector<std::complex<double>>(2 * N));
        	for (size_t i = 0; i < N; ++i) {
            	for (size_t j = 0; j < N; ++j) {
                	augmented[i][j] = A->getValueAt(i, j);
            	}
            	augmented[i][N + i] = 1.0;
        	}
    	} else {
        	throw std::runtime_error("Matrix type promotion failed or unsupported type.");
    	}

    	for (size_t k = 0; k < N; ++k) {
        	size_t maxRow = k;
        	double maxAbs = std::abs(augmented[k][k]);
        	for (size_t i = k + 1; i < N; ++i) {
            	double absVal = std::abs(augmented[i][k]);
            	if (absVal > maxAbs) {
                	maxAbs = absVal;
                	maxRow = i;
            	}
        	}

        	std::swap(augmented[k], augmented[maxRow]);

        	auto pivot = augmented[k][k];
        	for (size_t j = k; j < 2 * N; ++j) {
            	augmented[k][j] /= pivot;
        		if (std::abs(augmented[k][j]) < 1e-10) {
					augmented[k][j] = 0.0;
				}
			}

        	for (size_t i = 0; i < N; ++i) {
            	if (i != k) {
                	auto factor = augmented[i][k];
                	for (size_t j = k; j < 2 * N; ++j) {
                    	augmented[i][j] -= factor * augmented[k][j];
                		if (std::abs(augmented[i][j]) < 1e-10) {
							augmented[i][j] = 0.0;
						}
					}
            	}
        	}
    	}

    	if (matrixPtr) {
        	auto& A = *matrixPtr;
        	for (size_t i = 0; i < N; ++i) {
            	for (size_t j = 0; j < N; ++j) {
                	A->setValueAt(i, j, augmented[i][j + N].real());
            	}
        	}
    	} else if (matrixPtrComplex) {
        	auto& A = *matrixPtrComplex;
        	for (size_t i = 0; i < N; ++i) {
            	for (size_t j = 0; j < N; ++j) {
                	A->setValueAt(i, j, augmented[i][j + N]);
            	}
        	}
    	}
	}

	static void invert(const PyMatrix& A, PyMatrix& A1) {
        A1 = A;  
        GaussInvert(A1);  
    }

    static PyMatrix inverse(const PyMatrix& A) {
        PyMatrix result(A);  
        GaussInvert(result);  
        return result;  
    }

	py::object get_element(size_t row, size_t col) const {
        return std::visit([row, col](const auto& matrixPtr) -> py::object {
            using MatrixType = typename std::decay<decltype(*matrixPtr)>::type;
            return py::cast(matrixPtr->getValueAt(row, col));
        }, matrixVariant);
    }
	
	std::string toString() const {
		return std::visit([](auto&& arg) -> std::string {
			return arg->toString();
		}, matrixVariant);
	}
};	


void init_matrix(py::module_ &torus) {
	torus.def("unity", &PyMatrix::unity, "Creates an identity matrix of size NxN and type.");
	torus.def("multiply", [](const PyMatrix& A, const PyMatrix& B, PyMatrix& C) {
        PyMatrix::multiply(A, B, C);
    }, "Multiplies two matrices A and B and updates C with the result.");
	torus.def("multiplyZ", [](const PyMatrix& A, const PyMatrix& B, PyMatrix& C) {
        PyMatrix::multiplyZ(A, B, C);
    }, "Performs matrix multiplication, ignoring zero elements, and stores the result in the third matrix.");
	torus.def("GaussInvert", &PyMatrix::GaussInvert, "Sets the matrix passed in as the inverse of itself.");
	torus.def("invert", [](const PyMatrix& A, PyMatrix& A1) {
        PyMatrix::invert(A, A1);
    }, "Sets A1 to the inverse of A.");
    torus.def("inverse", [](const PyMatrix& A) {
        return PyMatrix::inverse(A);
    }, "Returns the inverse of A.");
	torus.def("as_ml_ml", [](const PyMatrix& A, const PyMatrix& B, const py::object& scalar, PyMatrix& C) {
    	A.as_ml_ml(A, B, scalar, C);
	}, "Performs matrix multiplication (only non-zero elements) followed by scalar multiplication and stores the result in C.");
	torus.def("apl_ml", [](PyMatrix& A, const py::object& scalar, PyMatrix& C) {
        PyMatrix::apl_ml(A, scalar, C);
    }, "Applies the apl_ml operation on matrices A and C with a scalar value.");
	torus.def("multiplyZ", [](PyMatrix& A, const py::object& scalar, PyMatrix& C) {
		PyMatrix::multiplyZ(A, scalar, C);
	}, "Applies the multiplyZ operation on matrices A and C with a scalar value.");
	py::class_<PyMatrix::RowProxy>(torus, "RowProxy")
	    .def("__getitem__", &PyMatrix::RowProxy::operator[], "Access a column value within the row")
		.def("__setitem__", [](PyMatrix::RowProxy& self, size_t colIndex, const py::object& value) {
            self.set(colIndex, value); 
        }, "Set an item in the row.");
	
	py::class_<PyMatrix>(torus, "Matrix")
		.def(py::init<size_t, size_t>())
		.def("set_column", &PyMatrix::setColumn, "Sets a column from a PyVector")
		.def("set_row", &PyMatrix::setRow, "Sets a row from a PyVector")
		.def(py::self += py::self)
		.def(py::self -= py::self)
		.def(py::self /= py::self)
		.def(py::self *= py::self)
		.def(py::self + py::self)
		.def(py::self * py::self)
		.def("transpose", &PyMatrix::operator!, "Returns the transpose of the matrix")
    	.def("mul_left", &PyMatrix::multiplyLeft, "Returns the transpose of a matrix multiplied by vector.")
		.def("__neg__", &PyMatrix::operator!, "Returns the transpose of the matrix")
		.def("__mul__", [](const PyMatrix &self, const PyVector &vec) {
            return self * vec;
        }, "Multiplies the matrix with a PyVector.")
			
		.def("__getitem__", [](PyMatrix &self, size_t rowIndex) -> PyMatrix::RowProxy {
            return PyMatrix::RowProxy(self, rowIndex);
        }, "Access a row proxy object")
		
		.def("__iadd__", [](PyMatrix &self, const py::object& scalar) {
            if (py::isinstance<py::int_>(scalar)) {
                self += scalar.cast<int>();
            } else if (py::isinstance<py::float_>(scalar)) {
                self += scalar.cast<double>();
			} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {	
				self += scalar.cast<std::complex<double>>();
            } else {
                throw std::runtime_error("Unsupported type for addition.");
            }
            return self;
        }, py::is_operator())
        .def("__isub__", [](PyMatrix &self, const py::object& scalar) {
            if (py::isinstance<py::int_>(scalar)) {
                self -= scalar.cast<int>();
            } else if (py::isinstance<py::float_>(scalar)) {
                self -= scalar.cast<double>();
			} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {
				self -= scalar.cast<std::complex<double>>();
            } else {
                throw std::runtime_error("Unsupported type for subtraction.");
            }
            return self;
        }, py::is_operator())
		.def("__imul__", [](PyMatrix &self, const py::object& scalar) {
            if (py::isinstance<py::int_>(scalar)) {
                self *= scalar.cast<int>();
            } else if (py::isinstance<py::float_>(scalar)) {
                self *= scalar.cast<double>();
			} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {
				self *= scalar.cast<std::complex<double>>();
            } else {
                throw std::runtime_error("Unsupported type for multiplication.");
            }
            return self;
        }, py::is_operator())
        .def("__itruediv__", [](PyMatrix &self, const py::object& scalar) {
            if (py::isinstance<py::int_>(scalar)) {
                self /= scalar.cast<int>();
            } else if (py::isinstance<py::float_>(scalar)) {
                self /= scalar.cast<double>();
			} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {
				self /= scalar.cast<std::complex<double>>();
            } else {
                throw std::runtime_error("Unsupported type for division.");
            }
            return self;
        }, py::is_operator())
		.def("__add__", [](const PyMatrix& self, const py::object& scalar) -> PyMatrix {
	        if (py::isinstance<py::int_>(scalar)) {
    	        return self + scalar.cast<int>();
        	} else if (py::isinstance<py::float_>(scalar)) {
            	return self + scalar.cast<double>();
        	} else if (py::isinstance(scalar, py::module_::import("numbers").attr("Complex"))) {
            	return self + scalar.cast<std::complex<double>>();
        	} else {
            	throw std::runtime_error("Unsupported scalar type for addition.");
        	}
    	}, "Add a scalar (int, float, or complex) to all elements of the matrix.")
		.def("column", &PyMatrix::getColumnAsPyVector, "Extracts a column and returns it as a PyVector")
		.def("row", &PyMatrix::getRowAsPyVector, "Extracts a row and returns it as a PyVector")

		.def("fill_column", &PyMatrix::fillColumn<int>, "Fills a column with an int value")
		.def("fill_column", &PyMatrix::fillColumn<double>, "Fills a column with a double value")
		.def("fill_column", &PyMatrix::fillColumn<std::complex<double>>, "Fills a column with a std::complex<double> value")

		.def("fill_row", &PyMatrix::fillRow<int>, "Fills a row with an int value")
		.def("fill_row", &PyMatrix::fillRow<double>, "Fills a row with a double value")	
		.def("fill_row", &PyMatrix::fillRow<std::complex<double>>, "Fills a row with a std::complex<double> value")

		.def("multiply_column", &PyMatrix::multiplyColumn<int>, "Multiplies a column by an int scalar")
		.def("multiply_column", &PyMatrix::multiplyColumn<double>, "Multiplies a column by a double scalar")
		.def("multiply_column", &PyMatrix::multiplyColumn<std::complex<double>>, "Multiplies a column by a std::complex<double> scalar")

		.def("multiply_row", &PyMatrix::multiplyRow<int>, "Multiplies a row by an int scalar")
		.def("multiply_row", &PyMatrix::multiplyRow<double>, "Multiplies a row by a double scalar")
		.def("multiply_row", &PyMatrix::multiplyRow<std::complex<double>>, "Multiplies a row by a std::complex<double> scalar")

		.def("set_column", &PyMatrix::setMatrixColumn, "Sets a column from a PyVector")
		.def("set_row", &PyMatrix::setMatrixRow, "Sets a row from a PyVector")
	
		.def("add_to_column", &PyMatrix::add_to_column, "Adds a PyVector to a specified column")
		.def("add_to_row", &PyMatrix::add_to_row, "Adds a PyVector to a specified row")

		.def("subtract_from_column", &PyMatrix::subtract_from_column, "Subtracts a PyVector from a specified column.")
		.def("subtract_from_row", &PyMatrix::subtract_from_row, "Subtracts a PyVector from a specified row.")

		.def("set_to_unity", &PyMatrix::setToUnity, "Sets PyMatrix to an idnetity matrix.")
		.def("trace", &PyMatrix::trace, "Calculates the trace of the matrix (sum of diagonal elements).")

		.def("__eq__", &PyMatrix::operator==, py::is_operator())
		.def("__repr__", &PyMatrix::toString);
}

