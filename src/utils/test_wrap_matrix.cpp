#include "../../Torus/src/utils/Matrix.h"
#include <vector>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <sstream>
#include <complex>
#include <type_traits>
#include <complex>
#include <iostream>
#include <memory>
#include <variant>


template<typename T, typename U>
struct PromotedType;

template<> struct PromotedType<int, int> { using type = int; };
template<> struct PromotedType<int, double> { using type = double; };
template<> struct PromotedType<double, int> { using type = double; };
template<> struct PromotedType<double, double> { using type = double; };
template<> struct PromotedType<int, std::complex<double>> { using type = std::complex<double>; };
template<> struct PromotedType<std::complex<double>, int> { using type = std::complex<double>; };
template<> struct PromotedType<double, std::complex<double>> { using type = std::complex<double>; };
template<> struct PromotedType<std::complex<double>, double> { using type = std::complex<double>; };
template<> struct PromotedType<std::complex<double>, std::complex<double>> { using type = std::complex<double>; };

template<typename T, typename U>
using PromotedType_t = typename PromotedType<T, U>::type;

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
    virtual PyVector getColumnAsPyVector(size_t colIndex) const = 0;
    virtual PyVector getRowAsPyVector(size_t rowIndex) const = 0;  // Add this
    virtual void setColumnFromPyVector(size_t colIndex, const PyVector& pyVec) = 0;
	virtual void setRowFromPyVector(size_t rowIndex, const PyVector& pyVec) = 0;
    virtual MatrixVariant variant() const = 0;
	virtual std::string toString() const = 0;
	virtual size_t getRows() const = 0;
    virtual size_t getCols() const = 0;
	virtual py::object getElement(size_t rowIndex, size_t colIndex) const = 0;
    virtual void setElement(size_t rowIndex, size_t colIndex, py::handle value) = 0;	
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
	MatrixImpl() : data() {}
//    MatrixImpl(size_t rows, size_t cols, T defaultValue = T()) : data(rows, std::vector<T>(cols, defaultValue)) {}
	MatrixImpl(size_t rows, size_t cols, T defaultValue = T())
        : data(rows, std::vector<T>(cols, defaultValue)), rows(rows), cols(cols) {}
//	size_t getRows() const { return rows; }
//    size_t getCols() const { return cols; }

//    MatrixImpl(const MatrixImpl& M) : data(M.data) {}
    
	size_t getRows() const override {
        return data.size();
    }

    size_t getCols() const override {
        return data.empty() ? 0 : data[0].size();
    }

    /*
		ERROR HANDLING METHODS
	*/
	static void range_error() {
        throw std::runtime_error("Range error.");
    }

    static void division_by_zero_error() {
        throw std::runtime_error("Division by zero error.");
    }

    static void error(const char* message) {
        throw std::runtime_error(message);
    }
	/*
		VARIANT MATRIX
	*/   
	MatrixVariant variant() const override {
        return std::make_shared<MatrixImpl<T>>(*this); // Create shared_ptr to copy
    }

    /*
		ELEMENT ACCESS OPERATORS 
	*/
/*	
	template<typename Type>
	void setElement(size_t row, size_t col, Type value) {
        if (row >= data.size() || col >= data[0].size()) {
            throw std::out_of_range("Row or column index out of bounds.");
        }
        data[row][col] = value;
    }
*/
/*
	template<typename Type>
	void setElement(size_t row, size_t col, Type value) {
    	if (row >= data.size() || col >= data[0].size()) {
        	throw std::out_of_range("Row or column index out of bounds.");
		}
    	if constexpr (std::is_same<T, Type>::value) {
        	data[row][col] = value;
    	}
		else if constexpr (std::is_same<T, double>::value) {
			data[row][col] = static_cast<double>(value);
		}
    	else if constexpr (std::is_same<T, std::complex<double>>::value && 
                       (std::is_integral<Type>::value || std::is_floating_point<Type>::value)) {
        	data[row][col] = static_cast<std::complex<double>>(value);
    	}
    	else {
        	throw std::logic_error("Unsupported type conversion in setElement.");
    	}
	}
*/
/*
	template<typename Type>
	void setElement(size_t row, size_t col, Type value) {
    	if (row >= data.size() || col >= data[0].size()) {
        	throw std::out_of_range("Row or column index out of bounds.");
    	}
    	if constexpr (std::is_convertible<Type, T>::value) {
        	data[row][col] = static_cast<T>(value);
    	} else {
        	throw std::logic_error("Unsupported type conversion in setElement.");
    	}
	}
*/
	
	void setColumnFromPyVector(size_t colIndex, const PyVector& pyVec) override {
        if (colIndex >= getCols()) {
            throw std::out_of_range("Column index out of bounds.");
        }
        if (pyVec.size() != data.size()) {
            throw std::runtime_error("Vector size does not match number of rows.");
        }
        for (size_t i = 0; i < pyVec.size(); ++i) {
            data[i][colIndex] = pyVec.__getitem__(i).cast<T>();
        }
    }


	std::vector<T>& operator[](const int rowIndex) {
        if (rowIndex < 0 || rowIndex >= data.size()) {
            throw std::runtime_error("Row index out of bounds.");
        }
        return data[rowIndex];
    }

    const std::vector<T>& operator[](const int rowIndex) const {
        if (rowIndex < 0 || rowIndex >= data.size()) {
            throw std::runtime_error("Row index out of bounds.");
        }
        return data[rowIndex];
    }

    T& operator()(const int rowIndex, const int colIndex) {
        if (rowIndex < 0 || rowIndex >= data.size() || colIndex < 0 || colIndex >= data[0].size()) {
            throw std::runtime_error("Index out of bounds.");
        }
        return data[rowIndex][colIndex];
    }

    const T& operator()(const int rowIndex, const int colIndex) const {
        if (rowIndex < 0 || rowIndex >= data.size() || colIndex < 0 || colIndex >= data[0].size()) {
            throw std::runtime_error("Index out of bounds.");
        }
        return data[rowIndex][colIndex];
    }

	/*
		COLUMN OPERATIONS 
	*/

	PyVector column(int colIndex) const {
        std::vector<T> columnData;
        for (const auto& row : data) {
            columnData.push_back(row[colIndex]);
        }
        return PyVector(py::cast(columnData));
    }

    void fill_column(const T value, int colIndex) {
        for (auto& row : data) {
            row[colIndex] = value;
        }
    }

    void multiply_column(const T value, int colIndex) {
        for (auto& row : data) {
            row[colIndex] *= value;
        }
    }

    void set_column(const PyVector& columnValues, int colIndex) {
        if (columnValues.size() != data.size()) throw std::out_of_range("Size mismatch.");
        for (size_t i = 0; i < data.size(); ++i) {
            data[i][colIndex] = columnValues.__getitem__(i).cast<T>();  // Simplified for example
        }
    }
	void set_column(const T value, int colIndex) {
		for (auto& row : data) {
			row[colIndex] = value;
		}
	}
	void add_to_column(const PyVector& col, unsigned long colIndex) {
		if (colIndex < 0 || colIndex >= data[0].size() || col.size() != data.size()) {
			throw std::runtime_error("Column index out of bounds or size mismatch.");
		}
		for (size_t i = 0; i < data.size(); ++i) {
			py::object pyItem = col.__getitem__(i);
			T item = pyItem.cast<T>();
			data[i][colIndex] += item;
		}
	}
	void subtract_from_column(const PyVector& col, unsigned long colIndex) {
		if (colIndex < 0 || colIndex >= data[0].size() || col.size() != data.size()) {
			throw std::runtime_error("Column index out of bounds or size mismatch.");
		}
		for (size_t i = 0; i < data.size(); ++i) {
			py::object pyItem = col.__getitem__(i);
			T item = pyItem.cast<T>();
			data[i][colIndex] -= item;
		}
	}
/*	void setColumnFromPyVector(size_t colIndex, const PyVector& pyVec) override {
		if (colIndex >= data[0].size()) {
			range_error();
		}
		if (pyVec.size() != data.size()) {
			throw std::runtime_error("tuple size does not match the number of rows.");
		}
		for (size_t i = 0; i < data.size(); ++i) {
			py::object pyItem = pyVec.__getitem__(i);
			T item = pyItem.cast<T>();
			data[i][colIndex] = item;
		}
	}
*/
	
	PyVector getColumnAsPyVector(size_t colIndex) const override {
        if (colIndex >= data[0].size()) {
            range_error(); // Use static error handling method
        }
        std::vector<T> columnData;
        for (const auto& row : data) {
            columnData.push_back(row[colIndex]);
        }
        return PyVector(py::cast(columnData)); // Convert vector to PyVector (tuple)
    }
	/*
		ROW OPERATIONS
	*/
    
	PyVector row(int rowIndex) const {
        return PyVector(py::cast(data[rowIndex]));
    }

    void fill_row(const T value, int rowIndex) {
        std::fill(data[rowIndex].begin(), data[rowIndex].end(), value);
    }

    void multiply_row(const T value, int rowIndex) {
        std::transform(data[rowIndex].begin(), data[rowIndex].end(), data[rowIndex].begin(),
                       [value](const T& elem) { return elem * value; });
    }

    void set_row(const PyVector& rowValues, int rowIndex) {
        if (rowValues.size() != data[0].size()) throw std::out_of_range("Size mismatch.");
        for (size_t i = 0; i < data[0].size(); ++i) {
            data[rowIndex][i] = rowValues.__getitem__(i).cast<T>();  // Simplified for example
        }
    }
	void set_row(const T value, int rowIndex) {
		std::fill(data[rowIndex].begin(), data[rowIndex].end(), value);
	}

	void add_to_row(const PyVector& row, unsigned long rowIndex) {
		if (rowIndex < 0 || rowIndex >= data[0].size() || row.size() != data.size()) {
			throw std::runtime_error("Row index out of bounds or size mismatch.");
		}
		for (size_t i = 0; i < data.size(); ++i) {
			py::object pyItem = row.__getitem__(i);
			T item = pyItem.cast<T>();
			data[rowIndex][i] += item;
		}
	}
	void subtract_from_row(const PyVector& row, unsigned long rowIndex) {
		if (rowIndex < 0 || rowIndex >= data[0].size() || row.size() != data.size()) {
			throw std::runtime_error("Row index out of bounds or size mismatch.");
		}
		for (size_t i = 0; i < data.size(); ++i) {
			py::object pyItem = row.__getitem__(i);
			T item = pyItem.cast<T>();
			data[rowIndex][i] -= item;
		}
	}
	
	void setRowFromPyVector(size_t rowIndex, const PyVector& pyVec) override {
		if (rowIndex >= data[0].size()) {
			range_error();
		}
		if (pyVec.size() != data.size()) {
			throw std::runtime_error("Tuple size does not match the number of rows.");
		}
		for (size_t i = 0; i < data.size(); ++i) {
			py::object pyItem = pyVec.__getitem__(i);
			T item = pyItem.cast<T>();
			data[rowIndex][i] = item;
		}
	}

	PyVector getRowAsPyVector(size_t rowIndex) const override {
        if (rowIndex >= data.size()) {
            throw std::out_of_range("Row index out of bounds.");
        }
        py::list rowList;
        for (const auto& elem : data[rowIndex]) {
            rowList.append(py::cast(elem));
        }
        return PyVector(rowList);
    }
	void addScalar(const T& scalar) {
        for (auto& row : data) {
            for (auto& elem : row) {
                elem += scalar;
            }
        }
    }

    void subtractScalar(const T& scalar) {
        for (auto& row : data) {
            for (auto& elem : row) {
                elem -= scalar;
            }
        }
    }

    void multiplyScalar(const T& scalar) {
        for (auto& row : data) {
            for (auto& elem : row) {
                elem *= scalar;
            }
        }
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

    template<typename T>
    void updateMatrixVariantType() {
        if (!std::holds_alternative<std::shared_ptr<MatrixImpl<T>>>(matrixVariant)) {
            matrixVariant = std::make_shared<MatrixImpl<T>>();
        }
    }

public:
	PyMatrix(size_t rows, size_t cols) : matrixVariant(std::make_shared<MatrixImpl<double>>(rows, cols)) {}
	const MatrixVariant& getMatrixVariant() const { return matrixVariant; }

/*
	template<typename TargetType>
	void promoteAndSetColumn(size_t colIndex, const PyVector& pyVec) {
    	auto targetMatrix = std::visit([&](auto&& arg) -> std::shared_ptr<MatrixImpl<TargetType>> {
        	using CurrentType = typename std::decay_t<decltype(arg)>::element_type::value_type;
        	if constexpr (std::is_same_v<TargetType, CurrentType>) {
            	return arg; // No promotion needed, types match
        	} else {
            	auto promotedMatrix = std::make_shared<MatrixImpl<TargetType>>();
            	return promotedMatrix;
        	}
    	}, matrixVariant);

    	matrixVariant = targetMatrix;

    	for (size_t i = 0; i < pyVec.size(); ++i) {
        	py::handle item = pyVec.__getitem__(i);
        	targetMatrix->setElement(i, colIndex, item.cast<TargetType>());
    	}
	}
*/
/*
	template<typename TargetType>
	void promoteAndSetColumn(size_t colIndex, const PyVector& pyVec) {
    	auto promoteMatrix = [&]() -> std::shared_ptr<BaseMatrix> {
        	if (auto currentPtr = std::get_if<std::shared_ptr<MatrixImpl<TargetType>>>(&matrixVariant)) {
            	return *currentPtr;
        	} else {
            	auto newMatrix = std::make_shared<MatrixImpl<TargetType>>();
            	matrixVariant = newMatrix; 
            	return newMatrix;
        	}
		};
    };
*/
/* stuff
	template<typename TargetType>
	void promoteAndSetColumn(size_t colIndex, const PyVector& pyVec) {
    	size_t currentRows = std::visit([](const auto& matrixPtr) -> size_t { return matrixPtr->getRows(); }, matrixVariant);
    	size_t currentCols = std::visit([](const auto& matrixPtr) -> size_t { return matrixPtr->getCols(); }, matrixVariant);

    	std::shared_ptr<BaseMatrix> promotedMatrix = [&]() -> std::shared_ptr<BaseMatrix> {
        	if (std::holds_alternative<std::shared_ptr<MatrixImpl<TargetType>>>(matrixVariant)) {
				return std::get<std::shared_ptr<MatrixImpl<TargetType>>>(matrixVariant);
			} else {
            	throw std::runtime_error("unknwon");
				auto newMatrix = std::make_shared<MatrixImpl<TargetType>>(currentRows, currentCols, TargetType{});
            	matrixVariant = newMatrix; 
            	return newMatrix;
			}
    	}();

    	auto targetMatrix = std::dynamic_pointer_cast<MatrixImpl<TargetType>>(promotedMatrix);

    	if (!targetMatrix) {
        	throw std::runtime_error("Failed to promote matrix to the target type.");
    	}
		for (size_t i = 0; i < pyVec.size(); ++i) {
			py::handle item = pyVec.__getitem__(i);
			TargetType value = item.cast<TargetType>();
			targetMatrix->setElement(i, colIndex, value);
		}
    	targetMatrix->setColumnFromPyVector(colIndex, pyVec);
	}
*/
/*
	template<typename TargetType>
    void promoteAndSetColumn(size_t colIndex, const PyVector& pyVec) {
        // Promote the matrix to the necessary type, if not already of that type.
        if (!std::holds_alternative<std::shared_ptr<MatrixImpl<TargetType>>>(matrixVariant)) {
            size_t rows = std::visit([](auto&& arg) { return arg->getRows(); }, matrixVariant);
            size_t cols = std::visit([](auto&& arg) { return arg->getCols(); }, matrixVariant);
            matrixVariant = std::make_shared<MatrixImpl<TargetType>>(rows, cols);
        }

        // Now that we're sure matrixVariant holds the correct type, set the column.
        auto matrixPtr = std::get<std::shared_ptr<MatrixImpl<TargetType>>>(matrixVariant);
        for (size_t i = 0; i < pyVec.size(); ++i) {
            TargetType value = pyVec.__getitem__(i).cast<TargetType>();
            matrixPtr->setElement(i, colIndex, value);
        }
    }
*/
/* works wit complex
	template<typename T>
    void promoteAndSetColumn(size_t colIndex, const PyVector& pyVec) {
        auto targetMatrix = std::visit([&](auto& arg) -> std::shared_ptr<BaseMatrix> {
            using MatrixType = std::decay_t<decltype(*arg)>;
            if constexpr (!std::is_same_v<MatrixType, MatrixImpl<T>>) {
                auto newMatrix = std::make_shared<MatrixImpl<T>>(arg->getRows(), arg->getCols());
                matrixVariant = newMatrix;
                return newMatrix;
            } else {
                return arg;
            }
        }, matrixVariant);

        targetMatrix->setColumnFromPyVector(colIndex, pyVec);
    }
*/

	template<typename T>
	void promoteAndSetColumn(size_t colIndex, const PyVector& pyVec) {
    	if (!std::holds_alternative<std::shared_ptr<MatrixImpl<T>>>(matrixVariant)) {
        	size_t rows = std::visit([](auto&& arg) { return arg->getRows(); }, matrixVariant);
        	size_t cols = std::visit([](auto&& arg) { return arg->getCols(); }, matrixVariant);
        	matrixVariant = std::make_shared<MatrixImpl<T>>(rows, cols);
    	}
/*
    	auto matrixPtr = std::get<std::shared_ptr<MatrixImpl<T>>>(matrixVariant);
    	for (size_t i = 0; i < pyVec.size(); ++i) {
        	T value = pyVec.__getitem__(i).cast<T>();
        	matrixPtr->setElement(i, colIndex, value);
    	}
*/
		std::shared_ptr<MatrixImpl<T>> matrixPtr = std::get<std::shared_ptr<MatrixImpl<T>>>(matrixVariant);
    	matrixPtr->setColumnFromPyVector(colIndex, pyVec);
	
	}

/*
	template<typename TargetType>
    void promoteAndSetColumn(size_t colIndex, const PyVector& pyVec) {
        std::visit([&](auto& arg) {
            using MatType = typename std::remove_cv_t<typename std::remove_reference_t<decltype(*arg)>::type>;
            if constexpr (!std::is_same_v<MatType, MatrixImpl<TargetType>>) {
                auto newMatrix = std::make_shared<MatrixImpl<TargetType>>(arg->getRows(), arg->getCols(), TargetType{});
                matrixVariant = newMatrix;
            }
            std::get<std::shared_ptr<MatrixImpl<TargetType>>>(matrixVariant)->setColumnFromPyVector(colIndex, pyVec);
        }, matrixVariant);


	    auto targetMatrix = std::static_pointer_cast<MatrixImpl<TargetType>>(promoteMatrix());

    	for (size_t i = 0; i < pyVec.size(); ++i) {
        	py::handle item = pyVec.__getitem__(i);
        	TargetType value = item.cast<TargetType>(); // Cast each item to TargetType
        	targetMatrix->setElement(i, colIndex, value);
    	}
	};
*/
/* stuff
	void setColumnFromPyVector(size_t colIndex, const PyVector& pyVec) {
    	py::handle firstElem = pyVec.__getitem__(0);
    	if (py::isinstance<py::int_>(firstElem)) {
        	promoteAndSetColumn<int>(colIndex, pyVec);
    	} else if (py::isinstance<py::float_>(firstElem)) {
        	promoteAndSetColumn<double>(colIndex, pyVec);
    	} else if (py::isinstance(firstElem, py::module_::import("numbers").attr("Complex"))) {
        	promoteAndSetColumn<std::complex<double>>(colIndex, pyVec);
    	} else {
        	throw std::runtime_error("Unsupported type in PyVector.");
    	}
	}
*/
/* works
	void setColumnFromPyVector(size_t colIndex, const PyVector& pyVec) {
        py::handle firstElem = pyVec.__getitem__(0);
        if (py::isinstance<py::int_>(firstElem)) {
            promoteAndSetColumn<int>(colIndex, pyVec);
        } else if (py::isinstance<py::float_>(firstElem)) {
            promoteAndSetColumn<double>(colIndex, pyVec);
		} else if (py::isinstance(firstElem, py::module_::import("numbers").attr("Complex"))) {   
			promoteAndSetColumn<std::complex<double>>(colIndex, pyVec);
        } else {
            throw std::runtime_error("Unsupported type in PyVector.");
        }
    }
*/
/*
	void setColumnFromPyVector(size_t colIndex, const PyVector& pyVec) {
	    py::handle firstElem = pyVec.__getitem__(0);
    	if (py::isinstance<py::int_>(firstElem)) {
        	promoteAndSetColumn<int>(colIndex, pyVec);
	    } else if (py::isinstance<py::float_>(firstElem)) {
    	    promoteAndSetColumn<double>(colIndex, pyVec);
		} else if (py::isinstance(firstElem, py::module_::import("numbers").attr("Complex"))) {     
			promoteAndSetColumn<std::complex<double>>(colIndex, pyVec);
    	} else {
        	throw std::runtime_error("Unsupported type in PyVector.");
    	}
	}
*/
	void setColumnFromPyVector(size_t colIndex, const PyVector& pyVec) {
    try {
        py::handle firstElem = pyVec.__getitem__(0);
        if (py::isinstance<py::int_>(firstElem)) {
            promoteAndSetColumn<int>(colIndex, pyVec);
        } else if (py::isinstance<py::float_>(firstElem)) {
            if (!std::holds_alternative<std::shared_ptr<MatrixImpl<double>>>(matrixVariant)) {
                size_t rows = std::visit([](auto&& arg) { return arg->getRows(); }, matrixVariant);
                size_t cols = std::visit([](auto&& arg) { return arg->getCols(); }, matrixVariant);
                matrixVariant = std::make_shared<MatrixImpl<double>>(rows, cols);
            }
            auto matrixPtr = std::get<std::shared_ptr<MatrixImpl<double>>>(matrixVariant);
            for (size_t i = 0; i < pyVec.size(); ++i) {
                double value = pyVec.__getitem__(i).cast<double>();
                matrixPtr->setElement(i, colIndex, value);
            }
		} else if (py::isinstance(firstElem, py::module_::import("numbers").attr("Complex"))) {
			promoteAndSetColumn<std::complex<double>>(colIndex, pyVec);
        } else {
            throw std::runtime_error("Unsupported type in PyVector.");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error in setColumnFromPyVector: " << e.what() << std::endl;
        throw; // Rethrow to ensure Python sees it.
    }
}

	 MatrixVariant& getMatrixVariant() {
        return matrixVariant;
    }

	py::object getElement(size_t rowIndex, size_t colIndex) {
        return std::visit([rowIndex, colIndex](auto&& arg) -> py::object {
            return arg->getElement(rowIndex, colIndex);
        }, matrixVariant);
    }

    void setElement(size_t rowIndex, size_t colIndex, py::handle value) {
        std::visit([rowIndex, colIndex, &value](auto&& arg) {
            arg->setElement(rowIndex, colIndex, value);
        }, matrixVariant);
    }

    std::string toString() const {
        return std::visit([](auto&& arg) -> std::string {
            return arg->toString();
        }, matrixVariant);
    }	
};


void init_matrix(py::module_ &torus) {
	py::class_<PyMatrix>(torus, "Matrix")
		/*
		.def(py::init<size_t, size_t, const std::string&>())

		.def("fill_column", [](PyMatrix& self, double value, int colIndex) {
			self.fill_column(value, colIndex);
		})
		.def("multiply_column", [](PyMatrix& self, double value, int colIndex) {
			self.multiply_column(value, colIndex);
		})
		.def("set_column", &PyMatrix::set_column_from_pyvector)
		.def("set_row", &PyMatrix::set_row_from_pyvector)
		.def("__getitem__", [](const PyMatrix& self, int rowIndex) -> PyVector {
			return self.getRow(rowIndex);
		})
		.def("__repr__", &PyMatrix::toString)
		.def("get_row", &PyMatrix::getRow)
		.def("get_column", &PyMatrix::getColumn);
		*/

		.def(py::init<size_t, size_t>())
        .def("set_column", [](PyMatrix &self, size_t colIndex, const PyVector& pyVec) {
        	self.setColumnFromPyVector(colIndex, pyVec);
    	})
		.def("__repr__", &PyMatrix::toString);
}

