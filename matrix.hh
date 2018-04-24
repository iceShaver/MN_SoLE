//
// Created by kamil on 4/6/18.
//

#ifndef MN_SYSTEMSOFLINEAREQUATIONS_MATRIX_HH
#define MN_SYSTEMSOFLINEAREQUATIONS_MATRIX_HH

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <typeinfo>
#include <functional>
#include "exception.hh"



using std::cout;
using std::endl;
using std::setw;
using std::transform;
using std::back_inserter;
using std::move;
using std::ostream;
using std::vector;
using std::unique_ptr;
using std::shared_ptr;
using std::plus;
using std::minus;
using std::multiplies;
using std::divides;
using std::initializer_list;
using std::conditional;
using std::bidirectional_iterator_tag;

//____________________________________MATRIX HEADER_____________________________________________
////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class Matrix {
public:
	template<bool isConstIterator = false>
	class IteratorTemplate;
	using ConstIterator = IteratorTemplate<true>;
	using Iterator = IteratorTemplate <false>;
	Matrix(size_t rows, size_t cols, const T &initValue = T());
	Matrix(const initializer_list<initializer_list<T>> &initList);
	Matrix(const Matrix &other);
	Matrix(Matrix &&other) noexcept;
	virtual ~Matrix() = default;
	T &operator()(size_t row, size_t col);
	T operator()(size_t row, size_t col) const;
	Matrix<T> operator+(const Matrix &other) const;
	Matrix<T> operator-(const Matrix &other) const;
	Matrix<T> operator*(const Matrix &other) const;
	vector<T> operator*(const vector<T> &vec) const;
	size_t rows() const noexcept;
	size_t cols() const noexcept;
	bool isSquare()const noexcept;
	Iterator begin();
	Iterator end();
	ConstIterator begin()const;
	ConstIterator end()const;
	friend ostream &operator<<(ostream &o, const Matrix &m) {
		o << "Matrix<" << typeid(T).name() << ">[" << m.rows_ << ", " << m.cols_ << "]\n";
		for (int row = 0; row < m.rows_; ++row) {
			for (int col = 0; col < m.cols_; ++col) {
				o << setw(6) << m.operator()(row, col) << ", ";
			}
			o << "\n";
		}
		return o;
	}

private:
	vector<T> data_;
	size_t rows_{};
	size_t cols_{};
	Matrix(size_t rows, size_t cols, vector<T> data);
	Matrix() = default;
};




//____________________________________MATRIX  ITERATOR HEADER_____________________________________________
//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
template<bool isConstIterator = true>
class Matrix<T>::IteratorTemplate {
	using TRef = typename conditional<isConstIterator, const T&, T&>::type;
	using VectorIterator = typename conditional<isConstIterator, typename vector<T>::const_iterator, typename vector<T>::iterator>::type;
public:
	IteratorTemplate(const IteratorTemplate& other);
	IteratorTemplate(IteratorTemplate&& other) noexcept;
	virtual ~IteratorTemplate() = default;
	virtual TRef operator*();
	IteratorTemplate& operator++();
	IteratorTemplate operator++(int);
	IteratorTemplate& operator--();
	IteratorTemplate operator--(int);
	IteratorTemplate& operator+=(size_t offset);
	IteratorTemplate operator+(size_t offset) const;
	IteratorTemplate& operator-=(size_t offset);
	IteratorTemplate operator-(size_t offset) const;
	bool operator==(const IteratorTemplate& other) const;
	bool operator!=(const IteratorTemplate& other) const;
	bool operator>(const IteratorTemplate& other) const;
	bool operator>=(const IteratorTemplate& other) const;
	bool operator<(const IteratorTemplate& other) const;
	bool operator<=(const IteratorTemplate& other) const;

	// types required for STL
	using difference_type = long long;
	using value_type = T;
	using pointer = const T*;
	using reference = TRef;
	using iterator_category = std::bidirectional_iterator_tag;
private:
	explicit IteratorTemplate(VectorIterator iter);
	VectorIterator vecIterator_;
};



//____________________________________________MATRIX IMPLEMENTATION_________________________________________________
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> Matrix<T>::Matrix(const initializer_list<initializer_list<T>> &initList) {
	const auto cols = initList.begin()->size();
	for (const auto &l : initList) {
		if (l.size() != cols) throw MatrixSizeMismatchException();
		data_.insert(data_.end(), l.begin(), l.end());
	}
	rows_ = initList.size();
	cols_ = data_.size() / rows_;
}
template <typename T>
Matrix<T>::Matrix(const size_t rows, const size_t cols, const T& initValue) :
	data_(rows * cols, initValue), rows_(rows), cols_(cols) {}

template <typename T>
Matrix<T>::Matrix(const size_t rows, const size_t cols, vector<T> data) :
	data_(move(data)), rows_(rows), cols_(cols) {}


template <typename T>
Matrix<T>::Matrix(const Matrix& other) :
	data_(other.data_), rows_(other.rows_), cols_(other.cols_) {}

template <typename T>
Matrix<T>::Matrix(Matrix&& other) noexcept :
	data_(move(other.data_)), rows_(other.rows_), cols_(other.cols_) {}

template<typename T>
T& Matrix<T>::operator()(const size_t row, const size_t col) {
	return data_[row * cols_ + col];
}

template<typename T>
T Matrix<T>::operator()(const size_t row, const size_t col) const {
	return data_[row * cols_ + col];
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &other) const {
	if (data_.size() != other.data_.size())
		throw MatrixSizeMismatchException();
	auto result = Matrix<T>();
	result.cols_ = cols_;
	result.rows_ = rows_;
	result.data_.reserve(data_.size());
	transform(begin(), end(), other.begin(), back_inserter(result.data_), plus<>());
	return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix &other) const {
	if (data_.size() != other.data_.size()) {
		throw MatrixSizeMismatchException();
	}
	auto result = Matrix<T>();
	result.cols_ = cols_;
	result.rows_ = rows_;
	result.data_.reserve(data_.size());
	transform(data_.begin(), data_.end(), other.data_.begin(), back_inserter(result.data_), minus<>());
	return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &other) const {
	if (cols_ != other.rows_)
		throw MatrixSizeMismatchException();
	auto result = Matrix<T>(rows_, cols_, vector<T>(data_.size()));
	for (int row = 0; row < rows_; ++row) {
		for (int col = 0; col < cols_; ++col) {
			for (int val = 0; val < cols_; ++val) {
				result(row, col) += this->operator()(row, val) * other(val, col);
			}
		}
	}
	return result;
}

template<typename T>
vector<T> Matrix<T>::operator*(const vector<T> &vec) const {
	if (rows_ != vec.size())
		throw MatrixSizeMismatchException();
	auto result = vector<T>(cols_, 0);
	for (int row = 0; row < rows_; ++row) {
		for (int col = 0; col < cols_; ++col) {
			result[row] += vec[col] * this->operator()(row, col);
		}
	}
	return result;
}

template <typename T>
size_t Matrix<T>::rows() const noexcept {
	return rows_;
}

template <typename T>
size_t Matrix<T>::cols() const  noexcept {
	return cols_;
}

template <typename T>
bool Matrix<T>::isSquare() const noexcept {
	return cols_ == rows_;
}

template <typename T>
typename Matrix<T>::Iterator
Matrix<T>::begin() {
	return Iterator(data_.begin());
}

template <typename T>
typename Matrix<T>::Iterator
Matrix<T>::end() {
	return Iterator(data_.end());
}

template <typename T>
typename Matrix<T>::ConstIterator
Matrix<T>::begin() const {
	return ConstIterator(data_.cbegin());
}

template <typename T>
typename Matrix<T>::ConstIterator
Matrix<T>::end() const {
	return ConstIterator(data_.cend());
}



//____________________________________MATRIX ITERATOR IMPLEMENTATION_____________________________________________
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename T>
template <bool isConstIterator>
Matrix<T>::IteratorTemplate<isConstIterator>::IteratorTemplate(const IteratorTemplate& other) :
	vecIterator_(other.vecIterator_) {}

template <typename T>
template <bool isConstIterator>
Matrix<T>::IteratorTemplate<isConstIterator>::IteratorTemplate(IteratorTemplate&& other) noexcept :
	vecIterator_(move(other.vecIterator_)) {}

template <typename T>
template <bool isConstIterator>
Matrix<T>::IteratorTemplate<isConstIterator>::IteratorTemplate(VectorIterator iter) :
	vecIterator_(iter) {}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>::TRef
Matrix<T>::IteratorTemplate<isConstIterator>::operator*() {
	return *vecIterator_;
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>&
Matrix<T>::IteratorTemplate<isConstIterator>::operator++() {
	auto result = *this;
	++vecIterator_;
	return result;
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>
Matrix<T>::IteratorTemplate<isConstIterator>::operator++(int) {
	++vecIterator_;
	return *this;
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>&
Matrix<T>::IteratorTemplate<isConstIterator>::operator--() {
	auto result = *this;
	--vecIterator_;
	return result;
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>
Matrix<T>::IteratorTemplate<isConstIterator>::operator--(int) {
	--vecIterator_;
	return *this;
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>&
Matrix<T>::IteratorTemplate<isConstIterator>::operator+=(size_t offset) {
	vecIterator_ += offset;
	return *this;
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>
Matrix<T>::IteratorTemplate<isConstIterator>::operator+(size_t offset) const {
	return IteratorTemplate{ vecIterator_ + offset };
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>&
Matrix<T>::IteratorTemplate<isConstIterator>::operator-=(size_t offset) {
	vecIterator_ -= offset;
	return *this;
}

template <typename T>
template <bool isConstIterator>
typename Matrix<T>::IteratorTemplate<isConstIterator>
Matrix<T>::IteratorTemplate<isConstIterator>::operator-(size_t offset) const {
	return IteratorTemplate{ vecIterator_ - offset };
}

template <typename T>
template <bool isConstIterator>
bool
Matrix<T>::IteratorTemplate<isConstIterator>::operator==(const IteratorTemplate& other) const {
	return vecIterator_ == other.vecIterator_;
}

template <typename T>
template <bool isConstIterator>
bool
Matrix<T>::IteratorTemplate<isConstIterator>::operator!=(const IteratorTemplate& other) const {
	return vecIterator_ != other.vecIterator_;
}

template <typename T>
template <bool isConstIterator>
bool
Matrix<T>::IteratorTemplate<isConstIterator>::operator>(const IteratorTemplate& other) const {
	return vecIterator_ > other.vecIterator_;
}

template <typename T>
template <bool isConstIterator>
bool
Matrix<T>::IteratorTemplate<isConstIterator>::operator>=(const IteratorTemplate& other) const {
	return vecIterator_ >= other.vecIterator_;
}

template <typename T>
template <bool isConstIterator>
bool
Matrix<T>::IteratorTemplate<isConstIterator>::operator<(const IteratorTemplate& other) const {
	return vecIterator_ < other.vecIterator_;
}

template <typename T>
template <bool isConstIterator>
bool
Matrix<T>::IteratorTemplate<isConstIterator>::operator<=(const IteratorTemplate& other) const {
	return vecIterator_ <= other.vecIterator_;
}

#endif //MN_SYSTEMSOFLINEAREQUATIONS_MATRIX_HH
