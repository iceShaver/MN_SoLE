#include "linear_system.hh"

namespace linear_system {

	unsigned getDigit(unsigned number, unsigned index) {
		return unsigned((std::to_string(number)[index] - '0'));
	}

	unsigned getN(unsigned index_no) {
		auto d = getDigit(index_no, 5); // last digit of index no
		auto c = getDigit(index_no, 4); // second last digit of index no
		return 9 * c * d;
	}

	Matrix<double> genA(size_t n, double a1, double a2, double a3) {
		auto result = Matrix<double>(n, n);
		for (int i = 0; i < result.rows(); ++i)
			for (int j = max(i - 2, 0); j < min(i + 3, int(result.cols())); ++j)
				switch (abs(i - j)) {
				case 0:
					result(i, j) = a1;
					break;
				case 1:
					result(i, j) = a2;
					break;
				default:
					result(i, j) = a3;
				}
		return result;
	}
	vector<double> genB(unsigned n, size_t f) {
		auto result = vector<double>(n);
		std::generate(result.begin(), result.end(), [f, n = 0]()mutable { return std::sin((n++) * (f + 1)); });
		return result;
	}
}