//
// Created by kamil on 3/30/18.
//
#include <iostream>
#include <fstream>
#include <chrono>
#include "matrix.hh"
#include "linear_system.hh"
#include "iostream_tools.hh"
#include "tools.hh"
#include <thread>

using std::cout;
using std::endl;
using std::generate;
using std::invoke;
using std::flush;
using std::ofstream;
using namespace iostream_tools;
using namespace tools;
using namespace std::chrono;
using namespace linear_system;

constexpr const unsigned INDEX_NO = 155555;
constexpr const auto STOP_NORM = 10e-9;


const auto N = getN(INDEX_NO);
const auto e = getDigit(INDEX_NO, 3);
const auto f = getDigit(INDEX_NO, 2);
const auto a1 = 5.0 + e;
const auto a2 = -1;
const auto a3 = -1;





int main(int argc, const char **argv) {
	//_____________________________________Zadanie A____________________________________________

	auto A = genA(N, a1, a2, a3);
	auto b = genB(N, f);

	//cout << A << endl;
	//cout << b << endl;

	//_____________________________________Zadanie B____________________________________________
	{	// Jacobi
		auto start = NOW();
		auto[x, iters] = solveJacobi(A, b, STOP_NORM);
		auto stop = NOW();
		//cout << x << '\n';
		cout << "Jacobi method\n" <<
			"Iters: " << iters <<
			"\nTime: " << duration_cast<DoubleMilliseconds>(stop - start)<< " ms" <<
			"\nNorm: " << euclideanNorm(A*x - b) << "\n\n";
	}

	{	// Gauss-Seidel
		auto start = NOW();
		auto[x, iters] = solveGaussSeidel(A, b, STOP_NORM);
		auto stop = NOW();
		//cout << x << '\n';
		cout << "Gauss-Seidel method\n" <<
			"Iters: " << iters <<
			"\nTime: " << duration_cast<DoubleMilliseconds>(stop - start) << " ms" <<
			"\nNorm: " << euclideanNorm(A*x - b) << "\n\n";
	}


	//_____________________________________Zadanie C____________________________________________

	// Dla podanych warto�ci metody iteracyjne nie zbiegaj� si�
	/*{
		auto M = gen_A(getN(INDEX_NO), 3, -1, -1);
		auto[x, iters] = solveJacobi(M, b, 10e-9);
		auto[x1, iters1] = solveGaussSeidel(M, b, 10e-9);
	}*/



	//_____________________________________Zadanie D____________________________________________
	{ // LU
		auto A = genA(N, 3, -1, -1);

		auto start = NOW();
		auto x = luFactorization(A, b);
		auto stop = NOW();
		//cout << x << '\n';
		cout << "LU factorization\n" <<
			"Time: " << duration_cast<DoubleMilliseconds>(stop - start) << " ms\n" <<
			"Norm: " << euclideanNorm(A*x - b) << "\n\n";

	}


	//_____________________________________Zadanie E____________________________________________
	{
		auto numbersOfTests = 7;

		// Prepare tests
		cout << string(5, '*') + "Time testing" + string(5, '*') + "\n\n";
		auto N = vector<size_t>(numbersOfTests);
		auto jacobiTimes = vector<DoubleMilliseconds>();
		auto gaussSeidelTimes = vector<DoubleMilliseconds>();
		auto luTimes = vector<DoubleMilliseconds>();
		N[0] = 100;
		N[1] = 500;
		generate(N.begin() + 2, N.end(), [n = 0]() mutable { return n += 1000; }); // [100, 500, 1000, 2000, 3000, ...]


		// test all Ns
		for (auto &n : N) {
			cout << string(4, '_') << "N = " << n << string(4, '_') << '\n';
			auto A = genA(n, a1, a2, a3);
			auto b = genB(n, f);

			// time Jacobi
			jacobiTimes.push_back(timeFunction<DoubleMilliseconds>([&]() { solveJacobi(A, b, STOP_NORM); }));
			cout << "Jacobi: " << jacobiTimes.back() << " ms\n";

			// time Gauss-Seidel
			gaussSeidelTimes.push_back(timeFunction<DoubleMilliseconds>([&]() { solveGaussSeidel(A, b, STOP_NORM); }));
			cout << "Gauss-Seidel: " << gaussSeidelTimes.back() << " ms\n";

			// time LU factorization
			luTimes.push_back(timeFunction<DoubleMilliseconds>([&]() { luFactorization(A, b); }));
			cout << "LU factorization: " << luTimes.back() << " ms\n\n";
		}
		cout << flush;

		// Write result to file
		auto fs = ofstream("zadanieE.txt");
		fs << "N\t" << "Jacobi\t" << "Gauss-Seidel\t" << "LU factorization\n";
		for (int i = 0; i < N.size(); ++i) {
			fs << N[i] << '\t' << jacobiTimes[i] << '\t' << gaussSeidelTimes[i] << '\t' << luTimes[i] << '\n';
		}
	}

	system("pause");
	return 0;
}
