//
// Created by kamil on 4/22/18.
//
#ifndef MN_SYSTEMSOFLINEAREQUATIONS_TOOLS_H
#define MN_SYSTEMSOFLINEAREQUATIONS_TOOLS_H
#include <vector>
#include <algorithm>
namespace tools {
	const auto NOW = std::chrono::high_resolution_clock::now;
	using DoubleMicroseconds = std::chrono::duration<double, std::micro>;
	using DoubleMilliseconds = std::chrono::duration<double, std::milli>;

	template <typename Time_t, typename F>
	auto timeFunction(F func) {
		const auto start = NOW();
		std::invoke(func);
		return std::chrono::duration_cast<Time_t>(NOW() - start);
	}

	template<typename T>
	std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>&b) {
		auto result = std::vector<T>();
		std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<>());
		return result;
	}

}
#endif //MN_SYSTEMSOFLINEAREQUATIONS_TOOLS_H