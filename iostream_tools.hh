//
// Created by kamil on 3/30/18.
//

#ifndef MN_SYSTEMSOFLINEAREQUATIONS_OSTREAM_TOOLS_H
#define MN_SYSTEMSOFLINEAREQUATIONS_OSTREAM_TOOLS_H

#include <iostream>
#include <vector>
#include <chrono>

namespace iostream_tools {

	template<typename T>
	std::ostream &operator<<(std::ostream &o, const std::vector<T> &v) {
		o << "vector<" << typeid(T).name() << ">[" << v.size() << "]\n";
		for (auto& val : v) o << val << '\n';
		return o << "\n";
	}

	template<typename _Rep, typename _Period>
	std::ostream& operator<<(std::ostream& o, const std::chrono::duration<_Rep, _Period>& duration) {
		return o << std::setprecision(3)<< duration.count();
	}
};


#endif //MN_SYSTEMSOFLINEAREQUATIONS_OSTREAM_TOOLS_H
