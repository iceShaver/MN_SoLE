//
// Created by kamil on 4/6/18.
//

#ifndef MN_SYSTEMSOFLINEAREQUATIONS_EXCEPTION_HH
#define MN_SYSTEMSOFLINEAREQUATIONS_EXCEPTION_HH
class Exception : std::exception {};
class MatrixSizeMismatchException : Exception{};
class NotASquareMatrixException : Exception{};
#endif //MN_SYSTEMSOFLINEAREQUATIONS_EXCEPTION_HH
