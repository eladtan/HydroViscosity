/*! \file utils.hpp
\brief A collection of utility functions
\author Elad Steinberg
*/
#ifndef UTILS_HPP
#define UTILS_HPP 1

#include <vector>
#include <stdexcept>
#include <string>

/*! \brief Calculates elementwise addition of two vectors using +=
\param lhs The vector that should be chagned
\param rhs The vector that we want to add
\return The new vector which is the elementwise addition of lhs and rhs
*/
template<typename T>
std::vector<T>& operator+=(std::vector<T>& lhs, const std::vector<T>& rhs) 
{
	if (lhs.size() != rhs.size())
		throw std::length_error("vectors must be same size to add");
	size_t N = lhs.size();
#pragma ivdep
	for (size_t i = 0; i < N; ++i)
		lhs[i] += rhs[i];
	return lhs;
}
/*! \brief Calculates elementwise addition of two vectors using +
\param lhs The left hand side of the addition
\param rhs The right hand side of the addition
\return A new vector which is the elementwise addition of lhs and rhs
*/
template<typename T>
std::vector<T> operator+ (const std::vector<T>& lhs, const std::vector<T>& rhs) 
{
	if (lhs.size() != rhs.size())
		throw std::length_error("vectors must be same size to add");
	size_t N = lhs.size();
	std::vector<T> res(N);
#pragma ivdep
	for (size_t i = 0; i < N; ++i)
		res[i] = lhs[i] + rhs[i];
	return res;
}
/*! \brief Calculates the addition of a scalar to all elements in a vector
\param s The scalar to be added
\param rhs The vector that the scalar is added to
\return A new vector which is the elementwise addition of the odl vector and the scalar
*/
template<typename T>
std::vector<T> operator+ (T const& s, const std::vector<T>& rhs)
{
	size_t N = rhs.size();
	std::vector<T> res(N);
#pragma ivdep
	for (size_t i = 0; i < N; ++i)
		res[i] = s + rhs[i];
	return res;
}
/*! \brief Calculates elementwise subtraction of two vectors using +
\param lhs The left hand side of the subtraction
\param rhs The right hand side of the subtraction
\return A new vector which is the elementwise subtraction of lhs and rhs
*/
template<typename T>
std::vector<T> operator- (const std::vector<T>& lhs, const std::vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		throw std::length_error("vectors must be same size to add");
	size_t N = lhs.size();
	std::vector<T> res(N);
#pragma ivdep
	for (size_t i = 0; i < N; ++i)
		res[i] = lhs[i] - rhs[i];
	return res;
}
/*! \brief Calculates elementwise multiplication between a vector and a scalar
\param s The scalar
\param rhs The vector
\return A new vector which is the elementwise multiplication with s
*/
template<typename T>
std::vector<T> operator* (T const& s, const std::vector<T>& rhs)
{
	size_t N = rhs.size();
	std::vector<T> res(N);
#pragma ivdep
	for (size_t i = 0; i < N; ++i)
		res[i] = s * rhs[i];
	return res;
}
/*! \brief Calculates elementwise multiplication between a vector and a scalar
\param s The scalar
\param lhs The vector
\return A new vector which is the elementwise multiplication with s
*/
template<typename T>
std::vector<T> operator* (const std::vector<T>& lhs, T const& s)
{
	return s * lhs;
}
/*! \brief Calculates elementwise multiplication between a vector and a vector
\param rhs The second vector
\param lhs The first vector
\return A new vector which is the elementwise multiplication of the two vectors
*/
template<typename T>
std::vector<T> operator* (const std::vector<T>& lhs, const std::vector<T>& rhs)
{
	size_t N = rhs.size();
	std::vector<T> res(N);
#pragma ivdep
	for (size_t i = 0; i < N; ++i)
		res[i] = lhs[i] * rhs[i];
	return res;
}
/*! \brief Calculates elementwise division of two vectors using /
\param lhs The numerator
\param rhs The denominator
\return A new vector which is the elementwise division of lhs and rhs
*/
template<typename T>
std::vector<T> operator/ (const std::vector<T>& lhs, const std::vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		throw std::length_error("vectors must be same size to add");
	size_t N = lhs.size();
	std::vector<T> res(N);
#pragma ivdep
	for (size_t i = 0; i < N; ++i)
		res[i] = lhs[i] / rhs[i];
	return res;
}

/*! \brief Writes a list of numbers to a file
\param v List of numbers
\param fname Name of the file
\param prec Precision
*/
void write_vector(std::vector<double> const& v, std::string const& fname, size_t prec = 6);
#endif // UTILS_HPP