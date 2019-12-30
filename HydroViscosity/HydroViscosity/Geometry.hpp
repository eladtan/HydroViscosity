/*! \file Geometry.hpp
\brief Abstract class for geometry and its instances
\author Elad Steinberg
*/
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1

#include <vector>

/*! \brief Abstract class for geometry
*/
class Geometry
{
public:
	/*! \brief Calculates the area of an interface
	\param edges The location of the interfaces
	\param index The index of the requested interface area
	\return The area of the interface
	*/
	virtual double CalcArea(std::vector<double> const& edges, size_t index)const = 0;
	/*! \brief Calculates the volume of a cell
	\param edges The location of the interfaces
	\param index The index of the requested cell volume
	\return The volume of the cell
	*/
	virtual double CalcVolume(std::vector<double> const& edges, size_t index)const = 0;
};

/*! \brief Class for planar geometry, all interface have area unity
*/
class Planar : public Geometry
{
public:
	double CalcArea(std::vector<double> const& edges, size_t index)const;

	double CalcVolume(std::vector<double> const& edges, size_t index)const;
};
/*! \brief Class for spherical geometry
*/
class Spherical : public Geometry
{
public:
	double CalcArea(std::vector<double> const& edges, size_t index)const;

	double CalcVolume(std::vector<double> const& edges, size_t index)const;
};

#endif // GEOMETRY_HPP