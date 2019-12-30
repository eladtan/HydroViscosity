/*! \file Boundary.hpp
\brief Abstract class for boundary conditions as well as instances of it.
\author Elad Steinberg
*/

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP 1

#include <vector>
/*! \brief Abstract class for boundary conditions
*/
class Boundary
{
public:
	/*! \brief Calculates the pressure inside the cell to the side of the domain
	\param left_side Flag for choosing left boundary or right boundary
	\param pressure The pressure inside the computational domain
	\return The pressure outside the computational domain
	*/
	virtual double GetSidePressure(bool left_side, std::vector<double> const& pressure)const = 0;
	/*! \brief Calculates the mass inside the cell to the side of the domain
	\param left_side Flag for choosing left boundary or right boundary
	\param mass The masses inside the computational domain
	\return The mass outside the computational domain
	*/
	virtual double GetSideMass(bool left_side, std::vector<double> const& mass)const = 0;
	/*! \brief Calculates the viscosity inside the cell to the side of the domain
	\param left_side Flag for choosing left boundary or right boundary
	\param viscosity The viscosity inside the computational domain
	\return The viscosity outside the computational domain
	*/
	virtual double GetSideViscosity(bool left_side, std::vector<double> const& viscosity)const = 0;
};

/*! \brief Class for rigid wall boundary conditions, the side cell is duplicated and it's velocity changes sign
*/
class RigidWall : public Boundary
{
public:
	double GetSidePressure(bool left_side, std::vector<double> const& pressure)const;

	double GetSideMass(bool left_side, std::vector<double> const& mass)const;

	double GetSideViscosity(bool left_side, std::vector<double> const& viscosity)const;
};

#endif // BOUNDARY_HPP