/*! \file Courant.hpp
\brief Class for calculating the time step based on the courant condition dt = cfl * dx / speed_of_sound
\author Elad Steinberg
*/
#ifndef COURANT_HPP
#define COURANT_HPP 1

#include <vector>
/*! \brief Class for calculating the time step based on the courant condition dt = cfl * dx / speed_of_sound
*/
class Courant
{
private:
	double cfl_;
public:
	/*! \brief Class constructor
	\param cfl The courant factor
	*/
	Courant(double cfl);
	/*! \brief Calculates the time step
	\param edges The location of the interfaces
	\param pressure The pressure inside the computational domain
	\param density The density inside the computational domain
	\param gamma The adiabatic index (assumes ideal gas law)
	\return The calculated time step
	*/
	double operator()(std::vector<double> const& edges, std::vector<double> const& pressure, std::vector<double> const& density, double gamma) const;
};

#endif // COURANT_HPP