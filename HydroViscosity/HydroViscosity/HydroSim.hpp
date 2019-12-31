/*! \file HydroSim.hpp
\brief Class for running a 1D hydro simulation. For the RM method the mass, pressure, viscosity and internal energy are cell centered while the velocity is edge centered. For the Godunov scheme, all the variables are cell centered.
\author Elad Steinberg
*/
#ifndef HYDROSIM_HPP
#define HYDROSIM_HPP 1

#include "Courant.hpp"
#include "Geometry.hpp"
#include "Boundary.hpp"
#include <string>
#include "Hllc.hpp"

//! Class for running a 1D hydro simulation
class HydroSim
{
private:
	Courant const& courant_;
	const double gamma_;
	Geometry const& geo_;
	Boundary const& boundary_left_; 
	Boundary const& boundary_right_;
	bool godunov_;
	double time_;
	size_t cycle_;
	const double viscosity_sigma_;
	std::vector<double> edges_, mass_, density_, velocity_, pressure_, energy_, acc_, volume_, viscosity_;
	const Hllc hllc_;
public:
	/*! \brief Class constructor
	\param edges The initial location of the interfaces
	\param density The initial values of the density inside the cells
	\param pressure The initial values of the pressure inside the cells
	\param velocity The initial values of the velocity of the interfaces/cells, depends on run type
	\param courant The class that calculates the time step size
	\param gamma The adiabatic index (assumes ideal gas law)
	\param geo The geometry of the problem
	\param left The left side boundary condition
	\param right The right side boundary condition
	\param godunov Flag if run should use artificial viscosity or Godunov scheme
	*/
	HydroSim(std::vector<double> const& edges, std::vector<double> const& density, std::vector<double> const& pressure, std::vector<double> const& velocity, 
		Courant const& courant, double gamma, Geometry const& geo, Boundary const& left,
		Boundary const& right, bool godunov = false);
	//! \brief Advances the simulation a single time step using artifical viscosity
	void TimeAdvanceViscosity();
	//! \brief Advances the simulation a single time step using a godunov scheme
	void TimeAdvanceGodunov();
	//! \brief Advances the simulation a single time using the scheme set in the constructor
	void TimeAdvance();
	/*! \brief Outputs the simulation data into ascii files
	\param prefix The prefix to add to all output files
	\param suffix The suffix to add to all output files
	*/
	void Output(std::string const& prefix, std::string const& suffix)const;
	/*!
	\brief Returns the time of the simulation
	\return The time of the simulation
	*/
	double GetTime()const;
	/*!
	\brief Returns the cycle number of the simulation
	\return The the cycle number of the simulation
	*/
	size_t GetCycle()const;
};

#endif // HYDROSIM_HPP

