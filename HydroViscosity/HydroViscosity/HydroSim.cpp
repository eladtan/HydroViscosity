#include "HydroSim.hpp"
#include "utils.hpp"
#include <cassert>
#include <cmath>

namespace
{
	void GetAllVolumes(Geometry const& geo, std::vector<double> const& edges, std::vector<double>& res)
	{
		size_t N = edges.size() - 1;
		res.resize(N);		
		for (size_t i = 0; i < N; ++i)
			res[i] = geo.CalcVolume(edges, i);
	}

	void CalcDensity(std::vector<double> const& mass, std::vector<double> const& volumes, std::vector<double>& density)
	{
		size_t N = density.size();
		for (size_t i = 0; i < N; ++i)
			density[i] = mass[i + 1] / volumes[i];
	}

	void CalcViscosity(std::vector<double> const& velocity, std::vector<double> const& density, Boundary const& left, Boundary const& right,
		double viscosity_sigma, std::vector<double> &viscosity)
	{
		size_t N = density.size();
		for (size_t i = 0; i < N; ++i)
			if (velocity[i + 1] < velocity[i])
				viscosity[i + 1] = viscosity_sigma * (velocity[i + 1] - velocity[i]) * (velocity[i + 1] - velocity[i]) * density[i];
		// Deal with boundaries
		viscosity[0] = left.GetSideViscosity(true, viscosity);
		viscosity.back() = right.GetSideViscosity(false, viscosity);
	}

	void UpdateInternalEnergy(std::vector<double> &energy, std::vector<double> const& pressure, std::vector<double> const& viscosity,
		std::vector<double> const& viscositynew, std::vector<double> const& volumesnew, std::vector<double> const& volume, std::vector<double> const& density,
		std::vector<double> const& mass, double gamma)
	{
		size_t N = energy.size();
		for (size_t i = 0; i < N; ++i)
			energy[i] = (energy[i] - 0.5 * (pressure[i + 1] + viscosity[i + 1] + viscositynew[i + 1]) * (volumesnew[i] - volume[i]) / mass[i + 1]) / (1.0 + 0.5 *
			(gamma - 1.0) * density[i] * (volumesnew[i] - volume[i]) / mass[i + 1]);
	}

	void CalcAcceleration(std::vector<double> const& pressure, std::vector<double> const& viscosity, std::vector<double> const& mass, std::vector<double> &acc)
	{
		size_t N = acc.size();
		for (size_t i = 0; i < N; ++i)
			acc[i] = -2.0 * (pressure[i + 1] - pressure[i] + viscosity[i + 1] - viscosity[i]) / (mass[i] + mass[i + 1]);
	}

	void de2p(std::vector<double> const& density, std::vector<double> const& energy, double gamma, std::vector<double>& pressure)
	{
		size_t N = density.size();
		for (size_t i = 0; i < N; ++i)
			pressure[i + 1] = (gamma - 1.0) * density[i] * energy[i];
	}

	void dp2c(std::vector<double> const& density, std::vector<double> const& pressure, double gamma,
		std::vector<double>& cs)
	{
		size_t N = density.size();
		cs.resize(N);
		for (size_t i = 0; i < N; ++i)
			cs[i] = std::sqrt(gamma * pressure[i + 1] / density[i]);
	}
}

HydroSim::HydroSim(std::vector<double> const& edges, std::vector<double> const& density, std::vector<double> const& pressure, std::vector<double> const& velocity,
	Courant const& courant, double gamma, Geometry const& geo,Boundary const& left,
	Boundary const& right) : edges_(edges), density_(density), velocity_(velocity), pressure_(pressure), courant_(courant), gamma_(gamma),
	geo_(geo), boundary_left_(left), boundary_right_(right), time_(0.0), cycle_(0), viscosity_sigma_(2.0), mass_(std::vector<double>()), 
	viscosity_(std::vector<double>()), volume_(std::vector<double>()), energy_(std::vector<double>())
{
	assert(edges_.size() > 1);
	size_t N = edges_.size() - 1;
	// Allocate the data
	mass_.resize(N + 2);
	energy_.resize(N);

	// Set the initial mass and energy
	for (size_t i = 0; i < N; ++i)
	{
		mass_[i + 1] = geo_.CalcVolume(edges_, i) * density[i];
		energy_[i] = pressure_[i] / ((gamma_ - 1.0) * density[i]);
	}
	mass_[0] = boundary_left_.GetSideMass(true, mass_);
	mass_.back() = boundary_right_.GetSideMass(false, mass_);

	// Increase size of pressure vector to include ghost cells
	pressure_.insert(pressure_.begin(), boundary_left_.GetSidePressure(true, pressure_));
	pressure_.push_back(boundary_right_.GetSidePressure(false, pressure_));

	// Calc artificial viscosity
	viscosity_.resize(N + 2);
	CalcViscosity(velocity_, density_, boundary_left_, boundary_right_, viscosity_sigma_, viscosity_);
		
	// Calculate initial acceleration
	acc_.resize(N + 1);
	CalcAcceleration(pressure_, viscosity_, mass_, acc_);

	// Calculate initial volume
	GetAllVolumes(geo_, edges_, volume_);
}

void HydroSim::TimeAdvanceViscosity()
{
	// Get the time step
	double dt = courant_(edges_, pressure_, density_, gamma_);
	// Kick interfaces
	velocity_ += (0.5 * dt) * acc_;
	// Move interfaces
	edges_ += dt * velocity_;
	// Calculate volumes and density
	std::vector<double> volumesnew;
	GetAllVolumes(geo_, edges_, volumesnew);
	CalcDensity(mass_, volumesnew, density_);
	// Calculate the viscosity
	std::vector<double> viscositynew(viscosity_.size());
	CalcViscosity(velocity_, density_, boundary_left_, boundary_right_, viscosity_sigma_, viscositynew);
	// Calculate new internal energy
	UpdateInternalEnergy(energy_, pressure_, viscosity_, viscositynew, volumesnew, volume_, density_, mass_, gamma_);
	// Update the pressure
	de2p(density_, energy_, gamma_, pressure_);
	pressure_[0] = boundary_left_.GetSidePressure(true, pressure_);
	pressure_.back() = boundary_right_.GetSidePressure(false, pressure_);
	// Calculate the acceleration
	CalcAcceleration(pressure_, viscositynew, mass_, acc_);
	// Kick interfaces again
	velocity_ += (0.5 * dt) * acc_;

	viscosity_ = viscositynew;
	volume_ = volumesnew;

	++cycle_;
	time_ += dt;
}

void HydroSim::TimeAdvanceGodunov()
{
	// Calculate the sound speed
	std::vector<double> cs;
	dp2c(density_, pressure_, gamma_, cs);
	// Calcualte the time step
	double dt = courant_(edges_, pressure_, density_, gamma_);
	// Solve the reimann problem for the main domain
	std::vector<double> Pstar, Ustar;
	hllc_.CalcPstarUstar(density_, pressure_, velocity_, cs, Pstar, Ustar);
	// Solve the reimann problem for the boundaries

	// Update extensives

	// Move grid

	// Update primitives

	++cycle_;
	time_ += dt;
}

void HydroSim::Output(std::string const& prefix, std::string const& suffix) const
{
	write_vector(density_, prefix + "density" + suffix + ".txt");
	write_vector(velocity_, prefix + "velocity" + suffix + ".txt");
	write_vector(pressure_, prefix + "pressure" + suffix + ".txt");
	write_vector(edges_, prefix + "x" + suffix + ".txt");
	std::vector<double> temp(1, time_);
	write_vector(temp, prefix + "time" + suffix + ".txt");
}

double HydroSim::GetTime() const
{
	return time_;
}

size_t HydroSim::GetCycle() const
{
	return cycle_;
}
