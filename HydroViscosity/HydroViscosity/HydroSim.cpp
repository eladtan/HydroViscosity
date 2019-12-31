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

	void GetAllAreas(Geometry const& geo, std::vector<double> const& edges, std::vector<double>& res)
	{
		size_t N = edges.size();
		res.resize(N);
		for (size_t i = 0; i < N; ++i)
			res[i] = geo.CalcArea(edges, i);
	}

	void CalcDensity(std::vector<double> const& mass, std::vector<double> const& volumes, std::vector<double>& density)
	{
		size_t N = density.size();
		for (size_t i = 0; i < N; ++i)
			density[i] = mass[i] / volumes[i];
	}

	void CalcViscosity(std::vector<double> const& velocity, std::vector<double> const& density, Boundary const& left, Boundary const& right,
		double viscosity_sigma, std::vector<double> &viscosity)
	{
		size_t N = density.size();
		for (size_t i = 0; i < N; ++i)
			if (velocity[i + 1] < velocity[i])
				viscosity[i] = viscosity_sigma * (velocity[i + 1] - velocity[i]) * (velocity[i + 1] - velocity[i]) * density[i];
	}

	void UpdateInternalEnergy(std::vector<double> &energy, std::vector<double> const& pressure, std::vector<double> const& viscosity,
		std::vector<double> const& viscositynew, std::vector<double> const& volumesnew, std::vector<double> const& volume, std::vector<double> const& density,
		std::vector<double> const& mass, double gamma)
	{
		size_t N = energy.size();
		for (size_t i = 0; i < N; ++i)
			energy[i] = (energy[i] - 0.5 * (pressure[i] + viscosity[i] + viscositynew[i]) * (volumesnew[i] - volume[i]) / mass[i]) / (1.0 + 0.5 *
			(gamma - 1.0) * density[i] * (volumesnew[i] - volume[i]) / mass[i]);
	}

	void CalcAcceleration(std::vector<double> const& pressure, std::vector<double> const& viscosity, std::vector<double> const& mass, std::vector<double> &acc,
		Boundary const& left, Boundary const& right)
	{
		size_t N = acc.size() - 1;
		for (size_t i = 1; i < N; ++i)
			acc[i] = -2.0 * (pressure[i] - pressure[i - 1] + viscosity[i] - viscosity[i - 1]) / (mass[i] + mass[i - 1]);
		// Deal with boundaries
		acc[0] = -2.0 * (pressure[0] - left.GetSidePressure(true, pressure) + viscosity[0] - left.GetSideViscosity(true, viscosity)) / 
			(mass[0] + left.GetSideMass(true, mass));
		acc.back() = -2.0 * (right.GetSidePressure(false, pressure) + right.GetSideViscosity(false, viscosity) - pressure.back() - viscosity.back()) / 
			(mass.back() + right.GetSideMass(false, mass));
	}

	void de2p(std::vector<double> const& density, std::vector<double> const& energy, double gamma, std::vector<double>& pressure)
	{
		size_t N = density.size();
		for (size_t i = 0; i < N; ++i)
			pressure[i] = (gamma - 1.0) * density[i] * energy[i];
	}

	void dp2c(std::vector<double> const& density, std::vector<double> const& pressure, double gamma,
		std::vector<double>& cs)
	{
		size_t N = density.size();
		cs.resize(N);
		for (size_t i = 0; i < N; ++i)
			cs[i] = std::sqrt(gamma * pressure[i + 1] / density[i]);
	}

	double dp2csingle(double density, double pressure, double gamma)
	{
		return std::sqrt(gamma * pressure / density);
	}

	void BoundaryRiemannSolve(Boundary const& boundary_left, Boundary const& boundary_right, std::vector<double> const& density, std::vector<double> const& pressure,
		std::vector<double> const& velocity, std::vector<double> const& cs, double gamma, Hllc const& hllc, double &pstarl, double &pstarr, double &ustarl, double &ustarr)
	{
		double dl = boundary_left.GetSideDensity(true, density);
		double pl = boundary_left.GetSidePressure(true, pressure);
		double vl = boundary_left.GetSideVelocity(true, velocity);
		double csl = dp2csingle(dl, pl, gamma);
		double dr = boundary_right.GetSideDensity(false, density);
		double pr = boundary_right.GetSidePressure(false, pressure);
		double vr = boundary_right.GetSideVelocity(false, velocity);
		double csr = dp2csingle(dr, pr, gamma);
		hllc.CalcPstarUstarSingle(dl, pl, vl, csl, density[0], pressure[0], velocity[0], cs[0], ustarl, pstarl);
		hllc.CalcPstarUstarSingle(density.back(), pressure.back(), velocity.back(), cs.back(), dr, pr, vr, csr, ustarr, pstarr);
	}

	void UpdateExtensive(std::vector<double>& momentum, std::vector<double>& energy, std::vector<double> const& pstar, std::vector<double> const& ustar,
		std::vector<double> const& areas, double pstarl, double pstarr, double ustarl, double ustarr, double dt)
	{
		size_t N = momentum.size() - 2;
		for (size_t i = 0; i < N; ++i)
		{
			momentum[i + 1] += dt * (pstar[i + 1] * areas[i + 2] - pstar[i] * areas[i + 1]);
			energy[i + 1] += dt * (pstar[i + 1] * areas[i + 2] * ustar[i + 1] - pstar[i] * areas[i + 1] * ustar[i]);
		}
		// Deal with boundaries
		momentum[0] += dt * (pstar[0] * areas[1] - pstarl * areas[0]);
		momentum.back() += dt * (pstarr * areas.back()  - pstar.back() * areas[N + 1]);
		energy[0] += dt * (pstar[0] * areas[1] * ustar[0] - pstarl * areas[0] * ustarl);
		energy.back() += dt * (pstarr * areas.back() * ustarr - pstar.back() * areas[N + 1] * ustar.back());
	}

	void UpdateMesh(std::vector<double>& edges, std::vector<double> const& ustar, double ul, double ur, double dt)
	{
		size_t N = edges.size() - 1;
		for (size_t i = 1; i < N; ++i)
			edges[i] += ustar[i - 1] * dt;
		edges[0] += ul * dt;
		edges.back() += ur * dt;
	}
}

HydroSim::HydroSim(std::vector<double> const& edges, std::vector<double> const& density, std::vector<double> const& pressure, std::vector<double> const& velocity,
	Courant const& courant, double gamma, Geometry const& geo,Boundary const& left,
	Boundary const& right, bool godunov) : edges_(edges), density_(density), velocity_(velocity), pressure_(pressure), courant_(courant), gamma_(gamma),
	geo_(geo), boundary_left_(left), boundary_right_(right), godunov_(godunov), time_(0.0), cycle_(0), viscosity_sigma_(2.0), mass_(std::vector<double>()), 
	viscosity_(std::vector<double>()), volume_(std::vector<double>()), energy_(std::vector<double>()), hllc_(Hllc())
{
	assert(edges_.size() > 1);
	size_t N = edges_.size() - 1;
	// Allocate the data
	mass_.resize(N);
	energy_.resize(N);

	// Set the initial mass and energy
	for (size_t i = 0; i < N; ++i)
	{
		mass_[i] = geo_.CalcVolume(edges_, i) * density[i];
		energy_[i] = pressure_[i] / ((gamma_ - 1.0) * density[i]);
	}
	// Calc artificial viscosity
	viscosity_.resize(N);
	CalcViscosity(velocity_, density_, boundary_left_, boundary_right_, viscosity_sigma_, viscosity_);
		
	// Calculate initial acceleration
	acc_.resize(N + 1);
	CalcAcceleration(pressure_, viscosity_, mass_, acc_, boundary_left_, boundary_right_);

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
	// Calculate the acceleration
	CalcAcceleration(pressure_, viscositynew, mass_, acc_, boundary_left_, boundary_right_);
	// Kick interfaces again
	velocity_ += (0.5 * dt) * acc_;

	viscosity_ = viscositynew;
	volume_ = volumesnew;

	++cycle_;
	time_ += dt;
}

void HydroSim::TimeAdvanceGodunov()
{
	std::vector<double> momentum = velocity_ * mass_;
	// Calculate the sound speed
	std::vector<double> cs;
	dp2c(density_, pressure_, gamma_, cs);
	// Calcualte the time step
	double dt = courant_(edges_, pressure_, density_, gamma_);
	// Solve the reimann problem for the main domain
	std::vector<double> Pstar, Ustar;
	hllc_.CalcPstarUstar(density_, pressure_, velocity_, cs, Pstar, Ustar);
	// Solve the reimann problem for the boundaries
	double pstarl, pstarr, ustarl, ustarr;
	BoundaryRiemannSolve(boundary_left_, boundary_right_, density_, pressure_, velocity_, cs, gamma_, hllc_, pstarl, pstarr, ustarl, ustarr);
	// Update extensives
	std::vector<double> areas;
	GetAllAreas(geo_, edges_, areas);
	UpdateExtensive(momentum, energy_, Pstar, Ustar, areas, pstarl, pstarr, ustarl, ustarr, dt);
	// Move grid
	UpdateMesh(edges_, Ustar, ustarl, ustarr, dt);
	// Update primitives
	GetAllVolumes(geo_, edges_, volume_);
	density_ = mass_ / volume_;
	velocity_ = momentum / mass_;
	de2p(density_, energy_ / mass_ - 0.5 * velocity_ * velocity_, gamma_, pressure_);

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
