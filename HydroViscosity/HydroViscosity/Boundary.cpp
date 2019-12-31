#include "Boundary.hpp"

double RigidWall::GetSidePressure(bool left_side, std::vector<double> const& pressure) const
{
	if (left_side)
		return pressure[1];
	else
		return pressure[pressure.size() - 2];
}

double RigidWall::GetSideMass(bool left_side, std::vector<double> const& mass) const
{
	if (left_side)
		return mass[1];
	else
		return mass[mass.size() - 2];
}

double RigidWall::GetSideViscosity(bool left_side, std::vector<double> const& viscosity) const
{
	if (left_side)
		return viscosity[1];
	else
		return viscosity[viscosity.size() - 2];
}

double RigidWall::GetSideVelocity(bool left_side, std::vector<double> const& velocity) const
{
	if (left_side)
		return velocity[1];
	else
		return velocity[velocity.size() - 2];
}

double RigidWall::GetSideDensity(bool left_side, std::vector<double> const& density) const
{
	if (left_side)
		return density[1];
	else
		return density[density.size() - 2];
}
