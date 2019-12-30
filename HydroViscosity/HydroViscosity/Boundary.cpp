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
