#include "Geometry.hpp"
#define _USE_MATH_DEFINES
#include <math.h>

double Planar::CalcArea(std::vector<double> const& edges, size_t index) const
{
	return 1.0;
}

double Planar::CalcVolume(std::vector<double> const& edges, size_t index) const
{
	return edges[index + 1] - edges[index];
}

double Spherical::CalcArea(std::vector<double> const& edges, size_t index) const
{
	return 4 * M_PI * edges[index] * edges[index];
}

double Spherical::CalcVolume(std::vector<double> const& edges, size_t index) const
{
	return 4 * M_PI * (edges[index + 1] * edges[index + 1] * edges[index + 1] - edges[index] * edges[index] * edges[index]) / 3.0;
}
