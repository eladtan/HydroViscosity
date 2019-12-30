#include "Courant.hpp"
#include <algorithm>
#include <cmath>

Courant::Courant(double cfl) :cfl_(cfl) {}

double Courant::operator()(std::vector<double> const& edges, std::vector<double> const& pressure, std::vector<double> const& density, double gamma) const
{
	double dt_inv = 0.0;
	size_t N = density.size();
	for (size_t i = 0; i < N; ++i)
	{
		double cs = std::sqrt(pressure[i + 1] * gamma / density[i]);
		dt_inv = std::max(dt_inv, cs / (edges[i + 1] - edges[i]));
	}
	return cfl_ / dt_inv;
}

double Courant::operator()(std::vector<double> const& edges, std::vector<double> const& cs) const
{
	double dt_inv = 0.0;
	size_t N = cs.size();
	for (size_t i = 0; i < N; ++i)
		dt_inv = std::max(dt_inv, cs[i] / (edges[i + 1] - edges[i]));
	return cfl_ / dt_inv;
}
