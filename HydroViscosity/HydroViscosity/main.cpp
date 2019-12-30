#include "HydroSim.hpp"

namespace
{
	void init_sod(std::vector<double>& density, std::vector<double>& pressure, std::vector<double>& velocity, std::vector<double> & edges, size_t resolution)
	{
		// allocate the data
		density.resize(resolution);
		pressure.resize(resolution);
		velocity.resize(resolution + 1, 0.0);
		edges.resize(resolution + 1);
		edges[0] = 0.0;

		double dx = 1.0 / resolution;
		for (size_t i = 0; i < resolution; ++i)
		{
			edges[i + 1] = (i + 1) * dx;
			if (edges[i + 1] < 0.5)
			{
				pressure[i] = 1.0;
				density[i] = 1.0;
			}
			else
			{
				pressure[i] = 0.1;
				density[i] = 0.125;
			}
		}
	}
}

int main(void)
{
	size_t N = 256;
	std::vector<double> density, pressure, velocity, edges;
	init_sod(density, pressure, velocity, edges, N);
	RigidWall left, right;
	Planar geo;
	double gamma = 1.4;
	Courant courant(0.3);
	HydroSim sim(edges, density, pressure, velocity, courant, gamma, geo, left, right);
	while (sim.GetTime() < 0.2)
	{
		sim.TimeAdvanceViscosity();
	}
	sim.Output("c:/sim_data/", "");
	return 0;
}