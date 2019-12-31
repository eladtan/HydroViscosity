#include "HydroSim.hpp"
#include <chrono>
#include <iostream>

namespace
{
	void init_sod(std::vector<double>& density, std::vector<double>& pressure, std::vector<double>& velocity, std::vector<double> & edges, 
		size_t resolution, bool godunov)
	{
		// allocate the data
		density.resize(resolution);
		pressure.resize(resolution);
		if(godunov)
			velocity.resize(resolution, 0.0);
		else
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
	size_t N = 128;
	std::vector<double> density, pressure, velocity, edges;
	init_sod(density, pressure, velocity, edges, N, true);
	RigidWall left, right;
	Planar geo;
	double gamma = 1.4;
	Courant courant(0.3);
	HydroSim sim(edges, density, pressure, velocity, courant, gamma, geo, left, right, true);
	auto start = std::chrono::high_resolution_clock::now();
	while (sim.GetTime() < 0.2)
		sim.TimeAdvance();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << duration.count() << std::endl;
	sim.Output("c:/sim_data/", "g");

	init_sod(density, pressure, velocity, edges, N, false);
	HydroSim sim2(edges, density, pressure, velocity, courant, gamma, geo, left, right);
	start = std::chrono::high_resolution_clock::now();
	while (sim2.GetTime() < 0.2)
		sim2.TimeAdvance();
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << duration.count() << std::endl;
	sim2.Output("c:/sim_data/", "rm");
	return 0;
}