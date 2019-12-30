#include "Hllc.hpp"
#include <algorithm>

namespace
{
	void GetWaveSpeedsVector(std::vector<double> const& density, std::vector<double> const& pressure, std::vector<double> const& velocity,
		std::vector<double> const& cs, std::vector<double> &sl, std::vector<double>& sr, std::vector<double>& sc)
	{
		size_t N = density.size() - 1;
		sl.resize(N);
		sr.resize(N);
		for (size_t i = 0; i < N; ++i)
		{
			sl[i] = std::min(velocity[i] - cs[i], velocity[i + 1] - cs[i + 1]);
			sr[i] = std::max(velocity[i] + cs[i], velocity[i + 1] + cs[i + 1]);
			sc[i] = (pressure[i + 1] - pressure[i] + density[i] * velocity[i] *
				(sl[i] - velocity[i]) - density[i + 1] * velocity[i + 1] * (sr[i] -
					velocity[i + 1])) / (density[i] * (sl[i] - velocity[i]) - density[i + 1]
						* (sr[i] - velocity[i + 1]));
		}
	}
}

void Hllc::CalcPstarUstar(std::vector<double> const& density, std::vector<double> const& pressure, std::vector<double> const& velocity,
	std::vector<double> const& cs, std::vector<double>& Pstar, std::vector<double>& Ustar)const
{
	std::vector<double> sl, sr;
	GetWaveSpeedsVector(density, pressure, velocity, cs, sl, sr, Ustar);
	Pstar.resize(Ustar.size());
	size_t N = Pstar.size();
	for (size_t i = 0; i < N; ++i)
		Pstar[i] = pressure[i] + density[i] * (sl[i] - velocity[i]) * (Ustar[i] - velocity[i]);
}

void Hllc::CalcPstarUstarSingle(double dl, double pl, double vl, double csl, double dr, double pr, double vr, double csr, double& Ustar, double& Pstar) const
{
	double sl = std::min(vl - csl, vr - csr);
	double sr = std::max(vl + csl, vr + csr);
	Ustar = (pr - pl + dl * vl * (sl - vl) - dr * vr * (sr - vr)) / (dl * (sl - vl) - dr * (sr - vr));
	Pstar = pl + dl * (sl - vl) * (Ustar - vl);
}
