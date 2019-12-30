#ifndef HLLC_HPP
#define HLLC_HPP 1
#include <vector>
class Hllc
{
public:
	void CalcPstarUstar(std::vector<double> const& density, std::vector<double> const& pressure, std::vector<double> const& velocity,
		std::vector<double> const& cs, std::vector<double>& Pstar, std::vector<double>& Ustar)const;

	void CalcPstarUstarSingle(double dl, double pl, double vl, double csl, double dr, double pr,
		double vr, double csr, double& Ustar, double& Pstar)const;
};

#endif // HLLC_HPP