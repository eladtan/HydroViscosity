#include "utils.hpp"
#include <iostream>
#include <fstream>

void write_vector(std::vector<double> const& v, std::string const& fname, size_t prec)
{
	std::ofstream f(fname.c_str());
	f.precision(prec);
	for (size_t i = 0; i < v.size(); ++i)
		f << v[i] << "\n";
	f.close();
}