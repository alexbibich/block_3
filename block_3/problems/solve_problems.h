#pragma once
#include <fstream>
#include <vector>
#include <cstdlib>

#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

class EulerWithMOC
{
	vector<double> QP_Euler_solver(vector<double> density, double p, double speed, bool direction = true)
	{
		double eps = pipe.wall.getCompressionRatio(); // Расчёт относительной шероховатости
		double Re = speed * pipe.wall.diameter / oil.viscosity.nominal_viscosity; // Расчёт числа Рейнольдса
		

		double dz = pipe.profile.heights[1] - pipe.profile.heights[0];
		double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
		size_t num_dots = pipe.profile.getPointCount();
		vector<double> press_prof(num_dots, p);

		for (size_t i = 1; i < num_dots; i++)
		{
			
			size_t index = direction ? i : num_dots - 1 - i;
			double lambda = pipe.resistance_function(Re, eps); // Расчёт коэффициента лямбда
			double diff = lambda * (1 / pipe.wall.diameter) * density[index] * pow(speed, 2) / 2 - dz / dx * density[index];
			if (direction)
				press_prof[index] = press_prof[index - 1] + dx * diff;
			else
				press_prof[index] = press_prof[index + 1] - dx * diff;
		}
		
		return press_prof;
	}

	EulerWithMOC(pipe_properties_t& pipe, oil_parameters_t& oil)
		:pipe{ pipe }, oil{ oil }
	{}

protected:
	pipe_properties_t& pipe;
	oil_parameters_t& oil;
};