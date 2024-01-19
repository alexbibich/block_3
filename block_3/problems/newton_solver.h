#pragma once

typedef std::function<double(const double& v)> residual_func_t;

/// @brief Класс с функцией невязок для решения методом Ньютона-Рафсона
class solver_Newton : public fixed_system_t<1>
{
public:
	solver_Newton(const residual_func_t& res_fun)
		: res_func{ res_fun }
	{
	}

	double residuals(const double& v) {
		double result = res_func(v);
		return result;
	}

	/// @brief Решение задачи РР методом Ньютона
	/// @param p0 Давление в начале трубопровода
	/// @param pl Давление в конце трубопровода
	/// @return Возвращает расход
	double solve()
	{
		fixed_solver_parameters_t<1, 0> parameters;
		// Создание структуры для записи результатов расчета
		fixed_solver_result_t<1> result;
		// Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
		// { 0, 0 } - Начальное приближение
		fixed_newton_raphson<1>::solve_dense(*this, { 1 }, parameters, &result);

		double speed = result.argument;

		return speed;
	}
protected:
	const residual_func_t& res_func;
};

/// @brief Солвер для решения задачи PP методом Ньютона-Рафсона
class PP_solver
{
public:
	PP_solver(const pipe_properties_t& pipe, const oil_parameters_t& oil, double p0, double pl)
		: pipe{ pipe }, oil{ oil }, p0{ p0 }, pl{ pl }
	{}

	

protected:
	const pipe_properties_t& pipe;
	const oil_parameters_t& oil;
	double p0;
	double pl;
};