#pragma once

/// @brief Уравнение трубы для задачи PQ
class Pipe_model_for_PQ_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    vector<double>& density;
    vector<double>& viscosity;
    pipe_properties_t& pipe;
    oil_parameters_t& oil;
    double flow;

public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    Pipe_model_for_PQ_t(pipe_properties_t& pipe, oil_parameters_t& oil, 
        vector<double>& density, vector<double>& viscosity, double flow)
        : pipe(pipe)
        , oil(oil)
        , viscosity(viscosity)
        , density(density)
        , flow(flow)
    {
    }

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {
        double rho = density[grid_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / viscosity[grid_index];
        double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);
        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
        grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

        double height_derivative = (pipe.profile.heights[grid_index] - pipe.profile.heights[grid_index - 1]) /
            (pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);

        return { ((-4) / pipe.wall.diameter) * tau_w - rho * M_G * height_derivative };
    }

};