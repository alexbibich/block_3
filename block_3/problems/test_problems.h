#pragma once

/// @brief Инициализация данных в структурах
/// @param pipe Ссылка на структуру трубы
/// @param oil Ссылка на структуру нефти
void init_cond(pipe_properties_t& pipe, oil_parameters_t& oil)
{
    double L = 8e+4;
    double x0 = 0;
    double xl = 8e4;
    double D = 0.720;
    double thickness = 0.010;
    double delta = 15e-6;
    double z0 = 50;
    double zl = 100;
    double ro = 860;
    double visc = 15e-6;
    double p_capacity = 10e6;
    size_t n = 1000;

    pipe.profile = PipeProfile::create(n, x0, xl, z0, zl, p_capacity);
    pipe.wall.wallThickness = thickness;
    pipe.wall.diameter = D - 2 * pipe.wall.wallThickness;
    pipe.wall.equivalent_roughness = delta;

    oil.density.nominal_density = ro;
    oil.viscosity.nominal_viscosity = visc;
}

TEST(BLOCK3, EulerWithCharacteristic)
{
    // Модель трубопровода
    pipe_properties_t pipe_prop;
    // Модель нефти
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    solve_euler_corrector

}