#pragma once

#include "transport_moc_solver.h"
#include "newton_solver.h"

/// @brief тип данных для лямбды-функции в методе эйлера
typedef std::function<double(size_t& index)> diff_function_t;
/// @brief тип данных для хранения слоёв
typedef profile_collection_t<3> layer_t;

/// @brief Функция для записи только профилей давления 
/// в разные моменты времени
/// @param press ссылка на профиль давления
/// @param dx шаг по координате
/// @param dt шаг по времени
/// @param step текущий шаг моделирования
/// @param filename название файла для записи
void write_press_profile_only(vector<double>& press, double& dx, double& time_moment, std::string filename = "output/press_prof.csv") {
    std::ofstream press_file;
    size_t profCount = press.size();
    if (time_moment == 0)
    {
        press_file.open(filename);
        press_file << "time,x,Давление" << std::endl;
    }
    else
        press_file.open(filename, std::ios::app);

    for (int i = 0; i < profCount; i++)
    {
        press_file << time_moment << "," << i * dx;
        press_file << "," << press[i] << std::endl;
    }

    press_file.close();
};

void uni_write(double& dx, double& time_moment, vector<vector<double>> params, std::string params_name, std::string filename = "output/all_prof.csv") {
    std::ofstream output_file;
    size_t profCount = params.size();
    if (time_moment == 0)
    {
        output_file.open(filename);
        output_file << params_name << std::endl;
    }
    else
        output_file.open(filename, std::ios::app);

    for (int i = 0; i < params[0].size(); i++)
    {
        output_file << time_moment << "," << i * dx;
        for (int index = 0; index < profCount; index++)
        {
            output_file << params[index][i];
            if (index == (profCount - 1))
                output_file << std::endl;
        }
    }

    output_file.close();
};


void write_profiles_problem(
    density_viscosity_layer& layer, double& dx, double& time_moment,
    std::string filename = "output/profiles_problems.csv")
{
    std::ofstream file;
    if (time_moment == 0)
    {
        file.open(filename);
        file << "time,x,Плотность,Вязкость" << std::endl;
    }
    else
        file.open(filename, std::ios::app);
    
    for (int i = 0; i < layer.density.size(); i++)
    {
        file << time_moment << "," << i * dx;
        file << "," << layer.density[i] << "," << layer.viscosity[i] << std::endl;
    }
    file.close();
}

/// @brief Класс для решения задач по квазистационару 
class Quasistationary : public ::testing::Test
{
public:

    /// @brief Алгоритм решения методом Эйлера
    /// @param press_prof Ссылка на профиль давления
    /// @param right_part Производная в точке, умноженная на шаг по координате 
    /// @param direction Направление расчёта давления
    void QP_Euler_solver(vector<double>& press_prof, const diff_function_t& right_part, size_t& direction)
    {
        size_t start_index = direction > 0 ? 1 : (press_prof.size()) - 2;
        size_t end_index = direction < 0 ? 0 : (press_prof.size());
        for (size_t index = start_index; index != end_index; index += direction)
        {
            size_t prev_index = index - direction;
            press_prof[index] = press_prof[prev_index] + direction * right_part(prev_index);
        }
    }

    /// @brief Задание функции производной и запуск алгоритма
    /// @param layer Ссылка на текущий слой
    /// @param press_prof Ссылка на профиль давления
    /// @param speed Скорость потока
    /// @param direction Направление расчёта давления
    void euler_solve(const density_viscosity_layer& layer, vector<double>& press_prof, double speed, size_t direction = 1)
    {
        // Профиль плотности, для учёта при рисчёте движения партий 
        const vector<double>& density = layer.density;
        // Профиль вязкости, для учёта при рисчёте движения партий 
        const vector<double>& viscosity = layer.density;

        // функция производной
        // возвращает значение производной в точке, умноженное на шаг по координате
        diff_function_t right_part =
            [this, speed, density, viscosity, direction](size_t& index)
            {
                double eps = pipe.wall.relativeRoughness(); // Расчёт относительной шероховатости
                double Re = speed * pipe.wall.diameter / viscosity[index]; // Расчёт числа Рейнольдса
                double lambda = pipe.resistance_function(Re, eps); // Расчёт коэффициента лямбда
                double dz = pipe.profile.heights[index] - pipe.profile.heights[index + direction]; // Расчёт перепада высот
                double dx = pipe.profile.coordinates[index] - pipe.profile.coordinates[index + direction]; // Расчёт шага по координате
                // Расчёт производной
                double diff = lambda * (1 / pipe.wall.diameter) * density[index] * pow(speed, 2) / 2 - dz / dx * density[index] * M_G;
                return dx * diff;
            };

        QP_Euler_solver(press_prof, right_part, direction);
    }

    /// @brief Решает задачу PP методом Ньютона поверх Эйлера
    /// @return Возвращает расход 
    double solve_newton_euler(const density_viscosity_layer& layer, double& p_n, double& p_L)
    {
        residual_func_t res_fun =
            [this, layer, p_n, p_L](const double& v)
            {
                // Количество точек профиля
                size_t dots_count = pipe.profile.getPointCount();
                // Профиль давлений
                vector<double> press_profile(dots_count, p_n);
                // Решение задачи QP методом Эйлера
                euler_solve(layer, press_profile, v);
                // Функция невязок
                return press_profile.back() - p_L;
            };

        solver_Newton Newton_solver(res_fun);
        double Q = Newton_solver.solve() * M_PI * pow(pipe.wall.diameter, 2) / 4;

        return Q;
    }


protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe; 
    /// @brief Параметры нефти
    oil_parameters_t oil;
    // Время моделирования
    size_t T = 700;
    // Давление в начале участка трубопровода
    double p0 = 6e6;
    // Расход потока
    double flow = 0.5;


    /// @brief Конструктор, инициализирующий параметры трубы и нефти 
    virtual void SetUp() override
    {
        double x0 = 0;
        double xl = 500;
        double D = 0.720;
        double thickness = 0.010;
        double delta = 15e-6;
        double z0 = 100;
        double zl = 50;
        double ro = 900;
        double visc = 15e-6;
        double p_capacity = 10e6;
        size_t n = 100;

        pipe.profile = PipeProfile::create(n, x0, xl, z0, zl, p_capacity);
        pipe.wall.wallThickness = thickness;
        pipe.wall.diameter = D - 2 * pipe.wall.wallThickness;
        pipe.wall.equivalent_roughness = delta;

        oil.density.nominal_density = ro;
        oil.viscosity.nominal_viscosity = visc;
    }
};

TEST_F(Quasistationary, EulerWithMOC_line_inter)
{
    // Зададим постоянный шаг для всех временных рядов
    double dt_par = 50;
    // Временные ряды для параметров
    vector<double> rho = { 900, 850, 870, 830, 860, 850 };
    vector<double> visc = { 15e-6, 12e-6, 14e-6, 11e-6, 10e-6, 13e-6 };
    vector<double> p_l = { 60.5e5, 61e5, 60.5e5, 59.8e5, 59e5, 60e5 };
    vector<double> Q = { 0.6, 0.7, 0.65, 0.58, 0.52, 0.47 };
    // Создадим сущность, хранящую временные ряды
    parameters_series_t parameters;
    parameters.input_dens_visc(dt_par, rho, visc);
    parameters.input_parameters(dt_par, { p_l, Q });

    // Буффер с проблемно-ориентированными слоями
    ring_buffer_t<density_viscosity_layer> buffer(2, pipe.profile.getPointCount());
    // Инициализация начальной плотности и вязкости в трубе
    auto& rho_initial = buffer.previous().density;
    auto& viscosity_initial = buffer.previous().viscosity;
    rho_initial = vector<double>(rho_initial.size(), oil.density.nominal_density);
    viscosity_initial = vector<double>(viscosity_initial.size(), oil.viscosity.nominal_viscosity);

    double modeling_time = 0;
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
    
    // Получаем начальный профиль давления в трубе
    // Возможно стоит включить профиль давления в проблемно-ориентированный слой
    vector<double> press_profile(pipe.profile.getPointCount(), p0);
    double speed = parameters.param_series[1][0][PAR_INDEX] / pipe.wall.getArea();
    euler_solve(buffer.previous(), press_profile, speed);
    
    // Создаём вектор изменения давлений относитльно начального профиля
    vector<double> diff_prof(press_profile.size(), 0);
    vector<double> initial_press_profile = press_profile;

    // Вывод в файлы начального состояния трубопровода
    write_profiles_problem(buffer.previous(), dx, modeling_time);
    write_press_profile_only(press_profile, dx, modeling_time);
    write_press_profile_only(diff_prof, dx, modeling_time, "output/diff_press_profile.csv");

    while (modeling_time <= T) 
    {

        density_viscosity_layer& prev = buffer.previous();
        density_viscosity_layer& next = buffer.current();

        // Вектор скорости потока по координате
        double Q_n = transport_moc_solver::interpolation(modeling_time, parameters.param_series[1]);
        vector<double> Q(pipe.profile.getPointCount(), Q_n);

        transport_moc_solver moc_solv(pipe, Q, prev, next);
        modeling_time += moc_solv.prepare_step(); // получим шаг dt для Cr = 1
        // Создаём массив краевых условий для метода характеристик
        array<double, 2> par_in = moc_solv.get_par_in(modeling_time, parameters.density_series, parameters.viscosity_series);
        moc_solv.step(par_in);
        
        // Решение задачи QP методом Эйлера
        double p_n = transport_moc_solver::interpolation(modeling_time, parameters.param_series[0]);
        press_profile[0] = p_n;
        euler_solve(next, press_profile, speed);
        
        // получение вектора изменения давлений относитльно начального профиля
        std::transform(initial_press_profile.begin(), initial_press_profile.end(), press_profile.begin(), diff_prof.begin(), 
            [](double initial, double current) {return initial - current; });

        // Вывод в файлы
        write_profiles_problem(next, dx, modeling_time);
        write_press_profile_only(press_profile, dx, modeling_time);
        write_press_profile_only(diff_prof, dx, modeling_time, "output/diff_press_profile.csv");

        buffer.advance(+1);
        
    }

}

TEST_F(Quasistationary, EulerWithMOC_step_inter)
{
    // Зададим постоянный шаг для всех временных рядов
    double dt_par = 50;
    // Временные ряды для параметров
    vector<double> rho = { 800, 850, 870, 830, 860, 850 };
    vector<double> visc = { 10e-6, 12e-6, 14e-6, 11e-6, 14e-6, 13e-6 };
    vector<double> p_l = { 60.5e5, 61e5, 60.5e5, 59.8e5, 59e5, 60e5 };
    vector<double> Q = { 0.6, 0.7, 0.65, 0.58, 0.52, 0.47 };
    // Создадим сущностьЮ хранящую временные ряды
    parameters_series_t parameters;
    parameters.input_dens_visc(dt_par, rho, visc);
    parameters.input_parameters(dt_par, { p_l, Q });

    // Буффер с проблемно-ориентированными слоями
    ring_buffer_t<density_viscosity_layer> buffer(2, pipe.profile.getPointCount());
    // Инициализация начальной плотности и вязкости в трубе
    auto& rho_initial = buffer.previous().density;
    auto& viscosity_initial = buffer.previous().viscosity;
    rho_initial = vector<double>(rho_initial.size(), oil.density.nominal_density);
    viscosity_initial = vector<double>(viscosity_initial.size(), oil.viscosity.nominal_viscosity);

    double modeling_time = 0;
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];

    // Получаем начальный профиль давления в трубе
    // Возможно стоит включить профиль давления в проблемно-ориентированный слой
    vector<double> press_profile(pipe.profile.getPointCount(), p0);
    double speed = flow / pipe.wall.getArea();
    euler_solve(buffer.previous(), press_profile, speed);

    // Создаём вектор изменения давлений относитльно начального профиля
    vector<double> diff_prof(press_profile.size(), 0);
    vector<double> initial_press_profile = press_profile;

    // Вывод в файлы начального состояния трубопровода
    write_profiles_problem(buffer.previous(), dx, modeling_time);
    write_press_profile_only(press_profile, dx, modeling_time);
    write_press_profile_only(diff_prof, dx, modeling_time, "output/diff_press_profile.csv");

    while (modeling_time <= T)
    {

        density_viscosity_layer& prev = buffer.previous();
        density_viscosity_layer& next = buffer.current();

        // Вектор скорости потока по координате
        double Q_n = transport_moc_solver::interpolation(modeling_time, parameters.param_series[1], "step");
        vector<double> Q(pipe.profile.getPointCount(), Q_n);

        transport_moc_solver moc_solv(pipe, Q, prev, next);
        modeling_time += moc_solv.prepare_step(); // получим шаг dt для Cr = 1
        // Создаём массив краевых условий для метода характеристик
        array<double, 2> par_in = moc_solv.get_par_in(modeling_time, parameters.density_series, parameters.viscosity_series, "step");
        moc_solv.step(par_in);

        // Решение задачи QP методом Эйлера
        double p_n = transport_moc_solver::interpolation(modeling_time, parameters.param_series[0], "step");
        press_profile[0] = p_n;
        euler_solve(next, press_profile, speed);

        // получение вектора изменения давлений относитльно начального профиля
        std::transform(initial_press_profile.begin(), initial_press_profile.end(), press_profile.begin(), diff_prof.begin(),
            [](double initial, double current) {return initial - current; });

        // Вывод в файлы
        write_profiles_problem(next, dx, modeling_time);
        write_press_profile_only(press_profile, dx, modeling_time);
        write_press_profile_only(diff_prof, dx, modeling_time, "output/diff_press_profile.csv");

        buffer.advance(+1);

    }

}

TEST_F(Quasistationary, NewtonWithMOC_step_inter)
{
    // Зададим постоянный шаг для всех временных рядов
    double dt_par = 50;
    // Временные ряды для параметров
    vector<double> rho = { 800, 850, 870, 830, 860, 850 };
    vector<double> visc = { 10e-6, 12e-6, 14e-6, 11e-6, 14e-6, 13e-6 };
    vector<double> p_n = { 60.5e5, 61e5, 60.5e5, 59.8e5, 59e5, 60e5 };
    vector<double> p_L = { 52.5e5, 54.6e5, 53.5e5, 55.8e5, 54e5, 52e5 };

    // Создадим сущность, хранящую временные ряды
    parameters_series_t parameters;
    parameters.input_dens_visc(dt_par, rho, visc);
    parameters.input_parameters(dt_par, { p_n, p_L });

    // Буффер с проблемно-ориентированными слоями
    ring_buffer_t<density_viscosity_layer> buffer(2, pipe.profile.getPointCount());
    // Инициализация начальной плотности и вязкости в трубе
    auto& rho_initial = buffer.previous().density;
    auto& viscosity_initial = buffer.previous().viscosity;
    rho_initial = vector<double>(rho_initial.size(), oil.density.nominal_density);
    viscosity_initial = vector<double>(viscosity_initial.size(), oil.viscosity.nominal_viscosity);

    double modeling_time = 0;
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
    double pn = transport_moc_solver::interpolation(modeling_time, parameters.param_series[0]);
    double pL = transport_moc_solver::interpolation(modeling_time, parameters.param_series[1]);
    vector<double> Q_new = vector<double>(solve_newton_euler(buffer.previous(), pn, pL));

    while (modeling_time <= T)
    {

        density_viscosity_layer& prev = buffer.previous();
        density_viscosity_layer& next = buffer.current();


        transport_moc_solver moc_solv(pipe, Q_new, prev, next);
        modeling_time += moc_solv.prepare_step(); // получим шаг dt для Cr = 1
        // Создаём массив краевых условий для метода характеристик
        array<double, 2> par_in = moc_solv.get_par_in(modeling_time, parameters.density_series, parameters.viscosity_series, "step");
        moc_solv.step(par_in);

        pn = transport_moc_solver::interpolation(modeling_time, parameters.param_series[0]);
        pL = transport_moc_solver::interpolation(modeling_time, parameters.param_series[1]);

        Q_new = vector<double>(solve_newton_euler(next, pn, pL));

        buffer.advance(+1);

    }

}