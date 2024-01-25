#pragma once

//--gtest_catch_exceptions=0 --gtest_filter=Quasistationary.EulerWithMOC_step_inter

#include "transport_moc_solver.h"
#include "newton_solver.h"
#include "euler_solver.h"

/// @brief тип данных для лямбды-функции в методе эйлера
typedef std::function<double(size_t& index)> diff_function_t;

double pars_time_line(std::string data)
{
    std::string ds{ data };
    std::tm tm{};
    std::istringstream iss(ds);
    iss >> std::get_time(&tm, "%d.%m.%Y %H:%M:%S"); // второй аргумент это маска
    if (iss.fail())
        std::cerr << "Parse failed\n";
    else
    {
        int day = tm.tm_mday;
        int minutes = tm.tm_min;
        int hours = tm.tm_hour;
        int seconds = tm.tm_sec;
        
        return seconds + minutes * 60 + hours * 3600 + (day - 1) * 3600 * 12;
    }
}

vector<vector<double>> pars_coord_heights(std::string filename = "coord_heights.csv")
{
    std::ifstream file("data\\" + filename);
    std::string line;
    char delimiter = ';';
    vector<double> coords;
    vector<double> heights;

    getline(file, line);
    while (getline(file, line))
    {
        std::stringstream stream(line);
        std::string coord;
        std::string height;
        
        getline(stream, coord, delimiter);
        getline(stream, height, delimiter);
        coords.push_back(stod(coord));
        heights.push_back(stod(height));
    }
}

// Начало моделирования - 01.08.2021 00:00:00
vector<vector<double>> parser_parameter(std::string filename)
{
    std::ifstream file("data\\" + filename);
    std::string line;
    char delimiter = ';';

    vector<double> vals;
    vector<double> times;

    while (getline(file, line))
    {
        std::replace(line.begin(), line.end(), ',', '.');
        std::stringstream stream(line);
        std::string val;
        std::string moment;

        getline(stream, moment, delimiter);

        if (moment.find(".08.2021") != std::string::npos)
        {
            getline(stream, val, delimiter);
            vals.push_back(stod(val));
            times.push_back(pars_time_line(moment));
        }

        
    }
    return { times, vals };
}

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
        size_t end_index = direction < 0 ? 0 : (press_prof.size()-1);
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
                double dz = pipe.profile.heights[index + direction] - pipe.profile.heights[index]; // Расчёт перепада высот
                double dx = pipe.profile.coordinates[index + direction] - pipe.profile.coordinates[index]; // Расчёт шага по координате
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
    size_t T = 1.2e5;
    // Давление в начале участка трубопровода
    double p0 = 6e6;
    // Расход потока
    double flow = 0.5;


    /// @brief Конструктор, инициализирующий параметры трубы и нефти 
    virtual void SetUp() override
    {
        double x0 = 0;
        double xl = 1e5;
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
    double dt_par = 5000;
    // Временные ряды для параметров
    vector<double> rho = { 860, 850, 870, 830, 860, 850 };
    vector<double> visc = { 15e-6, 12e-6, 14e-6, 11e-6, 10e-6, 13e-6 };
    vector<double> p_l = { 60e5, 61e5, 60.5e5, 59.8e5, 59e5, 60e5 };
    vector<double> Q = { 0.5, 0.7, 0.65, 0.58, 0.52, 0.47 };
    // Создадим сущность, хранящую временные ряды
    parameters_series_t parameters;
    parameters.input_dens_visc(dt_par, rho, visc);
    parameters.input_parameters(dt_par, { p_l, Q });

    // Буффер с проблемно-ориентированными слоями
    ring_buffer_t<density_viscosity_layer> buffer(2, pipe.profile.getPointCount());
    // Инициализация начальной плотности и вязкости в трубе
    auto& rho_initial = buffer.previous().density;
    auto& viscosity_initial = buffer.previous().viscosity;
    rho_initial = vector<double>(rho_initial.size(), interpolation(0, parameters.density_series, parameters.density_left_index));
    viscosity_initial = vector<double>(viscosity_initial.size(), interpolation(0, parameters.viscosity_series, parameters.viscosity_left_index));

    double modeling_time = 0;
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];

    // Получаем начальный профиль давления в трубе
    vector<double> press_profile(pipe.profile.getPointCount());
    // Профиль для начального профиля давления
    vector<double> initial_press_profile(pipe.profile.getPointCount());
    // Создаём вектор изменения давлений относитльно начального профиля
    vector<double> diff_prof(press_profile.size(), 0);


    while (modeling_time <= T)
    {

        // Вектор скорости потока по координате
        double Q_n = interpolation(modeling_time, parameters.param_series[1], parameters.params_left_index[1]);
        vector<double> Q(pipe.profile.getPointCount(), Q_n);

        density_viscosity_layer& prev = buffer.previous();
        density_viscosity_layer& next = buffer.current();

        transport_moc_solver moc_solv(pipe, Q, prev, next);


        if (modeling_time == 0)
        {
            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, prev.density, prev.viscosity, Q_n);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0]);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            initial_press_profile = press_profile;

            write_profiles_problem(prev, dx, modeling_time);
        }
        else
        {
            array<double, 2> par_in = moc_solv.get_par_in(modeling_time, parameters);
            moc_solv.step(par_in);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0]);

            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, next.density, next.viscosity, Q_n);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            // получение вектора изменения давлений относитльно начального профиля
            std::transform(initial_press_profile.begin(), initial_press_profile.end(), press_profile.begin(), diff_prof.begin(),
                [](double initial, double current) {return initial - current; });

            buffer.advance(+1);
        }

        // Вывод в файлы
        write_profiles_problem(prev, dx, modeling_time);
        write_press_profile_only(press_profile, dx, modeling_time);
        write_press_profile_only(diff_prof, dx, modeling_time, "output/diff_press_profile.csv");

        modeling_time += moc_solv.prepare_step();
    }
}

TEST_F(Quasistationary, EulerWithMOC_step_inter)
{
    // "step"

};

TEST_F(Quasistationary, NewtonWithMOC_step_inter)
{

}

TEST_F(Quasistationary, Test_euler)
{
    
}

TEST(ReadingCSV, TEST)
{
    std::string data = "02.08.2021 03:02:01";
    double moment = pars_time_line(data);

}