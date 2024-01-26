#pragma once

//--gtest_catch_exceptions=0 --gtest_filter=Quasistationary.EulerWithMOC_step_inter

#include "transport_moc_solver.h"
#include "newton_solver.h"
#include "euler_solver.h"

/// @brief тип данных для лямбды-функции в методе эйлера
typedef std::function<double(size_t& index)> diff_function_t;

/// @brief Функция для записи только профилей давления 
/// в разные моменты времени
/// @param press ссылка на профиль давления
/// @param dx шаг по координате
/// @param dt шаг по времени
/// @param step текущий шаг моделирования
/// @param filename название файла для записи
void write_press_profile_only(vector<double>& press, double dx, double time_moment, std::string filename = "output/press_prof.csv") {
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

void uni_write(vector<double>& x, double time_moment, vector<vector<double>> params, std::string params_name, std::string filename = "output/all_prof.csv") {
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
        output_file << time_moment << "," << x[i];
        for (int index = 0; index < profCount; index++)
        {
            output_file << ',' << params[index][i];
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
        
        return seconds + minutes * 60 + hours * 3600 + (day - 1) * 3600 * 24;
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
        coords.push_back(stod(coord) * 1000);
        heights.push_back(stod(height));
    }

    return { coords, heights };
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



/// @brief Класс для решения задач по квазистационару 
class Quasistationary : public ::testing::Test
{
public:


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
        vector<vector<double>> coord_height = pars_coord_heights();

        time_series_t c_h = parameters_series_t::build_parameters(coord_height[0], coord_height[1]);

        double x0 = coord_height[0].front();
        double xl = coord_height[0].back();
        double d = 1;
        double z0 = 100;
        double zl = 50;
        double p_capacity = 10e6;
        size_t n = 100;

        pipe.profile = PipeProfile::create(n, x0, xl, z0, zl, p_capacity);
        pipe.wall.diameter = d;

        vector<double> heights_inter(pipe.profile.getPointCount());
        size_t left_index_heights = 0;
        for (size_t index = 0; index < pipe.profile.getPointCount(); index++)
        {
            double height = interpolation(pipe.profile.coordinates[index], c_h, left_index_heights);
            heights_inter[index] = height;
        }

        pipe.profile.heights = heights_inter;

        uni_write(pipe.profile.coordinates, 0, { pipe.profile.heights }, "Время,Координаты,Высотки", "output\\test_pars_pipe.csv");

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
};

TEST_F(Quasistationary, EulerWithMOC_step_inter)
{
    // "step"

};


TEST_F(Quasistationary, TEST)
{
    std::string data = "02.08.2021 03:02:01";
    double moment = pars_time_line(data);

    vector<vector<double>> rho_pars = parser_parameter("rho_in.csv");
    vector<vector<double>> visc_pars = parser_parameter("visc_in.csv");
    vector<vector<double>> p_n_pars = parser_parameter("p_in.csv");
    vector<vector<double>> Q_pars = parser_parameter("Q_in.csv"); // При большой погрешности взять средний между началом и концом
    vector<vector<double>> p_l_pars = parser_parameter("p_out.csv");
    vector<vector<double>> Q_out_pars = parser_parameter("Q_out.csv");

    std::transform(visc_pars[1].begin(), visc_pars[1].end(), visc_pars[1].begin(), [](double visc) { return visc * 1e-6; });
    std::transform(p_n_pars[1].begin(), p_n_pars[1].end(), p_n_pars[1].begin(), [](double p) { return p * 1e6; });
    std::transform(Q_pars[1].begin(), Q_pars[1].end(), Q_pars[1].begin(), [](double p) { return p / 3600; });
    std::transform(p_l_pars[1].begin(), p_l_pars[1].end(), p_l_pars[1].begin(), [](double p) { return p * 1e6; });
    std::transform(Q_out_pars[1].begin(), Q_out_pars[1].end(), Q_out_pars[1].begin(), [](double p) { return p / 3600; });

    parameters_series_t parameters;
    parameters.input_dens_visc(
        rho_pars[0], rho_pars[1], 
        visc_pars[0], visc_pars[1]
    );
    parameters.input_parameters(
        { p_n_pars[0], Q_pars[0], p_l_pars[0], Q_out_pars[0]},
        { p_n_pars[1], Q_pars[1], p_l_pars[1], Q_out_pars[1]}
    );

    // Буффер с проблемно-ориентированными слоями
    ring_buffer_t<density_viscosity_layer> buffer(2, pipe.profile.getPointCount());
    // Инициализация начальной плотности и вязкости в трубе
    auto& rho_initial = buffer.previous().density;
    auto& viscosity_initial = buffer.previous().viscosity;
    rho_initial = vector<double>(rho_initial.size(), interpolation(0, parameters.density_series, parameters.density_left_index));
    viscosity_initial = vector<double>(viscosity_initial.size(), interpolation(0, parameters.viscosity_series, parameters.viscosity_left_index));

    double modeling_time = 0;
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];


    vector<double> press_profile(pipe.profile.getPointCount());
    // Профиль для начального профиля давления
    vector<double> initial_press_profile(pipe.profile.getPointCount());
    // Создаём вектор изменения давлений относитльно начального профиля
    vector<double> diff_prof(press_profile.size(), 0);
    vector<double> diff_prof_pout;
    vector<double> modeling_moments;


    while (modeling_time <= T)
    {

        // Вектор скорости потока по координате
        double Q_in = interpolation(modeling_time, parameters.param_series[1], parameters.params_left_index[1]);
        double Q_out = interpolation(modeling_time, parameters.param_series[3], parameters.params_left_index[3]);
        double Q_n = (Q_in + Q_out) / 2;
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

        double inter_pout = interpolation(modeling_time, parameters.param_series[2], parameters.params_left_index[2]);
        diff_prof_pout.push_back(fabs(inter_pout - press_profile.back()) / inter_pout * 100);
        // Вывод в файлы
        write_profiles_problem(prev, dx, modeling_time);
        write_press_profile_only(press_profile, dx, modeling_time);
        write_press_profile_only(diff_prof, dx, modeling_time, "output/diff_press_profile.csv");

        modeling_moments.push_back(modeling_time);
        modeling_time += moc_solv.prepare_step();
    }

    uni_write(modeling_moments, 0, { diff_prof_pout }, "Время,Время моделирования,Погрешность (%)", "output/diff_press_pout.csv");
}