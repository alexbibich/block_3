#pragma once

//--gtest_catch_exceptions=0 --gtest_filter=Quasistationary.EulerWithMOC_step_inter

#include "transport_moc_solver.h"
#include "newton_solver.h"
#include "euler_solver.h"


inline std::string get_test_string() {
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    auto test_string = std::string(test_info->test_case_name()) + "/" + std::string(test_info->name());
    return test_string;
}

inline std::string prepare_test_folder()
{
    std::string path = std::string("research_out/") + get_test_string() + "/";
    std::filesystem::create_directories(path);
    return path;
}

/// @brief Функция для записи только профилей давления 
/// в разные моменты времени
/// @param press ссылка на профиль давления
/// @param dx шаг по координате
/// @param dt шаг по времени
/// @param step текущий шаг моделирования
/// @param filename название файла для записи
void write_press_profile_only(vector<double>& press, double dx, double time_moment, std::string filename) {
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

template<typename argType=double>
void uni_write(vector<argType>& x, double time_moment, vector<vector<double>> params, std::string params_name, std::string filename) {
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
    std::string filename)
{
    std::ofstream file;
    if (time_moment == 0)
    {
        file.open(filename + "profiles_dens_visc.csv");
        file << "time,x,Плотность,Вязкость" << std::endl;
    }
    else
        file.open(filename + "profiles_dens_visc.csv", std::ios::app);

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

    file.close();
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
    //vector<std::string> dates;

    while (getline(file, line))
    {
        std::replace(line.begin(), line.end(), ',', '.');
        std::stringstream stream(line);
        std::string val;
        std::string moment;

        getline(stream, moment, delimiter);

        if (moment.find(".08.2021") != std::string::npos)
        {
            //dates.push_back(moment);
            getline(stream, val, delimiter);
            double value = stod(val);

            vals.push_back(value);
            times.push_back(pars_time_line(moment));

            if (value > 0)
            {
                vals.push_back(value);
                times.push_back(pars_time_line(moment));
            }
            
            
        }
        
    }

    file.close();


    //filename.erase(filename.length() - 4, 4);
    //uni_write(times, 0, { vals }, "time,время_" + filename + ',' + filename, "time_series/" + filename + "_series.csv");
    //uni_write<std::string>(dates, 0, { vals }, "time,время_" + filename + ',' + filename, "time_series/" + filename + "_series.csv");

    return { times, vals };
}



/// @brief Класс для решения задач по квазистационару 
class Quasistationary : public ::testing::Test
{
public:
    void write_press(double dx, double moment, vector<double>& press_profile, vector<double>& diff_prof, std::string path)
    {
        write_press_profile_only(press_profile, dx, moment, path + "press_prof.csv");
        write_press_profile_only(diff_prof, dx, moment, path + "diff_press_profile.csv");
    }

    void write_diff(vector<double>& moments, vector<double>& diff_prof, std::string name, std::string path)
    {
        std::transform(diff_prof.begin(), diff_prof.end(), diff_prof.begin(), [](double p) {return p / 1000;  });
        uni_write(moments, 0, { diff_prof },
            "Время,Время моделирования," + name, path + "diff_press_pout.csv");
    }

protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Параметры нефти
    oil_parameters_t oil;
    // Время моделирования
    size_t T = 3600 * 24 * 31;
    // Давление в начале участка трубопровода
    double p0 = 6e6;
    // Расход потока
    double flow = 0.5;

    vector<vector<double>> rho_pars;
    vector<vector<double>> visc_pars;
    vector<vector<double>> p_n_pars;
    vector<vector<double>> Q_pars;
    vector<vector<double>> p_l_pars;
    vector<vector<double>> Q_out_pars;


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

        uni_write(pipe.profile.coordinates, 0, { pipe.profile.heights }, "Время,Координаты,Высотки", "time_series/pipe_coord_heights.csv");

        rho_pars = parser_parameter("rho_in.csv");
        visc_pars = parser_parameter("visc_in.csv");
        p_n_pars = parser_parameter("p_in.csv");
        Q_pars = parser_parameter("Q_in.csv");
        p_l_pars = parser_parameter("p_out.csv");
        Q_out_pars = parser_parameter("Q_out.csv");

        std::transform(visc_pars[1].begin(), visc_pars[1].end(), visc_pars[1].begin(), [](double visc) { return visc * 1e-6; });
        std::transform(p_n_pars[1].begin(), p_n_pars[1].end(), p_n_pars[1].begin(), [](double p) { return p * 1e6; });
        std::transform(Q_pars[1].begin(), Q_pars[1].end(), Q_pars[1].begin(), [](double p) { return p / 3600; });
        std::transform(p_l_pars[1].begin(), p_l_pars[1].end(), p_l_pars[1].begin(), [](double p) { return p * 1e6; });
        std::transform(Q_out_pars[1].begin(), Q_out_pars[1].end(), Q_out_pars[1].begin(), [](double p) { return p / 3600; });
    }
};



TEST_F(Quasistationary, WithStationaryMean)
{
    std::string method = "step";
    std::string path = prepare_test_folder();

    double mean_density = std::accumulate(rho_pars[1].begin(), rho_pars[1].end(), 0.0) / rho_pars[1].size();
    double mean_visc = std::accumulate(visc_pars[1].begin(), visc_pars[1].end(), 0.0) / visc_pars[1].size();

    parameters_series_t parameters;
    parameters.input_parameters(
        { p_n_pars[0], Q_pars[0], p_l_pars[0], Q_out_pars[0] },
        { p_n_pars[1], Q_pars[1], p_l_pars[1], Q_out_pars[1] }
    );

    vector<double> rho_prof(pipe.profile.getPointCount(), mean_density);
    vector<double> visc_prof(pipe.profile.getPointCount(), mean_visc);

    double modeling_time = 0;
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
    double dt = 500;


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
        double Q_in = interpolation(modeling_time, parameters.param_series[1], parameters.params_left_index[1], method);
        double Q_out = interpolation(modeling_time, parameters.param_series[3], parameters.params_left_index[3], method);
        double Q_n = (Q_in + Q_out) / 2;
        vector<double> Q(pipe.profile.getPointCount(), Q_n);


        if (modeling_time == 0)
        {
            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, rho_prof, visc_prof, Q_n);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            initial_press_profile = press_profile;
        }
        else
        {

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, rho_prof, visc_prof, Q_n);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            // получение вектора изменения давлений относитльно начального профиля
            std::transform(initial_press_profile.begin(), initial_press_profile.end(), press_profile.begin(), diff_prof.begin(),
                [](double initial, double current) {return initial - current; });
        }

        double inter_pout = interpolation(modeling_time, parameters.param_series[2], parameters.params_left_index[2], method);
        diff_prof_pout.push_back(inter_pout - press_profile.back());
        // Вывод в файлы
        write_press(dx, modeling_time, press_profile, diff_prof, path);

        modeling_moments.push_back(modeling_time);
        modeling_time += dt;
    }

    write_diff(modeling_moments, diff_prof_pout, "Стац - средняя реология", path);
};

TEST_F(Quasistationary, WithStationaryCurrent)
{
    std::string method = "step";
    std::string path = prepare_test_folder();

    double mean_density = std::accumulate(rho_pars[1].begin(), rho_pars[1].end(), 0.0) / rho_pars[1].size();
    double mean_visc = std::accumulate(visc_pars[1].begin(), visc_pars[1].end(), 0.0) / visc_pars[1].size();

    parameters_series_t parameters;
    parameters.input_dens_visc(
        rho_pars[0], rho_pars[1],
        visc_pars[0], visc_pars[1]
    );
    parameters.input_parameters(
        { p_n_pars[0], Q_pars[0], p_l_pars[0], Q_out_pars[0] },
        { p_n_pars[1], Q_pars[1], p_l_pars[1], Q_out_pars[1] }
    );

    double modeling_time = 0;
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
    double dt = 500;


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
        double Q_in = interpolation(modeling_time, parameters.param_series[1], parameters.params_left_index[1], method);
        double Q_out = interpolation(modeling_time, parameters.param_series[3], parameters.params_left_index[3], method);
        double Q_n = (Q_in + Q_out) / 2;
        vector<double> Q(pipe.profile.getPointCount(), Q_n);

        vector<double> rho_prof(pipe.profile.getPointCount(), 
            interpolation(modeling_time, parameters.density_series, parameters.density_left_index, method));
        vector<double> visc_prof(pipe.profile.getPointCount(), 
            interpolation(modeling_time, parameters.viscosity_series, parameters.viscosity_left_index, method));

        if (modeling_time == 0)
        {
            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, rho_prof, visc_prof, Q_n);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            initial_press_profile = press_profile;
        }
        else
        {

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, rho_prof, visc_prof, Q_n);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            // получение вектора изменения давлений относитльно начального профиля
            std::transform(initial_press_profile.begin(), initial_press_profile.end(), press_profile.begin(), diff_prof.begin(),
                [](double initial, double current) {return initial - current; });
        }

        double inter_pout = interpolation(modeling_time, parameters.param_series[2], parameters.params_left_index[2], method);
        diff_prof_pout.push_back(inter_pout - press_profile.back());
        // Вывод в файлы
        write_press(dx, modeling_time, press_profile, diff_prof, path);

        modeling_moments.push_back(modeling_time);
        modeling_time += dt;
    }

    write_diff(modeling_moments, diff_prof_pout, "Стац - текущая реология", path);
};


TEST_F(Quasistationary, QuasiWithStepInter)
{
    std::string method = "step";
    std::string path = prepare_test_folder();

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
    rho_initial = vector<double>(rho_initial.size(), interpolation(0, parameters.density_series, parameters.density_left_index, method));
    viscosity_initial = vector<double>(viscosity_initial.size(), interpolation(0, parameters.viscosity_series, parameters.viscosity_left_index, method));

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
        double Q_in = interpolation(modeling_time, parameters.param_series[1], parameters.params_left_index[1], method);
        double Q_out = interpolation(modeling_time, parameters.param_series[3], parameters.params_left_index[3], method);
        double Q_n = (Q_in + Q_out) / 2;
        vector<double> Q(pipe.profile.getPointCount(), Q_n);

        density_viscosity_layer& prev = buffer.previous();
        density_viscosity_layer& next = buffer.current();

        transport_moc_solver moc_solv(pipe, Q, prev, next);


        if (modeling_time == 0)
        {
            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, prev.density, prev.viscosity, Q_n);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            initial_press_profile = press_profile;
        }
        else
        {
            array<double, 2> par_in = moc_solv.get_par_in(modeling_time, parameters, method);
            moc_solv.step(par_in);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, next.density, next.viscosity, Q_n);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            // получение вектора изменения давлений относитльно начального профиля
            std::transform(initial_press_profile.begin(), initial_press_profile.end(), press_profile.begin(), diff_prof.begin(),
                [](double initial, double current) {return initial - current; });

            buffer.advance(+1);
        }

        double inter_pout = interpolation(modeling_time, parameters.param_series[2], parameters.params_left_index[2], method);
        diff_prof_pout.push_back(inter_pout - press_profile.back());
        // Вывод в файлы
        write_profiles_problem(prev, dx, modeling_time, path);
        write_press(dx, modeling_time, press_profile, diff_prof, path);

        modeling_moments.push_back(modeling_time);
        modeling_time += moc_solv.prepare_step();
    }

    write_diff(modeling_moments, diff_prof_pout, "Квазистац - ступеньки", path);
};

TEST_F(Quasistationary, QuasiWithLineInter)
{
    std::string method = "line";
    std::string path = prepare_test_folder();

    parameters_series_t parameters;
    parameters.input_dens_visc(
        rho_pars[0], rho_pars[1],
        visc_pars[0], visc_pars[1]
    );
    parameters.input_parameters(
        { p_n_pars[0], Q_pars[0], p_l_pars[0], Q_out_pars[0] },
        { p_n_pars[1], Q_pars[1], p_l_pars[1], Q_out_pars[1] }
    );

    // Буффер с проблемно-ориентированными слоями
    ring_buffer_t<density_viscosity_layer> buffer(2, pipe.profile.getPointCount());
    // Инициализация начальной плотности и вязкости в трубе
    auto& rho_initial = buffer.previous().density;
    auto& viscosity_initial = buffer.previous().viscosity;
    rho_initial = vector<double>(rho_initial.size(), interpolation(0, parameters.density_series, parameters.density_left_index, method));
    viscosity_initial = vector<double>(viscosity_initial.size(), interpolation(0, parameters.viscosity_series, parameters.viscosity_left_index, method));

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
        double Q_in = interpolation(modeling_time, parameters.param_series[1], parameters.params_left_index[1], method);
        double Q_out = interpolation(modeling_time, parameters.param_series[3], parameters.params_left_index[3], method);
        double Q_n = (Q_in + Q_out) / 2;
        vector<double> Q(pipe.profile.getPointCount(), Q_n);

        density_viscosity_layer& prev = buffer.previous();
        density_viscosity_layer& next = buffer.current();

        transport_moc_solver moc_solv(pipe, Q, prev, next);


        if (modeling_time == 0)
        {
            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, prev.density, prev.viscosity, Q_n);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            initial_press_profile = press_profile;
        }
        else
        {
            array<double, 2> par_in = moc_solv.get_par_in(modeling_time, parameters, method);
            moc_solv.step(par_in);

            double p_n = interpolation(modeling_time, parameters.param_series[0], parameters.params_left_index[0], method);

            // Создаем расчетную модель трубы
            Pipe_model_for_PQ_t pipeModel(pipe, next.density, next.viscosity, Q_n);

            solve_euler_corrector<1>(pipeModel, 1, p_n, &press_profile);

            // получение вектора изменения давлений относитльно начального профиля
            std::transform(initial_press_profile.begin(), initial_press_profile.end(), press_profile.begin(), diff_prof.begin(),
                [](double initial, double current) {return initial - current; });

            buffer.advance(+1);
        }

        double inter_pout = interpolation(modeling_time, parameters.param_series[2], parameters.params_left_index[2], method);
        diff_prof_pout.push_back(inter_pout - press_profile.back());
        // Вывод в файлы
        write_profiles_problem(prev, dx, modeling_time, path);
        write_press(dx, modeling_time, press_profile, diff_prof, path);

        modeling_moments.push_back(modeling_time);
        modeling_time += moc_solv.prepare_step();
    }

    write_diff(modeling_moments, diff_prof_pout, "Квазистац - линейная", path);
};