#pragma once

#include "transport_moc_solver.h"

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
void write_press_profile_only(vector<double>& press, double& dx, double& dt, size_t step, std::string filename = "output/press_prof.csv") {
    std::ofstream press_file;
    size_t profCount = press.size();
    if (step == 0)
    {
        press_file.open(filename);
        press_file << "time,x,Давление" << std::endl;
    }
    else
        press_file.open(filename, std::ios::app);

    for (int i = 0; i < profCount; i++)
    {
        press_file << dt * step << "," << i * dx;
        press_file << "," << press[i] << std::endl;
    }

    press_file.close();
};

/// @brief Функция для записи профилей всех параметров
/// в разные моменты времени
/// @param layer ссылка на слой 
/// @param press_prof ссылка на профиль давления
/// @param dx шаг по координате
/// @param dt шаг по времени
/// @param step текщий шаг моделирования
/// @param filename название файла для записи
void write_profiles(
    layer_t& layer, double& dx, double& dt, size_t step,
    std::string filename = "output/profiles.csv")
{
    std::ofstream file;
    std::ofstream press_file;
    size_t profCount = layer.point_double[0].size();
    if (step == 0)
    {
        file.open(filename);
        file << "time,x,Плотность,Вязкость,Давление" << std::endl;
    }
    else
        file.open(filename, std::ios::app);
        press_file.open(filename, std::ios::app);
    
    for (int i = 0; i < profCount; i++)
    {
        file << dt * step << "," << i * dx;
        file << "," << layer.point_double[0][i] << "," << layer.point_double[1][i] << "," << layer.point_double[2][i] << std::endl;
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
    void euler_solve(layer_t& layer_prev, layer_t& layer_next, double& speed, size_t direction = 1)
    {
        // Профиль плотности, для учёта при рисчёте движения партий 
        vector<double>& density = layer_next.point_double[0];
        // Профиль вязкости, для учёта при рисчёте движения партий 
        vector<double>& viscosity = layer_next.point_double[1];
        // Перенос значений давления с предыдущего на текущий слой
        layer_next.point_double[2] = layer_prev.point_double[2];

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

        QP_Euler_solver(layer_next.point_double[2], right_part, direction);
    }

protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe; 
    /// @brief Параметры нефти
    oil_parameters_t oil;
    // Время моделирования
    size_t T = 1000;
    // Давление в начале участка трубопровода
    double p0 = 6e6;
    // Плотность вытесняющей партии
    double rho_in = 800;
    // Вязкость вытесняющей партии
    double visc_in = 10e-6;
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

TEST_F(Quasistationary, EulerWithMOC)
{
    //ring_buffer_t<layer_t> buffer(2, pipe.profile.getPointCount());

    //buffer.advance(+1);
    //// инициализация начальной плотности
    //buffer.previous().point_double[0] = vector<double>(pipe.profile.getPointCount(), oil.density.nominal_density); 
    //// инициализация начальной вязкости
    //buffer.previous().point_double[1] = vector<double>(pipe.profile.getPointCount(), oil.viscosity.nominal_viscosity);
    //// инициализация профиля давления
    //buffer.previous().point_double[2] = vector<double>(pipe.profile.getPointCount(), p0);
   
    //double speed = flow / pipe.wall.getArea(); // Расчёт скорости потока
    //double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
    //double dt = (dx) / speed;
    //size_t N = static_cast<size_t>(T / dt);
    //double input_parameters[] = { rho_in, visc_in };

    //for (size_t i = 1; i < N; i++)
    //{
    //    moc_solve(buffer.previous(), buffer.current(), input_parameters);
    //    euler_solve(buffer.previous(), buffer.current(), speed);

    //    write_profiles(buffer.current(), dx, dt, i);
    //    write_press_profile_only(buffer.current().point_double[2], dx, dt, i);

    //    buffer.advance(+1);
    //}

}

TEST_F(Quasistationary, Testing)
{
    double dt = 0.5;

    vector<double> Q = { 0, 1, 2, 3, 4, 5};
    vector<double> rho = { 3, 4, 5, 6, 7, 8 };
    vector<double> visc = { 6, 7, 8, 9, 10, 11 };

    vector<vector<double>> parameters_val{ Q, rho, visc };
    size_t count_input_series{ parameters_val.size() };

    input_parameters_t parameters(count_input_series);

    vector<vector<double>> moments(count_input_series, parameters.build_series(dt, count_input_series));

    vector<vector<vector<double>>> test =
    {
        {{dt}, Q},
        {{dt}, rho},
        {{dt}, visc}
    };

    parameters.input_parameters(dt, parameters_val);

    double test_time = 1.5;

    transport_moc_solver::interpolation(test_time, parameters.parameters_series[0]);

}
