#pragma once

typedef std::function<double(size_t& index)> diff_function_t;
typedef composite_layer_t<profile_collection_t<2>,
    moc_solver<2>::specific_layer> layer_t;

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

void write_profiles(
    layer_t& layer, vector<double>& press_prof, double& dx, double& dt, size_t step,
    std::string filename = "output/profiles.csv")
{
    std::ofstream file;
    size_t profCount = layer.vars.point_double[0].size();
    if (step == 0)
    {
        file.open(filename);
        file << "time,x,Плотность,Вязкость, Давление" << std::endl;
    } else
        file.open(filename, std::ios::app);
    

    for (int i = 0; i < profCount; i++)
    {
        file << dt * step << "," << i * dx;
        file << "," << layer.vars.point_double[0][i] << "," << layer.vars.point_double[1][i] * 100000000 << "," << press_prof[i] << std::endl;
    }
    file.close();
}

class EulerWithMOC : public ::testing::Test
{
public:

    void moc_solv(layer_t& prev, layer_t& next, double parametrs_in[])
    {
        double num_points = prev.vars.point_double[0].size();
        double num_profiles = prev.vars.point_double.size();
        for (size_t index = 0; index < num_points; index++)
        {
            for (size_t p = 0; p < num_profiles; p++)
            {
                if (index == 0)
                    next.vars.point_double[p][index] = parametrs_in[p];
                else
                    next.vars.point_double[p][index] = prev.vars.point_double[p][index];
            }
        }
    }

    void QP_Euler_solver(vector<double>& press_prof, const diff_function_t& right_part, int& direction)
    {
        size_t start_index = direction > 0 ? 1 : (press_prof.size()) - 2;
        size_t end_index = direction < 0 ? 0 : (press_prof.size()) - 1;
        for (size_t index = start_index; index != end_index; index += direction)
        {
            double p = press_prof[index - direction];
            press_prof[index] = press_prof[index - direction] + direction * right_part(index);
        }
    }

    void euler_solve(layer_t& layer, vector<double>& press_prof, double& speed, int direction = 1)
    {
        vector<double>& density = layer.vars.point_double[0];
        vector<double>& viscosity = layer.vars.point_double[1];
        diff_function_t right_part =
            [this, speed, density, viscosity, direction](size_t& index)
            {
                double eps = pipe.wall.relativeRoughness(); // Расчёт относительной шероховатости
                double Re = speed * pipe.wall.diameter / viscosity[index]; // Расчёт числа Рейнольдса
                double lambda = pipe.resistance_function(Re, eps);
                double dz = pipe.profile.heights[index] - pipe.profile.heights[index + direction];
                double dx = pipe.profile.coordinates[index] - pipe.profile.coordinates[index + direction];
                double diff = lambda * (1 / pipe.wall.diameter) * density[index] * pow(speed, 2) / 2 - dz / dx * density[index];
                return dx * diff;
            };

        QP_Euler_solver(press_prof, right_part, direction);
    }

protected:
    pipe_properties_t pipe;
    oil_parameters_t oil;

    virtual void SetUp() override
    {
        init_cond(pipe, oil);
    }
};

TEST_F(EulerWithMOC, EulerWithCharacteristic)
{
    // Время моделирования
    size_t T = 10000;

    double p0 = 6e6;
    double rho_in = 850;
    double visc_in = 13e-6;
    double input_par[] = { rho_in, visc_in };

    ring_buffer_t<layer_t> buffer(2, pipe.profile.getPointCount());

    buffer.advance(+1);
    // инициализация начальной плотности
    buffer.previous().vars.point_double[0] = vector<double>(pipe.profile.getPointCount(), oil.density.nominal_density); 
    // инициализация начальной вязкости
    buffer.current().vars.point_double[1] = vector<double>(pipe.profile.getPointCount(), oil.viscosity.nominal_viscosity); 
    // инициализация профиля давления
    vector<double> press_prof(pipe.profile.getPointCount(), p0);
    // инициализация профиля расхода
    vector<double> Q(pipe.profile.getPointCount(), 0.5);

    PipeQAdvection advection_model(pipe, Q);
    double speed = advection_model.getEquationsCoeffs(0, 0);
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
    double dt = (dx) / speed;
    size_t N = static_cast<size_t>(T / dt);
    
    euler_solve(buffer.previous(), press_prof, speed);
    write_profiles(buffer.previous(), press_prof, dx, dt, 0);

    for (size_t i = 1; i < N; i++)
    {
        moc_solv(buffer.previous(), buffer.current(), input_par);
        euler_solve(buffer.current(), press_prof, speed);
        write_profiles(buffer.current(), press_prof, dx, dt, i);
        buffer.advance(+1);
    }

}