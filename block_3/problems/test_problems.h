#pragma once
/// @brief ��� ������ ��� ������-������� � ������ ������
typedef std::function<double(size_t& index)> diff_function_t;
/// @brief ��� ������ ��� �������� ����
typedef composite_layer_t<profile_collection_t<2>,
    moc_solver<2>::specific_layer> layer_t;

/// @brief ������� ��� ������ ������ �������� �������� 
/// � ������ ������� �������
/// @param press ������ �� ������� ��������
/// @param dx ��� �� ����������
/// @param dt ��� �� �������
/// @param step ������� ��� �������������
/// @param filename �������� ����� ��� ������
void write_press_profile_only(vector<double>& press, double& dx, double& dt, size_t step, std::string filename = "output/press_prof.csv") {
    std::ofstream press_file;
    size_t profCount = press.size();
    if (step == 0)
    {
        press_file.open(filename);
        press_file << "time,x,��������" << std::endl;
    }
    else
        press_file.open(filename, std::ios::app);
  
    for (int i = 0; i < profCount; i++)
    {
        press_file << dt * step << "," << i * dx;
        press_file << "," << press[i] << std::endl;
    }

    press_file.close();
}
/// @brief ������� ��� ������ �������� ���� ����������
/// � ������ ������� �������
/// @param layer ������ �� ���� 
/// @param press_prof ������ �� ������� ��������
/// @param dx ��� �� ����������
/// @param dt ��� �� �������
/// @param step ������ ��� �������������
/// @param filename �������� ����� ��� ������
void write_profiles(
    layer_t& layer, vector<double>& press_prof, double& dx, double& dt, size_t step,
    std::string filename = "output/profiles.csv")
{
    std::ofstream file;
    std::ofstream press_file;
    size_t profCount = layer.vars.point_double[0].size();
    if (step == 0)
    {
        file.open(filename);
        file << "time,x,���������,��������,��������" << std::endl;
    }
    else
        file.open(filename, std::ios::app);
        press_file.open(filename, std::ios::app);
    
    for (int i = 0; i < profCount; i++)
    {
        file << dt * step << "," << i * dx;
        file << "," << layer.vars.point_double[0][i] << "," << layer.vars.point_double[1][i] << "," << press_prof[i] << std::endl;
    }
    file.close();
}

/// @brief ����� ��� ������� ����� �� ��������������� 
class Quasistationary : public ::testing::Test
{
public:
    
    /// @brief ������� ��� ������� ������� �������������
    /// @param prev ������ �� ���������� ����
    /// @param next ������ �� ������� ����
    /// @param rho_in ��������� ����������� ������
    /// @param visc_in �������� ����������� ������
    /// @param direction ����������� ������� ������
    void moc_solve(layer_t& prev, layer_t& next, double& rho_in, double& visc_in, int direction = 1)
    {
        size_t start_index = direction > 0 ? 1 : (prev.vars.point_double[0].size()) - 2;
        size_t end_index = direction < 0 ? 0 : (prev.vars.point_double[0].size());
        size_t num_profiles = prev.vars.point_double.size();
        double parametrs_in[] = { rho_in, visc_in };
        for (size_t index = start_index; index != end_index; index += direction)
        {
            for (size_t p = 0; p < num_profiles; p++)
            {
                if (index == start_index)
                    next.vars.point_double[p][index - direction] = parametrs_in[p];
                
                next.vars.point_double[p][index] = prev.vars.point_double[p][index - direction];
            }
        }
    }

    /// @brief �������� ������� ������� ������
    /// @param press_prof ������ �� ������� ��������
    /// @param right_part ����������� � �����, ���������� �� ��� �� ���������� 
    /// @param direction ����������� ������� ��������
    void QP_Euler_solver(vector<double>& press_prof, const diff_function_t& right_part, int& direction)
    {
        size_t start_index = direction > 0 ? 1 : (press_prof.size()) - 2;
        size_t end_index = direction < 0 ? 0 : (press_prof.size());
        for (size_t index = start_index; index != end_index; index += direction)
        {
            press_prof[index] = press_prof[index - direction] + direction * right_part(index);
        }
    }

    /// @brief ������� ������� ����������� � ������ ���������
    /// @param layer ������ �� ������� ����
    /// @param press_prof ������ �� ������� ��������
    /// @param speed �������� ������
    /// @param direction ����������� ������� ��������
    void euler_solve(layer_t& layer, vector<double>& press_prof, double& speed, int direction = 1)
    {
        // ������� ���������, ��� ����� ��� ������� �������� ������ 
        vector<double>& density = layer.vars.point_double[0];
        // ������� ��������, ��� ����� ��� ������� �������� ������ 
        vector<double>& viscosity = layer.vars.point_double[1];
        // ������� �����������
        // ���������� �������� ����������� � �����, ���������� �� ��� �� ����������
        diff_function_t right_part =
            [this, speed, density, viscosity, direction](size_t& index)
            {
                double eps = pipe.wall.relativeRoughness(); // ������ ������������� �������������
                double Re = speed * pipe.wall.diameter / viscosity[index]; // ������ ����� ����������
                double lambda = pipe.resistance_function(Re, eps); // ������ ������������ ������
                double dz = pipe.profile.heights[index - direction] - pipe.profile.heights[index]; // ������ �������� �����
                double dx = pipe.profile.coordinates[index - direction] - pipe.profile.coordinates[index]; // ������ ���� �� ����������
                // ������ �����������
                double diff = lambda * (1 / pipe.wall.diameter) * density[index] * pow(speed, 2) / 2 - dz / dx * density[index] * M_G;
                return dx * diff;
            };

        QP_Euler_solver(press_prof, right_part, direction);
    }

protected:
    /// @brief ��������� �����
    pipe_properties_t pipe; 
    /// @brief ��������� �����
    oil_parameters_t oil;
    // ����� �������������
    size_t T = 1000;
    // �������� � ������ ������� ������������
    double p0 = 6e6;
    // ��������� ����������� ������
    double rho_in = 800;
    // �������� ����������� ������
    double visc_in = 10e-6;
    // ������ ������
    double flow = 0.5;

    /// @brief �����������, ���������������� ��������� ����� � ����� 
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
    ring_buffer_t<layer_t> buffer(2, pipe.profile.getPointCount());

    buffer.advance(+1);
    // ������������� ��������� ���������
    buffer.previous().vars.point_double[0] = vector<double>(pipe.profile.getPointCount(), oil.density.nominal_density); 
    // ������������� ��������� ��������
    buffer.previous().vars.point_double[1] = vector<double>(pipe.profile.getPointCount(), oil.viscosity.nominal_viscosity);
    // ������������� ������� ��������
    vector<double> press_prof(pipe.profile.getPointCount(), p0);
   
    double speed = flow / pipe.wall.getArea(); // ������ �������� ������
    double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
    double dt = (dx) / speed;
    size_t N = static_cast<size_t>(T / dt);
    
    euler_solve(buffer.previous(), press_prof, speed); // ������ ������� �������� � ��������� ������ �������

    write_profiles(buffer.previous(), press_prof, dx, dt, 0);
    write_press_profile_only(press_prof, dx, dt, 0);

    for (size_t i = 1; i < N; i++)
    {
        moc_solve(buffer.previous(), buffer.current(), rho_in, visc_in);
        euler_solve(buffer.current(), press_prof, speed);

        write_profiles(buffer.current(), press_prof, dx, dt, i);
        write_press_profile_only(press_prof, dx, dt, i);

        buffer.advance(+1);
    }

}