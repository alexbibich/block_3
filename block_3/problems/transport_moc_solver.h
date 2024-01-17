#pragma once

#define TIME_INDEX 0
#define PAR_INDEX 1

typedef vector < array<double, 2> > time_series_t;
/// @brief тип данных для хранения слоёв
typedef profile_collection_t<3> layer_t;

class transport_moc_solver
{
public:

    transport_moc_solver(const pipe_properties_t& pipe,
        const vector<double>& vol_flow, const vector<time_series_t>& par)
        : pipe(pipe),
        Q(vol_flow),
        parameters{ par }
    {}

    /// @brief Алгоритм расчёта методом характеристик
    /// @param prev Ссылка на предыдущий профиль
    /// @param next Ссылка на текущий профиль
    /// @param par_in Параметр вытесняющей партии
    /// @param direction Направление течения потока
    void solve(vector<double>& prev, vector<double>& next, double& par_in, int& direction)
    {
        size_t start_index = direction > 0 ? 1 : (prev.size()) - 2;
        size_t end_index = direction < 0 ? 0 : (prev.size());
        next[start_index - direction] = par_in;
        for (size_t index = start_index; index != end_index; index += direction)
            next[index] = prev[index - direction];
    }

    /// @brief Функция для решения методом характеристик
    /// @param prev Ссылка на предыдущий слой
    /// @param next Ссылка на текущий слой
    /// @param rho_in Плотность вытесняющей партии
    /// @param visc_in Вязкость вытесняющей партии
    /// @param direction Направление течения потока
    void step(layer_t& prev, layer_t& next, double& time_moment)
    {
        size_t num_profiles = prev.point_double.size();
        vector<double> parameters_in = get_inter_parameters(time_moment);
        int direction = getEigenvals(0) > 0 ? 1 : -1;
        for (size_t p = 0; p < num_profiles; p++)
            solve(prev.point_double[p], next.point_double[p], parameters_in[p], direction);
    }

    vector<double> get_inter_parameters(double& time_moment)
    {
        vector<double> input_parameters{ interpolation(time_moment, parameters[0]), interpolation(time_moment, parameters[1]) };
        return input_parameters;
    }

    vector<double> getGrid()
    {
        return pipe.profile.coordinates;
    }

    double getEigenvals(size_t index)
    {
        double S = pipe.wall.getArea();
        return Q[index] / S;
    }

    double prepare_step()
    {
        vector<double> grid = getGrid();
        double max_eigen = 0;

        for (size_t index = 0; index < grid.size(); index++)
        {
            double eigen_index = getEigenvals(index);
            max_eigen = std::max(max_eigen, eigen_index);
        }

        double dx = grid[1] - grid[0];
        double courant_step = dx / max_eigen;

        return courant_step;
    }

    static double interpolation(double& time_moment, const time_series_t& parameter)
    {
        size_t left_index = 0;
        size_t right_index = parameter.size()-1;
        size_t center_index = right_index / 2;

        if (parameter[left_index][TIME_INDEX] <= time_moment)
        {
            if (parameter[right_index][TIME_INDEX] <= time_moment)
                return parameter[right_index][PAR_INDEX];
        }
        else
            return parameter[left_index][PAR_INDEX];

        while (left_index != center_index)
        {
            if (parameter[center_index][TIME_INDEX] >= time_moment)
                right_index = center_index;
            else
                left_index = center_index;

            center_index = (right_index + left_index) / 2;
        }
        double dt = parameter[right_index][TIME_INDEX] - parameter[left_index][TIME_INDEX];
        double u1 = (parameter[right_index][TIME_INDEX] - time_moment) / dt;
        double u2 = (time_moment - parameter[left_index][TIME_INDEX]) / dt;
        double res_par = parameter[left_index][PAR_INDEX] * u1 + parameter[right_index][PAR_INDEX] * u2;

        return res_par;
    }

protected:
    /// @brief модель трубы
    const pipe_properties_t& pipe;
    /// @brief Объемный расход
    const vector<double>& Q;
    /// @brief Временные ряды краевых условий
    const vector<time_series_t>& parameters;
};


struct input_parameters_t
{
    vector<time_series_t> parameters_series;
    size_t count_parameters;

    input_parameters_t(size_t& count_par)
        : count_parameters{ count_par }
    {}

    void input_parameters(vector<vector<double>>& moments, vector<vector<double>>& par)
    {
        for (size_t index = 0; index < count_parameters; index++)
        {
            time_series_t series = build_parameters(moments[index], par[index]);
            parameters_series.push_back(series);
        }
    }

    void input_parameters(vector<vector<vector<double>>>& parameters_moments)
    {
        for (size_t index = 0; index < count_parameters; index++)
        {
            vector<double> moments;
            if (parameters_moments[index][TIME_INDEX].size() == 1)
                moments = build_series(parameters_moments[index][TIME_INDEX][0], parameters_moments[index][PAR_INDEX].size());
            else
                moments = parameters_moments[index][TIME_INDEX];

            time_series_t series = build_parameters(moments, parameters_moments[index][PAR_INDEX]);
            parameters_series.push_back(series);
        }
    }

    void input_parameters(double& dt, vector<vector<double>>& par)
    {
        vector<double> moments = build_series(dt, par[0].size());
        for (size_t index = 0; index < count_parameters; index++)
        {
            time_series_t series = build_parameters(moments, par[index]);
            parameters_series.push_back(series);
        }
    }

    time_series_t build_parameters(vector<double>& moments, vector<double>& parameter)
    {
        time_series_t par(parameter.size());
        for (size_t index = 0; index < parameter.size(); index++)
        {
            par[index][0] = moments[index];
            par[index][1] = parameter[index];
        }

        return par;
    }

    vector<double> build_series(double& dt, size_t length)
    {
        vector<double> time_ser(length);
        for (size_t index = 0; index < length; index++)
            time_ser[index] = dt * index;
        return time_ser;
    }

    
  
};