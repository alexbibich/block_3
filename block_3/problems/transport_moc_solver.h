#pragma once

typedef vector < array<double, 2> > time_series_t;

class transport_moc_solver
{
public:

    transport_moc_solver(const pipe_properties_t& pipe,
        const vector<double>& vol_flow)
        : pipe(pipe)
        , Q(vol_flow)
    {}

    vector<double> getGrid()
    {
        return pipe.profile.coordinates;
    }

    double getEigenvals(size_t index)
    {
        double S = pipe.wall.getArea();
        return Q[index] / S;
    }

    double prepare_step(double time_step = std::numeric_limits<double>::quiet_NaN())
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

        if (std::isnan(time_step) || time_step > courant_step) {
            time_step = courant_step;
        }
        return time_step;
    }

protected:
    /// @brief модель трубы
    const pipe_properties_t& pipe;
    /// @brief Объемный расход
    const vector<double>& Q;
};

template <typename time_type>
struct input_parameters_t
{
    vector<time_series_t> parameters_series;
    size_t count_parameters;

    input_parameters_t(size_t& count_par)
        : count_parameters{ count_par }
    {}

    void input_parameters(vector<time_type>& moments, vector<vector<double>>& par)
    {
        for (size_t index = 0; index < count_parameters; index++)
        {
            time_series_t series = build_parameters(moments[index], par[index]);
            parameters_series.push_back(series);
        }
    }

    time_series_t build_parameters(double& dt, vector<double>& parameter)
    {
        time_series_t par(parameter.size());
        for (size_t index = 0; index < parameter.size(); index++)
        {
            par[index][0] = parameter[index];
            par[index][1] = dt * index;
        }

        return par;
    }

    time_series_t build_parameters(vector<double>& moments, vector<double>& parameter)
    {
        time_series_t par(parameter.size());
        for (size_t index = 0; index < parameter.size(); index++)
        {
            par[index][0] = parameter[index];
            par[index][1] = moments[index];
        }

        return par;
    }

    
};


//template<size_t count_parameters>
//struct parameters_series_t
//{
//    vector<time_series_t> parameters_series(count_parameters);
//
//    void input_parameters()
//    {
//
//    }
//};