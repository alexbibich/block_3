#pragma once

#define TIME_INDEX 0
#define PAR_INDEX 1

#include<testing/test_moc.h>

typedef vector < array<double, 2> > time_series_t;
/// @brief тип данных для хранения слоёв
typedef profile_collection_t<3> layer_t;

class transport_moc_solver
{
public:

    transport_moc_solver(const pipe_properties_t& pipe,
        const vector<double>& vol_flow,
        density_viscosity_layer& prev, density_viscosity_layer& next)
        : pipe(pipe),
        Q(vol_flow),
        prev{ prev }, next{ next }
    {}

    /// @brief Алгоритм расчёта методом характеристик
    /// @param prev Ссылка на предыдущий профиль
    /// @param next Ссылка на текущий профиль
    /// @param par_in Параметр вытесняющей партии
    /// @param direction Направление течения потока
    void step(array<double, 2> par_in)
    {
        int direction = getEigenvals(0) > 0 ? 1 : -1;
        size_t start_index = direction > 0 ? 1 : (prev.density.size()) - 2;
        size_t end_index = direction < 0 ? 0 : (prev.density.size());
        next.density[start_index - direction] = par_in[0];
        next.viscosity[start_index - direction] = par_in[1];
        for (size_t index = start_index; index != end_index; index += direction)
        {
            next.density[index] = prev.density[index - direction];
            next.viscosity[index] = prev.viscosity[index - direction];
        }
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

    array<double, 2> get_par_in(double& time_moment, time_series_t& density, time_series_t& viscosity)
    {
        double rho_in = interpolation(time_moment, density);
        double visc_in = interpolation(time_moment, viscosity);

        return { rho_in, visc_in };
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
    /// @brief Проблемно-ориентированные слои
    density_viscosity_layer& prev;
    density_viscosity_layer& next;
};


struct parameters_series_t
{
    vector<time_series_t> param_series;
    time_series_t density_series;
    time_series_t viscosity_series;

    void input_parameters(vector<vector<double>>& par_time, vector<vector<double>>& par)
    {
        for (size_t index = 0; index < par.size(); index++)
        {
            vector<double> moments;
            if (par_time[index].size() == 1)
                moments = build_series(par_time[index][0], par[index].size());
            else
                moments = par_time[index];

            time_series_t series = build_parameters(moments, par[index]);
            param_series.push_back(series);
        }
    }

    void input_dens_visc(double& dt, vector<double>& density_vals, vector<double>& visc_vals)
    {
        vector<double> moments = build_series(dt, density_vals.size());
        density_series = build_parameters(moments, density_vals);
        viscosity_series = build_parameters(moments, visc_vals);
    }

    void input_parameters(double& dt, vector<vector<double>> par)
    {
        vector<double> moments = build_series(dt, par[0].size());
        for (size_t index = 0; index < par.size(); index++)
        {
            time_series_t series = build_parameters(moments, par[index]);
            param_series.push_back(series);
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