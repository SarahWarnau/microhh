/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BOUNDARY_SURFACE_SOLAR_H
#define BOUNDARY_SURFACE_SOLAR_H

#include "boundary.h"
#include "stats.h"

template<typename> class Diff;

template<typename TF>
class Boundary_surface_solar : public Boundary<TF>
{
    public:
        Boundary_surface_solar(Master&, Grid<TF>&, Soil_grid<TF>&, Fields<TF>&, Input&);
        ~Boundary_surface_solar();

        void init(Input&, Thermo<TF>&, const Sim_mode);
        void create_cold_start(Netcdf_handle&);
        void create(Input&, Netcdf_handle&, Stats<TF>&, Column<TF>&, Cross<TF>&, Timeloop<TF>&);
        void set_values();

        // // Solar evaporator surface values
        // void get_surface_values_solar(
        //     TF* thl_fld_bot, //(potential temperature at the surface)
        //     TF* qt_fld_bot, //(specific humidity at the surface)
        //     TF* rnet_fld_bot, //(net radiation at the surface)
        //     TF* thl_fld, //(temperature of the atmosphere just above the surface)
        //     TF* qt_fld, //(specific humidity just above the surface)
        //     TF* rho_air_fld, //(air density previous timestep)
        //     TF* p_surf_fld, //(surface pressure)
        //     TF* ra_fld, //(aerodynamic resistance, keep constant for now)
        //     TF* rs_fld, //(surface resistance, constant, set to 0)
        //     const TF epsilon, //(emissivity of the surface)
        //     const TF cp, //(specific heat capacity of air)
        //     const TF rd, //(specific gas constant for dry air)
        //     const TF rv, //(specific gas constant for water vapor)
        //     const TF sigma, //(Stefan-Boltzmann constant)
        //     const TF lv, //(latent heat of vaporization)
        //     const int istart, const int iend,
        //     const int jstart, const int jend,
        //     const int icells, const int jcells
        // );

        const std::vector<TF>& get_z0m()  const { return z0m; };
        const std::vector<TF>& get_dudz() const { return dudz_mo; }
        const std::vector<TF>& get_dvdz() const { return dvdz_mo; }
        const std::vector<TF>& get_dbdz() const { return dbdz_mo; }

        void exec(Thermo<TF>&, Radiation<TF>&, Microphys<TF>&, Timeloop<TF>&);
        void exec_stats(Stats<TF>&);
        void exec_column(Column<TF>&);
        void exec_cross(Cross<TF>&, unsigned long);

        void load(const int, Thermo<TF>&);
        void save(const int, Thermo<TF>&);

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device(Thermo<TF>&);
        void forward_device(Thermo<TF>&);
        void backward_device(Thermo<TF>&);
        void clear_device(Thermo<TF>&);

        cuda_vector<TF>& get_z0m_g()  { return z0m_g; };
        cuda_vector<TF>& get_dudz_g() { return dudz_mo_g; };
        cuda_vector<TF>& get_dvdz_g() { return dvdz_mo_g; };
        cuda_vector<TF>& get_dbdz_g() { return dbdz_mo_g; };
        #endif

    protected:
        void process_input(Input&, Thermo<TF>&); // Process and check the surface input
        void init_surface(Input&, Thermo<TF>&); // Allocate and initialize the surface arrays
        void init_solver(); // Prepare the lookup table's for the surface layer solver
        void set_ustar(); // Set fixed ustar

    private:
        using Boundary<TF>::master;
        using Boundary<TF>::grid;
        using Boundary<TF>::soil_grid;
        using Boundary<TF>::fields;
        using Boundary<TF>::boundary_cyclic;
        using Boundary<TF>::swboundary;
        using Boundary<TF>::field3d_io;

        using Boundary<TF>::process_bcs;

        using Boundary<TF>::mbcbot;
        using Boundary<TF>::ubot;
        using Boundary<TF>::vbot;

        typedef std::map<std::string, Field3dBc<TF>> BcMap;
        using Boundary<TF>::sbc;

        TF ustarin;

        // Sarah: put parameter arrays here, for first tests
        std::vector<TF> rnetin;
        std::vector<TF> ra;
        std::vector<TF> rs;

        std::vector<float> zL_sl;
        std::vector<float> f_sl;
        std::vector<int> nobuk;

        std::vector<TF> z0m;
        std::vector<TF> z0h;

        std::vector<TF> ustar;
        std::vector<TF> obuk;

        std::vector<TF> dudz_mo;
        std::vector<TF> dvdz_mo;
        std::vector<TF> dbdz_mo;

        #ifdef USECUDA
        cuda_vector<TF> z0m_g;
        cuda_vector<TF> z0h_g;
        cuda_vector<TF> obuk_g;
        cuda_vector<TF> ustar_g;

        cuda_vector<TF> dudz_mo_g;
        cuda_vector<TF> dvdz_mo_g;
        cuda_vector<TF> dbdz_mo_g;

        int* nobuk_g = nullptr;

        float* zL_sl_g = nullptr;
        float* f_sl_g = nullptr;
        #endif

        Boundary_type thermobc;
        bool sw_constant_z0;

        bool sw_charnock;
        TF alpha_m;
        TF alpha_ch;
        TF alpha_h;

    protected:
        // Cross sections
        std::vector<std::string> cross_list;         // List of active cross variables

        void update_slave_bcs();
};
#endif
