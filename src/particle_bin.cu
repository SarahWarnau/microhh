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

#include "grid.h"
#include "fields.h"
#include "particle_bin.h"
#include "tools.h"

namespace
{
    template<typename TF> __global__
    void settle_particles_g(
            TF* const __restrict__ st,
            const TF* const __restrict__ s,
            const TF* const __restrict__ dzhi,
            const TF w_particles,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        // Simple upwind advection, like in subsidence.
        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jstride + k*kstride;
            st[ijk] -= w_particles * (s[ijk+kstride]-s[ijk])*dzhi[k+1];
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Particle_bin<TF>::exec(Stats<TF>& stats)
{
    if (!sw_particle)
        return;

    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU(gridi, gridj, gd.ktot);
    dim3 blockGPU(blocki, blockj, 1);

    for (auto& w : w_particle)
        settle_particles_g<TF><<<gridGPU, blockGPU>>>(
                fields.st.at(w.first)->fld_g,
                fields.sp.at(w.first)->fld_g,
                gd.dzhi_g,
                w.second,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

    cuda_check_error();
}
#endif

#ifdef FLOAT_SINGLE
template class Particle_bin<float>;
#else
template class Particle_bin<double>;
#endif
