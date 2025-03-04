/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
 * Copyright (c) 2022-2022 Stijn Heldens
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

#include <utility>
#include "cuda_buffer.h"

#ifdef USECUDA
#include <cuda_runtime_api.h>
#include "tools.h"
#endif //USECUDA

cuda_raw_buffer::cuda_raw_buffer(size_t size_in_bytes)
{
    reallocate(size_in_bytes);
}

cuda_raw_buffer::cuda_raw_buffer(cuda_raw_buffer&& that) noexcept
{
    *this = std::move(that);
}

cuda_raw_buffer& cuda_raw_buffer::operator=(cuda_raw_buffer&& that) noexcept
{
    std::swap(ptr_, that.ptr_);
    return *this;
}

#ifdef USECUDA
static std::string allocation_failed_message(size_t allocation_attempt_size)
{
    size_t free = 0, total = 0;
    int device = -1;

    // Intentionally do not check for errors. In the worst case, the variables will remain 0.
    cudaMemGetInfo(&free, &total);
    cudaGetDevice(&device);

    return std::string("memory allocation failed on device ") + std::to_string(device) +
        ": attempted to allocate " + std::to_string(allocation_attempt_size) + " bytes, but only " +
        std::to_string(free) + "bytes are available.";
}
#endif //USECUDA

void cuda_raw_buffer::reallocate(size_t size_in_bytes)
{
#ifdef USECUDA
    if (ptr_)
    {
        cuda_safe_call(cudaFree(ptr_));
        ptr_ = nullptr;
    }

    if (size_in_bytes > 0)
    {
        cudaError_t result = cudaMalloc(&ptr_, size_in_bytes);

        if (result == cudaErrorMemoryAllocation)
        {
            throw Tools_g::cuda_exception(result, allocation_failed_message(size_in_bytes));
        }

        cuda_safe_call(result);
    }
#else
    if (size_in_bytes > 0)
    {
        throw std::runtime_error("CUDA is not enabled, allocating GPU memory is not possible");
    }
#endif //USECUDA
}

void cuda_raw_copy(const void* src, void* dst, size_t nbytes)
{
#ifdef USECUDA
    cuda_safe_call(cudaMemcpy(dst, src, nbytes, cudaMemcpyDefault));
#else
    if (nbytes > 0)
    {
        throw std::runtime_error("CUDA is not enabled, copying GPU memory is not possible");
    }
#endif //USECUDA
}
