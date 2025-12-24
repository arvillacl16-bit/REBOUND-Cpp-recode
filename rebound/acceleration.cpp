/*
 * This file is part of the C++ translation of REBOUND.
 *
 * Original REBOUND (C) code by Hanno Rein and others.
 * This translation is licensed under the GNU General Public License v3 or later.
 *
 * REBOUND is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * REBOUND is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#include "rebound.hpp"
#include <vector>
#include <cmath>

namespace rebound {
  namespace _accel {
    std::vector<Vec3> temp_acc_flat; // flattened: n_threads * N

    void reset_temp_acc(size_t N, int n_threads) {
      if (temp_acc_flat.size() != n_threads * N) temp_acc_flat.assign(n_threads * N, Vec3{});
      else std::fill(temp_acc_flat.begin(), temp_acc_flat.end(), Vec3{});
    }

    void calc_accel_none(ParticleStore& particles) {
      size_t N = particles.size();
#pragma omp parallel for
      for (size_t i = 0; i < N; ++i) particles.accelerations[i] = {};
    }

    void calc_accel_basic(ParticleStore& particles, double softening2) {
      const size_t N = particles.size();

      // Zero particle accelerations
      #pragma omp parallel for
      for (size_t i = 0; i < N; ++i)
        particles.accelerations[i] = Vec3{};

      int n_threads = 1;
    #ifdef _OPENMP
      n_threads = omp_get_max_threads();
    #endif

      // Reset flattened temp accumulator
      reset_temp_acc(N, n_threads);

      #pragma omp parallel
      {
        int tid = 0;
  #ifdef _OPENMP
        tid = omp_get_thread_num();
  #endif

        #pragma omp for schedule(static)
        for (size_t i = 0; i < N; ++i) {
          if (particles.test_mass[i]) continue;

          for (size_t j = i + 1; j < N; ++j) {
            if (particles.test_mass[i] && particles.test_mass[j]) continue;

            Vec3 dx = particles.positions[j] - particles.positions[i];
            double r2 = dx.mag2() + softening2;
            double inv_r3 = 1.0 / (r2 * std::sqrt(r2));

            temp_acc_flat[tid * N + i] += dx * particles.mus[j] * inv_r3;
            temp_acc_flat[tid * N + j] -= dx * particles.mus[i] * inv_r3;
          }
        }
      }

      // Reduce thread-local sums
      #pragma omp parallel for
      for (size_t i = 0; i < N; ++i) {
        Vec3 acc{};
        for (int t = 0; t < n_threads; ++t)
          acc += temp_acc_flat[t * N + i];
        particles.accelerations[i] += acc;
      }
    }

    void calc_accel_jacobi(ParticleStore& particles, double softening2) {
      const size_t N = particles.size();

      // Zero particle accelerations
      #pragma omp parallel for
      for (size_t i = 0; i < N; ++i)
        particles.accelerations[i] = Vec3{};

      int n_threads = 1;
    #ifdef _OPENMP
      n_threads = omp_get_max_threads();
    #endif

      // Reset flattened temp accumulator
      reset_temp_acc(N, n_threads);

      #pragma omp parallel
      {
            int tid = 0;
    #ifdef _OPENMP
            tid = omp_get_thread_num();
    #endif

          #pragma omp for schedule(static)
      for (size_t i = 1; i < N; ++i) {
          for (size_t j = 0; j < i; ++j) {
            if (particles.test_mass[i] && particles.test_mass[j]) continue;

            Vec3 dx = particles.positions[j] - particles.positions[i];
            double r2 = dx.mag2() + softening2;
            double inv_r3 = 1.0 / (r2 * std::sqrt(r2));

            temp_acc_flat[tid * N + i] += dx * particles.mus[j] * inv_r3;
            temp_acc_flat[tid * N + j] -= dx * particles.mus[i] * inv_r3;
          }
        }
      }

      // Reduce thread-local sums
      #pragma omp parallel for
      for (size_t i = 0; i < N; ++i) {
        Vec3 acc{};
        for (int t = 0; t < n_threads; ++t)
          acc += temp_acc_flat[t * N + i];
        particles.accelerations[i] += acc;
      }
    }

    void calc_accel_mercurius(ParticleStore &particles, double softening2, Mercurius &settings) {
      calc_accel_basic(particles, softening2);
    }

    void calc_accel_compensated(ParticleStore &particles, double softening2) {
      size_t N = particles.size();
      if (particles.gravity_cs.size() < N) particles.gravity_cs.resize(N, Vec3{});

      auto &cs = particles.gravity_cs;
      auto &test_mass = particles.test_mass;
#pragma omp parallel for
      for (size_t i = 0; i < N; ++i) {
        particles.accelerations[i] = {};
        cs[i] = {};
      }

#pragma omp parallel for (guided)
      for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
          if (test_mass[j]) continue;
          if (i == j) continue;
          Vec3 dr = particles.positions[i] - particles.positions[j];
          double r2 = dr.mag2() + softening2;
          double inv_r3 = 1. / (r2 * std::sqrt(r2));
          Vec3 delta = -particles.mus[j] * inv_r3 * dr;
          Vec3 y = delta - cs[i];
          Vec3 t = particles.accelerations[i] + y;
          cs[i] = (t - particles.accelerations[i]) - y;
          particles.accelerations[i] = t;
        }
      }
    }
  }
}
