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
    void calc_accel_none(_ParticleStore& particles) {
      const size_t N = particles.size();
#pragma omp parallel for
      for (size_t i = 0; i < N; ++i) {
        particles.accelerations[i] = {};
      }
    }

    void calc_accel_basic(_ParticleStore& particles, double softening2) {
      const size_t N = particles.size();

      // Zero accelerations
#pragma omp parallel for
      for (size_t i = 0; i < N; ++i) particles.accelerations[i] = Vec3{};

      int n_threads = 1;
#ifdef _OPENMP
      n_threads = omp_get_max_threads();
#endif

      // Thread-local accumulators
      std::vector<std::vector<Vec3>> temp_acc(n_threads, std::vector<Vec3>(N, Vec3()));

#pragma omp parallel
      {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic)
      for (size_t i = 0; i < N; ++i) {
        if (particles.test_mass[i]) continue; // skip test mass as source?

        for (size_t j = i + 1; j < N; ++j) {
          if (particles.test_mass[i] && particles.test_mass[j]) continue; // ignore both test masses

          Vec3 dx = particles.positions[j] - particles.positions[i];
          double r2 = dx.mag2() + softening2;
          double inv_r3 = 1.0 / (r2 * std::sqrt(r2));

          Vec3 a_i = dx * particles.mus[j] * inv_r3;
          Vec3 a_j = -dx * particles.mus[i] * inv_r3;

          temp_acc[tid][i] += a_i;
          temp_acc[tid][j] += a_j;
          }
        }
      }

      for (int t = 0; t < n_threads; ++t)
        for (size_t i = 0; i < N; ++i)
          particles.accelerations[i] += temp_acc[t][i];
      }

    void calc_accel_jacobi(_ParticleStore& particles, double softening2) {
      const size_t N = particles.size();

#pragma omp parallel for
      for (size_t i = 0; i < N; ++i) particles.accelerations[i] = Vec3{};

      int n_threads = 1;
#ifdef _OPENMP
      n_threads = omp_get_max_threads();
#endif

      std::vector<std::vector<Vec3>> temp_acc(n_threads, std::vector<Vec3>(N, Vec3{}));

#pragma omp parallel
      {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic)
      for (size_t i = 1; i < N; ++i) {
        for (size_t j = 0; j < i; ++j) {
          if (particles.test_mass[i] && particles.test_mass[j]) continue; // ignore both test masses

          Vec3 dx = particles.positions[j] - particles.positions[i];
          double r2 = dx.mag2() + softening2;
          double inv_r3 = 1.0 / (r2 * std::sqrt(r2));

          Vec3 a_i = dx * particles.mus[j] * inv_r3;
          Vec3 a_j = -dx * particles.mus[i] * inv_r3;

          temp_acc[tid][i] += a_i;
          temp_acc[tid][j] += a_j;
          }
        }
      }

      for (int t = 0; t < n_threads; ++t)
        for (size_t i = 0; i < N; ++i)
          particles.accelerations[i] += temp_acc[t][i];
      }

    void calc_accel_compensated(_ParticleStore& particles, double softening2) {
      const size_t N = particles.size();

      if (particles.gravity_cs.size() < N)
        particles.gravity_cs.resize(N, Vec3{ 0.0, 0.0, 0.0 });

      auto& cs = particles.gravity_cs;

#pragma omp parallel for
      for (size_t i = 0; i < N; ++i) {
        particles.accelerations[i] = Vec3{ 0.0, 0.0, 0.0 };
        cs[i] = Vec3{ 0.0, 0.0, 0.0 };
        }

#pragma omp parallel for schedule(guided)
      for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
          if (i == j) continue;
          if (particles.test_mass[j]) continue; // skip contribution from test-mass particles

          Vec3 dr = particles.positions[i] - particles.positions[j];
          double r2 = dr.mag2() + softening2;
          double inv_r3 = 1.0 / (r2 * std::sqrt(r2));
          Vec3 delta = -particles.mus[j] * inv_r3 * dr;

          Vec3 y = delta - cs[i];
          Vec3 t = particles.accelerations[i] + y;
          cs[i] = (t - particles.accelerations[i]) - y;
          particles.accelerations[i] = t;
          }
        }
      }

    void calc_accel_mercurius(_ParticleStore& particles, double softening2, _MercuriusSettings settings) {
      calc_accel_basic(particles, softening2);
      }
    } // namespace _accel
  } // namespace rebound
