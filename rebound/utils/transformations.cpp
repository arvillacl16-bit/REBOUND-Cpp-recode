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

#include "transformations.hpp"
#include <iostream>

namespace rebound {
  namespace _transform {
    void inertial_to_jacobi_posvel(const ParticleStore& from, ParticleStore& to) {
      size_t N = from.size();
      double eta = from.mus[0];
      Vec3 s_pos = eta * from.positions[0];
      Vec3 s_vel = eta * from.velocities[0];

      for (size_t i = 1; i < N; ++i) {
        double ei = 1.0 / eta;
        double mu = from.mus[i];

        to.mus[i] = mu;
        to.positions[i] = from.positions[i] - ei * s_pos;
        to.velocities[i] = from.velocities[i] - ei * s_vel;

        if (!from.test_mass[i]) {
          eta += mu;
          double pme = eta * ei; // = eta / old_eta
          s_pos = s_pos * pme + mu * to.positions[i];
          s_vel = s_vel * pme + mu * to.velocities[i];
        }
      }
    }

    void inertial_to_jacobi_posvelacc(const ParticleStore& from, ParticleStore& to) {
      size_t N = from.size();
      double eta = from.mus[0];
      Vec3 s_pos = eta * from.positions[0];
      Vec3 s_vel = eta * from.velocities[0];
      Vec3 s_acc = eta * from.accelerations[0];

      for (size_t i = 1; i < N; ++i) {
        double ei = 1.0 / eta;
        double mu = from.mus[i];

        to.mus[i] = mu;
        to.positions[i] = from.positions[i] - ei * s_pos;
        to.velocities[i] = from.velocities[i] - ei * s_vel;
        to.accelerations[i] = from.test_mass[i] ? Vec3(0.0,0.0,0.0) : (from.accelerations[i] - ei * s_acc);

        if (!from.test_mass[i]) {
          eta += mu;
          double pme = eta * ei;
          s_pos = s_pos * pme + mu * to.positions[i];
          s_vel = s_vel * pme + mu * to.velocities[i];
          s_acc = s_acc * pme + mu * to.accelerations[i];
        }
      }
    }

    void inertial_to_jacobi_acc(const ParticleStore& from, ParticleStore& to) {
      size_t N = from.size();
      double eta = from.mus[0];
      Vec3 s_acc = eta * from.accelerations[0];

      for (size_t i = 1; i < N; ++i) {
        double ei = 1.0 / eta;
        double mu = from.mus[i];

        to.accelerations[i] = from.test_mass[i] ? Vec3(0.0,0.0,0.0) : (from.accelerations[i] - ei * s_acc);

        if (!from.test_mass[i]) {
          eta += mu;
          double pme = eta * ei;
          s_acc = s_acc * pme + mu * to.accelerations[i];
        }
      }
    }

    void jacobi_to_inertial_posvel(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      double eta = from.mus[0];
      Vec3 s_pos = from.positions[0] * eta;
      Vec3 s_vel = from.velocities[0] * eta;

      // Serial unwind from N-1 down to 1
      for (size_t i = N; i-- > 1;) {
        double ei = 1.0 / eta;
        double mu = from.mus[i];

        if (!from.test_mass[i]) {
          // subtract contribution of this massive body and recover previous s_pos
          s_pos = (s_pos - mu * from.positions[i]) * ei;
          s_vel = (s_vel - mu * from.velocities[i]) * ei;
          to.positions[i] = from.positions[i] + s_pos;
          to.velocities[i] = from.velocities[i] + s_vel;
          eta -= mu;
          s_pos *= eta;
          s_vel *= eta;
        } else {
          // test particle: simply copy
          to.positions[i] = from.positions[i];
          to.velocities[i] = from.velocities[i];
        }
      }

      if (eta != 0.0) {
        double mi = 1.0 / eta;
        to.positions[0] = s_pos * mi;
        to.velocities[0] = s_vel * mi;
      } else {
        to.positions[0] = {};
        to.velocities[0] = {};
      }
    }

    void jacobi_to_inertial_pos(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      double eta = from.mus[0];
      Vec3 s_pos = from.positions[0] * eta;

      for (size_t i = N; i-- > 1;) {
        double ei = 1.0 / eta;
        double mu = from.mus[i];

        if (!from.test_mass[i]) {
          s_pos = (s_pos - mu * from.positions[i]) * ei;
          to.positions[i] = from.positions[i] + s_pos;
          eta -= mu;
          s_pos *= eta;
        } else {
          to.positions[i] = from.positions[i];
        }
      }

      if (eta != 0.0) {
        double mi = 1.0 / eta;
        to.positions[0] = s_pos * mi;
      } else {
        to.positions[0] = Vec3(0.0,0.0,0.0);
      }
    }

    void jacobi_to_inertial_acc(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      double eta = from.mus[0];
      Vec3 s_acc = from.accelerations[0] * eta;

      for (size_t i = N; i-- > 1;) {
        double ei = 1.0 / eta;
        double mu = from.mus[i];

        if (!from.test_mass[i]) {
          s_acc = (s_acc - mu * from.accelerations[i]) * ei;
          to.accelerations[i] = from.accelerations[i] + s_acc;
          eta -= mu;
          s_acc *= eta;
        } else {
          to.accelerations[i] = from.accelerations[i];
        }
      }

      if (eta != 0.0) {
        double mi = 1.0 / eta;
        to.accelerations[0] = s_acc * mi;
      } else {
        to.accelerations[0] = Vec3(0.0,0.0,0.0);
      }
    }

    void inertial_to_democraticheliocentric_posvel(const ParticleStore& from, ParticleStore& to) {
      size_t N = from.size();

      double m0 = from.mus[0];
      double x0x = from.positions[0].x * from.mus[0];
      double x0y = from.positions[0].y * from.mus[0];
      double x0z = from.positions[0].z * from.mus[0];
      double v0x = from.velocities[0].x * from.mus[0];
      double v0y = from.velocities[0].y * from.mus[0];
      double v0z = from.velocities[0].z * from.mus[0];

      #pragma omp parallel for reduction(+: x0x, x0y, x0z, v0x, v0y, v0z, m0) schedule(static)
      for (size_t i = 1; i < N; ++i) {
        if (!from.test_mass[i]) {
          double mi = from.mus[i];
          x0x += from.positions[i].x * mi;
          x0y += from.positions[i].y * mi;
          x0z += from.positions[i].z * mi;
          v0x += from.velocities[i].x * mi;
          v0y += from.velocities[i].y * mi;
          v0z += from.velocities[i].z * mi;
          m0 += mi;
        }
      }

      // build Vec3 center-of-mass position & velocity
      Vec3 com_pos(x0x / m0, x0y / m0, x0z / m0);
      Vec3 com_vel(v0x / m0, v0y / m0, v0z / m0);

      to.positions[0] = com_pos;
      to.velocities[0] = com_vel;
      to.mus[0] = m0;

      // per-particle (embarrassingly parallel)
      #pragma omp parallel for schedule(static)
      for (size_t i = 1; i < N; ++i) {
        to.positions[i]  = from.positions[i] - com_pos;
        to.velocities[i] = from.velocities[i] - com_vel;
        to.mus[i]        = from.mus[i];
      }
    }

    void democraticheliocentric_to_inertial_pos(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      double mtot = from.mus[0];

      // accumulate x0 = sum(mi * positions_i / mtot) over real masses
      double x0x = 0.0, x0y = 0.0, x0z = 0.0;
      #pragma omp parallel for reduction(+: x0x, x0y, x0z) schedule(static)
      for (size_t i = 1; i < N; ++i) {
        if (!from.test_mass[i] && mtot != 0.0) {
          double factor = from.mus[i] / mtot;
          x0x += from.positions[i].x * factor;
          x0y += from.positions[i].y * factor;
          x0z += from.positions[i].z * factor;
          // set mass for to
        }
      }

      Vec3 x0(x0x, x0y, x0z);
      to.positions[0] = from.positions[0] - x0;

      #pragma omp parallel for schedule(static)
      for (size_t i = 1; i < N; ++i) {
        to.positions[i] = from.positions[i] + to.positions[0];
        to.mus[i] = from.mus[i];
      }
    }

    void democraticheliocentric_to_inertial_posvel(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      democraticheliocentric_to_inertial_pos(to, from);

      // velocities: per-particle addition
      #pragma omp parallel for schedule(static)
      for (size_t i = 1; i < N; ++i) {
        to.velocities[i] = from.velocities[i] + from.velocities[0];
      }

      // compute velocity correction v0 = sum(mi * vel_i / mtot)
      double mtot = from.mus[0];
      double v0x = 0.0, v0y = 0.0, v0z = 0.0;
      if (mtot != 0.0) {
        #pragma omp parallel for reduction(+: v0x, v0y, v0z) schedule(static)
        for (size_t i = 1; i < N; ++i) {
          if (!from.test_mass[i]) {
            double f = from.mus[i] / mtot;
            v0x += from.velocities[i].x * f;
            v0y += from.velocities[i].y * f;
            v0z += from.velocities[i].z * f;
          }
        }
      }
      Vec3 v0(v0x, v0y, v0z);
      to.velocities[0] = from.velocities[0] - v0;
    }

    // ------------------------------
    // WHDS Transformations (partial parallelization)
    // ------------------------------
    void inertial_to_whds_posvel(const ParticleStore& from, ParticleStore& to) {
      size_t N = from.size();

      // compute center-of-mass quantities
      double m0 = from.mus[0];
      double x0x = from.positions[0].x * from.mus[0];
      double x0y = from.positions[0].y * from.mus[0];
      double x0z = from.positions[0].z * from.mus[0];
      double v0x = from.velocities[0].x * from.mus[0];
      double v0y = from.velocities[0].y * from.mus[0];
      double v0z = from.velocities[0].z * from.mus[0];

      #pragma omp parallel for reduction(+: x0x, x0y, x0z, v0x, v0y, v0z, m0) schedule(static)
      for (size_t i = 1; i < N; ++i) {
        if (!from.test_mass[i]) {
          double mi = from.mus[i];
          x0x += from.positions[i].x * mi;
          x0y += from.positions[i].y * mi;
          x0z += from.positions[i].z * mi;
          v0x += from.velocities[i].x * mi;
          v0y += from.velocities[i].y * mi;
          v0z += from.velocities[i].z * mi;
          m0 += mi;
        }
      }

      Vec3 com_pos(x0x / m0, x0y / m0, x0z / m0);
      Vec3 com_vel(v0x / m0, v0y / m0, v0z / m0);

      to.positions[0] = com_pos;
      to.velocities[0] = com_vel;
      to.mus[0] = m0;

      // massive particles (1..n_real-1)
      size_t nreal = from.n_real();
      #pragma omp parallel for schedule(static)
      for (size_t i = 1; i < nreal; ++i) {
        double mi = from.mus[i];
        double mf = (m0 + mi) / m0; // safe because m0>0 for sane systems
        to.positions[i]  = from.positions[i] - from.positions[0];
        to.velocities[i] = mf * (from.velocities[i] - from.velocities[0]);
        to.mus[i] = mi;
      }

      // test particles (n_real .. N-1)
      #pragma omp parallel for schedule(static)
      for (size_t i = nreal; i < N; ++i) {
        to.positions[i]  = from.positions[i] - from.positions[0];
        to.velocities[i] = from.velocities[i] - com_vel;
        to.mus[i] = from.mus[i];
      }
    }

    void whds_to_inertial_pos(ParticleStore& to, const ParticleStore& from) {
      // same mapping as democraticheliocentric -> inertial pos
      democraticheliocentric_to_inertial_pos(to, from);
    }

    void whds_to_inertial_posvel(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      democraticheliocentric_to_inertial_posvel(to, from);

      double m0 = from.mus[0];

      #pragma omp parallel for schedule(static)
      for (size_t i = 1; i < N; ++i) {
        if (!from.test_mass[i]) {
          double mi = from.mus[i];
          double imf = (m0 + mi) / m0;
          to.velocities[i] = from.velocities[i] * imf + from.velocities[0];
        } else {
          to.velocities[i] = from.velocities[i] + from.velocities[0];
        }
      }

      // compute correction v0 term
      double v0x = 0.0, v0y = 0.0, v0z = 0.0;
      #pragma omp parallel for reduction(+: v0x, v0y, v0z) schedule(static)
      for (size_t i = 1; i < from.n_real(); ++i) {
        if (!from.test_mass[i]) {
          double mi = from.mus[i];
          double f = mi / (m0 + mi);
          v0x += from.velocities[i].x * f;
          v0y += from.velocities[i].y * f;
          v0z += from.velocities[i].z * f;
        }
      }
      Vec3 v0(v0x, v0y, v0z);
      to.velocities[0] = from.velocities[0] - v0;
    }

    // ------------------------------
    // Barycentric Transformations (parallelizable)
    // ------------------------------
    void barycentric_to_inertial_posvel(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();

      // initialize central scaled quantities
      to.positions[0] = from.mus[0] * from.positions[0];
      to.velocities[0] = from.mus[0] * from.velocities[0];
      to.mus[0] = from.mus[0];

      // accumulate s_pos and s_vel over real masses
      double s_pos_x = 0.0, s_pos_y = 0.0, s_pos_z = 0.0;
      double s_vel_x = 0.0, s_vel_y = 0.0, s_vel_z = 0.0;
      double s_m = 0.0;

      #pragma omp parallel for reduction(+: s_pos_x, s_pos_y, s_pos_z, s_vel_x, s_vel_y, s_vel_z, s_m) schedule(static)
      for (size_t i = 1; i < N; ++i) {
        to.positions[i] = from.positions[i] + from.positions[0];
        to.velocities[i] = from.velocities[i] + from.velocities[0];
        if (!from.test_mass[i]) {
          double m = from.mus[i];
          to.mus[i] = m;
          s_pos_x += to.positions[i].x * m;
          s_pos_y += to.positions[i].y * m;
          s_pos_z += to.positions[i].z * m;
          s_vel_x += to.velocities[i].x * m;
          s_vel_y += to.velocities[i].y * m;
          s_vel_z += to.velocities[i].z * m;
          s_m += m;
        } else {
          to.mus[i] = from.mus[i];
        }
      }

      Vec3 s_pos(s_pos_x, s_pos_y, s_pos_z);
      Vec3 s_vel(s_vel_x, s_vel_y, s_vel_z);

      to.positions[0] -= s_pos;
      to.velocities[0] -= s_vel;
      to.mus[0] -= s_m;

      if (to.mus[0] != 0.0) {
        double m0i = 1.0 / to.mus[0];
        to.positions[0] *= m0i;
        to.velocities[0] *= m0i;
      }
    }

    void barycentric_to_inertial_pos(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      to.positions[0] = from.mus[0] * from.positions[0];
      to.mus[0] = from.mus[0];

      double s_pos_x = 0.0, s_pos_y = 0.0, s_pos_z = 0.0;
      double s_m = 0.0;

      #pragma omp parallel for reduction(+: s_pos_x, s_pos_y, s_pos_z, s_m) schedule(static)
      for (size_t i = 1; i < N; ++i) {
        to.positions[i] = from.positions[i] + from.positions[0];
        if (!from.test_mass[i]) {
          double m = from.mus[i];
          to.mus[i] = m;
          s_pos_x += to.positions[i].x * m;
          s_pos_y += to.positions[i].y * m;
          s_pos_z += to.positions[i].z * m;
          s_m += m;
        } else {
          to.mus[i] = from.mus[i];
        }
      }

      Vec3 s_pos(s_pos_x, s_pos_y, s_pos_z);
      to.positions[0] -= s_pos;
      to.mus[0] -= s_m;

      if (to.mus[0] != 0.0) {
        double m0i = 1.0 / to.mus[0];
        to.positions[0] *= m0i;
      }
    }

    void barycentric_to_inertial_acc(ParticleStore& to, const ParticleStore& from) {
      size_t N = from.size();
      to.accelerations[0] = from.mus[0] * from.accelerations[0];
      to.mus[0] = from.mus[0];

      double s_acc_x = 0.0, s_acc_y = 0.0, s_acc_z = 0.0;
      double s_m = 0.0;

      #pragma omp parallel for reduction(+: s_acc_x, s_acc_y, s_acc_z, s_m) schedule(static)
      for (size_t i = 1; i < N; ++i) {
        to.accelerations[i] = from.accelerations[i] + from.accelerations[0];
        if (!from.test_mass[i]) {
          double m = from.mus[i];
          to.mus[i] = m;
          s_acc_x += to.accelerations[i].x * m;
          s_acc_y += to.accelerations[i].y * m;
          s_acc_z += to.accelerations[i].z * m;
          s_m += m;
        } else {
          to.mus[i] = from.mus[i];
        }
      }

      Vec3 s_acc(s_acc_x, s_acc_y, s_acc_z);
      to.accelerations[0] -= s_acc;
      to.mus[0] -= s_m;

      if (to.mus[0] != 0.0) {
        double m0i = 1.0 / to.mus[0];
        to.accelerations[0] *= m0i;
      }
    }

    void inertial_to_barycentric_posvel(const ParticleStore& from, ParticleStore& to) {
      size_t N = from.size();

      to.positions[0] = from.mus[0] * from.positions[0];
      to.velocities[0] = from.mus[0] * from.velocities[0];
      to.mus[0] = from.mus[0];

      double s_pos_x = 0.0, s_pos_y = 0.0, s_pos_z = 0.0;
      double s_vel_x = 0.0, s_vel_y = 0.0, s_vel_z = 0.0;
      double s_m = 0.0;

      #pragma omp parallel for reduction(+: s_pos_x, s_pos_y, s_pos_z, s_vel_x, s_vel_y, s_vel_z, s_m) schedule(static)
      for (size_t i = 1; i < N; ++i) {
        to.positions[i] = from.positions[i] - from.positions[0];
        to.velocities[i] = from.velocities[i] - from.velocities[0];
        if (!from.test_mass[i]) {
          double m = from.mus[i];
          to.mus[i] = m;
          s_pos_x += from.positions[i].x * m;
          s_pos_y += from.positions[i].y * m;
          s_pos_z += from.positions[i].z * m;
          s_vel_x += from.velocities[i].x * m;
          s_vel_y += from.velocities[i].y * m;
          s_vel_z += from.velocities[i].z * m;
          s_m += m;
        } else {
          to.mus[i] = from.mus[i];
        }
      }

      Vec3 s_pos(s_pos_x, s_pos_y, s_pos_z);
      Vec3 s_vel(s_vel_x, s_vel_y, s_vel_z);

      to.positions[0] += s_pos;
      to.velocities[0] += s_vel;
      to.mus[0] += s_m;

      if (to.mus[0] != 0.0) {
        double m0i = 1.0 / to.mus[0];
        to.positions[0] *= m0i;
        to.velocities[0] *= m0i;
      }
    }

  } // namespace _transform
} // namespace rebound
