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

#include "integrator.hpp"
#include "transformations.hpp"
#include <iostream>
#include <cassert>

namespace rebound {
  namespace _whfast {
#pragma region WHFAST Corrector Coefficients
    constexpr double corrector_a_1 = 0.41833001326703777398908601289259374469640768464934;
    constexpr double corrector_a_2 = 0.83666002653407554797817202578518748939281536929867;
    constexpr double corrector_a_3 = 1.2549900398011133219672580386777812340892230539480;
    constexpr double corrector_a_4 = 1.6733200530681510959563440515703749787856307385973;
    constexpr double corrector_a_5 = 2.0916500663351888699454300644629687234820384232467;
    constexpr double corrector_a_6 = 2.5099800796022266439345160773555624681784461078960;
    constexpr double corrector_a_7 = 2.9283100928692644179236020902481562128748537925454;
    constexpr double corrector_a_8 = 3.3466401061363021919126881031407499575712614771947;
    constexpr double corrector_b_31 = -0.024900596027799867499350357910273437184309981229127;
    constexpr double corrector_b_51 = -0.0083001986759332891664501193034244790614366604097090;
    constexpr double corrector_b_52 = 0.041500993379666445832250596517122395307183302048545;
    constexpr double corrector_b_71 = 0.0024926811426922105779030593952776964450539008582219;
    constexpr double corrector_b_72 = -0.018270923246702131478062356884535264841652263842597;
    constexpr double corrector_b_73 = 0.053964399093127498721765893493510877532452806339655;
    constexpr double corrector_b_111 = 0.00020361579647854651301632818774633716473696537436847;
    constexpr double corrector_b_112 = -0.0023487215292295354188307328851055489876255097419754;
    constexpr double corrector_b_113 = 0.012309078592019946317544564763237909911330686448336;
    constexpr double corrector_b_114 = -0.038121613681288650508647613260247372125243616270670;
    constexpr double corrector_b_115 = 0.072593394748842738674253180742744961827622366521517;
    constexpr double corrector_b_178 = 0.093056103771425958591541059067553547100903397724386;
    constexpr double corrector_b_177 = -0.065192863576377893658290760803725762027864651086787;
    constexpr double corrector_b_176 = 0.032422198864713580293681523029577130832258806467604;
    constexpr double corrector_b_175 = -0.012071760822342291062449751726959664253913904872527;
    constexpr double corrector_b_174 = 0.0033132577069380655655490196833451994080066801611459;
    constexpr double corrector_b_173 = -0.00063599983075817658983166881625078545864140848560259;
    constexpr double corrector_b_172 = 0.000076436355227935738363241846979413475106795392377415;
    constexpr double corrector_b_171 = -0.0000043347415473373580190650223498124944896789841432241;
    constexpr double corrector2_b = 0.03486083443891981449909050107438281205803;
#pragma endregion

    constexpr double inv_factorial[35] = { 1., 1., 1. / 2., 1. / 6., 1. / 24., 1. / 120.,
                                          1. / 720., 1. / 5040., 1. / 40320., 1. / 362880., 1. / 3628800., 1. / 39916800.,
                                          1. / 479001600., 1. / 6227020800., 1. / 87178291200., 1. / 1307674368000.,
                                          1. / 20922789888000., 1. / 355687428096000., 1. / 6402373705728000.,
                                          1. / 121645100408832000., 1. / 2432902008176640000.,
                                          1. / 51090942171709440000., 1. / 1124000727777607680000.,
                                          1. / 25852016738884976640000., 1. / 620448401733239439360000.,
                                          1. / 15511210043330985984000000., 1. / 403291461126605635584000000.,
                                          1. / 10888869450418352160768000000., 1. / 304888344611713860501504000000.,
                                          1. / 8841761993739701954543616000000., 1. / 265252859812191058636308480000000.,
                                          1. / 8222838654177922817725562880000000., 1. / 263130836933693530167218012160000000.,
                                          1. / 8683317618811886495518194401280000000., 1. / 295232799039604140847618609643520000000. };

    inline double fast_abs(double x) { return (x > 0.) ? x : -x; }

    inline void stumpff_cs(double* cs, double z) {
      unsigned int n = 0;
      while (fast_abs(z) > 0.1) {
        z = z / 4.;
        n++;
      }
      const int nmax = 15;
      double c_odd = inv_factorial[nmax];
      double c_even = inv_factorial[nmax - 1];
      for (int np = nmax - 2; np >= 5; np -= 2) {
        c_odd = inv_factorial[np] - z * c_odd;
        c_even = inv_factorial[np - 1] - z * c_even;
      }
      cs[5] = c_odd;
      cs[4] = c_even;
      cs[3] = inv_factorial[3] - z * cs[5];
      cs[2] = inv_factorial[2] - z * cs[4];
      cs[1] = inv_factorial[1] - z * cs[3];
      for (; n > 0; n--) {
        z = z * 4.;
        cs[5] = (cs[5] + cs[4] + cs[3] * cs[2]) * 0.0625;
        cs[4] = (1. + cs[1]) * cs[3] * 0.125;
        cs[3] = 1. / 6. - z * cs[5];
        cs[2] = 0.5 - z * cs[4];
        cs[1] = 1. - z * cs[3];
      }
    }

    inline void stumpff_cs3(double* cs, double z) {
      unsigned int n = 0;
      while (fast_abs(z) > 0.1) {
        z = z / 4.;
        n++;
      }

      const int nmax = 13;
      double c_odd = inv_factorial[nmax];
      double c_even = inv_factorial[nmax - 1];
      for (int np = nmax - 2; np >= 3; np -= 2) {
        c_odd = inv_factorial[np] - z * c_odd;
        c_even = inv_factorial[np - 1] - z * c_even;
      }

      cs[3] = c_odd;
      cs[2] = c_even;
      cs[1] = inv_factorial[1] - z * c_odd;
      cs[0] = inv_factorial[0] - z * c_even;
      for (; n > 0; n--) {
        cs[3] = (cs[2] + cs[0] * cs[3]) * 0.25;
        cs[2] = cs[1] * cs[1] * 0.5;
        cs[1] = cs[0] * cs[1];
        cs[0] = 2. * cs[0] * cs[0] - 1.;
      }
    }

    inline void stiefel_Gs(double* Gs, double beta, double X) {
      double X2 = X * X;
      stumpff_cs(Gs, beta * X2);
      Gs[1] *= X;
      Gs[2] *= X2;
      double _pow = X2 * X;
      Gs[3] *= _pow;
      _pow *= X;
      Gs[4] *= _pow;
      _pow *= X;
      Gs[5] *= _pow;
      return;
    }

    inline void stiefel_Gs3(double* Gs, double beta, double X) {
      double X2 = X * X;
      stumpff_cs3(Gs, beta * X2);
      Gs[1] *= X;
      Gs[2] *= X2;
      Gs[3] *= X2 * X;
      return;
    }

    constexpr unsigned int WHFAST_NMAX_QUART = 64;
    constexpr unsigned int WHFAST_NMAX_NEWT = 32;

    inline void kepler_solver(const WHFast&, ParticleStore& p_j, double M, size_t i, double dt) {
      Particle p1 = p_j[i];

      double r0 = p1.pos().mag();
      double r0i = 1. / r0;
      double v2 = p1.vel().mag2();
      double eta0 = p1.pos().dot(p1.vel());
      double beta = 2. * M / r0 - v2;
      double zeta0 = M - beta * r0;
      double X;
      double Gs[6];
      double invperiod = 0.;
      double X_per_period = nan("");

      if (beta > 0.) {
        double sqrt_beta = std::sqrt(beta);
        invperiod = sqrt_beta * beta / (2 * M_PI * M);
        X_per_period = 2 * M_PI / sqrt_beta;

        double dtr0i = dt * r0i;
        X = dtr0i * (1. - dtr0i * eta0 * 0.5 * r0i);
      } else
        X = 0.;

      unsigned int converged = 0;
      double oldX = X;

      stiefel_Gs3(Gs, beta, X);
      double eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
      double ri = 1. / (r0 + eta0Gs1zeta0Gs2);
      X = ri * (X * eta0Gs1zeta0Gs2 - eta0 * Gs[2] - zeta0 * Gs[3] + dt);

      if (fast_abs(X - oldX) > 0.01 * X_per_period) {
        X = beta * dt / M;
        double prevX[WHFAST_NMAX_QUART + 1];
        for (size_t nlag = 1; nlag < WHFAST_NMAX_QUART; nlag++) {
          stiefel_Gs3(Gs, beta, X);
          double f = r0 * X + eta0 * Gs[2] + zeta0 * Gs[3] - dt;
          double denom = eta0 * Gs[1] + zeta0 * Gs[2] + r0;
          X = (X * denom - 5. * f) / denom;
          for (size_t i = 1; i < nlag; i++) {
            if (X == prevX[i]) {
              converged = 1;
              nlag = WHFAST_NMAX_QUART;
              break;
            }
          }
          prevX[nlag] = X;
        }

        double eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
        ri = 1. / (r0 + eta0Gs1zeta0Gs2);
      } else {
        double oldX2 = nan("");
        for (size_t n_hg = 1; n_hg < WHFAST_NMAX_NEWT; n_hg++) {
          oldX2 = oldX;
          oldX = X;
          stiefel_Gs3(Gs, beta, X);
          double eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
          ri = 1. / (r0 + eta0Gs1zeta0Gs2);
          X = ri * (X * eta0Gs1zeta0Gs2 - eta0 * Gs[2] - zeta0 * Gs[3] + dt);

          if (X == oldX || X == oldX2) {
            converged = 1;
            break;
          }
        }
      }

      if (converged == 0) {
        double X_min, X_max;
        if (beta > 0.) {
          X_min = X_per_period * std::floor(dt * invperiod);
          X_max = X_min + X_per_period;
        } else {
          double h2 = r0 * r0 * v2 - eta0 * eta0;
          double q = h2 / M / (1. + std::sqrt(1. - h2 * beta / (M * M)));
          double vq = std::copysign(std::sqrt(h2) / q, dt);
          X_min = dt / (fast_abs(vq * dt) + r0);
          X_max = dt / q;
          if (dt < 0.) {
            double temp = X_min;
            X_min = X_max;
            X_max = temp;
          }
        }

        X = (X_max + X_min) / 2.;
        do {
          stiefel_Gs3(Gs, beta, X);
          double s = r0 * X + eta0 * Gs[2] + zeta0 * Gs[3] - dt;
          if (s >= 0.) {
            X_max = X;
          } else {
            X_min = X;
          }
        } while (fast_abs(X_max - X_min) > 1e-12);
      }

      if (isnan(ri)) {
        ri = 0.;
        Gs[1] = 0.;
        Gs[2] = 0.;
        Gs[3] = 0.;
      }

      double f = -M * Gs[2] * r0i;
      double g = dt - M * Gs[3];
      double fd = -M * Gs[1] * r0i * ri;
      double gd = -M * Gs[2] * ri;

      p_j[i].pos() += f * p1.pos() + g * p1.vel();
      p_j[i].vel() += fd * p1.pos() + gd * p1.vel();
    }

    void interaction_step(ParticleStore& particles, double dt, double softening2, WHFast& settings, GravityMethod method) {
      auto& p_j = settings.internals.p_jh;
      ParticleStore cpy = p_j;
      double m0 = particles.mus[0];
      size_t N = particles.size();
      switch (settings.coordinates) {
        case WHFast::Coordinates::JACOBI:
        {
          _transform::inertial_to_jacobi_acc(particles, p_j);
          double eta = m0;
          for (size_t i = 1; i < N; ++i) {
            Particle pji = p_j[i];
            if (!particles.test_mass[i])
              eta += pji.mu();
            p_j.velocities[i] += pji.acc() * dt;
            if (method != GravityMethod::JACOBI) {
              if (i > 1) {
                double rj2i = 1. / (pji.pos().mag2() + softening2);
                double rji = std::sqrt(rj2i);
                double rj3iM = rji * rj2i * eta;
                p_j.velocities[i] += dt * rj3iM * pji.pos();
              }
            }
          }
          break;
        }
        case WHFast::Coordinates::DEMOCRATIC_HELIOCENTRIC:
        {
#pragma omp parallel for
          for (size_t i = 1; i < N; ++i) {
            if (!particles.test_mass[i])
              p_j.velocities[i] += dt * p_j.accelerations[i];
          }
          break;
        }
        case WHFast::Coordinates::WHDS:
        {
#pragma omp parallel for
          for (size_t i = 1; i < N; ++i) {
            if (!particles.test_mass[i]) {
              double mi = particles.mus[i];
              p_j.velocities[i] += dt * (m0 + mi) * particles.accelerations[i] / m0;
            } else
              p_j.velocities[i] += dt * particles.accelerations[i];
          }
          break;
        }
        case WHFast::Coordinates::BARYCENTRIC:
        {
          for (size_t i = 1; i < N; ++i) {
            if (!particles.test_mass[i]) {
              double dr = p_j.positions[i].mag();
              double prefac = p_j.mus[0] / (dr * dr * dr);
              p_j.velocities[i] += dt * (prefac * p_j.positions[i] + p_j.accelerations[i]);
            }
          }
          break;
        }
      }

      LOG((p_j == cpy));
    }

    void jump_step(ParticleStore& particles, WHFast& settings, double dt) {
      ParticleStore& p_h = settings.internals.p_jh;
      ParticleStore cpy = p_h;
      size_t N = particles.size();
      double m0 = particles.mus[0];
      switch (settings.coordinates) {
        case WHFast::Coordinates::DEMOCRATIC_HELIOCENTRIC:
        {
          Vec3 p;
#pragma omp parallel for reduction(+ : p)
          for (size_t i = 1; i < N; ++i) {
            if (particles.test_mass[i])
              continue;
            p += particles.mus[i] * p_h.velocities[i];
          }

          for (size_t i = 1; i < N; ++i)
            p_h.positions[i] += dt * p / m0;
          break;
        }
        case WHFast::Coordinates::WHDS:
        {
          Vec3 p;
#pragma omp parallel for reduction(+ : p)
          for (size_t i = 1; i < N; ++i) {
            if (particles.test_mass[i])
              continue;
            double mu = particles.mus[i];
            p += mu * p_h.velocities[i] / (m0 + mu);
          }
#pragma omp parallel for
          for (size_t i = 1; i < N; ++i) {
            if (particles.test_mass[i]) p_h.positions[i] += p;
            else {
              double mu = particles.mus[i];
              p_h.positions[i] += dt * (p - mu / (m0 + mu) * p_h.velocities[i]);
            }
          }
          break;
        }
        default:
          break;
      }

      LOG((p_h == cpy));
    }

    void kepler_step(ParticleStore& particles, WHFast& settings, double dt) {
      double m0 = particles.mus[0];
      ParticleStore& p_j = settings.internals.p_jh;
      ParticleStore cpy = p_j;
      double eta = m0;
      size_t N = particles.size();
      switch (settings.coordinates) {
        case WHFast::Coordinates::JACOBI:
        {
          size_t N = particles.size();
          std::vector<double> etas(N);
          double running_eta = p_j.mus[0];
          for (size_t i = 1; i < N; ++i) {
            etas[i] = running_eta;
            if (!particles.test_mass[i])
              running_eta += p_j.mus[i];
          }
#pragma omp parallel for
          for (size_t i = 1; i < N; ++i)
            kepler_solver(settings, p_j, etas[i], i, dt);
          break;
        }
        case WHFast::Coordinates::DEMOCRATIC_HELIOCENTRIC:
        {
#pragma omp parallel for
          for (size_t i = 1; i < N; ++i)
            kepler_solver(settings, p_j, eta, i, dt);
          break;
        }
        case WHFast::Coordinates::WHDS:
        {
          for (size_t i = 1; i < N; ++i) {
            if (particles.test_mass[i])
              eta = m0;
            else
              eta = m0 + p_j.mus[i];
            kepler_solver(settings, p_j, eta, i, dt);
          }
          break;
        }
        case WHFast::Coordinates::BARYCENTRIC:
        {
          eta = p_j.mus[0];
          for (size_t i = 1; i < N; ++i)
            kepler_solver(settings, p_j, eta, i, dt);
          break;
        }
      }

      LOG((p_j == cpy));
    }

    void com_step(ParticleStore& p_j, double dt) { p_j.positions[0] += dt * p_j.velocities[0]; }

    void update_accel(ParticleStore& particles, GravityMethod method, double softening2) {
      switch (method) {
        case GravityMethod::BASIC:
          _accel::calc_accel_basic(particles, softening2);
          break;
        case GravityMethod::COMPENSATED:
          _accel::calc_accel_compensated(particles, softening2);
          break;
        case GravityMethod::JACOBI:
          _accel::calc_accel_jacobi(particles, softening2);
          break;
        case GravityMethod::NONE:
          _accel::calc_accel_none(particles);
          break;
        default:
          break;
      }
    }

    void corrector_Z(ParticleStore& particles, WHFast& settings, double a, double b) {
      switch (settings.coordinates) {
        case WHFast::Coordinates::JACOBI:
        {
          kepler_step(particles, settings, a);
          _transform::jacobi_to_inertial_pos(particles, settings.internals.p_jh);
          update_accel(particles, settings.gravity_method, settings.softening2);
          interaction_step(particles, -b, settings.softening2, settings, settings.gravity_method);
          kepler_step(particles, settings, -a * 2);
          _transform::jacobi_to_inertial_pos(particles, settings.internals.p_jh);
          update_accel(particles, settings.gravity_method, settings.softening2);
          interaction_step(particles, b, settings.softening2, settings, settings.gravity_method);
          kepler_step(particles, settings, a);
          break;
        }
        case WHFast::Coordinates::BARYCENTRIC:
        {
          kepler_step(particles, settings, a);
          _transform::barycentric_to_inertial_pos(particles, settings.internals.p_jh);
          update_accel(particles, settings.gravity_method, settings.softening2);
          interaction_step(particles, -b, settings.softening2, settings, settings.gravity_method);
          kepler_step(particles, settings, -a * 2);
          _transform::barycentric_to_inertial_pos(particles, settings.internals.p_jh);
          update_accel(particles, settings.gravity_method, settings.softening2);
          interaction_step(particles, b, settings.softening2, settings, settings.gravity_method);
          kepler_step(particles, settings, a);
          break;
        }
        case WHFast::Coordinates::WHDS:
        case WHFast::Coordinates::DEMOCRATIC_HELIOCENTRIC:
          std::cerr << "Coordinate system not supported" << '\n';
          break;
      }
    }

    void apply_corrector(ParticleStore& particles, WHFast& settings, double inv, double dt) {
      switch (settings.order) {
        case WHFast::Order::THIRD:
        {
          corrector_Z(particles, settings, corrector_a_1 * dt, -inv * corrector_b_31 * dt);
          corrector_Z(particles, settings, -corrector_a_1 * dt, inv * corrector_b_31 * dt);
          break;
        }
        case WHFast::Order::FIFTH:
        {
          corrector_Z(particles, settings, -corrector_a_2 * dt, -inv * corrector_b_51 * dt);
          corrector_Z(particles, settings, -corrector_a_1 * dt, -inv * corrector_b_52 * dt);
          corrector_Z(particles, settings, corrector_a_1 * dt, inv * corrector_b_52 * dt);
          corrector_Z(particles, settings, corrector_a_2 * dt, inv * corrector_b_51 * dt);
          break;
        }
        case WHFast::Order::SEVENTH:
        {
          corrector_Z(particles, settings, -corrector_a_3 * dt, -inv * corrector_b_71 * dt);
          corrector_Z(particles, settings, -corrector_a_2 * dt, -inv * corrector_b_72 * dt);
          corrector_Z(particles, settings, -corrector_a_1 * dt, -inv * corrector_b_73 * dt);
          corrector_Z(particles, settings, corrector_a_1 * dt, inv * corrector_b_73 * dt);
          corrector_Z(particles, settings, corrector_a_2 * dt, inv * corrector_b_72 * dt);
          corrector_Z(particles, settings, corrector_a_3 * dt, inv * corrector_b_71 * dt);
        }
        case WHFast::Order::ELEVENTH:
        {
          corrector_Z(particles, settings, -corrector_a_5 * dt, -inv * corrector_b_111 * dt);
          corrector_Z(particles, settings, -corrector_a_4 * dt, -inv * corrector_b_112 * dt);
          corrector_Z(particles, settings, -corrector_a_3 * dt, -inv * corrector_b_113 * dt);
          corrector_Z(particles, settings, -corrector_a_2 * dt, -inv * corrector_b_114 * dt);
          corrector_Z(particles, settings, -corrector_a_1 * dt, -inv * corrector_b_115 * dt);
          corrector_Z(particles, settings, corrector_a_1 * dt, inv * corrector_b_115 * dt);
          corrector_Z(particles, settings, corrector_a_2 * dt, inv * corrector_b_114 * dt);
          corrector_Z(particles, settings, corrector_a_3 * dt, inv * corrector_b_113 * dt);
          corrector_Z(particles, settings, corrector_a_4 * dt, inv * corrector_b_112 * dt);
          corrector_Z(particles, settings, corrector_a_5 * dt, inv * corrector_b_111 * dt);
          break;
        }
        case WHFast::Order::SEVENTEENTH:
        {
          corrector_Z(particles, settings, -corrector_a_8 * dt, -inv * corrector_b_171 * dt);
          corrector_Z(particles, settings, -corrector_a_7 * dt, -inv * corrector_b_172 * dt);
          corrector_Z(particles, settings, -corrector_a_6 * dt, -inv * corrector_b_173 * dt);
          corrector_Z(particles, settings, -corrector_a_5 * dt, -inv * corrector_b_174 * dt);
          corrector_Z(particles, settings, -corrector_a_4 * dt, -inv * corrector_b_175 * dt);
          corrector_Z(particles, settings, -corrector_a_3 * dt, -inv * corrector_b_176 * dt);
          corrector_Z(particles, settings, -corrector_a_2 * dt, -inv * corrector_b_177 * dt);
          corrector_Z(particles, settings, -corrector_a_1 * dt, -inv * corrector_b_178 * dt);
          corrector_Z(particles, settings, corrector_a_1 * dt, inv * corrector_b_178 * dt);
          corrector_Z(particles, settings, corrector_a_2 * dt, inv * corrector_b_177 * dt);
          corrector_Z(particles, settings, corrector_a_3 * dt, inv * corrector_b_176 * dt);
          corrector_Z(particles, settings, corrector_a_4 * dt, inv * corrector_b_175 * dt);
          corrector_Z(particles, settings, corrector_a_5 * dt, inv * corrector_b_174 * dt);
          corrector_Z(particles, settings, corrector_a_6 * dt, inv * corrector_b_173 * dt);
          corrector_Z(particles, settings, corrector_a_7 * dt, inv * corrector_b_172 * dt);
          corrector_Z(particles, settings, corrector_a_8 * dt, inv * corrector_b_171 * dt);
          break;
        }
        default:
          break;
      }
    }

    void operator_C(ParticleStore& particles, WHFast& settings, double a, double b) {
      kepler_step(particles, settings, a);
      _transform::jacobi_to_inertial_pos(particles, settings.internals.p_jh);
      update_accel(particles, settings.gravity_method, settings.softening2);
      interaction_step(particles, b, settings.softening2, settings, settings.gravity_method);
      kepler_step(particles, settings, -a);
    }

    void operator_Y(ParticleStore& particles, WHFast& settings, double a, double b) {
      operator_C(particles, settings, a, b);
      operator_C(particles, settings, -a, -b);
    }

    void operator_U(ParticleStore& particles, WHFast& settings, double a, double b) {
      kepler_step(particles, settings, a);
      operator_Y(particles, settings, a, b);
      operator_Y(particles, settings, a, -b);
      kepler_step(particles, settings, -a);
    }

    void apply_corrector2(ParticleStore& particles, WHFast& settings, double inv, double dt) {
      double a = 0.5 * inv * dt;
      double b = corrector2_b * inv * dt;
      operator_U(particles, settings, a, b);
      operator_U(particles, settings, -a, b);
    }

    void calculate_jerk(ParticleStore& particles, WHFast& settings) {
      size_t N = particles.size();
      ParticleStore& jerk = settings.internals.p_jh;
      Vec3 Rj{ 0, 0, 0 };
      double Mj = 0;
      Vec3 Aj{ 0, 0, 0 };
      for (size_t j = 0; j < N; ++j) {
        jerk.accelerations[j] = { 0, 0, 0 };
        for (size_t i = 0; i < j + 1; ++i) {
          if (j > 1) {
            double dQkrj = Mj;
            if (i < j)
              dQkrj = -particles.mus[j];
            Vec3 Qk = particles.positions[j] - Rj / Mj;
            Vec3 da = particles.accelerations[j] - Aj / Mj;

            double dr2 = Qk.mag2();

            double prefact2 = dQkrj / (dr2 * std::sqrt(dr2));

            jerk.accelerations[i] += prefact2 * da;

            double alphasum = da.dot(Qk);
            double prefact1 = 3 * alphasum * prefact2 / dr2;
            jerk.accelerations[i] -= prefact1 * Qk;
          }

          if (j != i && (i != 0 || j != 1)) {
            Vec3 d = particles.positions[j] - particles.positions[i];
            Vec3 da = particles.accelerations[j] - particles.accelerations[i];

            double dr2 = d.mag2();
            double alphasum = da.dot(d);
            double prefact2 = 1 / (dr2 * std::sqrt(dr2));
            double prefact2i = prefact2 * particles.mus[i];
            double prefact2j = prefact2 * particles.mus[j];
            jerk.accelerations[j] -= da * prefact2i;
            jerk.accelerations[i] += da * prefact2j;
            double prefact1 = 3 * alphasum * prefact2 / dr2;
            double prefact1i = prefact1 * particles.mus[i];
            double prefact1j = prefact1 * particles.mus[j];
            jerk.accelerations[j] += d * prefact1i;
            jerk.accelerations[i] += d * prefact1j;
          }
        }
        Aj += particles.accelerations[j] * particles.mus[j];
        Rj += particles.positions[j] * particles.mus[j];
        Mj += particles.mus[j];
      }
    }
  } // end _whfast

  bool rebound::WHFast::init(ParticleStore& particles) {
#ifdef _OPENMP
    if (coordinates != WHFast::Coordinates::DEMOCRATIC_HELIOCENTRIC && coordinates != WHFast::Coordinates::WHDS) {
      std::cerr << "OpenMP is only supported for Democratic Heliocentric and WHDS coordinates." << '\n';
      return true;
    }
#endif
    if (kernel != Kernel::DEFAULT && coordinates != Coordinates::JACOBI) {
      std::cerr << "Non-standard kernel requires Jacobi coordinates.\n";
      return true;
    }

    if (keep_unsynchronized && safe_mode) {
      std::cerr << "Cannot use keep_unsynchronized with safe_mode.\n";
      return true;
    }

    if (kernel == Kernel::MODIFIEDKICK || kernel == Kernel::LAZY)
      gravity_method = GravityMethod::JACOBI;

    size_t N = particles.size();
    if (internals.p_jh.positions.capacity() != N) {
      internals.p_jh.positions.reserve(N);
      internals.p_jh.velocities.reserve(N);
      internals.p_jh.accelerations.reserve(N);
      internals.p_jh.mus.reserve(N);
      internals.p_jh.test_mass.reserve(N);
      std::copy(particles.test_mass.begin(), particles.test_mass.end(), internals.p_jh.test_mass.begin());
      recalc_coords_this_timestep = true;
    }
    return false;
  }

  void WHFast::from_inertial(ParticleStore& particles) {
    switch (coordinates) {
      case Coordinates::JACOBI:
        _transform::inertial_to_jacobi_posvel(particles, internals.p_jh);
        break;
      case Coordinates::DEMOCRATIC_HELIOCENTRIC:
        _transform::inertial_to_democraticheliocentric_posvel(particles, internals.p_jh);
        break;
      case Coordinates::WHDS:
        _transform::inertial_to_whds_posvel(particles, internals.p_jh);
        break;
      case Coordinates::BARYCENTRIC:
        _transform::inertial_to_barycentric_posvel(particles, internals.p_jh);
        break;
    }
  }

  void WHFast::to_inertial(ParticleStore& particles) {
    switch (coordinates) {
      case Coordinates::JACOBI:
        _transform::jacobi_to_inertial_posvel(particles, internals.p_jh);
        break;
      case Coordinates::DEMOCRATIC_HELIOCENTRIC:
        _transform::democraticheliocentric_to_inertial_posvel(particles, internals.p_jh);
        break;
      case Coordinates::WHDS:
        _transform::whds_to_inertial_posvel(particles, internals.p_jh);
        break;
      case Coordinates::BARYCENTRIC:
        _transform::barycentric_to_inertial_posvel(particles, internals.p_jh);
        break;
    }
  }

  void WHFast::debug_operator_kepler(ParticleStore& particles, double dt) {
    if (init(particles))
      return;
    from_inertial(particles);
    _whfast::kepler_step(particles, *this, dt);
    _whfast::com_step(internals.p_jh, dt);
    to_inertial(particles);
  }

  void WHFast::debug_operator_interaction(ParticleStore& particles, double dt) {
    if (init(particles))
      return;
    from_inertial(particles);
    _whfast::update_accel(particles, gravity_method, softening2);
    _whfast::interaction_step(particles, dt, softening2, *this, gravity_method);
    to_inertial(particles);
  }

  void WHFast::step_p1(ParticleStore& particles, double dt) {
    if (init(particles))
      return;
    if (safe_mode || recalc_coords_this_timestep) {
      if (!internals.is_synchronized) {
        synchronize(particles, dt);
        if (!internals.recalc_coords_not_synchronized_warning) {
          std::cerr << "Recalculating coordinates but pos/vel were not synchronized" << '\n';
          internals.recalc_coords_not_synchronized_warning = true;
        }
      }
      from_inertial(particles);
      recalc_coords_this_timestep = false;
    }
    if (internals.is_synchronized) {
      _whfast::apply_corrector(particles, *this, 1., dt);
      if (use_corrector_2)
        _whfast::apply_corrector2(particles, *this, 1., dt);
      switch (kernel) {
        case Kernel::DEFAULT:
        case Kernel::MODIFIEDKICK:
        case Kernel::LAZY:
          _whfast::kepler_step(particles, *this, 0.5 * dt);
          _whfast::com_step(internals.p_jh, 0.5 * dt);
          break;
        case Kernel::COMPOSITION:
          _whfast::kepler_step(particles, *this, 0.625 * dt);
          _whfast::com_step(internals.p_jh, 0.625 * dt);
          break;
        default:
          std::cerr << "WHFast kernel not implemented." << '\n';
          return;
      }
    } else {
      _whfast::kepler_step(particles, *this, dt);
      _whfast::com_step(internals.p_jh, dt);
    }

    to_inertial(particles);
  }

  void WHFast::synchronize(ParticleStore& particles, double dt) {
    size_t N = particles.size();
    if (init(particles))
      return;
    if (!internals.is_synchronized) {
      ParticleStore sync_pj;
      if (keep_unsynchronized) {
        for (size_t i = 0; i < N; ++i) {
          sync_pj.positions.push_back(internals.p_jh.positions[i]);
          sync_pj.velocities.push_back(internals.p_jh.velocities[i]);
          sync_pj.accelerations.push_back(internals.p_jh.accelerations[i]);
          sync_pj.mus.push_back(internals.p_jh.mus[i]);
          sync_pj.test_mass.push_back(internals.p_jh.test_mass[i]);
          sync_pj.ids.push_back(internals.p_jh.ids[i]);
          sync_pj.versions.push_back(internals.p_jh.versions[i]);
        }
      }

      switch (kernel) {
        case Kernel::DEFAULT:
        case Kernel::MODIFIEDKICK:
        case Kernel::LAZY:
          _whfast::kepler_step(particles, *this, -0.5 * dt);
          _whfast::com_step(internals.p_jh, -0.5 * dt);
          break;
        case Kernel::COMPOSITION:
          _whfast::kepler_step(particles, *this, 0.375 * dt);
          _whfast::com_step(internals.p_jh, 0.375 * dt);
          break;
        default:
          std::cerr << "WHFast kernel not implemented." << '\n';
          return;
      }
      if (keep_unsynchronized) {
        for (size_t i = 0; i < N; ++i) {
          internals.p_jh.positions[i] = sync_pj.positions[i];
          internals.p_jh.velocities[i] = sync_pj.velocities[i];
          internals.p_jh.accelerations[i] = sync_pj.accelerations[i];
        }
      } else
        internals.is_synchronized = true;
    }
  }

  void WHFast::step_p2(ParticleStore& particles, double dt) {
    size_t N = particles.size();
    if (N == 0)
      return;
    switch (kernel) {
      case Kernel::DEFAULT:
        _whfast::interaction_step(particles, dt, softening2, *this, gravity_method);
        _whfast::jump_step(particles, *this, 0.5 * dt);
        break;
      case Kernel::MODIFIEDKICK:
        _whfast::calculate_jerk(particles, *this);
        for (size_t i = 0; i < N; ++i) {
          double prefact = dt * dt / 12.;
          particles.accelerations[i] += prefact * internals.p_jh.accelerations[i];
        }
        _whfast::interaction_step(particles, dt, softening2, *this, gravity_method);
        break;
      case Kernel::LAZY:
        if (internals.p_temp.size() != N) {
          internals.p_temp.positions.resize(N);
          internals.p_temp.velocities.resize(N);
          internals.p_temp.accelerations.resize(N);
          internals.p_temp.mus.resize(N);
          internals.p_temp.test_mass.resize(N);
          internals.p_temp.ids.resize(N);
          internals.p_temp.versions.resize(N);
        }
        _transform::inertial_to_jacobi_acc(particles, internals.p_temp);
        internals.p_jh = internals.p_temp;
        for (size_t i = 1; i < N; ++i) {
          double prefac1 = dt * dt / 12.;
          particles.positions[i] += prefac1 * internals.p_temp.accelerations[i];
        }

        _transform::jacobi_to_inertial_pos(particles, internals.p_jh);
        _whfast::update_accel(particles, gravity_method, softening2);
        _whfast::interaction_step(particles, dt, softening2, *this, gravity_method);
        for (size_t i = 1; i < N; ++i) {
          internals.p_jh.positions[i] = internals.p_temp.positions[i];
        }
        break;
      case Kernel::COMPOSITION:
        _whfast::interaction_step(particles, 1. / 6. * dt, softening2, *this, gravity_method);
        _whfast::kepler_step(particles, *this, 0.25 * dt);
        _whfast::com_step(internals.p_jh, 0.25 * dt);
        _transform::jacobi_to_inertial_pos(particles, internals.p_jh);
        _whfast::update_accel(particles, gravity_method, softening2);
        _whfast::interaction_step(particles, 1. / 6. * dt, softening2, *this, gravity_method);
        _whfast::kepler_step(particles, *this, 0.125 * dt);
        _whfast::com_step(internals.p_jh, 0.125 * dt);
        _transform::jacobi_to_inertial_pos(particles, internals.p_jh);
        _whfast::update_accel(particles, gravity_method, softening2);
        _whfast::interaction_step(particles, 1. / 6. * dt, softening2, *this, gravity_method);
        _whfast::kepler_step(particles, *this, 0.25 * dt);
        _whfast::com_step(internals.p_jh, 0.25 * dt);
        _transform::jacobi_to_inertial_pos(particles, internals.p_jh);
        _whfast::update_accel(particles, gravity_method, softening2);
        _whfast::interaction_step(particles, dt * 1. / 6., softening2, *this, gravity_method);
        break;
      default:
        return;
    }

    internals.is_synchronized = false;
    if (safe_mode)
      synchronize(particles, dt);
  }

  void WHFast::reset() {
    order = Order::NONE;
    use_corrector_2 = false;
    kernel = Kernel::DEFAULT;
    coordinates = Coordinates::JACOBI;
    internals.is_synchronized = true;
    keep_unsynchronized = false;
    safe_mode = true;
    recalc_coords_this_timestep = false;
    internals.recalc_coords_not_synchronized_warning = false;
    if (internals.p_jh.size() > 0) {
      internals.p_jh.positions.clear();
      internals.p_jh.velocities.clear();
      internals.p_jh.accelerations.clear();
      internals.p_jh.mus.clear();
      internals.p_jh.test_mass.clear();
      internals.p_jh.ids.clear();
      internals.p_jh.versions.clear();
    }
    if (internals.p_temp.size() > 0) {
      internals.p_temp.positions.clear();
      internals.p_temp.velocities.clear();
      internals.p_temp.accelerations.clear();
      internals.p_temp.mus.clear();
      internals.p_temp.test_mass.clear();
      internals.p_temp.ids.clear();
      internals.p_temp.versions.clear();
    }
  }

  void WHFast::step(ParticleStore& particles, double dt) {
    synchronize(particles, dt);
    step_p1(particles, dt);
    // particles.print_if_nan_or_inf();
    step_p2(particles, dt);
    // particles.print_if_nan_or_inf();
    synchronize(particles, dt);
  }
} // end rebound