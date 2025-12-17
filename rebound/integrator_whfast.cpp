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

    constexpr double inv_factorial[35] = {1., 1., 1./2., 1./6., 1./24., 1./120.,
      1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 
      1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 
      1./20922789888000., 1./355687428096000., 1./6402373705728000., 
      1./121645100408832000., 1./2432902008176640000., 
      1./51090942171709440000., 1./1124000727777607680000., 
      1./25852016738884976640000., 1./620448401733239439360000., 
      1./15511210043330985984000000., 1./403291461126605635584000000., 
      1./10888869450418352160768000000., 1./304888344611713860501504000000., 
      1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 
      1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 
      1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};

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
      double X2 = X*X;
      stumpff_cs3(Gs, beta*X2);
      Gs[1] *= X; 
      Gs[2] *= X2; 
      Gs[3] *= X2*X;
      return;
    }

    constexpr unsigned int WHFAST_NMAX_QUART = 64;
    constexpr unsigned int WHFAST_NMAX_NEWT = 32;

    inline void kepler_solver(const _WHFastSettings &settings, ParticleStore& p_j, double M, size_t i, double dt) {
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
      } else X = 0.;

      unsigned int converged = 0;
      double oldX = X;

      stiefel_Gs3(Gs, beta, X);
      double eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
      double ri = 1. / (r0 + eta0Gs1zeta0Gs2);
      X = ri * (X * eta0Gs1zeta0Gs2 - eta0 * Gs[2] - zeta0 * Gs[3] + dt);

      if (fast_abs(X - oldX) > 0.01 * X_per_period) {
        X = beta * dt / M;
        double prevX[WHFAST_NMAX_QUART + 1];
        for (int nlag = 1; nlag < WHFAST_NMAX_QUART; nlag++) {
          stiefel_Gs3(Gs, beta, X);
          double f = r0 * X + eta0 * Gs[2] + zeta0 * Gs[3] - dt;
          double denom = eta0 * Gs[1] + zeta0 * Gs[2] + r0;
          X = (X * denom - 5. * f) / denom;
          for (int i = 1; i < nlag; i++) {
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
        for (int n_hg = 1; n_hg < WHFAST_NMAX_NEWT; n_hg++) {
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

    void interaction_step(ParticleStore &particles, double dt, double softening2, _WHFastSettings &settings, GravityMethod method) {
      auto p_j = settings.internals.p_jh;
      size_t N = particles.size();
      switch (settings.coordinates) {
      case _WHFastSettings::Coordinates::JACOBI: {
        _transform::inertial_to_jacobi_acc(particles, *p_j);
        double eta = particles.mus[0];
        for (size_t i = 1; i < N; ++i) {
          Particle pji = (*p_j)[i];
          if (!particles.test_mass[i]) eta += pji.mu();
          p_j->velocities[i] += pji.acc() * dt;
          if (method != GravityMethod::JACOBI) {
            if (i > 1) {
              double rj2i = 1. / (pji.pos().mag2() + softening2);
              double rji = std::sqrt(rj2i);
              double rj3iM = rji * rj2i * eta;
              double prefac1 = dt * rj3iM;
              p_j->velocities[i] += pji.pos();
            }
          }
        }
        break;
      } 
      case _WHFastSettings::Coordinates::DEMOCRATIC_HELIOCENTRIC: {
        for (size_t i = 1; i < N; ++i) {
          if (!particles.test_mass[i]) p_j->velocities[i] += dt * p_j->accelerations[i];
        }
        break;
      }
      default:
        break;
      }
    }
  }
}