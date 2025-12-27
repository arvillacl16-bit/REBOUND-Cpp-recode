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

#pragma once
#include "vec.hpp"
#include <cstdint>
#include <vector>
#include "extension_preprocessor.hpp"

using boollist = std::vector<uint8_t>;

namespace rebound {
  struct ind_Particle;
  class Particle;

  template <typename T1, typename T2>
  struct pair {
    T1 first;
    T2 second;

    pair(const T1& first_, const T2& second_) : first(first_), second(second_) {}
  };

  template <typename T>
  using ExtParamMap = std::vector<pair<std::string, std::vector<T>>>;
  struct ParticleStore {
    std::vector<Vec3> positions;
    std::vector<Vec3> velocities;
    std::vector<Vec3> accelerations;
    std::vector<Vec3> gravity_cs;
    std::vector<double> mus;
    std::vector<double> radii;
    std::vector<uint32_t> ids;
    boollist test_mass;

    std::vector<uint32_t> versions;

    ExtParamMap<double> double_params;
    ExtParamMap<uint32_t> hash_params;
    ExtParamMap<int> int_params;
    ExtParamMap<void*> ptr_params;

    void reserve_double_ex_params(const std::vector<std::string> &names);
    void reserve_int_ex_params(const std::vector<std::string> &names);
    void reserve_hash_ex_params(const std::vector<std::string> &names);
    void reserve_ptr_ex_params(const std::vector<std::string> &names);

    void reserve_double_ex_param(const std::string &name);
    void reserve_int_ex_param(const std::string &name);
    void reserve_hash_ex_param(const std::string &name);
    void reserve_ptr_ex_param(const std::string &name);

    ParticleStore(size_t capacity = 0);

    size_t size() const;
    size_t n_real() const;

    Particle add_particle(const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass = false);
    void remove_particle(size_t idx);

    Particle operator[](size_t idx);
    const Particle operator[](size_t idx) const;
  };

  class Particle {
  private:
    size_t index;
    uint32_t version;
    ParticleStore& store;
    Particle(size_t idx, ParticleStore& store_) : index(idx), store(store_) {}
  public:
    Particle(ParticleStore& store_, const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass = false);

    Vec3& pos();
    Vec3& vel();
    Vec3& acc();
    double& mu();
    double& radius();
    uint32_t& id();
    void toggle_test_mass();

    const Vec3& pos() const;
    const Vec3& vel() const;
    const Vec3& acc() const;
    const double &mu() const;
    const double &radius() const;
    const uint32_t &id() const;
    const bool &get_test_mass() const;

    ind_Particle snap() const;

    friend ParticleStore;
  };

  struct Orbit {
    double a;
    double e;
    double inc;
    double Omega;
    double omega;
    double f;

    Orbit(double a_, double e_, double inc_, double Omega_, double omega_, double f_)
      : a(a_), e(e_), inc(inc_), Omega(Omega_), omega(omega_), f(f_) {}

    static Orbit from_particles(const Particle& primary, const Particle& secondary);

    static double M_to_E(double M, double e, double tol = 1e-10);
    static double E_to_f(double E, double e);
    static double M_to_f(double M, double e, double tol = 1e-10);
  };

  struct ind_Particle {
    Vec3 pos;
    Vec3 vel;
    double mu;
    double r;
    uint32_t id;
    bool test_mass;

    ind_Particle() : pos(), vel(), mu(0), r(0), id(0), test_mass(false) {}
    ind_Particle(const Vec3& p, const Vec3& v, double mu_, double r_, uint32_t id_, bool test_mass_ = false)
        : pos(p), vel(v), mu(mu_), r(r_), id(id_), test_mass(test_mass_) {}
    ind_Particle(const Particle& p) : pos(p.pos()), vel(p.vel()), mu(p.mu()), r(p.radius()), id(p.id()), test_mass(p.get_test_mass()) {}

    static ind_Particle from_orbit(const Orbit& orbit, const Particle& primary, const Particle& secondary, uint32_t id, bool test_mass = false);
  };
}