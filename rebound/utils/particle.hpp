/*
 * This file is part of the C++ translation of REBOUND.
 *
 * Original REBOUND (C) code by Hanno Rein and others.
 * This translation is licensed under the GNU General Public License v3 or
 * later.
 *
 * REBOUND is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * REBOUND is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "../repstl/pair"
#include "../repstl/string"
#include "../repstl/vector"
#include "vec.hpp"
#include <cassert>
#include <cstdint>

namespace rebound {
  struct ind_Particle;
  class Particle;
  class ConstParticle;

  template <typename T>
  using ExtParamMap =
    repstl::Vector<repstl::pair<repstl::String, repstl::Vector<T>>>;

  template <typename T>
  repstl::Vector<T>& get(ExtParamMap<T>& map, repstl::String param) {
    for (size_t i = 0; i < map.size(); ++i) {
      if (map[i].first == param)
        return map[i].second;
    }
    repstl::String msg("Cannot find parameter: ");
    msg += param;
    throw repstl::RuntimeError(msg.c_str());
  }

  struct ParticleStore {
    repstl::Vector<Vec3> positions;
    repstl::Vector<Vec3> velocities;
    repstl::Vector<Vec3> accelerations;
    repstl::Vector<Vec3> gravity_cs;
    repstl::Vector<double> mus;
    repstl::Vector<double> radii;
    repstl::Vector<uint32_t> ids;
    repstl::Vector<bool> test_mass;

    repstl::Vector<uint32_t> versions;

    ExtParamMap<double> double_params;
    ExtParamMap<uint32_t> hash_params;
    ExtParamMap<int> int_params;

    ParticleStore(size_t capacity = 0);

    size_t size() const;
    size_t n_real() const;

    Particle add_particle(const Vec3& position, const Vec3& velocity, double mu,
      double radius, uint32_t id, bool test_mass = false);
    void remove_particle(size_t idx);

    Particle operator[](size_t idx);
    ConstParticle operator[](size_t idx) const;

    void print_if_nan_or_inf() const;

    inline bool operator==(const ParticleStore& other) const {
      return positions == other.positions && velocities == other.velocities &&
        accelerations == other.accelerations &&
        gravity_cs == other.gravity_cs && mus == other.mus &&
        radii == other.radii && ids == other.ids &&
        test_mass == other.test_mass &&
        double_params == other.double_params &&
        int_params == other.int_params && hash_params == other.hash_params;
    }
  };

  class Particle {
  private:
    size_t index;
    uint32_t version;
    ParticleStore* store;
    Particle(size_t idx, ParticleStore* store_, uint32_t ver)
      : index(idx), store(store_), version(ver) {
    }

  public:
#define CHECK_VER assert(store->versions[index] == version)
    Vec3& pos() {
      CHECK_VER;
      return store->positions[index];
    }
    Vec3& vel() {
      CHECK_VER;
      return store->velocities[index];
    }
    Vec3& acc() {
      CHECK_VER;
      return store->accelerations[index];
    }
    double& mu() {
      CHECK_VER;
      return store->mus[index];
    }
    double& radius() {
      CHECK_VER;
      return store->radii[index];
    }
    void toggle_test_mass() {
      CHECK_VER;
      store->test_mass[index] = !store->test_mass[index];
    }

    const Vec3& pos() const {
      CHECK_VER;
      return store->positions[index];
    }
    const Vec3& vel() const {
      CHECK_VER;
      return store->velocities[index];
    }
    const Vec3& acc() const {
      CHECK_VER;
      return store->accelerations[index];
    }
    const double& mu() const {
      CHECK_VER;
      return store->mus[index];
    }
    const double& radius() const {
      CHECK_VER;
      return store->radii[index];
    }
    const uint32_t& id() const {
      CHECK_VER;
      return store->ids[index];
    }
    const bool& get_test_mass() const {
      CHECK_VER;
      return store->test_mass[index];
    }
#undef CHECK_VER

    friend ParticleStore;
  };

  class ConstParticle {
  private:
    size_t index;
    uint32_t version;
    const ParticleStore* store;
    ConstParticle(const size_t idx, const ParticleStore* store_,
      const uint32_t ver)
      : index(idx), store(store_), version(ver) {
    }

  public:
#define CHECK_VER assert(store->versions[index] == version)
    const Vec3& pos() const {
      CHECK_VER;
      return store->positions[index];
    }
    const Vec3& vel() const {
      CHECK_VER;
      return store->velocities[index];
    }
    const Vec3& acc() const {
      CHECK_VER;
      return store->accelerations[index];
    }
    const double& mu() const {
      CHECK_VER;
      return store->mus[index];
    }
    const double& radius() const {
      CHECK_VER;
      return store->radii[index];
    }
    const uint32_t& id() const {
      CHECK_VER;
      return store->ids[index];
    }
    const bool& get_test_mass() const {
      CHECK_VER;
      return store->test_mass[index];
    }
#undef CHECK_VER

    friend ParticleStore;
  };
} // namespace rebound