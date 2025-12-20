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
#include <chrono>
#include <unordered_set>

namespace chrono = std::chrono;

namespace rebound {
  struct Simulation::hashmap {
    std::unordered_map<uint32_t, size_t> hash_map;
  };

  Simulation::Simulation() : ptr_hash(new hashmap) {}
  Simulation::Simulation(const Simulation &other) : ptr_hash(new hashmap), integrator(other.integrator), coll_handler(other.coll_handler) {
    auto &parts = other.particles;
    for (size_t i = 0; i < parts.size(); ++i) add_particle(parts.positions[i], parts.velocities[i], parts.mus[i], parts.radii[i], parts.ids[i], parts.test_mass[i]);
  }

  Simulation::Simulation(Simulation &&other) : ptr_hash(other.ptr_hash) { other.ptr_hash = nullptr; }

  Simulation &Simulation::operator=(const Simulation &other) {
    ptr_hash->hash_map.clear();
    auto &parts = other.particles;
    for (size_t i = 0; i < parts.size(); ++i) add_particle(parts.positions[i], parts.velocities[i], parts.mus[i], parts.radii[i], parts.ids[i], parts.test_mass[i]);
    return *this;
  }

  Simulation &Simulation::operator=(Simulation &&other) { 
    ptr_hash = other.ptr_hash;
    other.ptr_hash = nullptr; 
    return *this;
  }

  Simulation::~Simulation() { delete ptr_hash; }

  Particle Simulation::add_particle(const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass) {
    ptr_hash->hash_map.emplace(id,curr_idx++);
    return particles.add_particle(position, velocity, mu, radius, id, test_mass);
  }

  void Simulation::remove_particle(size_t idx) {
    if (idx >= n()) throw std::out_of_range("Index outside of simulation");
    auto &hash_map = ptr_hash->hash_map;
    size_t last_idx = n() - 1; 
    uint32_t last_id = particles.ids[last_idx]; 
    if (idx != last_idx) {
      hash_map[last_id] = idx;
      hash_map.erase(particles.ids[idx]);
    } else hash_map.erase(last_id);
    particles.remove_particle(idx);
    --curr_idx;
  }

  void Simulation::remove_particle(uint32_t id) {
    auto &hash_map = ptr_hash->hash_map;
    if (auto it = hash_map.find(id); it != hash_map.end()) remove_particle(it->second);
  }

  double Simulation::time() const { return t; }

  void Simulation::integrate(double t_end, bool exact_finish_time) {
    double tleft = t_end - t;
    while (tleft > 0) {
      double dt_step = exact_finish_time ? std::min(tleft, dt) : dt;
      if (step(dt_step)) break;
      tleft -= dt_step;
    }
  }

  void Simulation::integrate_walltime(double t) {
    auto start = chrono::steady_clock::now();
    chrono::duration<double> max_duration(t);

    while (chrono::steady_clock::now() - start < max_duration) if (step(dt)) break;
  }

  bool Simulation::step(double dt_) {
    bool val = false;

    coll_handler->prev_pos = particles.positions;
    integrator->step(particles, dt_);
    t += dt_;

    val = coll_handler->detect_collision(particles);

    if (heartbeat(*this)) val = true;
    return val;
  }

  Vec3 Simulation::com_pos() {
    std::vector<Vec3> &positions = particles.positions;
    std::vector<double> &mus = particles.mus;
    if (positions.size() != n() || mus.size() != n()) return {};

    double total_mu = 0;
    Vec3 result = {};

#pragma omp parallel for reduction(+:total_mu)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) total_mu += mus[i];

#pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) result += positions[i] * (mus[i] / total_mu);

    return result;
  }

  Vec3 Simulation::com_vel() {
    std::vector<Vec3> &velocities = particles.velocities;
    std::vector<double> &mus = particles.mus;
    if (velocities.size() != n() || mus.size() != n()) return {};

    double total_mu = 0;
    Vec3 result = {};

#pragma omp parallel for reduction(+:total_mu)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) total_mu += mus[i];

#pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) result += velocities[i] * (mus[i] / total_mu);

    return result;
  }

  void Simulation::move_to_com() {
    Vec3 com_pos_ = com_pos(), com_vel_ = com_vel();
#pragma omp parallel for
    for (size_t i = 0; i < n(); ++i) {
      particles.positions[i] -= com_pos_;
      particles.velocities[i] -= com_vel_;
    }
  }
}