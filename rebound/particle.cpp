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

#include "particle.hpp"
#include <cassert>
#include <cmath>

namespace rebound {
  ParticleStore::ParticleStore(size_t capacity) {
    positions.reserve(capacity);
    velocities.reserve(capacity);
    accelerations.reserve(capacity);
    gravity_cs.reserve(capacity);
    mus.reserve(capacity);
    radii.reserve(capacity);
    ids.reserve(capacity);
    test_mass.reserve(capacity);
  }

  size_t ParticleStore::size() const { return positions.size(); }

  size_t ParticleStore::n_real() const {
    size_t n = 0;
    for (size_t i = 0; i < size(); i++) {
      if (!test_mass[i]) ++n;
    }
    return n;
  }

  Particle ParticleStore::add_particle(const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass_) {
    positions.push_back(position);
    velocities.push_back(velocity);
    accelerations.push_back(Vec3());
    gravity_cs.push_back(Vec3());
    mus.push_back(mu);
    radii.push_back(radius);
    ids.push_back(id);
    test_mass.push_back(test_mass_);
    versions.push_back(0);
    return Particle(size() - 1, *this);
  }

  void ParticleStore::remove_particle(size_t idx) {
    versions[idx]++;

    size_t last_idx = size() - 1;
    if (idx != last_idx) {
      positions[idx] = positions[last_idx];
      velocities[idx] = velocities[last_idx];
      accelerations[idx] = accelerations[last_idx];
      gravity_cs[idx] = gravity_cs[last_idx];
      mus[idx] = mus[last_idx];
      radii[idx] = radii[last_idx];
      ids[idx] = ids[last_idx];
      test_mass[idx] = test_mass[last_idx];
      versions[idx] = versions[last_idx];
    }
    positions.pop_back();
    velocities.pop_back();
    accelerations.pop_back();
    gravity_cs.pop_back();
    mus.pop_back();
    radii.pop_back();
    ids.pop_back();
    test_mass.pop_back();
    versions.pop_back();
  }

  Particle::Particle(ParticleStore& store_, const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass_)
    : index(store_.size()), store(store_), version(0) { store.add_particle(position, velocity, mu, radius, id, test_mass_); }

  Vec3& Particle::pos() {
    return store.positions[index];
  }

  Vec3& Particle::vel() {
    return store.velocities[index];
  }

  Vec3& Particle::acc() {
    return store.accelerations[index];
  }

  double& Particle::mu() {
    return store.mus[index];
  }

  double& Particle::radius() {
    return store.radii[index];
  }

  uint32_t& Particle::id() {
    return store.ids[index];
  }

  void Particle::toggle_test_mass() {
    store.test_mass[index] = !store.test_mass[index];
  }

  // const versions
  const Vec3& Particle::pos() const {
    return store.positions[index];
  }

  const Vec3& Particle::vel() const {
    return store.velocities[index];
  }

  const Vec3& Particle::acc() const {
    return store.accelerations[index];
  }

  const double &Particle::mu() const {
    return store.mus[index];
  }

  const double &Particle::radius() const {
    return store.radii[index];
  }

  const uint32_t &Particle::id() const {
    return store.ids[index];
  }

  const bool &Particle::get_test_mass() const {
    return reinterpret_cast<const bool&>(store.test_mass[index]);
  }

  ind_Particle Particle::snap() const {
    return ind_Particle(store.positions[index], store.velocities[index], store.mus[index], store.radii[index], store.ids[index], store.test_mass[index]);
  }

  double Orbit::M_to_E(double M, double e, double tol) {
    if (e == 0) return M;
    double E = M;
    for (int i = 0; i < 100; ++i) {
      double f = E - e * std::sin(E) - M;
      double f_prime = 1 - e * std::cos(E);
      double delta = -f / f_prime;
      E += delta;
      if (std::abs(delta) < tol) break;
    }
    return E;
  }

  double Orbit::E_to_f(double E, double e) {
    double cos_f = (std::cos(E) - e) / (1 - e * std::cos(E));
    double sin_f = (std::sqrt(1 - e * e) * std::sin(E)) / (1 - e * std::cos(E));
    return std::atan2(sin_f, cos_f);
  }

  double Orbit::M_to_f(double M, double e, double tol) {
    double E = M_to_E(M, e, tol);
    return E_to_f(E, e);
  }

  Orbit Orbit::from_particles(const Particle& primary, const Particle& secondary) {
    Vec3 r = secondary.pos() - primary.pos();
    Vec3 v = secondary.vel() - primary.vel();
    double mu = primary.mu() + secondary.mu();

    double r_mag = r.mag();
    double v2 = v.mag2();

    Vec3 h = r.cross(v);
    double h_mag = h.mag();

    Vec3 k(0.0, 0.0, 1.0);
    Vec3 n = k.cross(h);
    double n_mag = n.mag();

    Vec3 e_vec = (v.cross(h)) / mu - (r / r_mag);
    double e = e_vec.mag();

    double a = 1.0 / (2.0 / r_mag - v2 / mu);

    double inc = std::acos(h.z / h_mag);
    double Omega = 0.0;
    if (n_mag > 1e-16) Omega = std::atan2(n.y, n.x);
    double omega = 0.0;
    if (n_mag > 1e-16 && e > 1e-16) omega = std::acos(n.dot(e_vec) / (n_mag * e));
    if (e_vec.z < 0) omega = 2.0 * M_PI - omega;

    double f = 0.0;
    if (e > 1e-16)
        f = std::acos(e_vec.dot(r) / (e * r_mag));
    if (r.dot(v) < 0)
        f = 2.0 * M_PI - f;

    return Orbit(a, e, inc, Omega, omega, f);
  }

  ind_Particle ind_Particle::from_orbit(const Orbit& orbit, const Particle& primary, const Particle& secondary, uint32_t id, bool test_mass) {
    double mu = primary.mu() + secondary.mu();
    double a = orbit.a;
    double e = orbit.e;
    double inc = orbit.inc;
    double Omega = orbit.Omega;
    double omega = orbit.omega;
    double f = orbit.f;

    double r_mag = a * (1 - e * e) / (1 + e * std::cos(f));

    Vec3 r_perifocal(r_mag * std::cos(f), r_mag * std::sin(f), 0.0);
    double h = std::sqrt(mu * a * (1 - e * e));
    Vec3 v_perifocal(mu / h * std::sin(f), mu / h * (e + std::cos(f)), 0.0);

    double cos_Omega = std::cos(Omega);
    double sin_Omega = std::sin(Omega);
    double cos_inc = std::cos(inc);
    double sin_inc = std::sin(inc);
    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);

    Vec3 r_eci;
    r_eci.x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_inc) * r_perifocal.x +
              (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_inc) * r_perifocal.y;
    r_eci.y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_inc) * r_perifocal.x +
              (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_inc) * r_perifocal.y;
    r_eci.z = (sin_omega * sin_inc) * r_perifocal.x + (cos_omega * sin_inc) * r_perifocal.y;

    Vec3 v_eci;
    v_eci.x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_inc) * v_perifocal.x +
              (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_inc) * v_perifocal.y;
    v_eci.y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_inc) * v_perifocal.x +
              (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_inc) * v_perifocal.y;
    v_eci.z = (sin_omega * sin_inc) * v_perifocal.x + (cos_omega * sin_inc) * v_perifocal.y;
    Vec3 position = primary.pos() + r_eci;
    Vec3 velocity = primary.vel() + v_eci;
    return ind_Particle(position, velocity, mu, secondary.radius(), id, test_mass);
  }

  Particle ParticleStore::operator[](size_t idx) { return Particle(idx, *this); }
  const Particle ParticleStore::operator[](size_t idx) const { return Particle(idx, const_cast<ParticleStore&>(*this)); }

  namespace {
    template <typename T>
    void reserve_ex_param(ExtParamMap<T>& map, const std::string& name, size_t N) {
      auto it = std::find_if(map.begin(), map.end(), 
        [&name](const pair<std::string, std::vector<T>>& p) { return p.first == name; });
      if (it != map.end()) it->second.reserve(N);
    }
  }

  void ParticleStore::reserve_double_ex_params(const std::vector<std::string> &names) {
    size_t N = size();
    for (const auto& name : names) reserve_ex_param(double_params, name, N);
  }

  void ParticleStore::reserve_int_ex_params(const std::vector<std::string> &names) {
    size_t N = size();
    for (const auto& name : names) reserve_ex_param(int_params, name, N);
  }

  void ParticleStore::reserve_hash_ex_params(const std::vector<std::string> &names) {
    size_t N = size();
    for (const auto& name : names) reserve_ex_param(hash_params, name, N);
  }

  void ParticleStore::reserve_ptr_ex_params(const std::vector<std::string> &names) {
    size_t N = size();
    for (const auto& name : names) reserve_ex_param(ptr_params, name, N);
  }

  void ParticleStore::reserve_double_ex_param(const std::string &name) {
    size_t N = size();
    reserve_ex_param(double_params, name, N);
  }

  void ParticleStore::reserve_int_ex_param(const std::string &name) {
    size_t N = size();
    reserve_ex_param(int_params, name, N);
  }

  void ParticleStore::reserve_hash_ex_param(const std::string &name) {
    size_t N = size();
    reserve_ex_param(hash_params, name, N);
  }

  void ParticleStore::reserve_ptr_ex_param(const std::string &name) {
    size_t N = size();
    reserve_ex_param(ptr_params, name, N);
  }
}