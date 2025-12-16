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
#include <fstream>

namespace rebound {
  void Simulation::save(const std::string& filename) const {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) throw std::runtime_error("Failed to open file for saving: " + filename);

    int8_t method = static_cast<int8_t>(integrator.method);
    int8_t gravity_method = static_cast<int8_t>(integrator.gravity_method);
    outfile.write(reinterpret_cast<const char*>(&method), sizeof(int8_t));
    outfile.write(reinterpret_cast<const char*>(&gravity_method), sizeof(int8_t));
    outfile.write(reinterpret_cast<const char*>(&t), sizeof(double));
    uint64_t nparticles = particles.size();
    outfile.write(reinterpret_cast<const char*>(&nparticles), sizeof(uint64_t));
   
    outfile.write(reinterpret_cast<const char*> (particles.positions.data()), sizeof(Vec3) * nparticles);
    outfile.write(reinterpret_cast<const char*> (particles.velocities.data()), sizeof(Vec3) * nparticles);
    outfile.write(reinterpret_cast<const char*> (particles.accelerations.data()), sizeof(Vec3) * nparticles);
    outfile.write(reinterpret_cast<const char*> (particles.mus.data()), sizeof(double) * nparticles);
    outfile.write(reinterpret_cast<const char*> (particles.radii.data()), sizeof(double) * nparticles);
    outfile.write(reinterpret_cast<const char*> (particles.ids.data()), sizeof(uint32_t) * nparticles);
    outfile.write(reinterpret_cast<const char*> (particles.test_mass.data()), sizeof(uint8_t) * nparticles);

    outfile.close();
  }

  Simulation::Simulation(std::string filename) : integrator(Integrator(IntegratorMethod::NONE, 0.0)), t(0.0) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) throw std::runtime_error("Input file stream is not valid.");

    int8_t method_int;
    int8_t gravity_method_int;
    infile.read(reinterpret_cast<char*>(&method_int), sizeof(int8_t));
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    infile.read(reinterpret_cast<char*>(&gravity_method_int), sizeof(int8_t));
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    integrator.method = static_cast<IntegratorMethod>(method_int);
    integrator.gravity_method = static_cast<GravityMethod>(gravity_method_int);

    infile.read(reinterpret_cast<char*>(&t), sizeof(double));
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    uint64_t nparticles;
    infile.read(reinterpret_cast<char*>(&nparticles), sizeof(uint64_t));
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");

    particles = _ParticleStore(nparticles);

    particles.positions.resize(nparticles);
    particles.velocities.resize(nparticles);
    particles.accelerations.resize(nparticles);
    particles.mus.resize(nparticles);
    particles.radii.resize(nparticles);
    particles.ids.resize(nparticles);
    particles.test_mass.resize(nparticles);

    infile.read(reinterpret_cast<char*> (particles.positions.data()), sizeof(Vec3) * nparticles);
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    infile.read(reinterpret_cast<char*> (particles.velocities.data()), sizeof(Vec3) * nparticles);
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    infile.read(reinterpret_cast<char*> (particles.accelerations.data()), sizeof(Vec3) * nparticles);
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    infile.read(reinterpret_cast<char*> (particles.mus.data()), sizeof(double) * nparticles);
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    infile.read(reinterpret_cast<char*> (particles.radii.data()), sizeof(double) * nparticles);
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    infile.read(reinterpret_cast<char*> (particles.ids.data()), sizeof(uint32_t) * nparticles);
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");
    infile.read(reinterpret_cast<char*> (particles.test_mass.data()), sizeof(uint8_t) * nparticles);
    if (!infile) throw std::runtime_error("Failed reading particle data from file.");

    infile.close();
  }
}