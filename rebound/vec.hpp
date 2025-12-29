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

#define LOG(msg) std::cout << msg << std::endl

#if !(defined(__cplusplus))
#error "\n\
────────────────────────────────────────────────────────────────────\n\
  🚫  FATAL ERROR: C COMPILER DETECTED\n\
────────────────────────────────────────────────────────────────────\n\
  You are attempting to compile a *C++* physics engine using a\n\
  *C* compiler. This is not merely incorrect — it is unnecessary.\n\
\n\
  If you want a physics engine **in C**, you should be using the\n\
  real REBOUND: the mature, battle-tested, widely-cited N-body\n\
  integrator written in actual C.\n\
\n\
  What you are trying to compile here is a C++ rewrite that is\n\
  objectively worse, slower, younger, buggier, and generally\n\
  inferior to C REBOUND.\n\
\n\
  So not only is this file incompatible with C — there is also\n\
  absolutely no valid reason to compile it as C.\n\
\n\
  Please switch to:\n\
     • a C++ compiler (g++, clang++, MSVC)\n\
     • OR, better yet, use real REBOUND like a sane person.\n\
────────────────────────────────────────────────────────────────────\n\
"
#else
#include <cmath>

namespace rebound {
  struct Vec3 {
    double x, y, z;

    constexpr Vec3() : x(0), y(0), z(0) {}
    constexpr explicit Vec3(double x_) : x(x_), y(0), z(0) {}
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
    constexpr Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
    constexpr Vec3 operator*(double scalar) const { return {x * scalar, y * scalar, z * scalar}; }
    constexpr Vec3 operator/(double scalar) const { return {x / scalar, y / scalar, z / scalar}; }

    inline Vec3& operator+=(const Vec3& other) { x += other.x; y += other.y; z += other.z; return *this; }
    inline Vec3& operator-=(const Vec3& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
    inline Vec3& operator*=(double scalar) { x *= scalar; y *= scalar; z *= scalar; return *this; }
    inline Vec3& operator/=(double scalar) { x /= scalar; y /= scalar; z /= scalar; return *this; }

    constexpr Vec3 operator-() const { return {-x, -y, -z}; }

    constexpr double dot(const Vec3& other) const { return x * other.x + y * other.y + z * other.z; }
    constexpr double mag2() const { return dot(*this); }
    double mag() const { return std::sqrt(mag2()); }
    constexpr double distance2(const Vec3& other) const { return (*this - other).mag2(); }
    double distance(const Vec3& other) const { return (*this - other).mag(); }

    constexpr Vec3 cross(const Vec3& other) const { return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x}; }

    constexpr bool operator==(const Vec3& other) const { return x == other.x && y == other.y && z == other.z; }
    constexpr bool operator!=(const Vec3& other) const { return !operator==(other); }
  };

  constexpr Vec3 operator*(double scalar, const Vec3& vec) { return vec * scalar; }
}

#pragma omp declare reduction( \
  vec3_plus : rebound::Vec3 : omp_out += omp_in ) \
  initializer(omp_priv = rebound::Vec3{})
#endif