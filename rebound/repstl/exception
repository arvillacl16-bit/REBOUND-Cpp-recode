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

namespace rebound::repstl {
  class Exception {
  private:
    const char* message;
  public:
    constexpr Exception(const char* msg) noexcept : message(msg) {}
    constexpr const char* what() const noexcept { return message; }

    virtual ~Exception() = default;
  };

  class IndexOutOfBounds : public Exception {
  public:
    using Exception::Exception;
  };

  class RuntimeError : public Exception {
  public:
    using Exception::Exception;
  };

  class InvalidArgument : public Exception {
  public:
    using Exception::Exception;
  };
}