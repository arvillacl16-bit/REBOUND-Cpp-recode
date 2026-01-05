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
#include "exception.hpp"
#include <initializer_list>
#include <stddef.h>

namespace rebound::repstl {
  template <typename T>
  class Vector {
  private:
    T* data;
    size_t length;
    size_t capacity_;

    void resize_(size_t new_capacity) {
      T* new_data = new T[new_capacity];
      for (size_t i = 0; i < length; ++i) {
        new_data[i] = static_cast<T&&>(data[i]);
      }
      delete[] data;
      data = new_data;
      capacity_ = new_capacity;
    }

  public:
    using ThisVec = Vector<T>;
    using Iterator = T*;
    using ConstIterator = const T*;

    Vector() : data(nullptr), length(0), capacity_(0) {}
    Vector(const ThisVec& other) {
      length = other.length;
      capacity_ = other.capacity_;
      data = new T[capacity_];
      for (size_t i = 0; i < length; ++i) {
        data[i] = other.data[i];
      }
    }
    Vector(ThisVec&& other) noexcept
      : data(other.data), length(other.length), capacity_(other.capacity_) {
      other.data = nullptr;
      other.length = 0;
      other.capacity_ = 0;
    }
    ThisVec& operator=(const ThisVec& other) {
      if (this != &other) {
        delete[] data;
        length = other.length;
        capacity_ = other.capacity_;
        data = new T[capacity_];
        for (size_t i = 0; i < length; ++i) {
          data[i] = other.data[i];
        }
      }
      return *this;
    }
    ThisVec& operator=(ThisVec&& other) noexcept {
      if (this != &other) {
        delete[] data;
        data = other.data;
        length = other.length;
        capacity_ = other.capacity_;
        other.data = nullptr;
        other.length = 0;
        other.capacity_ = 0;
      }
      return *this;
    }

    ~Vector() { delete[] data; }

    Vector(std::initializer_list<T> init_list)
      : data(nullptr), length(0), capacity_(0) {
      reserve(init_list.size());
      for (const T& value : init_list) {
        push_back(value);
      }
    }

    size_t size() const noexcept { return length; }
    size_t capacity() const noexcept { return capacity_; }
    bool is_empty() const noexcept { return length == 0; }
    void clear() noexcept {
      delete[] data;
      data = nullptr;
      length = 0;
      capacity_ = 0;
    }

    Iterator begin() noexcept { return data; }
    Iterator end() noexcept { return data + length; }
    ConstIterator cbegin() const noexcept { return data; }
    ConstIterator cend() const noexcept { return data + length; }

    void push_back(const T& value) {
      if (length == capacity_) {
        resize_(capacity_ == 0 ? 1 : capacity_ * 2);
      }
      data[length++] = value;
    }

    T pop_back() {
      if (length == 0)
        return T();
      T value = static_cast<T&&>(data[length - 1]);
      data[length - 1].~T();
      --length;
      return value;
    }

    T& operator[](size_t index) { return data[index]; }
    const T& operator[](size_t index) const { return data[index]; }

    T& at(size_t index) {
      if (index >= length)
        throw IndexOutOfBounds("Index out of range");
      return data[index];
    }

    const T& at(size_t index) const {
      if (index >= length)
        throw IndexOutOfBounds("Index out of range");
      return data[index];
    }

    void reserve(size_t new_capacity) {
      if (new_capacity > capacity_) {
        resize_(new_capacity);
      }
    }

    void shrink_to_fit() {
      if (length < capacity_) {
        resize_(length);
      }
    }

    void resize(size_t new_size, const T& default_value = T()) {
      if (new_size > capacity_) {
        resize_(new_size);
      }
      for (size_t i = length; i < new_size; ++i) {
        data[i] = default_value;
      }
      length = new_size;
    }

    void assign(size_t new_size, const T& value) {
      if (new_size > capacity_) {
        resize_(new_size);
      }
      for (size_t i = 0; i < new_size; ++i) {
        data[i] = value;
      }
      length = new_size;
    }

    bool operator==(const ThisVec& other) const {
      if (length != other.length)
        return false;
      for (size_t i = 0; i < length; ++i)
        if (other.data[i] != data[i])
          return false;
      return true;
    }
  };
} // namespace rebound::repstl