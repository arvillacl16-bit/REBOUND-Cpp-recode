#pragma once
#include "exception.hpp"

namespace rebound::repstl {
  template <typename T>
  class Vector {
  private:
    T* data;
    size_t length;
    size_t capacity;

    void resize_(size_t new_capacity) {
      T* new_data = new T[new_capacity];
      for (size_t i = 0; i < length; ++i) {
        new_data[i] = (T&&)(data[i]);
      }
      delete[] data;
      data = new_data;
      capacity = new_capacity;
    }
  public:
    using ThisVec = Vector<T>;
    using Iterator = T*;
    using ConstIterator = const T*;

    Vector() : data(nullptr), length(0), capacity(0) {}
    Vector(const ThisVec& other) {
      length = other.length;
      capacity = other.capacity;
      data = new T[capacity];
      for (size_t i = 0; i < length; ++i) {
        data[i] = other.data[i];
      }
    }
    Vector(ThisVec&& other) noexcept : data(other.data), length(other.length), capacity(other.capacity) {
      other.data = nullptr;
      other.length = 0;
      other.capacity = 0;
    }
    ThisVec& operator=(const ThisVec& other) {
      if (this != &other) {
        delete[] data;
        length = other.length;
        capacity = other.capacity;
        data = new T[capacity];
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
        capacity = other.capacity;
        other.data = nullptr;
        other.length = 0;
        other.capacity = 0;
      }
      return *this;
    }

    ~Vector() { delete[] data; }

    size_t size() const noexcept { return length; }
    size_t capacity() const noexcept { return capacity; }
    bool is_empty() const noexcept { return length == 0; }
    void clear() noexcept {
      delete[] data;
      data = nullptr;
      length = 0;
      capacity = 0;
    }

    Iterator begin() noexcept { return data; }
    Iterator end() noexcept { return data + length; }
    ConstIterator cbegin() const noexcept { return data; }
    ConstIterator cend() const noexcept { return data + length; }

    void push_back(const T& value) {
      if (length == capacity) {
        resize_(capacity == 0 ? 1 : capacity * 2);
      }
      data[length++] = value;
    }

    T pop_back() {
      if (length == 0) return T();
      length--;
      auto val = data[length];
      data[length].~T();
      return val;
    }

    T& operator[](size_t index) { return data[index]; }
    const T& operator[](size_t index) const { return data[index]; }

    T& at(size_t index) {
      if (index >= length) throw IndexOutOfBounds("Index out of range");
      return data[index];
    }

    const T& at(size_t index) const {
      if (index >= length) throw IndexOutOfBounds("Index out of range");
      return data[index];
    }

    void reserve(size_t new_capacity) {
      if (new_capacity > capacity) {
        resize_(new_capacity);
      }
    }

    void shrink_to_fit() {
      if (length < capacity) {
        resize_(length);
      }
    }

    void resize(size_t new_size, const T& default_value = T()) {
      if (new_size > capacity) {
        resize_(new_size);
      }
      for (size_t i = length; i < new_size; ++i) {
        data[i] = default_value;
      }
      length = new_size;
    }

    void assign(size_t new_size, const T& value) {
      if (new_size > capacity) {
        resize_(new_size);
      }
      for (size_t i = 0; i < new_size; ++i) {
        data[i] = value;
      }
      length = new_size;
    }
  };
}