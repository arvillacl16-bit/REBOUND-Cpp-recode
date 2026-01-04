#pragma once
#include <stddef.h>
#include <cstring>
#include "exception.hpp"

namespace rebound::repstl {
  class String {
  private:
    char* data;
    size_t length;
  public:
    String() : data(nullptr), length(0) {}
    String(const char* str) {
      length = std::strlen(str);
      data = new char[length + 1];
      std::strcpy(data, str);
    }
    String(const String& other) {
      length = other.length;
      data = new char[length + 1];
      std::strcpy(data, other.data);
    }
    String(String&& other) noexcept : data(other.data), length(other.length) {
      other.data = nullptr;
      other.length = 0;
    }

    String& operator=(const String& other) {
      if (this != &other) {
        delete[] data;
        length = other.length;
        data = new char[length + 1];
        std::strcpy(data, other.data);
      }
      return *this;
    }

    String& operator=(String&& other) noexcept {
      if (this != &other) {
        delete[] data;
        data = other.data;
        length = other.length;
        other.data = nullptr;
        other.length = 0;
      }
      return *this;
    }

    ~String() { delete[] data; }

    const char* c_str() const { return data ? data : ""; }

    const char &operator[](size_t index) const { return data[index]; }
    char &operator[](size_t index) { return data[index]; }

    void clear() {
      delete[] data;
      data = nullptr;
      length = 0;
    }

    void push_back(char c) {
      char* new_data = new char[length + 2];
      if (data) {
        std::strcpy(new_data, data);
        delete[] data;
      }
      new_data[length] = c;
      new_data[length + 1] = '\0';
      data = new_data;
      length++;
    }

    void pop_back() {
      if (length == 0) return;
      char* new_data = new char[length];
      std::strncpy(new_data, data, length - 1);
      new_data[length - 1] = '\0';
      delete[] data;
      data = new_data;
      length--;
    }

    void del_char_at(size_t index) {
      if (index >= length) return;
      char* new_data = new char[length];
      std::strncpy(new_data, data, index);
      std::strcpy(new_data + index, data + index + 1);
      delete[] data;
      data = new_data;
      length--;
    }

    size_t len() const { return length; }

    void del_char(char c) {
      size_t new_length = 0;
      for (size_t i = 0; i < length; ++i) {
        if (data[i] != c) {
          data[new_length++] = data[i];
        }
      }
      data[new_length] = '\0';
      length = new_length;
    }

    const char &at(size_t index) const {
      if (index >= length) throw IndexOutOfBounds("Index out of range");
      return data[index];
    }

    char &at(size_t index) {
      if (index >= length) throw IndexOutOfBounds("Index out of range");
      return data[index];
    }

    bool operator==(const String& other) const {
      if (length != other.length) return false;
      return std::strcmp(data, other.data) == 0;
    }

    bool operator!=(const String& other) const { return !(*this == other); }
    bool operator<(const String& other) const { return std::strcmp(c_str(), other.c_str()) < 0; }
    bool operator>(const String& other) const { return std::strcmp(c_str(), other.c_str()) > 0; }
    bool operator<=(const String& other) const { return !(*this > other); }
    bool operator>=(const String& other) const { return !(*this < other); }

    bool is_empty() const { return length == 0; }

    char* begin() { return data; }
    char* end() { return data + length; }
    const char* cbegin() const { return data; }
    const char* cend() const { return data + length; }
  };
}