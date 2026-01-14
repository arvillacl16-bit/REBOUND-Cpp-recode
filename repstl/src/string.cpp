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

#include "../string"
#include "../exception"

namespace rebound::repstl {
  void String::grow(size_t cap_, bool allow_shrink) {
    if (cap_ == cap) return;
    if (cap_ > cap || allow_shrink) {
      size_t new_len = length > cap_ ? cap_ : length;
      char* new_data = new char[cap_ + 1];
      if (data) std::memcpy(new_data, data, new_len);
      new_data[new_len] = '\0';
      delete[] data;
      data = new_data;
      length = new_len;
      cap = cap_;
    }
  }

  void String::ensure_capacity(size_t min_cap) {
    if (cap >= min_cap) return;
    size_t new_cap = cap == 0 ? 8 : cap;
    while (new_cap < min_cap) new_cap = ceil(new_cap * GROW_FACTOR);
    grow(new_cap);
  }

  String::String(const char* str) {
    if (str) {
      length = std::strlen(str);
      cap = length;
      data = new char[cap + 1];
      std::strcpy(data, str);
    } else {
      data = nullptr;
      length = 0;
      cap = 0;
    }
  }

  String::String(const String& other) {
    if (other.data) {
      length = other.length;
      cap = other.cap;
      data = new char[cap + 1];
      std::strcpy(data, other.data);
    } else {
      data = nullptr;
      length = 0;
      cap = 0;
    }
  }

  String& String::operator=(const String& other) {
    if (this != &other) {
      length = other.length;
      cap = other.cap;
      char* new_data = new char[cap + 1];
      delete[] data;
      data = new_data;
      if (other.data) std::strcpy(data, other.data);
      else data[0] = '\0';
    }
    return *this;
  }

  String& String::operator=(String&& other) noexcept {
    if (this != &other) {
      delete[] data;
      data = other.data;
      cap = other.cap;
      length = other.length;
      other.data = nullptr;
      other.length = 0;
      other.cap = 0;
    }
    return *this;
  }

  void String::push_back(char c) {
    if (length == cap) {
      size_t new_cap = cap == 0 ? 8 : ceil(cap * GROW_FACTOR);
      grow(new_cap);
    }
    data[length] = c;
    data[length + 1] = '\0';
    ++length;
  }

  void String::pop_back() {
    if (length == 0) return;
    --length;
    data[length] = '\0';
  }

  void String::del_char_at(size_t index) {
    if (index >= length) return;
    std::memmove(data + index, data + index + 1, length - index);
    --length;
  }

  void String::del_char(char c) {
    if (!data) return;
    size_t new_length = 0;
    for (size_t i = 0; i < length; ++i) {
      if (data[i] != c) {
        data[new_length++] = data[i];
      }
    }
    data[new_length] = '\0';
    length = new_length;
  }

  const char& String::at(size_t index) const {
    if (index >= length) throw repstl::IndexOutOfBounds("Index is outside of string");
    return data[index];
  }

  char& String::at(size_t index) {
    if (index >= length) throw repstl::IndexOutOfBounds("Index is outside of string");
    return data[index];
  }

  String String::operator+(const String& other) const {
    if (data && !other.data) return *this;
    else if (!data && other.data) return other;
    else if (!data && !other.data) return "";
    String result;
    result.length = length + other.length;
    result.ensure_capacity(length + other.length);
    std::memcpy(result.data, data, length);
    std::memcpy(result.data + length, other.data, other.length + 1);
    return result;
  }

  String String::operator+(const char* str) const {
    if (data && !str) return *this;
    if (!data && str) return str;
    if (!data && !str) return "";
    String result;
    size_t len = std::strlen(str);
    result.length = len + length;
    result.ensure_capacity(length + len);
    std::strcpy(result.data, data);
    std::strcpy(result.data + length, str);
    return result;
  }

  String& String::operator+=(const String& other) {
    if (!other.data) return *this;
    ensure_capacity(length + other.length);
    std::strcpy(data + length, other.data);
    length += other.length;
    return *this;
  }

  String& String::operator+=(const char* str) {
    if (!str) return *this;
    size_t len = std::strlen(str);
    ensure_capacity(length + len);
    std::strcpy(data + length, str);
    length += len;
    return *this;
  }

  inline unsigned long long int_pow(unsigned long long base, unsigned long long exp) {
    unsigned long long res = 1;
    for (unsigned long long i = 0; i < exp; ++i) res *= base;
    return res;
  }

  String ull_to_string(unsigned long long val) {
    if (val == 0) return "0";
    String res;

    unsigned long long exp = 0;
    while (val >= int_pow(10, exp)) ++exp;
    while (exp > 0) {
      --exp;
      unsigned long long pow_of_10 = int_pow(10, exp);
      int digit = val / pow_of_10;
      val -= pow_of_10 * digit;
      res.push_back('0' + digit);
    }

    return res;
  }
}
