#pragma once

namespace rebound::repstl {
  template <typename T1, typename T2>
  struct pair {
    T1 first;
    T2 second;

    pair() {}
    pair(const T1& first_, const T2& second_) : first(first_), second(second_) {}

    bool operator==(const pair &other) { return first == other.first && second == other.second; }
    bool operator!=(const pair &other) { return !operator==(other); }
  };
}