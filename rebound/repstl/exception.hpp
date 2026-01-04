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