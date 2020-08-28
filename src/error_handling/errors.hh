#ifndef _ERRORS_H_
#define _ERRORS_H_

#include "exceptions.hh"

#include <sstream>

namespace Errors {

class Message : public Exceptions::Amanzi_exception {
 public:
  explicit Message() : message_() {};
  explicit Message(const char* message) : message_(message) {};
  explicit Message(const std::string& message) : message_(message) {};
  ~Message() noexcept {};

  const char* what() const noexcept override { return message_.c_str(); }

  void add_data(const char* data) { message_ += data; }
  void add_data(const std::string& data) { message_ += data; }

 public:
  std::string message_;
};

Message& operator<<(Message &message, const char* data);
Message& operator<<(Message &message, const std::string& data);
Message& operator<<(Message &message, double datum);
Message& operator<<(Message &message, int datum);
Message& operator<<(Message &message, std::size_t datum);

class CutTimeStep : public Exceptions::Amanzi_exception {};

}  // namespace Errors
#endif /* _ERRORS_H_ */
