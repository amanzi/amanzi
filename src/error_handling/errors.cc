/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "errors.hh"
#include <cstdio>

namespace Errors {

Message&
operator<<(Message& message, const char* data)
{
  message.add_data(data);
  return message;
}

Message&
operator<<(Message& message, const std::string& data)
{
  message.add_data(data);
  return message;
}

Message&
operator<<(Message& message, int datum)
{
  char number[24];
  snprintf(number, 23, "%d", datum);
  message.add_data(number);
  return message;
}

Message&
operator<<(Message& message, double datum)
{
  char number[24];
  snprintf(number, 23, "%g", datum);
  message.add_data(number);
  return message;
}

Message&
operator<<(Message& message, unsigned long datum)
{
  char number[24];
  snprintf(number, 23, "%lu", datum);
  message.add_data(number);
  return message;
}

Message&
operator<<(Message& message, unsigned datum)
{
  char number[24];
  snprintf(number, 23, "%u", datum);
  message.add_data(number);
  return message;
}


} // namespace Errors
