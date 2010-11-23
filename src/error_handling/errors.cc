#include "errors.hh"

namespace Errors
{

Message::~Message () throw() {  }

Message& operator<<(Message &message, const char* data) { message.add_data (data); return message; }


}

