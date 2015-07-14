  // dependency: {arg}
  {var}_key_ = plist_.get<std::string>("{argString} key",
          domain_name+"{arg}");
  dependencies_.insert({var}_key_);
