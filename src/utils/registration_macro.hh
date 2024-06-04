#define REGISTER(c) Utils::RegisteredFactory<Evaluator, c> c::reg_(c::eval_type)

#define REGISTER_MODEL(c)                                                                          \
  template <>                                                                                      \
  Utils::RegisteredFactory<Evaluator, EvaluatorModelCV<c, DefaultDevice>>                          \
  EvaluatorModelCV<c, DefaultDevice>::reg_(c<cView_type, View_type>::eval_type)


#define REGISTER_BY_MATERIAL(c)                                                                    \
  template <>                                                                                      \
  Utils::RegisteredFactory<Evaluator, EvaluatorModelCVByMaterial<c, DefaultDevice>>                \
  EvaluatorModelCVByMaterial<c, DefaultDevice>::reg_(c<cView_type, View_type>::eval_type +         \
                                                     " by region")


#define REGISTER_PK(c) RegisteredPKFactory<c> c::reg_(c::pk_type_)
