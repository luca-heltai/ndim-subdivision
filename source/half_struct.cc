#include <half_struct.h>

DEAL_II_NAMESPACE_OPEN;

template <>
const unsigned int
  HalfStruct<unsigned int>::invalid_index = static_cast<unsigned int>(-1);

template struct HalfStruct<int>;
template struct HalfStruct<long unsigned int>;

DEAL_II_NAMESPACE_CLOSE;