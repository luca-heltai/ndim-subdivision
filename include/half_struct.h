// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------
#ifndef dealii_half_struct_h
#define dealii_half_struct_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

DEAL_II_NAMESPACE_OPEN


template <class IndexType = unsigned int>
struct HalfStruct
{
  HalfStruct() = default;

  HalfStruct(const std::vector<IndexType> &owned_objects,
             const std::vector<IndexType> &owning_objects,
             const IndexType &             twin = invalid_index,
             const IndexType &             next = invalid_index)
    : owned_objects(owned_objects)
    , owning_objects(owning_objects)
    , twin(twin)
    , next(next)
  {}

  HalfStruct(const HalfStruct<IndexType> &other) = default;

  /**
   * Return the index of the ith owned object.
   */
  IndexType operator[](const IndexType &i)
  {
    AssertIndexRange(i, owned_objects.size());
    return owned_objects[i];
  }

  static const IndexType invalid_index;

  std::vector<IndexType> owned_objects;
  std::vector<IndexType> owning_objects;

  IndexType twin = invalid_index;
  IndexType next = invalid_index;
};

DEAL_II_NAMESPACE_CLOSE

#endif // HALF_STRUCT_H
