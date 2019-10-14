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

#include <deal.II/base/config.h>

#include <deal.II/base/iterator_range.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/tria.h>

#include <half_struct.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

#ifndef dealii_general_mesh_h
#  define dealii_general_mesh_h

template <int dim, int spacedim>
void
output_vertices(const Triangulation<dim, spacedim> &tria,
                const std::string &                 basename)
{
  std::ofstream out(basename + "_vertices.txt");
  for (const auto &v : tria.get_vertices())
    {
      out << v;
      for (unsigned int j = spacedim; j < 3; ++j)
        out << " 0.0" << std::endl;
    }
}


template <int dim, class IndexType = unsigned int>
class GeneralMeshConnectivity
{
public:
  template <int spacedim>
  GeneralMeshConnectivity(const Triangulation<dim, spacedim> &tria);

  template <int odim>
  IndexType
  new_object(const std::vector<IndexType> &sub_objects,
             const IndexType &twin = HalfStruct<IndexType>::invalid_index,
             const IndexType &next = HalfStruct<IndexType>::invalid_index);

  inline IndexType
  n_vertices() const
  {
    return connectivity[0].size();
  }


  inline IndexType
  n_faces() const
  {
    Assert(dim > 0, ExcMessage("Can't get faces in zero dimensions"));
    return connectivity[dim - 1].size();
  }


  inline IndexType
  n_cells() const
  {
    Assert(dim > 0, ExcMessage("Can't get faces in zero dimensions"));
    return connectivity[dim].size();
  }

  inline IndexType
  n_edges() const
  {
    Assert(dim > 0, ExcMessage("Can't get edges in zero dimensions"));
    return connectivity[1].size();
  }

  inline IndexType
  n_polygons() const
  {
    Assert(dim > 1, ExcMessage("Can't get polygons in < 2 dimensions"));
    return connectivity[2].size();
  }

  inline IndexType
  n_volumes() const
  {
    Assert(dim > 2, ExcMessage("Can't get volumes in < 3 dimensions"));
    return connectivity[3].size();
  }


  void
  output_connectivities(const std::string base_name) const;


  void
  pre_compute_subdivision_vertices()
  {
    // All points that participate in the given new cell point
    std::vector<std::vector<IndexType>> volume_point_indices;

    // All weights used to create a new cell point
    std::vector<std::vector<double>> volume_point_weights;

    // All points that participate in the given new face point
    std::vector<std::vector<IndexType>> polygon_point_indices;

    // All weights used to create a new face point
    std::vector<std::vector<double>> polygon_point_weights;

    if (dim == 2)
      {
        polygon_point_indices.resize(n_polygons());
        polygon_point_weights.resize(n_polygons());

        for (unsigned int i = 0; i < connectivity[2].size(); ++i)
          {
            const auto &edges = connectivity[2][i].owned_objects;
            const auto  ne    = edges.size();
            for (unsigned int j = 0; j < ne; ++j)
              {
                polygon_point_weights[i].push_back(1.0 / ne);
                polygon_point_indices[i].push_back(edges[j][0]);
              }
          }
      }
  }

private:
  /**
   * Store all connectivities.
   *
   * index = 0: half_vertices
   * index = 1: half_edges
   * index = 2: half_faces
   * index = 3: half_cells
   *
   * etc.
   */
  std::array<std::vector<HalfStruct<IndexType>>, dim + 1> connectivity;
};



template <int dim, class IndexType>
template <int spacedim>
GeneralMeshConnectivity<dim, IndexType>::GeneralMeshConnectivity(
  const Triangulation<dim, spacedim> &tria)
{
  connectivity[0].resize(tria.n_vertices());
  unsigned int nv = 0;
  for (auto &v : connectivity[0])
    v.owned_objects.push_back(nv++);

  if (dim == 1)
    for (const auto &cell : tria.active_cell_iterators())
      new_object<1>({cell->vertex_index(0), cell->vertex_index(1)});
  if (dim == 2)
    for (const auto &cell : tria.active_cell_iterators())
      {
        auto id0 =
          new_object<1>({cell->vertex_index(0), cell->vertex_index(1)});
        auto id1 =
          new_object<1>({cell->vertex_index(1), cell->vertex_index(3)});
        auto id2 =
          new_object<1>({cell->vertex_index(3), cell->vertex_index(2)});
        auto id3 =
          new_object<1>({cell->vertex_index(2), cell->vertex_index(0)});
        new_object<2>({id0, id1, id2, id3});
      }
}



template <int dim, class IndexType>
void
GeneralMeshConnectivity<dim, IndexType>::output_connectivities(
  const std::string base_name) const
{
  {
    std::ofstream out(base_name + "_edges.txt");
    for (auto edge : connectivity[1])
      out << edge[0] << " " << edge[1] << std::endl;
  }
  if (dim > 1)
    {
      std::ofstream out_face_edges(base_name + "_face_edges.txt");
      std::ofstream out_face_vertices(base_name + "_face_vertices.txt");
      std::ofstream out_face_starts(base_name + "_face_start.txt");
      unsigned int  s = 0;
      for (auto face : connectivity[2])
        {
          out_face_starts << s << " ";
          s += face.owned_objects.size();
          for (auto edge_id : face.owned_objects)
            {
              auto edge = connectivity[1][edge_id];
              out_face_vertices << edge[0] << " ";
              out_face_edges << edge_id << " ";
            }
        }
      out_face_vertices << std::endl;
      out_face_edges << std::endl;
      out_face_starts << std::endl;
    }
  if (dim > 2)
    {
      std::ofstream out_cell_faces(base_name + "_cell_faces.txt");
      std::ofstream out_cell_starts(base_name + "_cell_start.txt");
      unsigned int  s = 0;
      for (auto cell : connectivity[3])
        {
          out_cell_starts << s << " ";
          s += cell.owned_objects.size();
          for (auto face_id : cell.owned_objects)
            {
              auto face = connectivity[2][face_id];
              out_cell_faces << face_id << " ";
            }
        }
      out_cell_faces << std::endl;
      out_cell_starts << std::endl;
    }
}

template <int dim, class IndexType>
template <int odim>
IndexType
GeneralMeshConnectivity<dim, IndexType>::new_object(
  const std::vector<IndexType> &sub_objects,
  const IndexType &             twin,
  const IndexType &             next)
{
  IndexType new_id = connectivity[odim].size();
  connectivity[odim].emplace_back(HalfStruct<IndexType>(sub_objects, {}));
  if (odim > 0)
    for (const auto &sub : sub_objects)
      connectivity[odim - 1][sub].owning_objects.push_back(new_id);
  return new_id;
}

DEAL_II_NAMESPACE_CLOSE

#endif // HALF_STRUCT_H
