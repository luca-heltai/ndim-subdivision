#include <deal.II/grid/grid_generator.h>

#include <general_mesh.h>
#include <gtest/gtest.h>
#include <half_struct.h>

using namespace dealii;

TEST(HalfStruct, Circle)
{
  Triangulation<2> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(2);

  GeneralMeshConnectivity<2> connectivity(tria);
  output_vertices(tria, "../../blender/circle");
  connectivity.output_connectivities("../../blender/circle");

  ASSERT_EQ(connectivity.n_vertices(), tria.n_vertices());
  ASSERT_EQ(connectivity.n_cells(), tria.n_active_cells());
  ASSERT_EQ(connectivity.n_cells(), connectivity.n_polygons());
}

TEST(HalfStruct, Ball)
{
  Triangulation<3> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(1);

  GeneralMeshConnectivity<3> connectivity(tria);
  output_vertices(tria, "../../blender/ball");
  connectivity.output_connectivities("../../blender/ball");

  ASSERT_EQ(connectivity.n_vertices(), tria.n_vertices());
  ASSERT_EQ(connectivity.n_cells(), tria.n_active_cells());
  ASSERT_EQ(connectivity.n_cells(), connectivity.n_volumes());
}
