#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <half_struct.h>

#include <fstream>

using namespace dealii;
int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                            argv,
                                            numbers::invalid_unsigned_int);

  ParameterHandler  prm;
  const std::string parameter_file = "param.prm";
  prm.parse_input(parameter_file);


  std::ofstream logfile;
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      logfile.open("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_file(10);
    }

  deallog << "OK" << std::endl;
}
