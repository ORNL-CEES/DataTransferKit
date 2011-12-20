#ifndef simple_example_damper_hpp
#define simple_example_damper_hpp

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"

//---------------------------------------------------------------------------//

class Damper
{
  private:

    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    std::vector<double> wave_data;
    std::vector<double> damping;
    std::vector<double> grid;

  public:

    Damper(Teuchos::RCP<const Teuchos::Comm<int> > _comm,
	   double x_min,
	   double x_max,
	   int num_x);

    ~Damper();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > get_comm()
    {
	return comm;
    }

    // Get a reference to the local damping data.
    std::vector<double>& get_damping()
    {
	return damping;
    }

    // Get a const reference to the local grid.
    const std::vector<double>& get_grid()
    {
	return grid;
    }

    // Get the wave data to apply damping to from an external source.
    std::vector<double>& set_wave_data()
    {
	return wave_data;
    }

    // Apply damping to the local problem.
    void solve();
};

//---------------------------------------------------------------------------//

#endif // simple_example_damper_hpp

//---------------------------------------------------------------------------//
