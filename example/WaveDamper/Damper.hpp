//---------------------------------------------------------------------------//
/*!
 * \file Damper.hpp
 * \author Stuart R. Slattery
 * \brief Damper code declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EXAMPLE_DAMPER_HPP
#define DTK_EXAMPLE_DAMPER_HPP

#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

//---------------------------------------------------------------------------//

class Damper
{
  private:

    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    std::vector<double> data;
    std::vector<double> damping;
    std::vector<double> grid;

  public:

    Damper( Teuchos::RCP<const Teuchos::Comm<int> > _comm,
	    double x_min, double x_max, int num_x);

    ~Damper();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > get_comm()
    {
	return comm;
    }

    // Get a reference to the local damping data.
    const std::vector<double>& get_damping()
    {
	return damping;
    }

    // Get a reference to the local grid.
    const std::vector<double>& get_grid()
    {
	return grid;
    }

    // Get a reference to the memory space for external data to be applied to.
    std::vector<double>& get_external_data()
    {
	return wave_data;
    }

    // Apply damping to the local problem.
    void solve();
};

//---------------------------------------------------------------------------//

#endif // DTK_EXAMPLE_DAMPER_HPP

//---------------------------------------------------------------------------//
// end Damper.hpp
//---------------------------------------------------------------------------//

