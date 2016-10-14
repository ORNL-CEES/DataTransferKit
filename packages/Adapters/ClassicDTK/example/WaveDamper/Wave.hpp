//---------------------------------------------------------------------------//
/*!
 * \file Wave.hpp
 * \author Stuart R. Slattery
 * \brief Wave code declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EXAMPLE_WAVE_HPP
#define DTK_EXAMPLE_WAVE_HPP

#include <vector>

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

//---------------------------------------------------------------------------//

class Wave
{
  private:
    Teuchos::RCP<const Teuchos::Comm<int>> comm;
    Teuchos::RCP<std::vector<double>> grid;
    Teuchos::RCP<std::vector<double>> data;
    Teuchos::RCP<std::vector<double>> damping;

  public:
    Wave( Teuchos::RCP<const Teuchos::Comm<int>> _comm, double x_min,
          double x_max, int num_x );

    ~Wave();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> get_comm() const { return comm; }

    // Get a const reference to the local grid.
    const Teuchos::RCP<std::vector<double>> &get_grid() const { return grid; }

    // Get a reference to the local data.
    const Teuchos::RCP<std::vector<double>> &get_data() const { return data; }

    // Get a reference to the local data space storing the damping
    // coefficients.
    Teuchos::RCP<std::vector<double>> &get_damping() { return damping; }

    // Solve the local problem.
    void solve();
};

//---------------------------------------------------------------------------//

#endif // DTK_EXAMPLE_WAVE_HPP

//---------------------------------------------------------------------------//
// end Wave.hpp
//---------------------------------------------------------------------------//
