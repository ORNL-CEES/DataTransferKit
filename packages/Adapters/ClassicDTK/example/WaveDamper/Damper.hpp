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
    Teuchos::RCP<std::vector<double> > data;
    Teuchos::RCP<std::vector<double> > damping;
    Teuchos::RCP<std::vector<double> > grid;

  public:

    Damper( Teuchos::RCP<const Teuchos::Comm<int> > _comm,
            double x_min, double x_max, int num_x);

    ~Damper();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > get_comm() const
    { return comm; }

    // Get a reference to the local damping data.
    Teuchos::RCP<std::vector<double> > get_damping() const
    { return damping; }

    // Get a reference to the local grid.
    Teuchos::RCP<std::vector<double> > get_grid() const
    { return grid; }

    // Get a reference to the memory space for external data to be applied to.
    Teuchos::RCP<std::vector<double> > get_external_data()
    { return data; }

    // Apply damping to the local problem.
    void solve();
};

//---------------------------------------------------------------------------//

#endif // DTK_EXAMPLE_DAMPER_HPP

//---------------------------------------------------------------------------//
// end Damper.hpp
//---------------------------------------------------------------------------//

