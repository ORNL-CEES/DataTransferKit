//---------------------------------------------------------------------------//
/*!
 * \file Damper.cpp
 * \author Stuart R. Slattery
 * \brief Damper code definition.
 */
//---------------------------------------------------------------------------//

#include <cmath>

#include "Damper.hpp"

//---------------------------------------------------------------------------//
Damper::Damper(Teuchos::RCP<const Teuchos::Comm<int> > _comm,
               double x_min, double x_max, int num_x)
    : comm(_comm)
    , data( new std::vector<double>(num_x) )
    , damping( new std::vector<double>(num_x) )
    , grid( new std::vector<double>(num_x) )
{
    // Create the grid.
    double x_size = (x_max - x_min) / (num_x-1);

    std::vector<double>::iterator grid_iterator;
    int i = 0;

    for (grid_iterator = grid->begin();
         grid_iterator != grid->end();
         ++grid_iterator, ++i)
    {
        *grid_iterator = i*x_size + x_min;
    }
}

//---------------------------------------------------------------------------//
Damper::~Damper()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Apply damping to the local problem.
void Damper::solve()
{
    std::vector<double>::iterator damping_iterator;
    std::vector<double>::const_iterator data_iterator;
    for (damping_iterator = damping->begin(),
            data_iterator = data->begin();
         damping_iterator != damping->end();
         ++damping_iterator, ++data_iterator)
    {
        *damping_iterator = *data_iterator / 2;
    }
}

//---------------------------------------------------------------------------//
// end Damper.cpp
//---------------------------------------------------------------------------//
