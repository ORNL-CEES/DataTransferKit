//---------------------------------------------------------------------------//
/*!
 * \file Wave.cpp
 * \author Stuart R. Slattery
 * \brief Wave code definition.
 */
//---------------------------------------------------------------------------//

#include "Wave.hpp"

#include <algorithm>
#include <cmath>

//---------------------------------------------------------------------------//
Wave::Wave( Teuchos::RCP<const Teuchos::Comm<int>> _comm, double x_min,
            double x_max, int num_x )
    : comm( _comm )
    , grid( new std::vector<double>( num_x ) )
    , data( new std::vector<double>( num_x ) )
    , damping( new std::vector<double>( num_x ) )
{
    // Create the grid.
    double x_size = ( x_max - x_min ) / ( num_x - 1 );

    std::vector<double>::iterator grid_iterator;
    int i = 0;

    for ( grid_iterator = grid->begin(); grid_iterator != grid->end();
          ++grid_iterator, ++i )
    {
        *grid_iterator = i * x_size + x_min;
    }

    // Set initial conditions.
    std::fill( damping->begin(), damping->end(), 0.0 );
    std::vector<double>::iterator data_iterator;
    for ( data_iterator = data->begin(), grid_iterator = grid->begin();
          data_iterator != data->end(); ++data_iterator, ++grid_iterator )
    {
        *data_iterator = std::cos( *grid_iterator );
    }
}

//---------------------------------------------------------------------------//
Wave::~Wave() { /* ... */}

//---------------------------------------------------------------------------//
void Wave::solve()
{
    // Apply the dampened component.
    std::vector<double>::iterator data_iterator;
    std::vector<double>::const_iterator damping_iterator;
    for ( data_iterator = data->begin(), damping_iterator = damping->begin();
          data_iterator != data->end(); ++data_iterator, ++damping_iterator )
    {
        *data_iterator -= *damping_iterator;
    }
}

//---------------------------------------------------------------------------//
// end Wave.cpp
//---------------------------------------------------------------------------//
