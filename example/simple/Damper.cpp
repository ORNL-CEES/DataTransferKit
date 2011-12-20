#include "Wave.hpp"

//---------------------------------------------------------------------------//
Damper::Damper()
{ /* ... */ }

//---------------------------------------------------------------------------//
Damper::~Damper()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Apply damping to the local problem.
void Damper::solve()
{
    std::vector<double>::iterator damping_iterator;
    std::vector<double>::const_iterator wave_data_iterator;
    for (damping_iterator = damping.begin(),
       wave_data_iterator = wave_data.begin();
	 damping_iterator != damping.end();
	 ++damping_iterator, ++wave_data_iterator)
    {
	*damping_iterator = *wave_data_iterator / 2;
    }
}

//---------------------------------------------------------------------------//
