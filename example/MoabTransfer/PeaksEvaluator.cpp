//---------------------------------------------------------------------------//
/*!
 * \file PeaksEvaluator.cpp
 * \author Stuart R. Slattery
 * \brief Peaks function evaluator definition.
 */
//---------------------------------------------------------------------------//

#include <cmath>

#include "PeaksEvaluator.hpp"

ArrayField PeaksEvaluator::evaluate(  
    const Teuchos::ArrayRCP<global_ordinal_type>& elements,
    const Teuchos::ArrayRCP<double>& coords )
{
    ArrayField evaluations( elements.size(), 1 );

    // If the element is in the local partition, average the tagged values of
    // the source field for that element.
    for ( int n = 0; n < elements.size(); ++n )
    {
	// Insert MOAB code here....
	double x = coords[n];
	double y = coords[ elements.size() + n ];
	evaluations[n] = 3*(1-x)*(1-x)*std::exp(-x*x-(y+1)*(y+1))
			 -10*(x/5 - x*x*x - y*y*y*y*y)*std::exp(-x*x-y*y)
			 -(1/3)*std::exp(-(x+1)*(x+1)-y*y);
    }

    return evaluations;
}

//---------------------------------------------------------------------------//
// end PeaksEvaluator.cpp
//---------------------------------------------------------------------------//
