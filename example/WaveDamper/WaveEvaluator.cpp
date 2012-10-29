//---------------------------------------------------------------------------//
/*!
 * \file WaveEvaluator.cpp
 * \author Stuart R. Slattery
 * \brief Wave code evaluators for DTK.
 */
//---------------------------------------------------------------------------//

#include "WaveEvaluator.hpp"

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
WaveEvaluator::WaveEvaluator( const RCP_Wave& wave )
    : d_wave( wave )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*! 
  \brief Destructor.
*/
WaveEvaluator::~WaveEvaluator()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Wave function evaluator.
 */
WaveEvaluator::field_type
WaveEvaluator::evaluate(
    const Teuchos::ArrayRCP<int>& elements,
    const Teuchos::ArrayRCP<double>& coords )
{
    // Get the process rank.
    int my_rank = d_wave->get_comm()->getRank();
    
    // Setup a field to apply the data to.
    int field_size = elements.size();
    Teuchos::ArrayRCP<double> eval_data( field_size );
    DataTransferKit::FieldContainer<double> evaluations( eval_data, 1 );

    // Get the grid.
    Teuchos::RCP<std::vector<double> > grid = d_wave->get_grid();
    int num_grid_elements = grid->size() - 1;

    // Get the local element ids.
    int min_local_element_id = my_rank*num_grid_elements;
    int max_local_element_id = (my_rank+1)*num_grid_elements-1;

    // Get the wave data.
    Teuchos::RCP<std::vector<double> > data = d_wave->get_data();

    // Interpolate the local solution onto the given coordinates.
    Teuchos::ArrayRCP<int>::iterator element_iterator;
    int eval_id, element_id;
    for ( element_iterator = elements.begin();
	  element_iterator != elements.end();
	  ++element_iterator )
    {
	element_id = *element_iterator - min_local_element_id;
	eval_id = std::distance( elements.begin(), element_iterator );

	// If this is a valid element id for this process, interpolate with a
	// linear basis.
	if ( *element_iterator >= min_local_element_id &&
	     *element_iterator <= max_local_element_id )
	{
	    eval_data[eval_id] = (*data)[element_id]
				    + ( coords[eval_id] - (*grid)[element_id] )
				    * ( (*data)[element_id+1] - (*data)[element_id] )
				    / ( (*grid)[element_id+1] - (*grid)[element_id] );
	}

	// Otherwise put a zero for this element/coordinate pair per the
	// FieldEvaluator documentation.
	else
	{
	    eval_data[eval_id] = 0.0;
	}
    }

    // Return the field.
    return evaluations;
}

//---------------------------------------------------------------------------//
// end WaveEvaluator.cpp
//---------------------------------------------------------------------------//

