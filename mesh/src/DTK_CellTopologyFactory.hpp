//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file CellTopologyFactory.hpp
 * \author Stuart Slattery
 * \brief Factory method declaration for cell topology data.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_CELLTOPOLOGYFACTORY_HPP
#define FOOD_CELLTOPOLOGYFACTORY_HPP

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

namespace FOOD
{

class CellTopologyFactory
{
  public:

    // Consructor.
    CellTopologyFactory();

    // Destructor.
    ~CellTopologyFactory();

    // Factory method.
    Teuchos::RCP<shards::CellTopology> create( const int entity_topology,
					       const int num_entity_nodes );
};

} // end namespace FOOD

#endif // end FOOD_CELLTOPOLOGYFACTORY_HPP

//---------------------------------------------------------------------------//
// end CellTopologyFactory.hpp
//---------------------------------------------------------------------------//
