//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file DTK_TargetManager.hpp
 * \author Stuart R. Slattery
 * \brief Target manager declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TARGETMANAGER_HPP
#define DTK_TARGETMANAGER_HPP

#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_MeshManager.hpp>

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class TargetManager
 * \brief Manager for targets.
 */
//---------------------------------------------------------------------------//
template<class CoordinateField, class TargetField>
class TargetManager
{
  public:

    //@{
    //! Typedefs.
    typedef CoordinateField                           coordinate_field_type;
    typedef FieldManager<CoordinateField>             CoordManagerType;
    typedef Teuchos::RCP<CoordManagerType>            RCP_CoordManager;
    typedef TargetField                               target_field_type;
    typedef FieldManager<TargetField>                 TargetManagerType;
    typedef Teuchos::RCP<TargetManagerType>           RCP_TargetManager;
    //@}

    // Constructor.
    TargetManager( const RCP_CoordManager& target_coord_manager,
		   const RCP_TargetManager& target_space_manager );

    // Destructor.
    ~TargetManager();

    //! Get the coordinate field manager.
    const RCP_CoordManager& coords() const
    { return d_coord_manager; }

    //! Get the target field manager.
    const RCP_TargetManager& target() const
    { return d_target_manager; }

  private:

    // Coordinate field manager.
    RCP_CoordManager d_coord_manager;

    // Target field manager.
    RCP_TargetManager d_target_manager;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_TargetManager.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_TARGETMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_TargetManager.hpp
//---------------------------------------------------------------------------//
