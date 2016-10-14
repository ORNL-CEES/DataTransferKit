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
 * \brief DTK_ConsistentInterpolationOperator.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATIONOPERATOR_HPP
#define DTK_CONSISTENTINTERPOLATIONOPERATOR_HPP

#include "DTK_MapOperator.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class ConsistentInterpolationOperator
  \brief Map operator interface.

  A map operator maps a field in one entity set to another entity set.
*/
//---------------------------------------------------------------------------//
class ConsistentInterpolationOperator : virtual public MapOperator
{
  public:
    //! Root class tyepdef.
    typedef MapOperator Base;
    typedef typename Base::Root Root;
    typedef typename Root::scalar_type Scalar;
    typedef typename Root::local_ordinal_type LO;
    typedef typename Root::global_ordinal_type GO;
    typedef typename Base::TpetraMultiVector TpetraMultiVector;
    typedef typename Base::TpetraMap TpetraMap;

    /*!
     * \brief Constructor.
     */
    ConsistentInterpolationOperator(
        const Teuchos::RCP<const TpetraMap> &domain_map,
        const Teuchos::RCP<const TpetraMap> &range_map,
        const Teuchos::ParameterList &parameters );

    /*!
     * \brief Return the ids of the range entities that were not mapped during
     * the last setup phase (i.e. those that are guaranteed to not receive
     * data from the transfer).
     *
     * \return A view of the ids.
     */
    Teuchos::ArrayView<const EntityId> getMissedRangeEntityIds() const;

  protected:
    /*
     * \brief Setup the map operator from a domain entity set and a range
     * entity set.
     *
     * \param domain_map Parallel map for domain vectors this map should be
     * compatible with.
     *
     * \param domain_function The function that contains the data that will be
     * sent to the range. Must always be nonnull but the pointers it contains
     * may be null of no entities are on-process.
     *
     * \param range_map Parallel map for range vectors this map should be
     * compatible with.
     *
     * \param range_space The function that will receive the data from the
     * domain. Must always be nonnull but the pointers it contains to entity
     * data may be null of no entities are on-process.
     *
     * \param parameters Parameters for the setup.
     */
    void setupImpl( const Teuchos::RCP<FunctionSpace> &domain_space,
                    const Teuchos::RCP<FunctionSpace> &range_space ) override;

    /*!
     * \brief Apply the operator.
     */
    void applyImpl(
        const TpetraMultiVector &X, TpetraMultiVector &Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        double alpha = Teuchos::ScalarTraits<double>::one(),
        double beta = Teuchos::ScalarTraits<double>::zero() ) const override;

    /*
     * \brief Transpose apply option.
     */
    bool hasTransposeApplyImpl() const override;

  private:
    // Range entity topological dimension. Default is 0 (vertex).
    int d_range_entity_dim;

    // Boolean for keeping the original range data when range entities are not
    // mapped.
    bool d_keep_missed_sol;

    // Search sublist.
    Teuchos::ParameterList d_search_list;

    // The coupling matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO>> d_coupling_matrix;

    // An array of range entity ids that were not mapped during the last call
    // to setup.
    Teuchos::Array<EntityId> d_missed_range_entity_ids;

    // The missed range entity update vector.
    Teuchos::RCP<Tpetra::Vector<Scalar, LO, GO>> d_keep_range_vec;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CONSISTENTINTERPOLATIONOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolationOperator.hpp
//---------------------------------------------------------------------------//
