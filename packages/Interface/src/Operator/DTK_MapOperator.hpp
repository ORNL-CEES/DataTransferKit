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
 * \brief DTK_MapOperator.hpp
 * \author Stuart R. Slattery
 * \brief Map operator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAPOPERATOR_HPP
#define DTK_MAPOPERATOR_HPP

#include "DTK_FunctionSpace.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MapOperator
  \brief Map operator interface.

  A map operator maps a field in one entity set to another entity set.
*/
//---------------------------------------------------------------------------//
template<class Scalar>
class MapOperator : public Tpetra::Operator<Scalar,int,std::size_t>
{
  public:

    //! Root class typedef.
    typedef Tpetra::Operator<Scalar,int,std::size_t> Root;

    //! Map typedef.
    typedef Tpetra::Map<int,std::size_t,typename Root::node_type> TpetraMap;

    //! MultiVector typedef.
    typedef Tpetra::MultiVector<Scalar,int,std::size_t,typename Root::node_type>
    TpetraMultiVector;

    /*!
     * \brief Constructor.
     */
    MapOperator( const Teuchos::RCP<const TpetraMap>& domain_map,
		 const Teuchos::RCP<const TpetraMap>& range_map );

    /*!
     * \brief Destructor.
     */
    virtual ~MapOperator();

    /*
     * \brief Setup the map operator from a domain entity set and a range
     * entity set.
     * \param domain_map Parallel map for domain vectors this map should be
     * compatible with.
     * \param domain_function The function that contains the data that will be
     * sent to the range. Must always be nonnull but the pointers it contains
     * may be null of no entities are on-process.
     * \param range_map Parallel map for range vectors this map should be
     * compatible with.
     * \param range_space The function that will receive the data from the
     * domain. Must always be nonnull but the pointers it contains
     * may be null of no entities are on-process.
     * \param parameters Parameters for the setup.
     */
    virtual void 
    setup( const Teuchos::RCP<FunctionSpace>& domain_space,
	   const Teuchos::RCP<FunctionSpace>& range_space,
	   const Teuchos::RCP<Teuchos::ParameterList>& parameters );

    //@{
    //! Tpetra::Operator interface.
    Teuchos::RCP<const TpetraMap> getDomainMap() const override;
    Teuchos::RCP<const TpetraMap> getRangeMap() const override;
    virtual void apply( const TpetraMultiVector& X,
			TpetraMultiVector &Y,
			Teuchos::ETransp mode = Teuchos::NO_TRANS,
			Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
			Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override;
    virtual bool hasTransposeApply() const override { return false; }
    //@}

  private:

    //! Domain map.
    Teuchos::RCP<const TpetraMap> d_domain_map;

    //! Range map.
    Teuchos::RCP<const TpetraMap> d_range_map;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_MapOperator_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_MAPOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MapOperator.hpp
//---------------------------------------------------------------------------//
