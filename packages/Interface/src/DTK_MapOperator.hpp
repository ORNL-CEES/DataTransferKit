//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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

#include "DTK_AbstractBuilder.hpp"
#include "DTK_AbstractBuildableObject.hpp"
#include "DTK_FieldGroup.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_LinearOpBase.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MapOperator
  \brief Map operator interface.

  A map operator maps a field in one entity set to another entity set.
*/
//---------------------------------------------------------------------------//
class MapOperator : public AbstractBuildableObject<MapOperator>
{
  public:

    /*!
     * \brief Constructor.
     */
    MapOperator();

    /*!
     * \brief Destructor.
     */
    virtual ~MapOperator();

    //@{
    //! Identification functions.
    /*!
     * \brief Return a string indicating the derived map operator type.
     * \return A string key for the derived operator type.
     */
    virtual std::string name() const;
    //@}

    //@{
    //! Setup functions.
    /*
     * \brief Setup the map operator from a source entity set and a target
     * entity set.
     * \param source_group The field group that contains the data that will be
     * sent to the target.
     * mapped. 
     * \param target_set The field group that will receive the data from the
     * source.
     * \param parameters Parameters for the setup.
     */
    virtual void setup( const Teuchos::RCP<FieldGroup>& source_group,
			const Teuchos::RCP<FieldGroup>& target_group,
			const Teuchos::ParameterList>& parameters );
    //@}

    //@{
    //! Apply functions.
    /*
     * \brief Apply the map operator to data defined on the entities.
     * \param source_field_name Field defined on the source entities that will
     * be sent to the target.
     * \param target_field_name Field defined on the target entities that will
     * be received from the source.
     */
    virtual void apply( const std::string& source_field_name,
			const std::string& target_field_name ) const;
    //@}

  protected:

    //@{
    //! Vector construction functions.
    /* 
     * \brief Construct the source and target data vectors from the source and
     * target field.
     */ 
    virtual void buildSourceAndTargetVectors( 
	Teuchos::RCP<const Thyra::MultiVectorBase<double> >& source_data,
	Teuchos::RCP<Thyra::MultiVectorBase<double> >& target_data ) const;
    //@}

  protected:

    //! Source field group.
    Teuchos::RCP<FieldGroup> b_source_group;

    //! Target field group.
    Teuchos::RCP<FieldGroup> b_target_group;

    //! Mass matrix inverse.
    Teuchos::RCP<Thyra::LinearOpBase<double> > b_mass_matrix_inv;

    //! Coupling matrix.
    Teuchos::RCP<Thyra::LinearOpBase<double> > b_coupling_matrix;

    //! Forcing vector.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > b_forcing_vector;
};

//---------------------------------------------------------------------------//
// AbstractBuildableObjectPolicy implementation.
//---------------------------------------------------------------------------//
template<>
class AbstractBuildableObjectPolicy<MapOperator>
{
  public:

    typedef MapOperator object_type;

    static std::string objectType( const MapOperator& map_operator )
    {
	return map_operator.name();
    }

    static Teuchos::RCP<DataTransferKit::AbstractBuilder<MapOperator> > 
    getBuilder()
    {
	return MapOperator::getBuilder();
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MAPOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MapOperator.hpp
//---------------------------------------------------------------------------//
