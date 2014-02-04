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
 * \file DTK_CloudSearch.hpp
 * \author Stuart R. Slattery
 * \brief Spatial searching for point clouds.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLOUDSEARCH_HPP
#define DTK_CLOUDSEARCH_HPP

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <nanoflann.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Non-templated search tree base class.
//---------------------------------------------------------------------------//
class SearchTree
{
  public:

    // Default constructor.
    SearchTree()
    { /* ... */ }

    // Destructor.
    virtual ~SearchTree()
    { /* ... */ }

    // Perform an n-nearest neighbor search.
    virtual Teuchos::Array<unsigned> nnSearch( 
	const Teuchos::ArrayView<const double>& point,
	const unsigned num_neighbors ) const = 0;

    // Perform a nearest neighbor search within a specified radius.
    virtual Teuchos::Array<unsigned> radiusSearch( 
	const Teuchos::ArrayView<const double>& point, 
	const double radius ) const = 0;
};

//---------------------------------------------------------------------------//
// Point cloud structure.
//---------------------------------------------------------------------------//
template<int DIM>
class Cloud
{
  public:
    
    //! Default constructor.
    Cloud()
    { /* ... */ }

    //! Constructor.
    Cloud( const Teuchos::ArrayView<const double>& points )
	: d_points( points )
    { /* ... */ }

    //! Destructor.
    ~Cloud()
    { /* ... */ }

    //! Number of cloud points.
    inline std::size_t kdtree_get_point_count() const
    { return d_points.size() / DIM; }

    // Distance between points.
    double kdtree_distance( 
	const double *p1, const std::size_t idx_p2, std::size_t size) const;

    // Get the point coordinate at the given dimension.
    double kdtree_get_pt( const std::size_t idx, int dim ) const;

    //! Default bounding box calculation.
    template <class BBOX>
    bool kdtree_get_bbox( BBOX& bb ) const 
    { return false; }

  private:

    // Cloud points.
    Teuchos::ArrayView<const double> d_points;
};

//---------------------------------------------------------------------------//
/*!
 * \class CloudSearch
 * \brief Spatial searching for point clouds.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class CloudSearch : public SearchTree
{
  public:

    //! Tree typedef.
    typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double,Cloud<DIM> >,
      Cloud<DIM>,
      DIM,
      unsigned> TreeType;

  public:

    // Default constructor.
    CloudSearch( const Teuchos::ArrayView<const double>& cloud_centers,
		 const unsigned max_leaf_size );

    // Destructor.
    ~CloudSearch()
    { /* ... */ }

    // Perform an n-nearest neighbor search.
    Teuchos::Array<unsigned> nnSearch( 
	const Teuchos::ArrayView<const double>& point,
	const unsigned num_neighbors ) const;

    // Perform a nearest neighbor search within a specified radius.
    Teuchos::Array<unsigned> radiusSearch( 
	const Teuchos::ArrayView<const double>& point, 
	const double radius ) const;

 private:

    // Cloud.
    Cloud<DIM> d_cloud;

    // kD-tree.
    Teuchos::RCP<TreeType> d_tree;
};

//---------------------------------------------------------------------------//
// Base class creation method.
//---------------------------------------------------------------------------//
Teuchos::RCP<SearchTree> createSearchTree( 
    const unsigned dim,
    const Teuchos::ArrayView<const double>& cloud_centers,
    const unsigned leaf_size )
{
    Teuchos::RCP<SearchTree> tree;

    switch ( dim )
    {
	case 1:
	{
	    tree = Teuchos::rcp( 
		new CloudSearch<1>(cloud_centers, leaf_size) );
	}
	break;

	case 2:
	{
	    tree = Teuchos::rcp( 
		new CloudSearch<2>(cloud_centers, leaf_size) );
	}
	break;

	case 3:
	{
	    tree = Teuchos::rcp( 
		new CloudSearch<3>(cloud_centers, leaf_size) );
	}
	break;
    };

    return tree;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_CloudSearch_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_CLOUDSEARCH_HPP

//---------------------------------------------------------------------------//
// end DTK_CloudSearch.hpp
//---------------------------------------------------------------------------//

