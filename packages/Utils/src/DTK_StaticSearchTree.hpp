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
 * \file DTK_StaticSearchTree.hpp
 * \author Stuart R. Slattery
 * \brief Spatial searching for point clouds.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STATICSEARCHTREE_HPP
#define DTK_STATICSEARCHTREE_HPP

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <DTK_nanoflann.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Non-templated static search tree base class.
//---------------------------------------------------------------------------//
class StaticSearchTree
{
  public:

    // Default constructor.
    StaticSearchTree()
    { /* ... */ }

    // Destructor.
    virtual ~StaticSearchTree() = default;

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
class PointCloud
{
  public:
    
    //! Default constructor.
    PointCloud()
    { /* ... */ }

    //! Constructor.
    PointCloud( const Teuchos::ArrayView<const double>& points )
	: d_points( points )
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

    // PointCloud points.
    Teuchos::ArrayView<const double> d_points;
};

//---------------------------------------------------------------------------//
/*!
 * \class NanoflannTree
 * \brief Spatial searching for point clouds.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class NanoflannTree : public StaticSearchTree
{
  public:

    //! Tree typedef.
    typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double,PointCloud<DIM> >,
      PointCloud<DIM>,
      DIM,
      unsigned> TreeType;

  public:

    // Default constructor.
    NanoflannTree( const Teuchos::ArrayView<const double>& points,
		   const unsigned max_leaf_size );

    // Perform an n-nearest neighbor search.
    Teuchos::Array<unsigned> nnSearch( 
	const Teuchos::ArrayView<const double>& point,
	const unsigned num_neighbors ) const;

    // Perform a nearest neighbor search within a specified radius.
    Teuchos::Array<unsigned> radiusSearch( 
	const Teuchos::ArrayView<const double>& point, 
	const double radius ) const;

  private:

    // PointCloud.
    PointCloud<DIM> d_cloud;

    // kD-tree.
    Teuchos::RCP<TreeType> d_tree;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_StaticSearchTree_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_STATICSEARCHTREE_HPP

//---------------------------------------------------------------------------//
// end DTK_StaticSearchTree.hpp
//---------------------------------------------------------------------------//

