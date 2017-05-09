#ifndef DTK_BVHQUERY_DEF_HPP
#define DTK_BVHQUERY_DEF_HPP

#include "DTK_ConfigDefs.hpp"
#include "details/DTK_DetailsTreeTraversal.hpp"

namespace DataTransferKit
{

template <typename NO>
int BVHQuery<NO>::query( BVH<NO> const bvh, Details::Nearest const &predicates,
                         Kokkos::View<int *, typename NO::device_type> out )
{
    using Tag = typename Details::Nearest::Tag;
    return Details::query_dispatch( bvh, predicates, out, Tag{} );
}

template <typename NO>
int BVHQuery<NO>::query( BVH<NO> const bvh, Details::Within const &predicates,
                         Kokkos::View<int *, typename NO::device_type> out )
{
    using Tag = typename Details::Within::Tag;
    return Details::query_dispatch( bvh, predicates, out, Tag{} );
}
}

// Explicit instantiation macro
#define DTK_BVHQUERY_INSTANT( NODE ) template class BVHQuery<NODE>;

#endif
