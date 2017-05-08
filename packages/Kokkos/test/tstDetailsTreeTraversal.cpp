#include <details/DTK_DetailsTreeTraversal.hpp>

#include <Teuchos_UnitTestHarness.hpp>

TEUCHOS_UNIT_TEST( LinearBVH, sorted_list )
{
    // list is empty at construction
    DataTransferKit::Details::SortedList list( 2 );
    TEST_ASSERT( list.empty() );
    TEST_ASSERT( !list.full() );
    // insert a first element
    list.emplace( 0, 255.0 );
    TEST_ASSERT( !list.empty() );
    TEST_EQUALITY( list.back().second, 255.0 );
    TEST_EQUALITY( list.size(), 1 );
    TEST_ASSERT( !list.full() );
    // insert a second element that is closer
    // first one is at the back
    list.emplace( 0, 1.41 );
    TEST_EQUALITY( list.size(), 2 );
    TEST_ASSERT( list.full() );
    TEST_EQUALITY( list.back().second, 255.0 );
    // insert another element
    list.emplace( 0, 3.14 );
    TEST_EQUALITY( list.back().second, 3.14 );
    TEST_EQUALITY( list.size(), 2 );
    TEST_ASSERT( list.full() );
    // and so on
    list.emplace( 0, -1.0 );
    TEST_EQUALITY( list.back().second, 1.41 );
    TEST_EQUALITY( list.size(), 2 );
    // attempting to insert an element that is further than the one at the back
    TEST_THROW( list.emplace( 0, 3.14 ), std::runtime_error );
    TEST_EQUALITY( list.back().second, 1.41 );
}

TEUCHOS_UNIT_TEST( LinearBVH, priority_queue )
{
    // queue is empty at construction
    DataTransferKit::Details::PriorityQueue queue;
    TEST_ASSERT( queue.empty() );
    // insert element
    queue.emplace( nullptr, 1.41 );
    TEST_ASSERT( !queue.empty() );
    TEST_EQUALITY( queue.top().second, 1.41 );
    // smaller distance stays on top of the priority queue
    queue.emplace( nullptr, 3.14 );
    TEST_EQUALITY( queue.top().second, 1.41 );
    // remove highest priority element
    queue.pop();
    TEST_EQUALITY( queue.top().second, 3.14 );
    // insert element with higher priority and check it shows up on top
    queue.emplace( nullptr, 1.41 );
    TEST_EQUALITY( queue.top().second, 1.41 );
    // empty the queue
    queue.pop();
    queue.pop();
    TEST_ASSERT( queue.empty() );
}
