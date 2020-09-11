/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <DTK_ConfigDefs.hpp>
#include <DTK_Version.hpp>

#include <Teuchos_UnitTestHarness.hpp>

#include <iostream>
#include <sstream>

TEUCHOS_UNIT_TEST( DataTransferKitRuntimeAPI, return_version )
{
    auto const dtk_version = DataTransferKit::version();
    TEST_ASSERT( !dtk_version.empty() );
    std::cout << "DataTransferKit version " << dtk_version << std::endl;

    auto const dtk_commit_hash = DataTransferKit::gitCommitHash();
    TEST_ASSERT( !dtk_commit_hash.empty() );
    std::cout << "DataTransferKit commit hash " << dtk_commit_hash << std::endl;
}

namespace dummy
{
struct Foo
{
    Foo( std::ostream &os ) { os << DTK_MARK_REGION( "hello world" ); }
};
void bar( std::ostream &os ) { os << DTK_MARK_REGION( "it works" ); }
} // namespace dummy

TEUCHOS_UNIT_TEST( DataTransferKitMacros, mark_parallel_region )
{
    std::stringstream ss;
    dummy::Foo foo( ss );
    TEST_EQUALITY( ss.str(), "DTK_hello world" );

    ss.clear();
    ss.str( "" );
    dummy::bar( ss );
    TEST_EQUALITY( ss.str(), "DTK_it works" );
}
