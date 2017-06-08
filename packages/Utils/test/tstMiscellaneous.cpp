/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include <DTK_Version.hpp>

#include <Teuchos_UnitTestHarness.hpp>

#include <iostream>

TEUCHOS_UNIT_TEST( DataTransferKitRuntimeAPI, return_version )
{
    auto const dtk_version = DataTransferKit::version();
    TEUCHOS_ASSERT( !dtk_version.empty() );
    std::cout << "DataTransferKit version " << dtk_version << std::endl;

    auto const dtk_commit_hash = DataTransferKit::git_commit_hash();
    TEUCHOS_ASSERT( !dtk_commit_hash.empty() );
    std::cout << "DataTransferKit commit hash " << dtk_commit_hash << std::endl;
}
