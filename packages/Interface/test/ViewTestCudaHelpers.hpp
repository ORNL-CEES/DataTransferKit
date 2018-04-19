/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file   ViewTestCudaHelpers.hpp
 * \author Stuart Slattery
 * \brief  tstView.cpp cuda test helpers
 */
//---------------------------------------------------------------------------//

#ifndef DTK_VIEWTESTCUDAHELPERS_HPP
#define DTK_VIEWTESTCUDAHELPERS_HPP

#include "DTK_View.hpp"

namespace ViewTestCudaHelpers
{
//---------------------------------------------------------------------------//
/*
 * \brief Native CUDA fill function to test the view class.
 *
 * This function provides an implementation of the view fill with native
 * cuda. Passing the view directly to the device and using the raw pointer in
 * the view are tested.
 */
template <class Scalar>
void fillViewCuda( DataTransferKit::View<Scalar> dtk_view,
                   const std::vector<unsigned> &dims );

} // namespace ViewTestCudaHelpers

#include "ViewTestCudaHelpers_def.hpp"

#endif // end DTK_VIEWTESTCUDAHELPERS_HPP

//---------------------------------------------------------------------------//
// end ViewTestCudaHelpers.hpp
//---------------------------------------------------------------------------//
