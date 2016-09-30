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
 * \file   DTK_NonlinearProblemTraits.hpp
 * \author Stuart Slattery
 * \brief  Traits/policy class for nonlinear problems.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_NONLINEARPROBLEMTRAITS_HPP
#define DTK_NONLINEARPROBLEMTRAITS_HPP

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename UndefinedNonlinearProblem>
struct UndefinedNonlinearProblemTraits
{
    static inline UndefinedNonlinearProblem notDefined()
    {
        return UndefinedNonlinearProblem::
            this_type_is_missing_a_specialization();
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class NonlinearProblemTraits
 * \brief Traits class for nonlinear problems.
 */
//---------------------------------------------------------------------------//
template<typename NonlinearProblem>
class NonlinearProblemTraits
{
  public:

    //! Typedef for NonlinearProblem.
    typedef NonlinearProblem nonlinear_problem_type;

    //! Typedef for multidimensional array type.
    typedef typename NonlinearProblem::md_array_type MDArray;

    //! Typedef for multidimensional array scalar type.
    typedef typename MDArray::scalar_type Scalar;

    //! Update the state of the nonlinear problem given a new solution vector.
    static inline void updateState(
        NonlinearProblem& problem, const MDArray& u )
    {
        UndefinedNonlinearProblemTraits<NonlinearProblem>::notDefined();
    }

    //! Compute the nonlinear residual given a new solution vector. F must be
    //! allocated.
    static inline void evaluateResidual(
        const NonlinearProblem& problem, const MDArray& u, MDArray& F )
    {
        UndefinedNonlinearProblemTraits<NonlinearProblem>::notDefined();
    }

    //! Compute the Jacobian matrix given a new solution vector. J must be
    //! allocated.
    static inline void evaluateJacobian(
        const NonlinearProblem& problem, const MDArray& u, MDArray& J )
    {
        UndefinedNonlinearProblemTraits<NonlinearProblem>::notDefined();
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_NONLINEARPROBLEMTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_NonlinearProblemTraits.hpp
//---------------------------------------------------------------------------//

