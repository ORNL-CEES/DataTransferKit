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
 * \brief DTK_DerivedObjectRegistry_impl.hpp
 * \author Stuart R. Slattery
 * \brief Derived class registration with base class.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DERIVEDOBJECTREGISTRY_IMPL_HPP
#define DTK_DERIVEDOBJECTREGISTRY_IMPL_HPP

#include <Teuchos_CompileTimeAssert.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template<class Base,
	 class Derived1,class Derived2,class Derived3,
	 class Derived4,class Derived5,class Derived6,
	 class Derived7,class Derived8,class Derived9> 
void
DerivedObjectRegistry<Base,
		      Derived1,Derived2,Derived3,
		      Derived4,Derived5,Derived6,
		      Derived7,Derived8,Derived9>::registerDerivedClasses()
{
    DerivedObjectRegistrationPolicy<Derived1>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived2>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived3>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived4>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived5>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived6>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived7>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived8>::registerDerivedClassWithBaseClass();
    DerivedObjectRegistrationPolicy<Derived9>::registerDerivedClassWithBaseClass();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_DERIVEDOBJECTREGISTRY_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_DerivedObjectRegistry_impl.hpp
//---------------------------------------------------------------------------//
