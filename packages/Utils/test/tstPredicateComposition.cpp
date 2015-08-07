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
 * \file tstPredicateComposition.cpp
 * \author Stuart R. Slattery
 * \brief Abstract iterator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <functional>

#include <DTK_PredicateComposition.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

//---------------------------------------------------------------------------//
// Helper predicates.
//---------------------------------------------------------------------------//
std::function<bool(int&)> even_func = [](int& n){ return ((n%2) == 0); };
std::function<bool(int&)> odd_func = [](int& n){ return ((n%2) == 1); };
std::function<bool(int&)> div3_func = [](int& n){ return ((n%3) == 0); };
std::function<bool(int&)> div4_func = [](int& n){ return ((n%4) == 0); };
std::function<bool(int&)> div5_func = [](int& n){ return ((n%5) == 0); };

//---------------------------------------------------------------------------//
// Numbers.
//---------------------------------------------------------------------------//
int zero = 0;
int one = 1;
int two = 2;
int three = 3;
int four = 4;
int five = 5;
int seven = 7;
int nine = 9;
int twelve = 12;
int twentyeight = 28;
int thirty = 30;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Predefined predicate tests.
TEUCHOS_UNIT_TEST( PredicateComposition, predefined_predicate_test )
{
    using namespace DataTransferKit;

    // Check the even predicate.
    TEST_ASSERT( even_func(zero) );
    TEST_ASSERT( !even_func(one) );
    TEST_ASSERT( even_func(two) );
    TEST_ASSERT( !even_func(three) );
    TEST_ASSERT( even_func(four) );
    TEST_ASSERT( !even_func(five) );
    TEST_ASSERT( !even_func(seven) );
    TEST_ASSERT( !even_func(nine) );
    TEST_ASSERT( even_func(twelve) );
    TEST_ASSERT( even_func(twentyeight) );
    TEST_ASSERT( even_func(thirty) );

    // Check the odd predicate.
    TEST_ASSERT( !odd_func(zero) );
    TEST_ASSERT( odd_func(one) );
    TEST_ASSERT( !odd_func(two) );
    TEST_ASSERT( odd_func(three) );
    TEST_ASSERT( !odd_func(four) );
    TEST_ASSERT( odd_func(five) );
    TEST_ASSERT( odd_func(seven) );
    TEST_ASSERT( odd_func(nine) );
    TEST_ASSERT( !odd_func(twelve) );
    TEST_ASSERT( !odd_func(twentyeight) );
    TEST_ASSERT( !odd_func(thirty) );

    // Check the div3 predicate.
    TEST_ASSERT( div3_func(zero) );
    TEST_ASSERT( !div3_func(one) );
    TEST_ASSERT( !div3_func(two) );
    TEST_ASSERT( div3_func(three) );
    TEST_ASSERT( !div3_func(four) );
    TEST_ASSERT( !div3_func(five) );
    TEST_ASSERT( !div3_func(seven) );
    TEST_ASSERT( div3_func(nine) );
    TEST_ASSERT( div3_func(twelve) );
    TEST_ASSERT( !div3_func(twentyeight) );
    TEST_ASSERT( div3_func(thirty) );

    // Check the div4 predicate.
    TEST_ASSERT( div4_func(zero) );
    TEST_ASSERT( !div4_func(one) );
    TEST_ASSERT( !div4_func(two) );
    TEST_ASSERT( !div4_func(three) );
    TEST_ASSERT( div4_func(four) );
    TEST_ASSERT( !div4_func(five) );
    TEST_ASSERT( !div4_func(seven) );
    TEST_ASSERT( !div4_func(nine) );
    TEST_ASSERT( div4_func(twelve) );
    TEST_ASSERT( div4_func(twentyeight) );
    TEST_ASSERT( !div4_func(thirty) );

    // Check the div5 predicate.
    TEST_ASSERT( div5_func(zero) );
    TEST_ASSERT( !div5_func(one) );
    TEST_ASSERT( !div5_func(two) );
    TEST_ASSERT( !div5_func(three) );
    TEST_ASSERT( !div5_func(four) );
    TEST_ASSERT( div5_func(five) );
    TEST_ASSERT( !div5_func(seven) );
    TEST_ASSERT( !div5_func(nine) );
    TEST_ASSERT( !div5_func(twelve) );
    TEST_ASSERT( !div5_func(twentyeight) );
    TEST_ASSERT( div5_func(thirty) );
}

//---------------------------------------------------------------------------//
// Predicate composition tests.
TEUCHOS_UNIT_TEST( PredicateComposition, predicate_composition_test )
{
    using namespace DataTransferKit;

    std::function<bool(int&)> even_and_even =
	PredicateComposition::And( even_func, even_func );
    TEST_ASSERT( even_and_even(zero) );
    TEST_ASSERT( !even_and_even(one) );
    TEST_ASSERT( even_and_even(two) );
    TEST_ASSERT( !even_and_even(three) );
    TEST_ASSERT( even_and_even(four) );
    TEST_ASSERT( !even_and_even(five) );
    TEST_ASSERT( !even_and_even(seven) );
    TEST_ASSERT( !even_and_even(nine) );
    TEST_ASSERT( even_and_even(twelve) );
    TEST_ASSERT( even_and_even(twentyeight) );
    TEST_ASSERT( even_and_even(thirty) );

    std::function<bool(int&)> even_or_even =
	PredicateComposition::Or( even_func, even_func );
    TEST_ASSERT( even_or_even(zero) );
    TEST_ASSERT( !even_or_even(one) );
    TEST_ASSERT( even_or_even(two) );
    TEST_ASSERT( !even_or_even(three) );
    TEST_ASSERT( even_or_even(four) );
    TEST_ASSERT( !even_or_even(five) );
    TEST_ASSERT( !even_or_even(seven) );
    TEST_ASSERT( !even_or_even(nine) );
    TEST_ASSERT( even_or_even(twelve) );
    TEST_ASSERT( even_or_even(twentyeight) );
    TEST_ASSERT( even_or_even(thirty) );

    std::function<bool(int&)> even_andnot_even =
	PredicateComposition::AndNot( even_func, even_func );
    TEST_ASSERT( !even_andnot_even(zero) );
    TEST_ASSERT( !even_andnot_even(one) );
    TEST_ASSERT( !even_andnot_even(two) );
    TEST_ASSERT( !even_andnot_even(three) );
    TEST_ASSERT( !even_andnot_even(four) );
    TEST_ASSERT( !even_andnot_even(five) );
    TEST_ASSERT( !even_andnot_even(seven) );
    TEST_ASSERT( !even_andnot_even(nine) );
    TEST_ASSERT( !even_andnot_even(twelve) );
    TEST_ASSERT( !even_andnot_even(twentyeight) );
    TEST_ASSERT( !even_andnot_even(thirty) );

    std::function<bool(int&)> even_and_odd =
	PredicateComposition::And( even_func, odd_func );
    TEST_ASSERT( !even_and_odd(zero) );
    TEST_ASSERT( !even_and_odd(one) );
    TEST_ASSERT( !even_and_odd(two) );
    TEST_ASSERT( !even_and_odd(three) );
    TEST_ASSERT( !even_and_odd(four) );
    TEST_ASSERT( !even_and_odd(five) );
    TEST_ASSERT( !even_and_odd(seven) );
    TEST_ASSERT( !even_and_odd(nine) );
    TEST_ASSERT( !even_and_odd(twelve) );
    TEST_ASSERT( !even_and_odd(twentyeight) );
    TEST_ASSERT( !even_and_odd(thirty) );

    std::function<bool(int&)> even_or_odd =
	PredicateComposition::Or( even_func, odd_func );
    TEST_ASSERT( even_or_odd(zero) );
    TEST_ASSERT( even_or_odd(one) );
    TEST_ASSERT( even_or_odd(two) );
    TEST_ASSERT( even_or_odd(three) );
    TEST_ASSERT( even_or_odd(four) );
    TEST_ASSERT( even_or_odd(five) );
    TEST_ASSERT( even_or_odd(seven) );
    TEST_ASSERT( even_or_odd(nine) );
    TEST_ASSERT( even_or_odd(twelve) );
    TEST_ASSERT( even_or_odd(twentyeight) );
    TEST_ASSERT( even_or_odd(thirty) );

    std::function<bool(int&)> even_andnot_odd =
	PredicateComposition::AndNot( even_func, odd_func );
    TEST_ASSERT( even_andnot_odd(zero) );
    TEST_ASSERT( !even_andnot_odd(one) );
    TEST_ASSERT( even_andnot_odd(two) );
    TEST_ASSERT( !even_andnot_odd(three) );
    TEST_ASSERT( even_andnot_odd(four) );
    TEST_ASSERT( !even_andnot_odd(five) );
    TEST_ASSERT( !even_andnot_odd(seven) );
    TEST_ASSERT( !even_andnot_odd(nine) );
    TEST_ASSERT( even_andnot_odd(twelve) );
    TEST_ASSERT( even_andnot_odd(twentyeight) );
    TEST_ASSERT( even_andnot_odd(thirty) );

    std::function<bool(int&)> div3_and_div4 =
	PredicateComposition::And( div3_func, div4_func );
    TEST_ASSERT( div3_and_div4(twelve) );
    TEST_ASSERT( !div3_and_div4(three) );
    TEST_ASSERT( !div3_and_div4(four) );

    std::function<bool(int&)> div4_and_div3 =
	PredicateComposition::And( div4_func, div3_func );
    TEST_ASSERT( div4_and_div3(twelve) );
    TEST_ASSERT( !div4_and_div3(three) );
    TEST_ASSERT( !div4_and_div3(four) );

    std::function<bool(int&)> div3_or_div4 =
	PredicateComposition::Or( div3_func, div4_func );
    TEST_ASSERT( div3_or_div4(twelve) );
    TEST_ASSERT( div3_or_div4(three) );
    TEST_ASSERT( div3_or_div4(four) );

    std::function<bool(int&)> div4_or_div3 =
	PredicateComposition::Or( div4_func, div3_func );
    TEST_ASSERT( div4_or_div3(twelve) );
    TEST_ASSERT( div4_or_div3(three) );
    TEST_ASSERT( div4_or_div3(four) );

    std::function<bool(int&)> div3_andnot_div4 =
	PredicateComposition::AndNot( div3_func, div4_func );
    TEST_ASSERT( !div3_andnot_div4(twelve) );
    TEST_ASSERT( div3_andnot_div4(three) );
    TEST_ASSERT( !div3_andnot_div4(four) );

    std::function<bool(int&)> div4_andnot_div3 =
	PredicateComposition::AndNot( div4_func, div3_func );
    TEST_ASSERT( !div4_andnot_div3(twelve) );
    TEST_ASSERT( !div4_andnot_div3(three) );
    TEST_ASSERT( div4_andnot_div3(four) );

    std::function<bool(int&)> div3_ornot_div4 =
	PredicateComposition::OrNot( div3_func, div4_func );
    TEST_ASSERT( div3_ornot_div4(twelve) );
    TEST_ASSERT( div3_ornot_div4(three) );
    TEST_ASSERT( !div3_ornot_div4(four) );
    TEST_ASSERT( div3_ornot_div4(five) );

    std::function<bool(int&)> div4_ornot_div3 =
	PredicateComposition::OrNot( div4_func, div3_func );
    TEST_ASSERT( div4_ornot_div3(twelve) );
    TEST_ASSERT( !div4_ornot_div3(three) );
    TEST_ASSERT( div4_ornot_div3(four) );
    TEST_ASSERT( div4_ornot_div3(five) );
}

//---------------------------------------------------------------------------//
// end tstPredicateComposition.cpp
//---------------------------------------------------------------------------//

