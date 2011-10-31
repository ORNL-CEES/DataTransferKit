//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field/test/field_test.hh
 * \author Stuart Slattyer
 * \date   Tue May 24 16:18:08 2011
 * \brief  Testing harness for field.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: pkg_Test.hh,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef field_test_hh
#define field_test_hh

#include <iostream>
#include <string>

namespace field_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set field_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (field_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain nemesis tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

void unit_test(const bool pass, int line, char *file);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern bool passed;

} // end namespace field_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS      field_test::fail(__LINE__);
#define FAILURE      field_test::fail(__LINE__, __FILE__);
#define PASSMSG(a)   field_test::pass_msg(a);
#define FAILMSG(a)   field_test::fail_msg(a);
#define UNIT_TEST(x) field_test::unit_test(x, __LINE__, __FILE__)
    
#endif // field_test_hh

//---------------------------------------------------------------------------//
//     end of field/field_test.hh
//---------------------------------------------------------------------------//
