//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/core_test.hh
 * \author gqe
 * \date   Tue May 24 16:18:08 2011
 * \brief  Testing harness for core.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: pkg_Test.hh,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_test_hh
#define core_test_hh

#include <iostream>
#include <string>

namespace core_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set core_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (core_test::passed will have its default value of true)

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

} // end namespace core_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS      core_test::fail(__LINE__);
#define FAILURE      core_test::fail(__LINE__, __FILE__);
#define PASSMSG(a)   core_test::pass_msg(a);
#define FAILMSG(a)   core_test::fail_msg(a);
#define UNIT_TEST(x) core_test::unit_test(x, __LINE__, __FILE__)
    
#endif // core_test_hh

//---------------------------------------------------------------------------//
//     end of core/core_test.hh
//---------------------------------------------------------------------------//
