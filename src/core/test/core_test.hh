//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/core_test.hh
 * \author Stuart R. Slattery
 * \date   Tue May 24 16:18:08 2011
 * \brief  Testing harness for core.
 */
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

bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
              << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool pass_msg(const std::string &passmsg)
{
    std::cout << "Test: passed" << std::endl;
    std::cout << "     " << passmsg << std::endl;
    return true;
}

//---------------------------------------------------------------------------//

bool fail_msg(const std::string &failmsg)
{
    std::cout << "Test: failed" << std::endl;
    std::cout << "     " << failmsg << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

void unit_test(const bool pass, int line, char *file)
{
    if ( ! pass )
    {
        fail(line, file);
    }
}

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

bool passed = true;

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
