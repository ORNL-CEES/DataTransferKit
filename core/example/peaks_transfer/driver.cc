//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/driver.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:59 2011
 * \brief  Driver for super simple mesh based data transfer example.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Physics_A.hh"
#include "Physics_B.hh"
#include "TE_Physics_A.hh"
#include "TE_Physics_B.hh"
#include "../../src/core/Transfer_Evaluator.hh"
#include "../../src/core/Data_Transfer_Manager.hh"

//---------------------------------------------------------------------------//
// Main function driver.
int main()
{
    // Initialize Physics A.
    physics_A::Physics_A* a = new physics_A::Physics_A(-3.0, 3.0,
						       -3.0, 3.0,
						       250, 250);

    // Initialize Physics B.
    physics_B::Physics_B* b = new physics_B::Physics_B(-3.0, 3.0,
						       -3.0, 3.0,
						       75, 75);

    // Initialize Physics A Transfer_Evaluator.
    dtransfer::Transfer_Evaluator* te_a = new dtransfer::TE_Physics_A(a);

    // Initialize Physics B Transfer_Evaluator.
    dtransfer::Transfer_Evaluator* te_b = new dtransfer::TE_Physics_B(b);

    // Initialize the Data_Transfer_Manager.
    dtransfer::Data_Transfer_Manager manager(te_a, te_b);

    // Build a mapping to transfer the state vector from physics A to the
    // physics B source term.
    manager.map_A2B();

    // Physics A solve.
    a->solve();

    // Transfer physics A solution to physics B.
    manager.transfer_A2B();

    // Physics A plot the state vector.
    a->plot_state();

    // Physics B plot the source term.
    b->plot_source();

    return 0;
}

//---------------------------------------------------------------------------//
//                 end of driver.cc
//---------------------------------------------------------------------------//
