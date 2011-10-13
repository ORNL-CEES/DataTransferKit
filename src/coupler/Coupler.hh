//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Coupler.hh
 * \author Stuart R. Slattery
 * \date   Fri Jun 10 08:56:27 2011
 * \brief  Coupler manager class definition
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Coupler_hh
#define coupler_Coupler_hh

#include "Point_Map.hh"
#include "Messenger.hh"
#include "LG_Indexer.hh"
#include "mesh_type/Point.hh"
#include "neutronics/Neutronics.hh"
#include "fields/Field_View.hh"
#include "harness/DBC.hh"
#include "comm/global.hh"
#include "utils/SP.hh"

#include <string>
#include <vector>

namespace coupler
{

//===========================================================================//
/*!
 * \class Coupler
 * \brief Manages solution coupling for the neutronics package.
 *
 * The Coupler class has two template parameters:
 * \param Master_App  This is the application that determines the points where
 *                    data is required.  A typical example would be a CFD
 *                    code.
 * \param Slave_App   This is the application that returns data on a set of
 *                    points.  A typical example would be a neutronics
 *                    application. 
 * A typical use of the coupler class would be to transfer power densities
 * from a neutronics code to a CFD code.  The CFD code is the Master
 * Application, and determines where the powers should be evaluated.  The
 * neutronics code is the Slave_App that provides the data requested by the
 * CFD code at the provided points.
 *
 * \sa Coupler.cc for detailed descriptions.
 */
/*! 
 * \example coupler/test/tstCoupler.cc
 *
 * Test of Coupler.
 */
//===========================================================================//

template<class Master_App, class Slave_App>
class Coupler 
{
  public:
    //@{
    //! Useful typedefs.
    typedef int                                     OrdinateType;
    typedef double                                  DataType;
    typedef nemesis::Communicator_t                 Communicator_t;
    typedef Master_App                              Master_Application;
    typedef Slave_App                               Slave_Application;
    typedef denovo::SP<Master_Application>          SP_Master_Application;
    typedef denovo::SP<Slave_Application>           SP_Slave_Application;
    typedef LG_Indexer<OrdinateType>                LG_Indexer_t;
    typedef denovo::SP<LG_Indexer_t>                SP_LG_Indexer;
    typedef Point_Map<OrdinateType, DataType>       Point_Map_t;
    typedef denovo::SP<Point_Map_t>                 SP_Point_Map;
    typedef Messenger<DataType>                     Messenger_t;
    typedef denovo::SP<Messenger_t>                 SP_Messenger;
    typedef denovo::Point<OrdinateType, DataType>   Point_t;
    typedef std::vector<Point_t>                    Vec_Point;
    typedef std::vector<DataType>                   Vec_DataType;
    typedef std::vector<OrdinateType>               Vec_OrdinateType;
    typedef std::vector<char>                       Buffer;
    typedef fields::Field_View                      View_Field;
    //@}

  private:
    // Coupler communicator (union of neutronics and external app)
    Communicator_t d_world_comm;

    // Master application communicator
    Communicator_t d_master_comm;

    // Slave application communicator
    Communicator_t d_slave_comm;

    // Block-constant neutronics communicator.
    Communicator_t d_neutronics_block_comm;

    // Master application object
    SP_Master_Application d_master_app;

    // Slave application object.
    SP_Neutronics d_neutronics;

    // Master app local to global indexer.
    SP_LG_Indexer d_indexer_master;

    // Slave app local to global indexer.
    SP_LG_Indexer d_indexer_slave;
    
    // Neutronics messenger object.
    SP_Messenger d_messenger_neut;

    // External messenger object.
    SP_Messenger d_messenger_ext;

    // Neutronics local map object.
    SP_Point_Map d_neutronics_map;

    // External local map object
    SP_Point_Map d_external_map;

    // Local point vector.
    Vec_Point d_points;

    // Local volume vector.
    Vec_Dbl d_volumes;

  public:
    // Constructor.
    Coupler(Communicator_t comm_world, Communicator_t comm_neutronics, 
            Communicator_t comm_external, denovo::SP<T> neutronics_app_ptr,
            denovo::SP<X> external_app_ptr);

    // Register a neutronics object with the Coupler.
    void register_neutronics(SP_Neutronics neutronics);

    // Register a vector of interleaved point coordinates with the Coupler
    // and their associated handles.
    void register_points(const Vec_Dbl &points, const Vec_Int &handles);

    // Get the vector of points registered with the Coupler.
    const Vec_Point get_points() { return d_points; }

    // Build a local map on all processes in the communicator union.
    void build_map();
    void build_map(const std::string &STAR_geom_file);

    // Transfer power data from neutronics to an external application.
    void transfer_power();
    void transfer_power(const std::string &STAR_power_file);

    // Transfer temperature from CFD to neutronics.
    void transfer_temperature();
    void transfer_temperature(const std::string& STAR_temp_file);

    // Get the vector of local powers registered with the Coupler.
    void get_power(Vec_Dbl &powers);

  private:
    // Handle comparator.
    static bool handle_compare(Point_Dbl p1, Point_Dbl p2);
};

//---------------------------------------------------------------------------//
// DUMMY STAR_CCM OBJECT
//---------------------------------------------------------------------------//
class STAR_CCM {    };

//---------------------------------------------------------------------------//
// COUPLER TEMPLATE SPECIALIZATION ON STAR-CCM+/DENOVO
//---------------------------------------------------------------------------//
template<>
class Coupler<STAR_CCM, neutronics::Neutronics>
{
  public:
    //@{
    //! Useful typedefs.
    typedef int                                     OrdinateType;
    typedef double                                  DataType;
    typedef nemesis::Communicator_t                 Communicator_t;
    typedef Master_App                              Master_Application;
    typedef Slave_App                               Slave_Application;
    typedef denovo::SP<Master_Application>          SP_Master_Application;
    typedef denovo::SP<Slave_Application>           SP_Slave_Application;
    typedef LG_Indexer<OrdinateType>                LG_Indexer_t;
    typedef denovo::SP<LG_Indexer_t>                SP_LG_Indexer;
    typedef Point_Map<OrdinateType, DataType>       Point_Map_t;
    typedef denovo::SP<Point_Map_t>                 SP_Point_Map;
    typedef Messenger<DataType>                     Messenger_t;
    typedef denovo::SP<Messenger_t>                 SP_Messenger;
    typedef denovo::Point<OrdinateType, DataType>   Point_t;
    typedef std::vector<Point_t>                    Vec_Point;
    typedef std::vector<DataType>                   Vec_DataType;
    typedef std::vector<OrdinateType>               Vec_OrdinateType;
    typedef std::vector<char>                       Buffer;
    typedef fields::Field_View                      View_Field;
    //@}

  private:
    // Coupler communicator (union of neutronics and external app)
    Communicator_t d_world_comm;

    // Master application communicator
    Communicator_t d_master_comm;

    // Slave application communicator
    Communicator_t d_slave_comm;

    // Master application object
    SP_Master_Application d_master_app;

    // Slave application object.
    SP_Neutronics d_neutronics;

    // Master app local to global indexer.
    SP_LG_Indexer d_indexer_master;

    // Slave app local to global indexer.
    SP_LG_Indexer d_indexer_slave;

    // Master local map object
    SP_Point_Map d_master_map;

    // Slave local map object.
    SP_Point_Map d_slave_map;
    
    // Master messenger object.
    SP_Messenger d_messenger_master;

    // Slave messenger object.
    SP_Messenger d_messenger_slave;

    // Local point vector.
    Vec_Point d_points;

    // Local volume vector.
    Vec_Dbl d_volumes;

  public:
    // Constructor.
    Coupler(Communicator_t comm_world, Communicator_t comm_neutronics, 
            Communicator_t comm_external, denovo::SP<T> neutronics_app_ptr,
            denovo::SP<X> external_app_ptr);

    // Get the vector of points registered with the Coupler.
    const Vec_Point get_points() { return d_points; }

    // Build a local map on all processes in the communicator union.
    void build_map(const std::string &STAR_geom_file);

    // Transfer power data from neutronics to an external application.
    void transfer_power(const std::string &STAR_power_file);

    // Transfer temperature from CFD to neutronics.
    void transfer_temperature(const std::string& STAR_temp_file);

    // Get the vector of local powers registered with the Coupler.
    void get_power(Vec_DataType &powers);

  private:
    // Handle comparator.
    static bool handle_compare(Point_Dbl p1, Point_Dbl p2);
};


} // end namespace coupler

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "Coupler.i.hh"

#endif // coupler_Coupler_hh

//---------------------------------------------------------------------------//
//              end of coupler/Coupler.hh
//---------------------------------------------------------------------------//
