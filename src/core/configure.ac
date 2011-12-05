dnl-------------------------------------------------------------------------dnl
dnl core/configure.ac
dnl Stuart Slattery
dnl Tue May 24 16:18:08 2011
dnl $Id: configure.package.ac,v 1.7 2009/09/09 14:59:15 9te Exp $
dnl-------------------------------------------------------------------------dnl
dnl Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
dnl-------------------------------------------------------------------------dnl
dnl configure.ac for core
dnl-------------------------------------------------------------------------dnl

AC_INIT

dnl As an argument to this macro provide the name of a file that
dnl exists in the local package directory.  
AC_CONFIG_SRCDIR([aclocal.m4])  

AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(core/config.h:config.h.in)
package="core"
AC_COMPONENT_BANNER

dnl Nemesis version
AC_NEMESIS_PREREQ

dnl ----------------------------------------
dnl Vendor setups: 
dnl ----------------------------------------

dnl Uncomment any vendor setups that are required for this package.
dnl For example, if this package has parallel programming directives
dnl you need to uncomment AC_MPI_SETUP.

AC_MPI_SETUP(pkg)
dnl AC_SPRNG_SETUP(pkg)
AC_TRILINOS_SETUP(pkg)
dnl AC_LAPACK_SETUP(pkg)
dnl AC_METIS_SETUP(pkg)
dnl AC_GSL_SETUP(pkg)
dnl AC_SUPERLU_SETUP(pkg)
dnl AC_HDF5_SETUP(pkg)
dnl AC_SILO_SETUP(pkg)

dnl Turn on MC uncollided flux components
dnl AC_DENOVO_MC(mc, kba_mc)

dnl This macro should be uncommented if the current package is a Nemesis
dnl client (i.e.: This package exists in Denovo, AMP, etc.)

dnl AC_NEMESIS_VENDOR_SETUP

dnl ----------------------------------------
dnl Install options:
dnl ----------------------------------------

AC_INSTALL_HEADERS
AC_INSTALL_LIB
dnl AC_INSTALL_EXECUTABLE

dnl ----------------------------------------
dnl Package dependencies:
dnl ----------------------------------------

dnl List local package dependencies here as arguments to
dnl AC_NEEDS_LIBS.  If this is a Nemesis client (i.e.: a part of
dnl Capsaicin, Milagro, etc.), the dependencies on Nemesis libraries
dnl should be listed as arguments to AC_NEEDS_LIBS_NEMESIS.

AC_NEEDS_LIBS(harness \
              comm \
              release \
              utils)

dnl ----------------------------------------
dnl Setup for unit tests:
dnl ----------------------------------------

dnl List local package dependencies for unit tests here as arguments
dnl to AC_NEEDS_LIBS_TEST.  Unit test dependencies are usually
dnl different than package dependencies.  If this is a Nemesis client,
dnl the unit test dependencies on Nemesis libraries should be listed as
dnl arguments to AC_NEEDS_LIBS_TEST_NEMESIS.

dnl AC_NEEDS_LIBS_TEST(release database)
dnl AC_NEEDS_LIBS_TEST_NEMESIS(comm harness)

dnl List the unit tests that should be run and the number of
dnl processors that should be used when running the test. If this is a
dnl parallel package the number of processors should be a an integer.
dnl If this is a scalar package, the number of processors required
dnl should be listed as "scalar".

AC_RUNTESTS(tstLG_Indexer \
            tstTransfer_Map \
	    tstMessage_Buffer \
	    tstInterfaces \
	    tstTransfer_Data_Field \
            , 1 2 4)

dnl ----------------------------------------
dnl Finish Nemesis setup:
dnl ----------------------------------------

AC_NEMESIS_ENV
AC_NEMESIS_AUTODOC(coupler-0_0_1)

dnl Define FORTRAN name-mangling
AC_FORTRAN_SETUP

dnl Create Makefiles from config templates.

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])

dnl Create test/Makefile.options, if it exists.

AC_CHECK_FILE($srcdir/test/Makefile.options.in, \
              AC_CONFIG_FILES([test/Makefile.options:test/Makefile.options.in]))

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of core/configure.ac
dnl-------------------------------------------------------------------------dnl
