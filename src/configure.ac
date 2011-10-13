dnl-------------------------------------------------------------------------dnl 
dnl denovo/src configure.in
dnl Thomas M. Evans
dnl Mon Mar 29 14:53:37 1999
dnl $Id: configure.ac,v 1.20 2009/07/23 18:03:28 9te Exp $
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for denovo/src
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([aclocal.m4])
AC_CONFIG_AUX_DIR(../config)

# nemesis version
AC_NEMESIS_PREREQ

AC_ALL_VENDORS_SETUP
AC_NEMESIS_ARGS

AC_DENOVO_ARGS

# The order of packages in AC_CONFIG_SUBDIRS **is** important.
# Each AC_CONFIG_SUBDIRS corresponds to a level: lowest level
# components at the top; highest level at the bottom:

##---------------------------------------------------------------------------##
## NEMESIS PACKAGES
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(harness) # 1
AC_CONFIG_SUBDIRS(comm)    # 2

##---------------------------------------------------------------------------##
## RELEASE PACKAGE
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(release)

##---------------------------------------------------------------------------##
## DENOVO LEVEL 1
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(utils fields)

##---------------------------------------------------------------------------##
## DENOVO LEVEL 2
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(database meshing mesh_type krp_papi)

##---------------------------------------------------------------------------##
## DENOVO LEVEL 3
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(angle material kba_mesh)

if test "${enable_mc}" = yes; then
    AC_CONFIG_SUBDIRS(mc)
fi

##---------------------------------------------------------------------------##
## DENOVO LEVEL 4
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(source kba_equations geometry)

##---------------------------------------------------------------------------##
## DENOVO LEVEL 5
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(kba_uncflux)

if test "${enable_shift}" = yes; then
   AC_CONFIG_SUBDIRS(mc_physics)
fi

if test "${enable_mc}" = yes; then
    AC_CONFIG_SUBDIRS(kba_mc)
fi

##---------------------------------------------------------------------------##
## DENOVO LEVEL 6
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(kba)

if test "${enable_shift}" = yes; then
   AC_CONFIG_SUBDIRS(mc_sources)
fi

##---------------------------------------------------------------------------##
## DENOVO LEVEL 7
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(kba_solvers)

if test "${enable_shift}" = yes; then
   AC_CONFIG_SUBDIRS(mc_transport)
fi

##---------------------------------------------------------------------------##
## DENOVO LEVEL 8
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(kba_utils)

if test "${enable_shift}" = yes; then
   AC_CONFIG_SUBDIRS(mc_utils)
fi

##---------------------------------------------------------------------------##
## PYTHON AND DRIVERS
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(xkba)

if test "${enable_python}" = yes; then
   AC_CONFIG_SUBDIRS(pykba)

   if test "${enable_shift}" = yes; then
      AC_CONFIG_SUBDIRS(pyrtk)
   fi
fi

if test "${with_silo}" != no; then
   AC_CONFIG_SUBDIRS(hpckba)

   if test "${enable_neutronics}" != no; then
      AC_CONFIG_SUBDIRS(neutronics coupler)
   fi
fi


##---------------------------------------------------------------------------##
## Make the src Makefile

AC_SRC_MAKEFILE

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of denovo/src configure.in
dnl-------------------------------------------------------------------------dnl
