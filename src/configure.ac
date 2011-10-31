dnl-------------------------------------------------------------------------dnl 
dnl coupler/src configure.in
dnl Thomas M. Evans
dnl Mon Mar 29 14:53:37 1999
dnl $Id: configure.ac,v 1.20 2009/07/23 18:03:28 9te Exp $
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for coupler/src
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([aclocal.m4])
AC_CONFIG_AUX_DIR(../config)

# nemesis version
AC_NEMESIS_PREREQ

AC_ALL_VENDORS_SETUP
AC_NEMESIS_ARGS

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
## COUPLER LEVEL 1
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(utils)

##---------------------------------------------------------------------------##
## COUPLER LEVEL 2
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(core)

##---------------------------------------------------------------------------##
## Make the src Makefile

AC_SRC_MAKEFILE

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of coupler/src configure.in
dnl-------------------------------------------------------------------------dnl
