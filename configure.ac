dnl-------------------------------------------------------------------------dnl
dnl denovo configure.ac
dnl $Id: configure.ac,v 1.6 2009/07/28 21:33:34 9te Exp $
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for denovo
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([README])
AC_CONFIG_AUX_DIR(config)

# nemesis version
AC_NEMESIS_PREREQ

AC_ALL_VENDORS_SETUP
AC_NEMESIS_ARGS
AC_DENOVO_ARGS

# other arguments
AC_ARG_WITH(timing-diag,
    [  --with-timing-diag      set nemesis timing level; 0-3 (default: 1)])

# The order of packages in AC_CONFIG_SUBDIRS **is** important.
# The packages may depend on other packages to the left of them
# in the list.

AC_CONFIG_SUBDIRS(src)

# make the head makefile
AC_HEAD_MAKEFILE

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of denovo configure.ac
dnl-------------------------------------------------------------------------dnl
