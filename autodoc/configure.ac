dnl-------------------------------------------------------------------------dnl
dnl coupler/autodoc configure.in
dnl Stuart Slattery
dnl Fri Sep  3 13:04:05 1999
dnl $Id: configure.ac,v 1.20 2009/10/09 20:35:20 9te Exp $
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for autodoc
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([mainpage.dcc.in])
AC_CONFIG_AUX_DIR(../config)

# clubimc version
AC_NEMESIS_PREREQ
package="coupler"
AC_COMPONENT_BANNER
AC_NEMESIS_ENV

dnl
dnl doxygen setup
dnl

# check for and as sign the path to doxygen
AC_PATH_PROG(DOXYGEN_PATH, doxygen, null)   

#
# Run the generic package-level autodoc configuration:
#
AC_PACKAGE_AUTODOC

#
# Do clubimc-specific configuration:
#

# dotpath (for coupler levelization diagram)
dotpath="${srcdir}/../doc/design"
AC_SUBST(dotpath)

# find the release number
number='null'
AC_MSG_CHECKING("coupler release number")
number=coupler-3_0_0
AC_MSG_RESULT($number)
AC_SUBST(number)

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of coupler/autodoc configure.in
dnl-------------------------------------------------------------------------dnl
