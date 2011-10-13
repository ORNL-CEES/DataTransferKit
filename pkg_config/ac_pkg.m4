dnl-------------------------------------------------------------------------dnl
dnl COUPLER CONFIGURE MACROS
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DENOVO_ARGS], [dnl

   #
   # PYTHON DRIVER SETUP
   #
   AC_ARG_ENABLE(python,
      [  --enable-python         set python driver (off by default)])

   # turn on python if requested: the default is off
   AC_MSG_CHECKING("python front end requested")
   
   if test "${enable_python:=no}" != no; then
       enable_python='yes'
   fi

   AC_MSG_RESULT("${enable_python}")

   #
   # MC SETUP (FOR DENOVO-SN RELATED FEATURES)
   #
   AC_ARG_ENABLE(mc,
      [  --enable-mc             set mc on (off by default)])
   
   if test "${enable_mc:=no}" != no; then
       enable_mc='yes'
   fi

   #
   # ENABLE SHIFT HYBRID MC/SN CODE
   #
   AC_ARG_ENABLE([shift],
     [AC_HELP_STRING([--enable-shift],
                     [turn on Shift hybrid MC/SN code])])

   if test "${enable_shift:=no}" != no; then
       enable_shift='yes'
       # if shift is on make sure that mc components are on
       enable_mc='yes' 
   fi

   # turn on mc if requested: the default is off
   AC_MSG_CHECKING("Denovo Sn Monte Carlo packages requested")
   AC_MSG_RESULT("${enable_mc}")

   # turn on mc if requested: the default is off
   AC_MSG_CHECKING("Build Shift hybrid mc/Sn code")
   AC_MSG_RESULT("${enable_shift}")

   #
   # ENABLE NEUTRONICS PACKAGE
   #
   AC_ARG_ENABLE([neutronics],
     [AC_HELP_STRING([--enable-neutronics],
                     [turn on Neutronics package])])

   if test "${enable_neutronics:=no}" != no; then
       enable_neutronics='yes'
   fi

   # turn on neutronics if requested: the default is off
   AC_MSG_CHECKING("Build Neutronics package")
   AC_MSG_RESULT("${enable_neutronics}")

   #
   # KRP_PAPI
   #
   AC_ARG_ENABLE([krp-papi],
     [AC_HELP_STRING([--enable-krp-papi], 
                     [turn on K. Roaches PAPI HW counter implementation])])
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_pkg.m4
dnl-------------------------------------------------------------------------dnl
