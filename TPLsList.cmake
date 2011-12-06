#
# Extra add-on TPLs
#
# For a description of the fields, see:
#
#   Trilinos/cmake/TrilinosTPLs.cmake
#

SET(Coupler_TPLS_FINDMODS_CLASSIFICATIONS
  nemesis   TPLs_cmake/        EX
  Trilinos  TPLs_cmake/        EX
  )

# NOTE: Above, the paths to the FindTPL<TPLNAME> modules (with an implicit
# *.cmake extension) are relative to the Trilinos/cmake/TPLs directory.
