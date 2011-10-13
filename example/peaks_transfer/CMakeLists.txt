# CMakeLists.txt for dtransfer/example/peaks_transfer

# Headers and sources
SET(HEADERS
  Post_Proc.hh
  Physics_A.hh
  Physics_B.hh
  TE_Physics_A.hh
  TE_Physics_B.hh
  )

SET(SOURCES
  Post_Proc.cc
  Physics_A.cc
  Physics_B.cc
  TE_Physics_A.cc
  TE_Physics_B.cc
  driver.cc
  )

# Dependency list
SET(LIBS "")

# MOAB Dependency
INCLUDE_DIRECTORIES( ${MOAB_DIR}/include )
LINK_DIRECTORIES( ${MOAB_DIR}/lib )
SET(LIBS ${LIBS} MOAB core)


# Example Executable
ADD_EXECUTABLE( peaks_transfer ${SOURCES} )

# Link
TARGET_LINK_LIBRARIES( peaks_transfer ${LIBS} )