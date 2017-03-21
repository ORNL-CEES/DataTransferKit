SET_DEFAULT(DataTransferKit_REPOSITORY_MASTER_EMAIL_ADDRESS coupler-infrastructure@casl-dev.ornl.gov)
##---------------------------------------------------------------------------##
## Disable Fortran_API if Trilinos_ENABLE_Fortran is FALSE
##---------------------------------------------------------------------------##
IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    IF ("${${PROJECT_NAME}_ENABLE_DataTransferKitFortran_API}" STREQUAL "")
        MESSAGE("-- " "NOTE: Setting ${PROJECT_NAME}_ENABLE_DataTransferKitFortran_API=OFF because"
                "${PROJECT_NAME}_ENABLE_Fortran = '${${PROJECT_NAME}_ENABLE_Fortran}'")
        SET(${PROJECT_NAME}_ENABLE_DataTransferKitFortran_API OFF)
    ENDIF()
ENDIF()
