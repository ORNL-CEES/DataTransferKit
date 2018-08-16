if(NOT CLANG_FORMAT_EXECUTABLE)
    find_program(CLANG_FORMAT_EXECUTABLE
        NAMES
        clang-format-6.0
        clang-format-mp-6.0
        clang-format
    )
    if(CLANG_FORMAT_EXECUTABLE)
        message("-- Found clang-format: ${CLANG_FORMAT_EXECUTABLE}")
    else()
        message(FATAL_ERROR "-- clang-format not found")
    endif()
else()
    message("-- Using clang-format: ${CLANG_FORMAT_EXECUTABLE}")
    if(NOT EXISTS ${CLANG_FORMAT_EXECUTABLE})
        message(FATAL_ERROR "-- clang-format path is invalid")
    endif()
endif()

# Check that the version of clang-format is the correct one
execute_process(
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -version
    OUTPUT_VARIABLE CLANG_FORMAT_VERSION
)
if(NOT CLANG_FORMAT_VERSION MATCHES "6.0")
    message(FATAL_ERROR "You must use clang-format version 6.0")
endif()

# Add a custom target that applies the C++ code formatting style to the source
add_custom_target(format-cpp
    COMMAND
        CLANG_FORMAT_EXE=${CLANG_FORMAT_EXECUTABLE}
        ${${PACKAGE_NAME}_SOURCE_DIR}/scripts/check_format_cpp.sh --apply-patch
    WORKING_DIRECTORY ${${PACKAGE_NAME}_SOURCE_DIR}
)

# Add a test that checks the code is formatted properly
add_test(
    NAME check_format_cpp
    COMMAND ${${PACKAGE_NAME}_SOURCE_DIR}/scripts/check_format_cpp.sh
    WORKING_DIRECTORY ${${PACKAGE_NAME}_SOURCE_DIR}
)
set_tests_properties( check_format_cpp PROPERTIES
    ENVIRONMENT CLANG_FORMAT_EXE=${CLANG_FORMAT_EXECUTABLE}
)
