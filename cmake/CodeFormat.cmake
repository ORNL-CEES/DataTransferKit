if(NOT CLANG_FORMAT_EXECUTABLE)
    find_program(CLANG_FORMAT_EXECUTABLE
        NAMES
        clang-format-4.0
        clang-format-mp-4.0
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

# Check that the vesion of clang-format is 4.0
execute_process(
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -version
    OUTPUT_VARIABLE CLANG_FORMAT_VERSION
)
if(NOT CLANG_FORMAT_VERSION MATCHES "4.0")
    message(FATAL_ERROR "You must use clang-format version 4.0")
endif()
# Download diff-clang-format.py from ORNL-CEES/Cap
file(DOWNLOAD
    https://raw.githubusercontent.com/ORNL-CEES/Cap/master/diff-clang-format.py
    ${CMAKE_BINARY_DIR}/diff-clang-format.py
    STATUS status
)
list(GET status 0 error_code)
if(error_code)
    list(GET status 1 error_string)
    message(WARNING "Failed downloading diff-clang-format.py from GitHub"
            " (${error_string})")
    message("-- " "NOTE: Disabling C++ code formatting because "
            "diff-clang-format-cpp.py is missing")
    set(skip TRUE)
endif()

# Do not bother continuing if not able to fetch diff-clang-format.py
if(NOT skip)

# Download docopt command line argument parser
file(DOWNLOAD
    https://raw.githubusercontent.com/docopt/docopt/0.6.2/docopt.py
    ${CMAKE_BINARY_DIR}/docopt.py
)
# Add a custom target that applies the C++ code formatting style to the source
add_custom_target(format-cpp
    ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/diff-clang-format.py
        --file-extension='.hpp'
        --file-extension='.cpp'
        --binary=${CLANG_FORMAT_EXECUTABLE}
        --style=file
        --config=${${PACKAGE_NAME}_SOURCE_DIR}/.clang-format
        --apply-patch
        ${${PACKAGE_NAME}_SOURCE_DIR}/packages
)
# Add a test that checks the code is formatted properly
file(WRITE
    ${${PACKAGE_NAME}_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_format_cpp.sh
    "#!/usr/bin/env bash\n"
    "\n"
    "${PYTHON_EXECUTABLE} "
    "${CMAKE_BINARY_DIR}/diff-clang-format.py "
    "--file-extension='.hpp' --file-extension='.cpp' "
    "--binary=${CLANG_FORMAT_EXECUTABLE} "
    "--style=file "
    "--config=${${PACKAGE_NAME}_SOURCE_DIR}/.clang-format "
    "${${PACKAGE_NAME}_SOURCE_DIR}/packages"
)
file(COPY
    ${${PACKAGE_NAME}_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_format_cpp.sh
    DESTINATION
        ${${PACKAGE_NAME}_BINARY_DIR}
    FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)
add_test(
    NAME check_format_cpp
    COMMAND ${${PACKAGE_NAME}_BINARY_DIR}/check_format_cpp.sh
)

endif() # skip when download fails
