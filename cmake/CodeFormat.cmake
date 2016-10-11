find_program(CLANG_FORMAT_EXECUTABLE
    NAMES
    clang-format-3.9
    clang-format-mp-3.9
)
if(CLANG_FORMAT_EXECUTABLE)
    message("-- Found clang-format: ${CLANG_FORMAT_EXECUTABLE}")
else()
    message(FATAL_ERROR "-- clang-format not found")
endif()
file(DOWNLOAD
    https://raw.githubusercontent.com/ORNL-CEES/Cap/master/diff-clang-format.py
    ${CMAKE_BINARY_DIR}/diff-clang-format.py
)
file(DOWNLOAD
    https://raw.githubusercontent.com/docopt/docopt/0.6.2/docopt.py
    ${CMAKE_BINARY_DIR}/docopt.py
)
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
file(WRITE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_format_cpp.sh
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
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_format_cpp.sh
    DESTINATION
        ${CMAKE_BINARY_DIR}
    FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)
add_test(
    NAME check_format_cpp
    COMMAND ${CMAKE_BINARY_DIR}/check_format_cpp.sh
)
