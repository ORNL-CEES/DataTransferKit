FILE(WRITE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_no_trailing_whitespaces_or_tabs.sh
    "#!/usr/bin/env bash\n"
    "\n"
    "grep "
    "--exclude={*.doc,*.eps,*.gif,*.jpg,*.pdf,*.png} "
    "--exclude={.project,.cproject} "
    "--exclude={*Makefile*,*makefile*} "
    "--exclude=.gitmodules "
    "--exclude-dir={.git,data} "
    "--regexp '[[:blank:]]$' "
    "--regexp $'\\t' "
    "--line-number "
    "`git ls-files`"
    "\n"
    "test \$? -eq 1"
)
FILE(COPY
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_no_trailing_whitespaces_or_tabs.sh
    DESTINATION
        ${CMAKE_BINARY_DIR}
    FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)
ADD_TEST(
    NAME check_no_trailing_whitespaces_or_tabs
    COMMAND ${CMAKE_BINARY_DIR}/check_no_trailing_whitespaces_or_tabs.sh
    WORKING_DIRECTORY ${${PACKAGE_NAME}_SOURCE_DIR}
)
