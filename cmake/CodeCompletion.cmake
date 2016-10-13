set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

configure_file(
    ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/.ycm_extra_conf.py.in
    ${${PACKAGE_NAME}_SOURCE_DIR}/.ycm_extra_conf.py
    @ONLY
)
