set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Enable/Disable output of compile commands during generation." FORCE)

configure_file(
    ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/.ycm_extra_conf.py.in
    ${${PACKAGE_NAME}_SOURCE_DIR}/.ycm_extra_conf.py
    @ONLY
)
