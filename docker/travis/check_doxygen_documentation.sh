#!/usr/bin/env bash

sed docs/doxygen/Doxyfile.in \
  -e "s|@DTK_REVISION@|${DTK_REVISION}|g" \
  -e "s|@CMAKE_CURRENT_SOURCE_DIR@|${PWD}/docs/doxygen|g" \
  > /tmp/Doxyfile

doxygen /tmp/Doxyfile 2>&1 | tee /tmp/doxygen.log

grep -q 'warning:\|error:' /tmp/doxygen.log

test $? -eq 1
