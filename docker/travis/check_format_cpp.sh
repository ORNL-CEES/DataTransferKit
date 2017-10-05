#!/usr/bin/env bash

unformatted_files=0
for file in $(git ls-files packages | grep "\(.hpp\|.cpp\)"); do
    diff -u \
        <(cat $file) \
        --label a/$file \
        <(clang-format-5.0 $file) \
        --label b/$file
    if [ $? -eq 1 ]; then
        let "unformatted_files++"
    fi
done
exit $unformatted_files
