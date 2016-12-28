#!/bin/bash

cwd=`pwd`
if [[ "$cwd" != *"/scripts" ]]; then
    echo "./setup_hooks.sh must be called from scripts dir"
    exit 1
fi

cp hooks/pre-commit ../.git/hooks/pre-commit

## If there are formatting errors, print the offending lines and fail.
# Find clang-format that is used with git
git_clang_name=""
for name in "git-clang-format-mp-3.9" "git-clang-format-3.9" "git-clang-format" ; do
    which $name > /dev/null 2>&1
    if [[ $? == 0 ]]; then
        git_clang_name="$name"
        break;
    fi
done

# Fix the git-clang-format name in the hook (even if it's empty)
sed -i "" "s/GIT_CLANG_NAME/$git_clang_name/" ../.git/hooks/pre-commit

if [[ $git_clang_name == "" ]]; then
    echo "Did not find git-clang-format"
    exit 0
fi

# We actually do have some version of clang-format
clang_format=""
for name in "clang-format-mp-3.9" "clang-format-3.9" "clang-format"; do
    which $name > /dev/null 2>&1
    if [[ $? == 0 ]]; then
        clang_format="$name"
        break
    fi
done
if [[ $clang_format == "" ]]; then
    # Found git-clang-format but not clang-format, weird.
    echo "Found git-clang-format but not clang-format, weird"
    exit 1
fi

# Check clang-format version (only work with 3.9)
clang_version=`$clang_format --version | awk '{print $3}' | cut -f 1,2 -d .`
if [[ $clang_version != "3.9" ]]; then
    echo "Clang-format version is bad: $clang_version"
    git_clang_name=""
    clang_format=""
fi

git config clangFormat.style file
git config clangFormat.extension .hpp,.cpp
git config clangFormat.binary $clang_format
