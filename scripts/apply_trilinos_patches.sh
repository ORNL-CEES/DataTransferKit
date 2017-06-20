#!/usr/bin/env bash

hashes=('545269e1838282d98c4a113e17e6322d9c87fe8b')

cd ${TRILINOS_DIR}
for i in "${hashes[@]}"
do
  wget https://github.com/trilinos/trilinos/commit/${i}.patch
  git apply --whitespace=nowarn ${i}.patch
  rm ${i}.patch
done
