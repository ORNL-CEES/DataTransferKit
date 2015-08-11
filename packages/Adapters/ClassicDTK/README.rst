Adapters for Classic DTK

This package contains adapters for versions 1.x of DTK. With these adapters
previous users who have implemented DTK versions 1.x for use with their codes
can use those same interfaces while the implementations under the data
transfers have been replaced with version 2.x and higher code. These adapters
then effectively give version 1.x functionality using version 2.x code. This
functionality is verified by maintaining the version 1.x unit tests in the
test directory of this package.

It is highly recommended that any new users of DTK do not use the interfaces
in this directory but instead use the new interfaces and implementations.

It is also highly recommended that users of DTK 1.x switch to the new version
of the code. See the implementations in SharedDomainMap and VolumeSourceMap
for how the new code maps to the old.
