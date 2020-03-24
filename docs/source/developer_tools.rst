Developer Tools
===============

Using SWIG to generate Fortran
------------------------------
Change the declaration for ``DTK_UserApplicationHandle`` by replacing
``struct _DTK_UserApplicationHandle *`` to ``void *`` in ``DTK_C_API.h`` and run

.. code:: bash

    $ swig -fortran dtk.i

using Fortran-enabled SWIG `fork <https://github.com/sethrj/swig>`_ (``fortran``
branch). This would generate two files: ``dtk_wrap.cxx`` and
``DataTransferKit.f90``. Move them to ``DTK_Fortran_wrap.cpp`` and
``DTK_Fortran_API.F90``, respectively.

.. note:

    The Fortran file must have the uppercase extension: F90. The only reason for
    that is that preprocessing with Doxygen would not honor
    DOXYGEN_SHOULD_SKIP_THIS otherwise.
