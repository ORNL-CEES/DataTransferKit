Coding Style Guidelines
=======================

DataTransferKit developers follow a set of style guidelines and use the clang
format tool to create a consistent appearance to all source code committed to
the repository.

ClangFormat
-----------

ClangFormat (version 3.9) is used to check the C++ code formatting style in
DTK.  A pull request that does not comply will be rejected. Configure with
``-D DataTransferKit_ENABLE_ClangFormat=ON`` and do ``make format-cpp`` to
apply the formatting style before your commit.  Alternatively, run
``ctest -V -R check_format_cpp`` display the diff without applying the
changes.

Style Guide
-----------

The following conventions are used in the code.

Names of classes, structs, and enumerations are camel case and capitalized:

.. code-block:: c++

   class ExampleClassName {};

Function names are camel case and not capitalized:

.. code-block:: c++

   void exampleFunctionName() const;

Variable names are lower case and have underscores to separate words:

.. code-block:: c++

   double example_double_var;
   std::vector<int> example_vec_var;

If a variable is also class data that is private, protected, or public it is
prefixed with a ``d_``:

.. code-block:: c++

   class ExampleClass
   {
     private:
       double d_class_double_var;
       std::vector<int> d_example_vec_var;
   };

The clang format tool described above enforces spacing, line breaks, and other
general file formatting requirements. Header files are suffixed with ``.hpp``
and non-templated implementation files are suffixed with ``.cpp``. Header
guards are needed for all header files following the convention of
``DTK_CLASSNAME_HPP``. For example:

.. code-block:: c++

  #ifndef DTK_EXAMPLECLASS_HPP
  #define DTK_EXAMPLECLASS_HPP

  class ExampleClass
  {
      // Class definition...
  };
  
  #endif

