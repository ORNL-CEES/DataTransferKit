Coding Style Guidelines
=======================

DataTransferKit developers follow a set of style guidelines and use the clang
format tool to create a consistent appearance to all source code committed to
the repository.

ClangFormat
-----------

ClangFormat (version 6.0) is used to check the C++ code formatting style in
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

If a variable is class data prefixed with a ``_``:

.. code-block:: c++

   class ExampleClass
   {
     public:
       int _a_public_var;
     private:
       double _class_double_var;
       std::vector<int> _example_vec_var;
     protected:
       std::string _a_protected_string;
   };

Previously, the convention for class data was to prefix the variable name with
``d_`` so this will be seen throughout the code. We will be transitioning to
the ``_`` prefix convention in future work and slowly transition existing
code.

When a function as both input and output arguments, the inputs should come
first:

.. code-block:: c++

  void myFunction(int const input_1, int const input_2, int &output)

The ``clang-format`` tool described above enforces spacing, line breaks, and
other general file formatting requirements. Header files are suffixed with
``.hpp`` and non-templated implementation files are suffixed with
``.cpp``. Header guards are needed for all header files following the
convention of ``DTK_CLASSNAME_HPP``. For example:

.. code-block:: c++

  #ifndef DTK_EXAMPLECLASS_HPP
  #define DTK_EXAMPLECLASS_HPP

  class ExampleClass
  {
      // Class definition...
  };

  #endif

