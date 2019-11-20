%module DataTransferKit

%apply int { MPI_Comm };
%typemap(ftype) MPI_Comm
   "integer"
%typemap(fin, noblock=1) MPI_Comm {
    $1 = int($input, C_INT)
}
%typemap(fout, noblock=1) MPI_Comm {
    $result = int($1)
}

%typemap(in, noblock=1) MPI_Comm {
    $1 = MPI_Comm_f2c(%static_cast(*$input, MPI_Fint));
}
%typemap(out, noblock=1) MPI_Comm {
    $result = %static_cast(MPI_Comm_c2f($1), int);
}

%{
#include "DTK_C_API.h"
%}

%rename DTK_initializeCmd DTK_initialize_cmd;
%rename DTK_isInitialized DTK_is_initialized;

%rename DTK_isValidUserApplication DTK_is_valid_user_application;
%rename DTK_createUserApplication DTK_create_user_application;
%rename DTK_destroyUserApplication DTK_destroy_user_application;

#if 0
%rename DTK_createMap DTK_create_map;
%rename DTK_isValidMap DTK_is_valid_map;
%rename DTK_applyMap DTK_apply_map;
%rename DTK_destroyMap DTK_destroy_map;
#else
%ignore _DTK_MapHandle;
%ignore DTK_createMap;
%ignore DTK_isValidMap;
%ignore DTK_applyMap;
%ignore DTK_destroyMap;
#endif

// The following typemaps are required so that the treatment of `void (*f)()` is generic.
// The main reason is that in C (but not in C++), empty parenthesis have a special meaning
//   http://port70.net/~nsz/c/c99/n1256.html#6.7.5.3p14
// However, it is hard for SWIG to determine that. The new argument check
// funcionality in SWIG/Fortran changed the behavior. This restores it.
%apply void* { void (*f)() } ;
%typemap(ctype) void (*f)() "$1_ltype"
%typemap(imtype, in="type(C_FUNPTR), value")            void (*f)() "type(C_FUNPTR)"
%typemap(ftype, in="type(C_FUNPTR), intent(in), value") void (*f)() "type(C_FUNPTR)"
%typemap(bindc, in="type(C_FUNPTR), value")             void (*f)() "type(C_FUNPTR)"
%rename DTK_setUserFunction DTK_set_user_function;

%include <std_string.i>

%rename DTK_string_version DTK_version;
%rename DTK_string_git_commit_hash DTK_git_commit_hash;
%inline %{
  std::string DTK_string_version() {
    return std::string(DTK_version());
  }
  std::string DTK_string_git_commit_hash() {
    return std::string(DTK_gitCommitHash());
  }
%}
%ignore DTK_version;
%ignore DTK_gitCommitHash;
%ignore DTK_error(int);

%include "DTK_CellTypes.h"
%include "DTK_C_API.h"
