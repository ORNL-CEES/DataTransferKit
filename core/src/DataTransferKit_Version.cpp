
#include "DataTransferKit_Version.hpp"
#include "Trilinos_version.h"

std::string DataTransferKit::DataTransferKit_Version()
{ 
  return("DataTransferKit built against Trilinos " TRILINOS_VERSION_STRING); 
}
