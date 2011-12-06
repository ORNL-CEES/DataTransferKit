INCLUDE(TribitsTplDeclareLibraries)

TRIBITS_TPL_DECLARE_LIBRARIES( Nemesis
  REQUIRED_HEADERS DBC.hh
  REQUIRE_LIBS_NAMES nemesis_harness nemesis_comm
)