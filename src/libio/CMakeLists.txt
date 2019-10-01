target_sources(openrdm
               INTERFACE
	       ${CMAKE_CURRENT_LIST_DIR}/memRW.cc
	       ${CMAKE_CURRENT_LIST_DIR}/memRW.h
	       ${CMAKE_CURRENT_LIST_DIR}/diskRW.h
	       ${CMAKE_CURRENT_LIST_DIR}/rdmInterface.h
)

target_include_directories(openrdm INTERFACE ${CMAKE_CURRENT_LIST_DIR})
