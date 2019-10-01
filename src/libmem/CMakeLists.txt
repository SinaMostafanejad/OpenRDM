target_sources(openrdm
 	       INTERFACE
	       ${CMAKE_CURRENT_LIST_DIR}/libMem.cc
	       ${CMAKE_CURRENT_LIST_DIR}/libMem.h
)

target_include_directories(openrdm INTERFACE ${CMAKE_CURRENT_LIST_DIR})
