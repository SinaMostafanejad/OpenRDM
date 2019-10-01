target_sources(openrdm
	       INTERFACE
               ${CMAKE_CURRENT_LIST_DIR}/mcpdft.cc
               ${CMAKE_CURRENT_LIST_DIR}/mcpdft.h
               ${CMAKE_CURRENT_LIST_DIR}/energy.cc
               ${CMAKE_CURRENT_LIST_DIR}/energy.h
)

target_include_directories(openrdm INTERFACE ${CMAKE_CURRENT_LIST_DIR})
