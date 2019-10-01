target_sources(openrdm
               INTERFACE
               ${CMAKE_CURRENT_LIST_DIR}/functional.cc
               ${CMAKE_CURRENT_LIST_DIR}/functional.h
)

target_include_directories(openrdm INTERFACE ${CMAKE_CURRENT_LIST_DIR})

