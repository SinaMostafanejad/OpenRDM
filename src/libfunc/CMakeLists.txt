foreach(idx functional.cc
            functional.h)
target_sources(openrdm
               INTERFACE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/${idx}>
	       $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)
endforeach()

target_include_directories(openrdm
                           INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
                           $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)
