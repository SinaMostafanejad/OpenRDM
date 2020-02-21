foreach(idx mcpdft.cc
            mcpdft.h
            energy.cc
            energy.h
	    utility.cc)
target_sources(openrdm
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/${idx}>
               $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)
endforeach()

target_include_directories(openrdm
                           PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
                           $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)

target_include_directories(openrdm PUBLIC ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(openrdm PUBLIC ${ARMADILLO_LIBRARIES})

