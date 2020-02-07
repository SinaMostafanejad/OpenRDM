add_library(openrdm INTERFACE)
add_library(mcpdft::openrdm ALIAS openrdm)

foreach(idx libmcpdft
	    libmem
            libio
            libfunc)
       add_subdirectory(${idx})
endforeach()

target_include_directories(openrdm INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libmcpdft>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libmem>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libio>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libfunc>
  #  $<INSTALL_INTERFACE:${OPENRDM_STAGING_PREFIX}/openrdm/include
)

if(${WITH_OPENMP})
  find_package(OpenMP REQUIRED)
  target_link_libraries(openrdm INTERFACE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(openrdm INTERFACE coverage_config)

install(TARGETS openrdm
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	INCLUDES ${CMAKE_INSTALL_INCLUDEDIR}
)

