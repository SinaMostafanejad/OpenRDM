if(${WITH_PSI4})
   if(${psi4_FOUND})
      foreach(idx plugin.cc                          
                  mcpdft_solver.cc     
                  read_opdm.cc
                  read_tpdm.cc
                  lsda_gga_functionals.cc
                  rs_functionals.cc
                  read_range_separated_integrals.cc)
              list(APPEND psi4_i_list_ ${idx})
      endforeach()
      add_psi4_plugin(mcpdft
               	      MODULE
                      ${CMAKE_CURRENT_LIST_DIR}/${psi4_i_list_}
      )
      unset(psi4_i_list_)
   else()
      add_library(mcpdft STATIC "")
      foreach(idx plugin.cc                          
                  mcpdft_solver.cc     
                  read_opdm.cc
                  read_tpdm.cc
                  lsda_gga_functionals.cc
                  rs_functionals.cc
                  read_range_separated_integrals.cc
		  pymodule.py
		  __init__.py)
      target_sources(mcpdft
                     PUBLIC
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/${idx}>
                     $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
      )
      endforeach()
   endif()
   set_target_properties(mcpdft
      PROPERTIES
         POSITION_INDEPENDENT_CODE ON
         ARCHIVE_OUTPUT_NAME "mcpdft"
         DEBUG_POSTFIX "_d"
	 #RELEASE_POSTFIX "_s"
         SKIP_BUILD_RPATH OFF
         BUILD_WITH_INSTALL_RPATH OFF
         INSTALL_RPATH_USE_LINK_PATH ON
   )
   target_include_directories(mcpdft
                              PUBLIC
                              $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
                              $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
   )
   target_link_libraries(openrdm INTERFACE mcpdft)
endif()
