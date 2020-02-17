if(WITH_PSI4)
   if(psi4_FOUND)
      get_property(loc_ TARGET psi4::core PROPERTY LOCATION)
      message(STATUS "${Cyan}Psi4 (${psi4_VERSION}) is found in ${loc_} ${ColourReset}")
      unset(loc_)
      foreach(idx plugin.cc                          
                  mcpdft_solver.cc     
                  read_opdm.cc
                  read_tpdm.cc
                  lsda_gga_functionals.cc
                  rs_functionals.cc
                  read_range_separated_integrals.cc)
              list(APPEND psi4_i_list_ ${idx})
      endforeach()
      add_psi4_plugin(mcpdft_iface
                      STATIC
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/${psi4_i_list_}>
                     $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
      )
   else()
      add_library(mcpdft_iface STATIC "")
      foreach(idx plugin.cc                          
                  mcpdft_solver.cc     
                  read_opdm.cc
                  read_tpdm.cc
                  lsda_gga_functionals.cc
                  rs_functionals.cc
                  read_range_separated_integrals.cc)
      target_sources(mcpdft_iface
                     PUBLIC
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/${idx}>
                     $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
      )
      endforeach()
   endif()
   target_link_libraries(mcpdft_iface PUBLIC psi4::core)
   target_include_directories(mcpdft_iface
                              PUBLIC
                              ${PYTHON_INCLUDE_DIRS}
   )
   target_link_libraries(mcpdft_iface
                         PUBLIC
                         ${PYTHON_LIBRARIES}
   )
   set_target_properties(mcpdft_iface
      PROPERTIES
         POSITION_INDEPENDENT_CODE 1
         ARCHIVE_OUTPUT_NAME "mcpdft_iface"
         DEBUG_POSTFIX "_sd"
         RELEASE_POSTFIX "_s"
         SKIP_BUILD_RPATH OFF
         BUILD_WITH_INSTALL_RPATH OFF
         INSTALL_RPATH_USE_LINK_PATH ON
   )
   target_include_directories(mcpdft_iface
                              PUBLIC
                              $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
                              $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
   )
   target_link_libraries(openrdm PUBLIC mcpdft_iface)
endif()

