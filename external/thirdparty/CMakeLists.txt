if(WITH_ARMADILLO)
   list(APPEND THIRDLIBLIST "armadillo")
endif()

foreach(idx IN LISTS THIRDLIBLIST)
        add_subdirectory(${idx})
endforeach()
