if(WITH_LIBXC)
   list(APPEND EXTLIBLIST "libxc")
endif()

if(WITH_HDF5)
   list(APPEND EXTLIBLIST "hdf5")
endif()

if(WITH_PYBIND)
   list(APPEND EXTLIBLIST "pybind11")
endif()


foreach(idx IN LISTS EXTLIBLIST)
        add_subdirectory(${idx})
endforeach()
