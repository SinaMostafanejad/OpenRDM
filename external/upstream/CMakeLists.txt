foreach(idx IN LISTS EXTLIBLIST)
        add_subdirectory(${idx})
endforeach()
