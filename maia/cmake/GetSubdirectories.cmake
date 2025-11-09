# Include guard
if (__get_subdirectories)
  return()
endif ()
set(__get_subdirectories YES)

# Get a sorted list of subdirectories in directory _path and return it in _var
function(get_subdirectories _var _path)
  file(GLOB _nodes RELATIVE ${_path} ${_path}/*)
  set(dirlist "")
  foreach (_node ${_nodes})
    if (IS_DIRECTORY ${_path}/${_node})
      set(dirlist ${dirlist} ${_node})
    endif ()
  endforeach ()
  list(SORT dirlist)
  set(${_var} ${dirlist} PARENT_SCOPE)
endfunction()
