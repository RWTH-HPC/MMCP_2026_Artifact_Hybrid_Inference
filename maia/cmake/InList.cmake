# Include guard
if (__in_list)
  return()
endif ()
set(__in_list YES)

# Check if _value id among the passed arguments and save YES or NO in _var
function(in_list _var _value)
  list(FIND ARGN ${_value} _index)
  if (${_index} GREATER -1)
    set(${_var} "YES" PARENT_SCOPE)
  else ()
    set(${_var} "NO" PARENT_SCOPE)
  endif ()
endfunction()
