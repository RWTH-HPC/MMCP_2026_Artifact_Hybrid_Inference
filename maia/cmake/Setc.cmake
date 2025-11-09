# Include guard
if (__setc)
  return()
endif ()
set(__setc YES)

# Join all arguments into one string and save it to _var
function(setc _var)
  join(_var "" ${ARGN})
  set(${_var} ${_var} PARENT_SCOPE)
endfunction()

