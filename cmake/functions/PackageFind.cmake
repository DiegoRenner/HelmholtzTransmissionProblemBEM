# This file includes the macro to easily find a package using
# either machine or hunter
MACRO(numpde_find_package package_name)
  # First we try to find the package
  find_package(${package_name})

  if (${package_name}"_NOT_FOUND")
    hunter_add_package(${package_name})
    
    find_package(${package_name} REQUIRED)
  endif()
ENDMACRO()
