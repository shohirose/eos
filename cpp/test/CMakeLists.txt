function(add_unit_test name)
  add_executable(${name}
    ${name}.cpp
    )
  target_link_libraries(${name}
    PRIVATE
      eos
      gtest_main
      GSL::gsl
    )
  add_test(
    NAME ${name}
    COMMAND ${name}
    )
endfunction()

add_unit_test(cubic_eos_test)