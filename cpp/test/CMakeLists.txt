enable_testing()

function(add_unit_test name)
  add_executable(${name}
    ${name}.cpp
    )
  target_link_libraries(${name}
    PRIVATE
      eos
      GTest::gtest_main
      GSL::gsl
    )
endfunction()

add_unit_test(cubic_eos_test)

include(GoogleTest)
gtest_discover_tests(cubic_eos_test)