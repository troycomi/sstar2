# add in macros for testing
macro(package_add_test TESTNAME FILES LIBRARIES)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${FILES})
    # link the Google test infrastructure, mocking library,
    target_link_libraries(${TESTNAME} ${LIBRARIES})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    # gtest_discover_tests replaces gtest_add_tests,
    gtest_discover_tests(${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
        )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    target_include_directories(${TESTNAME} PUBLIC include)
    add_test(NAME, ${TESTNAME} COMMAND ${TESTNAME})
endmacro()

package_add_test(vcf_file_test test_vcf_file.cc vcf_file)
package_add_test(vcf_entry_test test_vcf_entry.cc vcf_file)
package_add_test(population_data_test test_population_data.cc population_data)
package_add_test(window_test test_window.cc window)
package_add_test(window_generator_test test_window_generator.cc window_generator)
package_add_test(sstar_test test_sstar.cc sstar)
package_add_test(validator_test test_validator.cc validator)
