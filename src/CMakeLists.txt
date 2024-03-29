add_library(vcf_file STATIC vcf_file.cc ${SStar_SOURCE_DIR}/include/sstar2/vcf_file.h)
target_include_directories(vcf_file PUBLIC ../include)

add_library(population_data population_data.cc ${SStar_SOURCE_DIR}/include/sstar2/vcf_file.h)
target_include_directories(population_data PUBLIC ../include)

add_library(validator validator.cc
    ${SStar_SOURCE_DIR}/include/sstar2/validator.h)
target_include_directories(validator PUBLIC ../include)
target_link_libraries(validator
    vcf_file)

add_library(window window.cc
    ${SStar_SOURCE_DIR}/include/sstar2/window.h)
target_include_directories(window PUBLIC ../include)
target_link_libraries(window
    vcf_file validator)

add_library(window_generator window_generator.cc
    ${SStar_SOURCE_DIR}/include/sstar2/window_generator.h)
target_include_directories(window_generator PUBLIC ../include)
target_link_libraries(window_generator
    vcf_file population_data validator window)

add_library(sstar sstar.cc ${SStar_SOURCE_DIR}/include/sstar2/sstar.h)
target_include_directories(sstar PUBLIC ../include)
target_link_libraries(sstar
    window_generator population_data vcf_file)

add_executable(sstar2 main.cc)
# sstar window_generator population_data vcf_file
target_include_directories(sstar2 PUBLIC ../include)
target_link_libraries(sstar2
    sstar window_generator validator CLI11::CLI11)
