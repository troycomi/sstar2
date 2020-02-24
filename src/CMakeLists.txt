add_library(vcf_file STATIC vcf_file.cc ${SStar_SOURCE_DIR}/include/sstar2/vcf_file.h)
target_include_directories(vcf_file PUBLIC ../include)

add_library(population_data population_data.cc ${SStar_SOURCE_DIR}/include/sstar2/vcf_file.h)
target_include_directories(population_data PUBLIC ../include)

add_library(window_generator window_generator.cc
    ${SStar_SOURCE_DIR}/include/sstar2/window_generator.h)
target_include_directories(window_generator PUBLIC ../include)
target_link_libraries(window_generator
    vcf_file population_data)

add_library(sstar sstar.cc ${SStar_SOURCE_DIR}/include/sstar2/sstar.h)
target_include_directories(sstar PUBLIC ../include)
target_link_libraries(sstar
    window_generator population_data vcf_file)

add_executable(sstar2 main.cc ${SStar_SOURCE_DIR}/include/CLI11/CLI11.hpp)
# sstar window_generator population_data vcf_file
target_include_directories(sstar2 PUBLIC ../include)
target_link_libraries(sstar2
    sstar window_generator)
