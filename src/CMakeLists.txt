
project(pnp_hybrid)


add_library(sim_lib SHARED sim_corres.cpp utility.cpp
	solve_5pts.cpp solve_translate.cpp translate_factor.cpp translate_pro_factor.cpp
	opt_solver.cpp rotation_only.cpp)

target_link_libraries(sim_lib ${CERES_LIBRARIES} ${OpenCV_LIBS} opengv)

add_library(opengv_helper SHARED experiment_helpers.cpp time_measurement.cpp random_generators.cpp)
target_link_libraries(opengv_helper opengv)

add_executable(main_sim_mt  main_monte_carlo.cpp)
target_link_libraries(main_sim_mt sim_lib )

add_executable(main_opt main_opt.cpp)
target_link_libraries(main_opt sim_lib)

add_executable(main_opengv_pnp main_opengv_pnp.cpp)
target_link_libraries(main_opengv_pnp sim_lib opengv_helper)

add_executable(main_mt_noise main_mt_noise.cpp)
target_link_libraries(main_mt_noise sim_lib opengv_helper)

add_executable(main_mt_points main_mt_points.cpp)
target_link_libraries(main_mt_points sim_lib opengv_helper)

# add_executable(test_rotation_only test_rotation_only.cpp)
# target_link_libraries(test_rotation_only sim_lib )

############## test iphone12's data #################
add_executable(test_iphone_data test_iphone_data.cpp)
target_link_libraries(test_iphone_data sim_lib)