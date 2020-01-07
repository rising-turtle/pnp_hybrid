
project(pnp_hybrid)

add_executable(main_sim  main.cpp sim_corres.cpp utility.cpp
	solve_5pts.cpp solve_translate.cpp translate_factor.cpp )
target_link_libraries(main_sim ${CERES_LIBRARIES} ${OpenCV_LIBS})