include("Include.jl")

# load the data dictionary -
# data_dictionary = maximize_urea_production(0,0,0)
data_dictionary = maximize_urea_production_open(0,0,0)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = calculate_optimal_flux_distribution(data_dictionary)

#display(flux_array);
println("Flux array", flux_array)

println("Maximimum urea production = ",flux_array[10],"mmol/gdw-hr")
