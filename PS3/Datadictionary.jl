function maximize_urea_production_open(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = maximize_urea_production(time_start,time_stop,time_step)

	# lets open up the side products of ec:1.14.13.39 -

	# 2: Update the reaction bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
	default_flux_bounds_array[15,1] = 0
	default_flux_bounds_array[15,2] = 10

	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
	default_flux_bounds_array[16,1] = 0
	default_flux_bounds_array[16,2] = 10

	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
	default_flux_bounds_array[17,1] = 0
	default_flux_bounds_array[17,2] = 10

	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
	default_flux_bounds_array[18,1] = 0
	default_flux_bounds_array[18,2] = 10

	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
	default_flux_bounds_array[19,1] = 0
	default_flux_bounds_array[19,2] = 10

	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
	default_flux_bounds_array[20,1] = 0
	default_flux_bounds_array[20,2] = 10

	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	default_flux_bounds_array[21,1] = 0
	default_flux_bounds_array[21,2] = 10

	# repackage -
	data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array

	# return the updated dictionary -
	return data_dictionary
end

function maximize_urea_production(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# 1: set the objective function -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[10] = -1

	# 2: Update the reaction bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# let all exchanges go from 0,10 mmol/gDW-hr
	range_of_exchange_reactions = collect(7:21)
	for reaction_index in range_of_exchange_reactions
		default_flux_bounds_array[reaction_index,1] = 0.0
		default_flux_bounds_array[reaction_index,2] = 10.0
	end

	# don't allow water exchange -
	default_flux_bounds_array[20,2] = 0
	default_flux_bounds_array[21,2] = 0
	#set up (s/s+Km) value array from Park et al. Supplementary doc 2

	s= [

		0.923*0.99*1	;	# 1 M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		1	            ;	# 2 M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		0.142*1	      	; 	# 3 M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		1*1	    		; 	# 4 M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		0.986*1      	; 	# 5 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		1             	; 	#6  2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c
	];

	# we have some specific values for v1 -> v5
	E = (0.01)*(1/1000)	# mmol/gDW
	metabolic_vmax_array = [
		203*(3600)*E*s[1]	;	# v1 ec:6.3.4.5 mmol/gDW-hr
		34.5*(3600)*E*s[2]	;	# v2 ec:4.3.2.1 mmol/gDW-hr
		249*(3600)*E*s[3]	;	# v3 ec:3.5.3.1 mmol/gDW-hr
		88.1*(3600)*E*s[4]	;	# v4 ec:2.1.3.3 mmol/gDW-hr
		13.7*(3600)*E*s[5]	;	# v5 ec:1.14.13.39 mmol/gDW-hr
		13.7*(3600)*E*s[6]	;	# v6 ec:1.14.13.39 mmol/gDW-hr
	]
	range_of_cycle_reactions = collect(1:6)
	for reaction_index in range_of_cycle_reactions
		default_flux_bounds_array[reaction_index,1] = 0.0
		default_flux_bounds_array[reaction_index,2] = metabolic_vmax_array[reaction_index]
	end

	# repackage -
		data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array
		data_dictionary["objective_coefficient_array"] = objective_coefficient_array

		# return the updated dictionary -
		return data_dictionary
	end






function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("network.dat");

	# What is the system dimension? -
	(number_of_species,number_of_reactions) = size(stoichiometric_matrix)

	E = 0.01*(1/1000); #Steady state concentration of enzyme (mmol/gDW)

	# Setup kcat values array
	k_cat  = [

	203.0	;	# 1 M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
	34.5	;	# 2 M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
	249.0	;	# 3 M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
	88.1	;	# 4 M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
	13.7    ;	# 5 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
	13.7    ;   # 6 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c

	];   # 1/s


		#(s/s+Km) from Park et al.

		s= [

			0.923*0.99*1	;	# 1 M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
			1	            ;	# 2 M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
			0.142*1	     	; 	# 3 M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
			1*1	    		; 	# 4 M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
			0.986*1      	; 	# 5 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
			1             	; 	#6  2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c
		];

		metabolic_vmax_array = v_max = 3600/10^3*E.*k_cat.*s;  #mmol/gdW.hr



	# Setup default flux bounds array -
	default_bounds_array = [
		0	metabolic_vmax_array[1]	;	# Vmax [mmol/gdw-hr] 1	M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		0	metabolic_vmax_array[2]	;	# Vmax [mmol/gdw-hr] 2	M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		0	metabolic_vmax_array[3]	;	# Vmax [mmol/gdw-hr] 3	M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		0	metabolic_vmax_array[4]	;	# Vmax [mmol/gdw-hr] 4	M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		0	metabolic_vmax_array[5]	;	# Vmax [mmol/gdw-hr] 5	2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		0	metabolic_vmax_array[6]	;	# Vmax [mmol/gdw-hr] 6	2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c

		0	2.2148976	;	# Vmax [mmol/gdw-hr] 7	[] --> M_Carbamoyl_phosphate_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 8	[] --> M_L-Aspartate_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 9	M_Fumarate_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 10	M_Urea_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 11	[] --> M_ATP_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 12	M_AMP_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 13	M_Diphosphate_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 14	M_Orthophosphate_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 15	[] --> M_NADPH_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 16	[] --> M_H_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 17	[] -->  M_Oxygen_c
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	];


	# Setup default species bounds array -
	species_bounds_array = [
		0.0	0.0	;	# 1    M_L-Aspartate_c
		0.0	0.0	;	# 2    M_N-(L-Arginino)succinate_c
		0.0	0.0	;	# 3    M_Fumarate_c
		0.0	0.0	;	# 4    M_L-Arginine_c
		0.0	0.0	;	# 5    M_Urea_c
		0.0	0.0	;	# 6    M_L-Ornithine_c
		0.0	0.0	;	# 7    M_Carbamoyl_phosphate_c
		0.0	0.0	;	# 8    M_L-Citrulline_c
		0.0	0.0	;	# 9    M_Diphosphate_c
		0.0	0.0	;	# 10   M_ATP_c
		0.0	0.0	;	# 11   M_AMP_c
		0.0	0.0	;	# 12   M_H2O_c
		0.0	0.0	;	# 13   M_Orthophosphate_c
		0.0	0.0	;	# 14   M_Oxygen_c
		0.0	0.0	;	# 15   M_NADPH_c
		0.0	0.0	;	# 16   M_H_c
		0.0	0.0	;	# 17   M_NADP_c
		0.0	0.0	;	# 18   M_Nitric_oxide_c
	];



	#set up (s/s+Km) value array from Park et al. Supplementary doc 2

	s= [

		0.923*0.99*1	;	# 1 M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		1	            ;	# 2 M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		0.142*1	      	; 	# 3 M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		0.0117*1	    ; 	# 4 M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		 0.986*1     	; 	# 5 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		1             	; 	#6  2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c
	];







	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# Setup the objective coefficient array -
	objective_coefficient_array = [

		0.0	;	# 1 v1::M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		0.0	;	# 2 v2::M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		0.0	;	# 3 v3::M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		0.0	;	# 4 v4::M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		0.0	;	# 5 v5::2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		0.0	;	# 6 v5_reverse::2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c
		0.0	;	# 7 b1::[] --> M_Carbamoyl_phosphate_c
		0.0	;	# 8 b2::[] --> M_L-Aspartate_c
		0.0	;	# 9 b3::M_Fumarate_c --> []
		0.0	;	# 10 b4::M_Urea_c --> []
		0.0	;	# 11 b5::[] --> M_ATP_c
		0.0	;	# 12 b6::M_AMP_c --> []
		0.0	;	# 13 b7::M_Diphosphate_c --> []
		0.0	;	# 14 b8::M_Orthophosphate_c --> []
		0.0	;	# 15 b9::[] --> M_H2O_c
		0.0	;	# 16 b10::[] --> M_Oxygen_c
		0.0	;	# 17 b11::[] -->  M_H_c
		0.0	;	# 18 b12::M_Nitric_oxide_c --> []
		0.0	;	# 19 b13::M_NADP_c --> []
		0.0	;	# 20 b14::M_H2O_c --> []
		0.0	;	# 21 b15::[] --> M_NADPH_c
	];



	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"M_L-Aspartate_c"	;	          				        # 1
		"M_N-(L-Arginino)succinate_c"	;  						# 2
  	    "M_Fumarate_c"		;						            # 3
		"M_L-Arginine_c";	  									# 4
	    "M_Urea_c"	;						                	# 5
	    "M_L-Ornithine_c"	;									# 6
	 	"M_Carbamoyl_phosphate_c"	;							# 7
		"M_L-Citrulline_c"	;									# 8
		"M_Diphosphate_c"	;									# 9
		"M_ATP_c"	;											# 10
		"M_AMP_c"	;											# 11
	    "M_H2O_c"		;										# 12
		"M_Orthophosphate_c"	;								# 13
		"M_Oxygen_c"	;										# 14
		"M_NADPH_c"	;											# 15
  	    "M_H_c"		;											# 16
		"M_NADP_c"	;											# 17
		"M_Nitric_oxide_c"	;									# 18
		];

	data_dictionary = Dict{AbstractString,Any}()
		data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
		data_dictionary["objective_coefficient_array"] = objective_coefficient_array
		data_dictionary["default_flux_bounds_array"] = default_bounds_array;

		data_dictionary["species_bounds_array"] = species_bounds_array

		data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
		data_dictionary["is_minimum_flag"] = is_minimum_flag
		data_dictionary["number_of_species"] = number_of_species
		data_dictionary["number_of_reactions"] = number_of_reactions

		data_dictionary["metabolic_vmax_array"] = metabolic_vmax_array

		return data_dictionary
end
