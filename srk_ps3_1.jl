#How to add the package
#add https://github.com/varnerlab/CoreEcoliModelKit.git

#KEGG Diagram
#https://www.genome.jp/kegg-bin/show_pathway?map00220+C00624

include("flux.jl")
include("convert.jl")

#-----------------------
#PART A
#------------------------

#Stoichiometic Matrix without considering extra metabolizes
stoichiometric_matrix_original=[ #v1 v2 v3 v4 v5,1 v5,-1 b1 b2 b3 b4
-1 0  0  0  0  0  0 1  0  0; #Aspartate
1  -1 0  0  0  0  0 0  0  0; #Arginosuccitate
0  1  0  0  0  0  0 0 -1  0; #Fumarate
0  1  -1 0  1  -1 0 0  0  0; #Arginine
0  0  1  0  0  0  0 0  0 -1; #Urea
0  0  1  -1 0  0  0 0  0  0; #Orinithine
0  0  0  -1 0  0  1 0  0  0; #Carbamoyl Phosphate
-1 0  0  1  -1 1  0 0  0  0; #Citruline
];

#---------------------------------------------------------------------------
#Part B
#-------------------------------------------------------------------------

#Define atom matrix for the original S matrix stoichiometric_matrix_original
atom_original = [
#Columns: Asp Arginosucc Fum Argini Ur  Orin CP Cit
4 10 4 6  1 5  1 6; #C
7 18 4 14 4 12 4 13; #H
1 4  0 4  2 2  1 3; #N
4 6  4 2  1 2  5 3; #O
0 0  0 0  0 0  1 0; #P
0 0  0 0  0 0  0 0; #S
];

#Define the balance matrix
E_original = atom_original * stoichiometric_matrix_original
#Print the E_orginal in the command line and notice there are many non-zero values

stoichiometric_matrix_balanced = [ #v1 v2 v3 v4 v5,1 v5,-1 b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b17 b18 b19 b20 b21
-1 0  0  0  0    0    0 1  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Aspartate
1  -1 0  0  0    0    0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Arginosuccitate
0  1  0  0  0    0    0 0 -1  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Fumarate
0  1  -1 0  1    -1   0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Arginine
0  0  1  0  0    0    0 0  0 -1 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Urea
0  0  1  -1 0    0    0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Orinithine
0  0  0  -1 0    0    1 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Carbamoyl Phosphate
-1 0  0  1  -1   1    0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #Citruline
-1 0  0  0  0    0    0 0  0  0 1  0   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #ATP
1  0  0  0  0    0    0 0  0  0 0 -1   0 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #AMP
1  0  0  0  0    0    0 0  0  0 0  0  -1 0  0  0  0  0  0  0  0  0  0  0  0  0 0; #PPi
0  0  -1 0  -2   2    0 0  0  0 0  0   0 1  0  0  0  1  0  0  0  0  0  0 -1  0 0; #H20
0  0  0  1  0    0    0 0  0  0 0  0   0 0 -1  0  0  0  0  0  0  0  0  0  0  0 0; #phosphate
0  0  0  0  1.5 -1.5  0 0  0  0 0  0   0 0  0  0  0  0 -1  0  1  0  0  0  0  0 0; #NADPH
0  0  0  0 -1.5  1.5  0 0  0  0 0  0   0 0  0  0  1  0  0  0  0  0  0 -1  0  0 0; #NADP+
0  0  0  0 -1    1    0 0  0  0 0  0   0 0  0  1  0  0  0  0  0  0 -1  0  0  0 0; #Nitric oxide
0  0  0  0  1.5 -1.5  0 0  0  0 0  0   0 0  0  0  0  0  0 -1  0  1  0  0  0  0 0; #H+
0  0  0  0  2   -2    0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0 -1 1 #O2
];

#Define atom matrix for the original S matrix stoichiometric_matrix_original
atom_balanced = [
#Columns: Asp Arginosucc Fum Argini Ur  Orin CP Cit ATP AMP PPi H20 Phosphate NADPH NADP+ Nitric_Oxide H+ O2
4 10 4 6  1 5  1 6  10 10 0 0 0 21 21 0 0 0; #C
7 18 4 14 4 12 4 13 16 14 4 2 3 30 29 0 1 0; #H
1 4  0 4  2 2  1 3  5  5  0 0 0  7  7 1 0 0; #N
4 6  4 2  1 2  5 3  13 7  7 1 4 17 17 1 0 2; #O
0 0  0 0  0 0  1 0  3  1  2 0 1  3  3 0 0 0; #P
0 0  0 0  0 0  0 0  0  0  0 0 0  0  0 0 0 0; #S
];

E_balanced = atom_balanced*stoichiometric_matrix_balanced;
#If printed in the command line, E_balances show that v1-v5,-1 are all balanced

#---------------------------------------------------------------------------
#Part C
#-------------------------------------------------------------------------

#Define Inputs to the flux array

#Bound definitions in separate file
include("upper_bound_math.jl")

#Lower Bound
flux_lower_bound = [
v_lower_bound_irrev; #v1
v_lower_bound_irrev; #v2
v_lower_bound_irrev; #v3
v_lower_bound_irrev; #v4
v_lower_bound_irrev; #v5,1
v_lower_bound_irrev; #v5,-1
b_lower_bound_irrev; #b1
b_lower_bound_irrev; #b2
b_lower_bound_irrev; #b3
b_lower_bound_irrev; #b4
b_lower_bound_irrev; #b5
b_lower_bound_irrev; #b6
b_lower_bound_irrev; #b7
b_lower_bound_irrev; #b8
b_lower_bound_irrev; #b9
b_lower_bound_irrev; #b10
b_lower_bound_irrev; #b11
b_lower_bound_irrev; #b12
b_lower_bound_irrev; #b13
b_lower_bound_irrev; #b14
b_lower_bound_irrev; #b15
b_lower_bound_irrev; #b16
b_lower_bound_irrev; #b17
b_lower_bound_irrev; #b18
b_lower_bound_irrev; #b19
b_lower_bound_irrev; #b20
b_lower_bound_irrev #b21
];

flux_upper_bound = [
v1_upper_bound; #v1
v2_upper_bound; #v2
v3_upper_bound; #v3
v4_upper_bound; #v4
v5pos1_upper_bound; #v5,1
v5neg1_upper_bound; #v5,-1
b_upper_bound; #b1
b_upper_bound; #b2
b_upper_bound; #b3
b_upper_bound; #b4
b_upper_bound; #b5
b_upper_bound; #b6
b_upper_bound; #b7
b_upper_bound; #b8
b_upper_bound; #b9
b_upper_bound; #b10
b_upper_bound; #b11
b_upper_bound; #b12
b_upper_bound; #b13
b_upper_bound; #b14
b_upper_bound; #b15
b_upper_bound; #b16
b_upper_bound; #b17
b_upper_bound; #b18
b_upper_bound; #b19
b_upper_bound; #b20
b_upper_bound; #b21
];

#Compile Bounds
default_bounds_array = [flux_lower_bound flux_upper_bound];

#Species Bounds Array
species_bounds_array = [
0.0 0.0; #Aspartate
0.0 0.0; #Arginosuccitate
0.0 0.0; #Fumarate
0.0 0.0; #Arginine
0.0 0.0; #Urea
0.0 0.0; #Orinithine
0.0 0.0; #Carbamoyl Phosphate
0.0 0.0; #Citruline
0.0 0.0; #ATP
0.0 0.0; #AMP
0.0 0.0; #PPi
0.0 0.0; #H20
0.0 0.0; #phosphate
0.0 0.0; #NADPH
0.0 0.0; #NADP+
0.0 0.0; #Nitric oxide
0.0 0.0; #H+
0.0 0.0; #O2
];

#Objective Coefficient Array
objective_coefficient_array =[
0.0; #v1
0.0; #v2
0.0; #v3
0.0; #v4
0.0; #v5,1
0.0; #v5,-1
0.0; #b1
0.0; #b2
0.0; #b3
-1.0; #b4
0.0; #b5
0.0; #b6
0.0; #b7
0.0; #b8
0.0; #b9
0.0; #b10
0.0; #b11
0.0; #b12
0.0; #b13
0.0; #b14
0.0; #b15
0.0; #b16
0.0; #b17
0.0; #b18
0.0; #b19
0.0; #b20
0.0; #b21
];

#Run the flux function
flux_answer = calculate_optimal_flux_distribution(stoichiometric_matrix_balanced,default_bounds_array,species_bounds_array,objective_coefficient_array)
#Pull out the answers
objective_value= flux_answer[1] #umol/gDW per second
objective_value_converted = objective_value/1000*3600 #mmol/gDW per hour
calculated_flux_array= flux_answer[2] #umol/gDW per second
calculated_flux_array_converted = calculated_flux_array/1000*3600 #mmol/gDW per hour
dual_value_array = flux_answer[3]
uptake_array = flux_answer[4]
