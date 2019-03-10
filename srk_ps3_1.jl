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

stoichiometric_matrix_balanced = [ #v1 v2 v3 v4 v5,1 v5,-1 b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b17 b18 b19
-1 0  0  0  0  0  0 1  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Aspartate
1  -1 0  0  0  0  0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Arginosuccitate
0  1  0  0  0  0  0 0 -1  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Fumarate
0  1  -1 0  1  -1 0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Arginine
0  0  1  0  0  0  0 0  0 -1 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Urea
0  0  1  -1 0  0  0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Orinithine
0  0  0  -1 0  0  1 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Carbamoyl Phosphate
-1 0  0  1  -1 1  0 0  0  0 0  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #Citruline
-1 0  0  0  0  0  0 0  0  0 1  0   0 0  0  0  0  0  0  0  0  0  0  0  0; #ATP
1  0  0  0  0  0  0 0  0  0 0 -1   0 0  0  0  0  0  0  0  0  0  0  0  0; #ADP
1  0  0  0  0  0  0 0  0  0 0  0  -1 0  0  0  0  0  0  0  0  0  0  0  0; #PPi
0  0  -1 0  -1 1  0 0  0  0 0  0   0 1  0  0  0  1  0  0  0  0  0  0 -1; #H20
0  0  0  1  0  0  0 0  0  0 0  0   0 0 -1  0  0  0  0  0  0  0  0  0  0; #phosphate
0  0  0  0  1 -1  0 0  0  0 0  0   0 0  0  0  0  0 -1  0  1  0  0  0  0; #NADPH
0  0  0  0 -1  1  0 0  0  0 0  0   0 0  0  0  1  0  0  0  0  0  0 -1  0; #NADP+
0  0  0  0 -1  1  0 0  0  0 0  0   0 0  0  1  0  0  0  0  0  0 -1  0  0; #Nitric oxide
0  0  0  0  1 -1  0 0  0  0 0  0   0 0  0  0  0  0  0 -1  0  1  0  0  0; #H+
];

#Define atom matrix for the original S matrix stoichiometric_matrix_original
atom_balanced = [
#Columns: Asp Arginosucc Fum Argini Ur  Orin CP Cit ATP ADP PPi H20 Phosphate NADPH NADP+ Nitric_Oxide H+*
4 10 4 6  1 5  1 6  10 10 0 0 0 21 21 0 0; #C
7 18 4 14 4 12 4 13 16 14 4 2 3 30 29 0 1; #H
1 4  0 4  2 2  1 3  5  5  0 0 0  7  7 1 0; #N
4 6  4 2  1 2  5 3  13 7  7 1 4 17 17 1 0; #O
0 0  0 0  0 0  1 0  3  1  2 0 1  3  3 0 0; #P
0 0  0 0  0 0  0 0  0  0  0 0 0  0  0 0 0; #S
];

E_balanced = atom_balanced*stoichiometric_matrix_balanced;

#---------------------------------------------------------------------------
#Part C
#-------------------------------------------------------------------------

#Define Inputs to the flux array

#Lower Bound
flux_lower_bound = [0;0;0;0;0;0;0;0;0;0];

#Upper bound (math in separate file)
include("upper_bound_math.jl")
flux_upper_bound = 0;

#default_bounds_array = [flux_lower_bound flux_upper_bound];

#Species Bounds Array
species_bounds_array = [0;0;0;0;0;0;0;0;0;0];

#Objective Coefficient Array
objective_coefficient_array = [0;0;0;0;0;0;0;0;0;1];
