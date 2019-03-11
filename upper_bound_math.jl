include("convert.jl") #allows use of "convert_factor" function to convert from uM to umol/gDW

k_cat = [203;34.5;249;88.1;13.7]; #kcat for v1 to v5 (0 as placeholder for b's),s^-1
E = 0.01; #umol/gDW, steady state enzyme concentration, given in PS3
theta = 1; #regulatory control, assume is 1

#metabolite_concentration converted to umol/gDW, taken from Park Paper
metabolite_concentration_Aspartate = convert_factor(1.49E-2 * 10^6); #Aspartate
metabolite_concentration_Citruline = convert_factor(2.7E-3 * 10^6); #Citruline
metabolite_concentration_ATP = convert_factor(2.25E-3 * 10^6); #ATP
metabolite_concentration_Fumarate = convert_factor(2.88E-4 * 10^6); #Fumarate, for E. coli (don't need it)
metabolite_concentration_Arginine = convert_factor(2.55e-4 * 10^6); #Arginine
metabolite_concentration_Orinithine = convert_factor(4490);#Ornithine
metabolite_concentration_NADPH = convert_factor(100E-2 * 10^6); #Aspartate


#Km  -- uM converted to umol/gDW
Km_v1_Aspartate = convert_factor(0.18e-3 * 10^6); #v1, 6.3.4.5, Aspartate (BRENDA Mus musculus)
Km_v1_Citruline = convert_factor(0.056e-3 * 10^6); #v1, 6.3.4.5, Citruline, (BRENDA, H. sapien, wild type)
Km_v1_ATP = convert_factor(0.051e-3 * 10^6); #v1, 6.3.4.5, ATP
Km_v2_Argino = convert_factor(0.1e-3 * 10^6); #v2, 4.3.2.1, Arginosuccitate, (BRENDA H. sapien)
Km_v3_Arginine = convert_factor(2.05e-3 * 10^6); #v3, 3.5.3.1, Arginine, (BRENDA, Rattus norvegicus)
Km_v4_Ornithine = convert_factor(0.8e-3 * 10^6); #v4, 2.1.3.3, Ornithine, (BRENDA, Pseudomonas aeruginosa)
Km_v4_Carbamoyl = convert_factor(0.13e-3 * 10^6); #v4, 2.1.3.3, Carbamoyl Phosphate, (BRENDA, Homo sapiens, wild type)
Km_v5_neg1_Arginine = convert_factor(3.50E-6 * 10^6); #v5,-1, 1.14.14.39, Arginine
Km_v5_neg1_NADPH = convert_factor(3.50E-6 * 10^6); #v5,-1, 1.14.14.39, NADPH
Km_v5_pos1 = 0; #Assume saturation

saturation_term_v1 = 0; #v1, 6.3.4.5
saturation_term_v2 = 0; #v2, 4.3.2.1
saturation_term_v3 = 0;#v3, 4.3.2.1
saturation_term_v4 = 0;#(metabolite_concentration(2)/(Km(5)+metabolite_concentration(2)))*(metabolite_concentration(2)/(Km(5)+metabolite_concentration(2))); #v4, 2.1.3.3
saturation_term_v5pos1 = 0;#(metabolite_concentration(2)/(Km(3)+metabolite_concentration(2))); #v5,1 1.14.14.39
saturation_term_v5neg1 = 0;#(metabolite_concentration(2)/(Km(3)+metabolite_concentration(2))); #v5,-1 1.14.14.39


v1_upper_bound = k_cat[1]*E*theta*saturation_term_v1;#v1, 6.3.4.5,
v2_upper_bound = k_cat[2]*E*theta*saturation_term_v2;#v2, 4.3.2.1
v3_upper_bound = k_cat[3]*E*theta*saturation_term_v3;#v3, 3.5.3.1
v4_upper_bound = k_cat[4]*E*theta*saturation_term_v4; #v4, 2.1.3.3
v5pos1_upper_bound = k_cat[5]*E*theta*saturation_term_v5pos1;#v5,1 forwards, 1.14.14.39
v5neg1_upper_bound = k_cat[5]*E*theta*saturation_term_v2;#v5,-1, backwards 1.14.14.39
