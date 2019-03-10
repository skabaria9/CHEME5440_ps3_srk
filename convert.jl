#Converting from mM to mmol/gDW
vol_cell = 3.7e-12;                                       #L/cell (HeLa) BIND: 105879
fraction_cell_water = 0.798;                               # RBC, BIND: 101723
mass_of_single_cell = 2.3e-9;                              # g/cell BIND:103720
dryweight_cell = (1-fraction_cell_water)*mass_of_single_cell;
convert_factor(uM) = uM*(vol_cell)/dryweight_cell       #umol/gDW
