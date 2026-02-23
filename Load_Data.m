% Helper: read a single-column numeric range from Excel as a column vector
readcol = @(file,sheet,range) table2array( ...
    readtable(file,'Sheet',sheet,'Range',range,'ReadVariableNames',false) );

% Exps. 1 & 2, P+M, low P
J_01  = (readcol('Oxygen_Data.xlsx','PYR+MAL (50uM pyr)','H2:H327')-3).*60./674;
eJ_01 = (readcol('Oxygen_Data.xlsx','PYR+MAL (50uM pyr)','I2:I327')).*60./674;
tJ_01 = (2*(0:325) - 2)';
J_02  = (readcol('Oxygen_Data.xlsx','PYR+MAL (50uM pyr)','Q2:Q327')-3).*60./674;
eJ_02 = (readcol('Oxygen_Data.xlsx','PYR+MAL (50uM pyr)','R2:R327')).*60./674;
tJ_02 = (2*(0:325) - 0)';
N_01  = readcol('NADH_Data.xlsx','PYR+MAL (50uM pyr)','M2:M150');
eN_01 = readcol('NADH_Data.xlsx','PYR+MAL (50uM pyr)','N2:N150');
tN_01 = readcol('NADH_Data.xlsx','PYR+MAL (50uM pyr)','A2:A150');
N_02  = readcol('NADH_Data.xlsx','PYR+MAL (50uM pyr)','AA2:AA142');
eN_02 = readcol('NADH_Data.xlsx','PYR+MAL (50uM pyr)','AB2:AB142');
tN_02 = readcol('NADH_Data.xlsx','PYR+MAL (50uM pyr)','A2:A142');

% Exps. 3 & 4, P+M, med P
J_03  = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr)','J2:J374')-3).*60./674;
eJ_03 = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr)','K2:K374')).*60./674;
tJ_03 = (2*(0:372) + 2)';
J_04  = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr)','BA2:BA300')-3).*60./674;
eJ_04 = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr)','BB2:BB300')).*60./674;
tJ_04 = (2*(0:298) + 6)';
N_03  = readcol('NADH_Data.xlsx','PYR+MAL (med pyr)','H2:H100');
eN_03 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr)','I2:I100');
tN_03 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr)','A2:A100');
N_04  = readcol('NADH_Data.xlsx','PYR+MAL (med pyr)','Q2:Q95');
eN_04 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr)','R2:R95');
tN_04 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr)','A2:A95');

% Exps. 7 & 8, P+M, med P, high ADP
J_07  = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr high ADP)','H2:H373')-3).*60./674;
eJ_07 = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr high ADP)','I2:I373')).*60./674;
tJ_07 = (2*(0:371) + 0)';
J_08  = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr high ADP)','Q2:Q332')-3).*60./674;
eJ_08 = (readcol('Oxygen_Data.xlsx','PYR+MAL (med pyr high ADP)','R2:R332')).*60./674;
tJ_08 = (2*(0:330) + 2)';
N_07  = readcol('NADH_Data.xlsx','PYR+MAL (med pyr high ADP)','H2:H100');
eN_07 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr high ADP)','I2:I100');
tN_07 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr high ADP)','A2:A100');
N_08  = readcol('NADH_Data.xlsx','PYR+MAL (med pyr high ADP)','Q2:Q100');
eN_08 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr high ADP)','R2:R100');
tN_08 = readcol('NADH_Data.xlsx','PYR+MAL (med pyr high ADP)','A2:A100');


% Exps. 11 & 12, PYR only, high P
J_11  = (readcol('Oxygen_Data.xlsx','PYR (high pyr)','H2:H324')-3).*60./674;
eJ_11 = (readcol('Oxygen_Data.xlsx','PYR (high pyr)','I2:I324')).*60./674;
tJ_11 = (2*(0:322) - 4)';
J_12  = (readcol('Oxygen_Data.xlsx','PYR (high pyr)','Q2:Q324')-3).*60./674;
eJ_12 = (readcol('Oxygen_Data.xlsx','PYR (high pyr)','R2:R324')).*60./674;
tJ_12 = (2*(0:322) - 4)';
N_11  = readcol('NADH_Data.xlsx','PYR (high pyr)','H2:H100');
eN_11 = readcol('NADH_Data.xlsx','PYR (high pyr)','I2:I100');
tN_11 = readcol('NADH_Data.xlsx','PYR (high pyr)','A2:A100');
N_12  = readcol('NADH_Data.xlsx','PYR (high pyr)','Q2:Q100');
eN_12 = readcol('NADH_Data.xlsx','PYR (high pyr)','R2:R100');
tN_12 = readcol('NADH_Data.xlsx','PYR (high pyr)','A2:A100');

% Exps. 17, 18, S only
N_17  = readcol('NADH_Data.xlsx','SUCC','H2:H100');
eN_17 = readcol('NADH_Data.xlsx','SUCC','I2:I100');
tN_17 = readcol('NADH_Data.xlsx','SUCC','A2:A100');
N_18  = readcol('NADH_Data.xlsx','SUCC','Q2:Q100');
eN_18 = readcol('NADH_Data.xlsx','SUCC','R2:R100');
tN_18 = readcol('NADH_Data.xlsx','SUCC','A2:A100');
J_17  = (readcol('Oxygen_Data.xlsx','SUCC','I2:I349')-3).*60./674;
eJ_17 = (readcol('Oxygen_Data.xlsx','SUCC','J2:J349')).*60./674;
tJ_17 = (2*(0:347) - 3)';
J_18  = (readcol('Oxygen_Data.xlsx','SUCC','S2:S349')-3).*60./674;
eJ_18 = (readcol('Oxygen_Data.xlsx','SUCC','S2:S349')).*60./674;
tJ_18 = (2*(0:347) - 3)';

% Exps. 19, 20, S + M
N_19  = readcol('NADH_Data.xlsx','SUCC+MAL (high MAL)','H2:H100');
eN_19 = readcol('NADH_Data.xlsx','SUCC+MAL (high MAL)','I2:I100');
tN_19 = readcol('NADH_Data.xlsx','SUCC+MAL (high MAL)','A2:A100');
N_20  = readcol('NADH_Data.xlsx','SUCC+MAL (high MAL)','Q2:Q100');
eN_20 = readcol('NADH_Data.xlsx','SUCC+MAL (high MAL)','R2:R100');
tN_20 = readcol('NADH_Data.xlsx','SUCC+MAL (high MAL)','A2:A100');
J_19  = (readcol('Oxygen_Data.xlsx','SUCC+MAL (high MAL)','H2:H336')-3).*60./674;
eJ_19 = (readcol('Oxygen_Data.xlsx','SUCC+MAL (high MAL)','I2:I336')).*60./674;
tJ_19 = (2*(0:334) + 10)';
J_20  = (readcol('Oxygen_Data.xlsx','SUCC+MAL (high MAL)','Q2:Q336')-3).*60./674;
eJ_20 = (readcol('Oxygen_Data.xlsx','SUCC+MAL (high MAL)','R2:R336')).*60./674;
tJ_20 = (2*(0:334) + 10)';

% Exps. 21, 22, M only
J_21  = (readcol('Oxygen_Data.xlsx','MAL (5mM)','H2:H443')-3).*60./674;
eJ_21 = (readcol('Oxygen_Data.xlsx','MAL (5mM)','I2:I443')).*60./674;
tJ_21 = (2*(0:441) - 4)';
J_22  = (readcol('Oxygen_Data.xlsx','MAL (5mM)','Q2:Q443')-3).*60./674;
eJ_22 = (readcol('Oxygen_Data.xlsx','MAL (5mM)','R2:R443')).*60./674;
tJ_22 = (2*(0:441) - 6)';
N_21  = readcol('NADH_Data.xlsx','MAL (5mM)','H2:H100');
eN_21 = readcol('NADH_Data.xlsx','MAL (5mM)','I2:I100');
tN_21 = readcol('NADH_Data.xlsx','MAL (5mM)','A2:A100');
N_22  = readcol('NADH_Data.xlsx','MAL (5mM)','Q2:Q100');
eN_22 = readcol('NADH_Data.xlsx','MAL (5mM)','R2:R100');
tN_22 = readcol('NADH_Data.xlsx','MAL (5mM)','A2:A100');

% Exps. 23, 24, M only (0.25 mM)
J_23  = (readcol('Oxygen_Data.xlsx','MAL (0.25mM)','H2:H443')-3).*60./674;
tJ_23 = (2*(0:314) - 10)';
J_24  = (readcol('Oxygen_Data.xlsx','MAL (0.25mM)','Q2:Q443')-3).*60./674;
tJ_24 = (2*(0:314) - 12)';

% Exps. 25, 26, A + M
J_25  = (readcol('Oxygen_Data.xlsx','AKG+MAL','I2:I315')-3).*60./674;
eJ_25 = (readcol('Oxygen_Data.xlsx','AKG+MAL','J2:J315')).*60./674;
tJ_25 = (2*(0:313) - 0)';
J_26  = (readcol('Oxygen_Data.xlsx','AKG+MAL','S2:S252')-3).*60./674;
eJ_26 = (readcol('Oxygen_Data.xlsx','AKG+MAL','T2:T252')).*60./674;
tJ_26 = (2*(0:250) + 1)';
N_25  = readcol('NADH_Data.xlsx','AKG+MAL','H2:H100');
eN_25 = readcol('NADH_Data.xlsx','AKG+MAL','I2:I100');
tN_25 = readcol('NADH_Data.xlsx','AKG+MAL','A2:A100');
N_26  = readcol('NADH_Data.xlsx','AKG+MAL','Q2:Q100');
eN_26 = readcol('NADH_Data.xlsx','AKG+MAL','R2:R100');
tN_26 = readcol('NADH_Data.xlsx','AKG+MAL','A2:A100');

% Exps. 35, 36, S + G (high GLU)
J_35  = (readcol('Oxygen_Data.xlsx','SUCC+GLU (high glu)','B2:B283')-3).*60./674;
tJ_35 = (2*(0:281) - 10)';
J_36  = (readcol('Oxygen_Data.xlsx','SUCC+GLU (high glu)','K2:K283')-3).*60./674;
tJ_36 = (2*(0:281) - 10)';

% Exps. 37, 38, S + G (low GLU)
J_37  = (readcol('Oxygen_Data.xlsx','SUCC+GLU (low glu)','B2:B283')-3).*60./674;
tJ_37 = (2*(0:281) - 10)';
J_38  = (readcol('Oxygen_Data.xlsx','SUCC+GLU (low glu)','K2:K283')-3).*60./674;
tJ_38 = (2*(0:281) - 10)';

% loading metabolite measurement data
tm    = readcol("Metabolite_Data.xlsx",'H_P+M','A9:A14');
succ3 = readcol("Metabolite_Data.xlsx",'L_P+M','F16:F21');
es3   = readcol("Metabolite_Data.xlsx",'L_P+M','F30:F35');
pyr3  = readcol("Metabolite_Data.xlsx",'L_P+M','B16:B21');
ep3   = readcol("Metabolite_Data.xlsx",'L_P+M','B30:B35');
mal3  = readcol("Metabolite_Data.xlsx",'L_P+M','C16:C21');
em3   = readcol("Metabolite_Data.xlsx",'L_P+M','C30:C35');
akg3  = readcol("Metabolite_Data.xlsx",'L_P+M','D16:D21');
ea3   = readcol("Metabolite_Data.xlsx",'L_P+M','D30:D35');
oaa3  = readcol("Metabolite_Data.xlsx",'L_P+M','E16:E21');
eo3   = readcol("Metabolite_Data.xlsx",'L_P+M','E30:E35');
amp3  = readcol("Metabolite_Data.xlsx",'L_P+M','G16:G21');
eamp3 = readcol("Metabolite_Data.xlsx",'L_P+M','G30:G35');
adp3  = readcol("Metabolite_Data.xlsx",'L_P+M','H16:H21');
eadp3 = readcol("Metabolite_Data.xlsx",'L_P+M','H30:H35');
atp3  = readcol("Metabolite_Data.xlsx",'L_P+M','I16:I21');
eatp3 = readcol("Metabolite_Data.xlsx",'L_P+M','I30:I35');
succ4 = readcol("Metabolite_Data.xlsx",'H_P+M','F16:F21');
es4   = readcol("Metabolite_Data.xlsx",'H_P+M','F30:F35');
pyr4  = readcol("Metabolite_Data.xlsx",'H_P+M','B16:B21');
ep4   = readcol("Metabolite_Data.xlsx",'H_P+M','B30:B35');
mal4  = readcol("Metabolite_Data.xlsx",'H_P+M','C16:C21');

em4   = readcol("Metabolite_Data.xlsx",'H_P+M','C30:C35');
akg4  = readcol("Metabolite_Data.xlsx",'H_P+M','D16:D21');
ea4   = readcol("Metabolite_Data.xlsx",'H_P+M','D30:D35');
oaa4  = readcol("Metabolite_Data.xlsx",'H_P+M','E16:E21');
eo4   = readcol("Metabolite_Data.xlsx",'H_P+M','E30:E35');
amp4  = readcol("Metabolite_Data.xlsx",'H_P+M','G16:G21');
eamp4 = readcol("Metabolite_Data.xlsx",'H_P+M','G30:G35');
adp4  = readcol("Metabolite_Data.xlsx",'H_P+M','H16:H21');
eadp4 = readcol("Metabolite_Data.xlsx",'H_P+M','H30:H35');
atp4  = readcol("Metabolite_Data.xlsx",'H_P+M','I16:I21');
eatp4 = readcol("Metabolite_Data.xlsx",'H_P+M','I30:I35');

%AKG + Mal - Low and High Pi - exps 25 & 26
tm1   = readcol("Metabolite_Data.xlsx",'H_P+M','A9:A13');
t25   = readcol("Metabolite_Data.xlsx",'H_A+M','A3:A7');
succ25 = readcol("Metabolite_Data.xlsx",'L_A+M','F3:F7');
pyr25  = readcol("Metabolite_Data.xlsx",'L_A+M','B3:B7');
mal25  = readcol("Metabolite_Data.xlsx",'L_A+M','C3:C7');
akg25  = readcol("Metabolite_Data.xlsx",'L_A+M','D3:D7');
oaa25  = readcol("Metabolite_Data.xlsx",'L_A+M','E3:E7');
amp25  = readcol("Metabolite_Data.xlsx",'L_A+M','G3:G7');
adp25  = readcol("Metabolite_Data.xlsx",'L_A+M','H3:H7');
atp25  = readcol("Metabolite_Data.xlsx",'L_A+M','I3:I7');
succ26 = readcol("Metabolite_Data.xlsx",'H_A+M','F3:F7');
pyr26  = readcol("Metabolite_Data.xlsx",'H_A+M','B3:B7');
mal26  = readcol("Metabolite_Data.xlsx",'H_A+M','C3:C7');
akg26  = readcol("Metabolite_Data.xlsx",'H_A+M','D3:D7');
oaa26  = readcol("Metabolite_Data.xlsx",'H_A+M','E3:E7');
amp26  = readcol("Metabolite_Data.xlsx",'H_A+M','G3:G7');
adp26  = readcol("Metabolite_Data.xlsx",'H_A+M','H3:H7');
atp26  = readcol("Metabolite_Data.xlsx",'H_A+M','I3:I7');

% loading substrate conditions
phosphate = readcol('Substrates.xlsx','Sheet1','C2:C39');
adp       = readcol('Substrates.xlsx','Sheet1','D2:D39');
pyr       = readcol('Substrates.xlsx','Sheet1','E2:E39');
mal       = readcol('Substrates.xlsx','Sheet1','F2:F39');
suc       = readcol('Substrates.xlsx','Sheet1','G2:G39');
glu       = readcol('Substrates.xlsx','Sheet1','H2:H39');
akg       = readcol('Substrates.xlsx','Sheet1','I2:I39');
cit       = readcol('Substrates.xlsx','Sheet1','J2:J39');
asp       = readcol('Substrates.xlsx','Sheet1','K2:K39');