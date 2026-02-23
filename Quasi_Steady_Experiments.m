%% Setting up simulation

% Set temperature
T = 37;
% Set proton cytoplasm sizes
BX(1) = 0.10;       % matrix
BX(2) = 100;     % cytoplasm
BX(3) = 100;     % intermembrane space DAB
% Proton cytoplasm binding constant(s)
K_BX(1) = 1e-7; % DAB
K_BX(2) = 10^(-7.2);
K_BX(3) = 10^(-7.2);
% % Setting physical parameters for computing the graphs
% W_m = 0.664*1.09;   % mitochondrial water space (ml water per ml mito) [VB2002]
% W_x = 0.7*W_m;      % Matrix water space, (l mito water) (l mito)^{-1}
% W_i = 0.3*W_m;      % IM water space, (l IM water) (l mito)^{-1}
% W_c = 1/((4*0.21e-3)/1.09);           % buffer, (l buffer water) (l mito)^{-1}
rho_m = 3.6697e-6;  % (l mito) (mg protein)^{-1}
VRegion_matrix = 0.001; % [=] l mito (l cuvette)^{-1} 

% defining total pools
NADtot = 3.0e-3; % NAD + NADH
NADPtot = 0.60e-3; % NADP + NADPH
Qtot = 2.0e-3; % coenzyme Q, corresponds to 10e-3 mol/(L lipid) with lambdaQ = 10; %%% CHANGED by DAB 8/27/2021
Atot = 10e-3; % ADP+ATP
Ctot = 0.632e-3; % cytochrome C (new value DAB 4/15/2021)

%% Setting initial conditions
%% Setting initial conditions
x0 = [];

% MATRIX VARIABLES
x0(1) = 1e-9;           % pyruvate, matrix
x0(2) = 2.0e-3 - 1e-9;  % coenzyme A
x0(3) = 1e-9;           % NAD
x0(4) =  2.14e-2;       % CO2 matrix
x0(5) =  1e-9;          % acetyl coA
x0(6) =  NADtot - 1e-9; % NADH
x0(7) =  1e-9;          % succinly coA
x0(8) =  1e-9;          % citrate, matrix
x0(9) =  1e-9;          % succinate, matrix
x0(10) = 1e-9;          % ATP, matrix
x0(11) = Atot - 1e-9;   % ADP, matrix
x0(12) = 1e-9;          % oxaloacetate_matrix
x0(13) = 0;             % AMP, matrix
x0(14) = 2.9916e-4;     % Pi, matrix
x0(15) = 1e-9;          % isocitrate
x0(16) = 1e-9;          % oxoglutarate
x0(17) = 30e-9;         % Ca matrix
x0(18) = 2e-3 - 1e-9;   % GDP, matrix
x0(19) = 1e-9;          % GTP, matrix
x0(20) = 1e-9;          % CoQ, oxidised
x0(21) = Qtot - 1e-9;   % CoQ, reduced
x0(22) = 1e-9;          % fumarate, matrix
x0(23) = 1e-9;          % malate, matrix
x0(24) = 1e-9;          % aspartate, matrix
x0(25) = 1e-9;          % glutamate, matrix
x0(26) = 1e-12;         % ammonia, matrix
x0(30) = 1e-12;         % mito H2O2
x0(31) = 1e-12;         % mito SO
x0(34) = 0.19e-3;       % O2, matrix     
x0(54) = 10^-(7.4);     % H, matrix
x0(55) = 1.2e-3;        % Mg, matrix
x0(56) = 0.070;         % K, matrix

% IM VARIABLES
x0(32) = Ctot - 1e-12;  % cytochrome c, ox
x0(33) = 1e-12;         % ctochrome c, red
x0(35) = 1e-12;         % ATP, IM
x0(36) =  1e-3;         % Pi, IM
x0(37) =  1e-12;        % ADP, IM
x0(38) =  1e-12;        % pyruvate, IM
x0(39) =  1e-12;        % glutamate, IM
x0(40) =  1e-12;        % aspartate, IM
x0(41) =  1e-12;        % citrate, IM
x0(42) =  1e-12;        % malate, IM
x0(43) =  1e-12;        % oxoglutarate, IM
x0(44) =  1e-12;        % succinate, IM
x0(45) =  0.19e-3;      % O2, IM
x0(60) = 10^-(7.2);     % H, IM
x0(61) = 1.0e-9;        % Mg, IM
x0(62) = 0.130;         % K, IM

% CYTOPLASM VARIABLES
x0(27) = 1e-12;         % ADP, cytoplasm
x0(28) = 1e-12;         % ATP, cytoplasm
x0(29) = 1e-12;         % Pi, cytoplasm
x0(46) =  1e-12;        % pyruvate, cytoplasm
x0(47) =  1e-12;        % citrate, cytoplasm
x0(48) =  1e-12;        % oxoglutarate, cytoplasm
x0(49) =  1e-12;        % succinate, cytoplasm
x0(50) =  1e-12;        % glutamate, cytoplasm
x0(51) =  1e-12;        % aspartate, cytoplasm
x0(52) =  1e-12;        % malate, cytoplasm
x0(53) =  0.19e-3;      % O2, cytoplasm
x0(57) = 10^-(7.2);     % H, buffer
x0(58) = 1.0e-9;        % Mg, buffer
x0(59) = 0.130;         % K, buffer

% OTHER VARIABLES
x0(63) = 0; % DPsi_im_to_matrix
x0(64) = 0; % DPsi_cytoplasm_to_im
x0(65) = 1.0; % initial PDH activity
x0(66) = NADPtot; % initial NADP_matrix
x0(67) = 0; % initial NADPH_matrix
x0(68) = 0; % initial AMP_cytoplasm
x0(69) = 0; % K+ leak activity 

%%

% Numerical tolerances
options1 = odeset('NonNegative', 1:62, 'abstol', 1e-7, 'reltol', 1e-8);
options2 = odeset('MaxStep', 1);

%% Run quasi-steady simulations to compare to data of Bazil et al. (2016)

% low Pi (1.0 mM) run
for i = 1:28
    i
  x_ATPase = 0;       % ATPase
  Pflag = 0;
  Kflag = 2;
  % State 1 simulation
  xsim0 = x0;
  xsim0(29) = 1.0e-3 ; % setting Pi, cytoplasm 
  [tsim1,xsim1] = ode15s(@dXdT, [-60 0], xsim0, options1,  T, BX, K_BX, Pflag, x_ATPase ); 
  
  % State 3 simulation
  xsim0 = xsim1(end,:);
  x_ATPase = (i-1)*1.0e-6;       % ATPase
  xsim0(46) = 2.5e-3 ; % setting pyruvate, cytoplasm 
  xsim0(52) = 0.5e-3 ; % setting malate, cytoplasm 
  xsim0(27) = 5.0e-3 ; % setting ATP_c 

  [tsim3,xsim3] = ode15s(@dXdT, [0 150], xsim0, options1,  T, BX, K_BX, Pflag, x_ATPase );    
   
  tsim_qs{i} = [tsim1; tsim3(2:end); ];
  xsim_qs{i} = [xsim1; xsim3(2:end,:); ];

  NAD  = xsim_qs{i}(end,3);
  NADH = max(1e-12,NADtot - NAD);
  NADPH = xsim_qs{i}(end,67)';
  Nr_qs1(i) = (NADH + NADPH)/(NADtot + NADPtot);
  cytC_qs1(i) = xsim_qs{i}(end,33)/Ctot;
  ADP_qs1(i) = xsim_qs{i}(end,28);
  dPsi_qs1(i) = xsim_qs{i}(end,63);

end

Vmito = 0.001; % [=] l mito (l cuvette)^{-1} 
clear J_ETC4_im_to_matrix
clear J_ROS
% computing fluxes
for i = 1:28
  x_ATPase = (i-1)*1.0e-6; 
  for ii = 1:length(tsim_qs{i})
    [fco,Jco] = dXdT(0,xsim_qs{i}(ii,:)', T,BX,K_BX,Pflag, x_ATPase );
    J_ETC4_im_to_matrix{i}(ii) = Jco(17)/2; % J_O2 in mol / sec / L mito
    J_ROS{i}(ii) = Jco(44);
  end
  % Accounting for electrode response time
  [t,J_electrode] = ode15s(@dXdT_electrode, tsim_qs{i}, 0, options2, tsim_qs{i}, J_ETC4_im_to_matrix{i}+0*J_ROS{i} );  
  J_o2_qs1(i) =  J_electrode(end)*60/674*(1e9)*Vmito;
end

% high Pi (5.0 mM) run
for i = 1:28
    i
  x_ATPase = 0;       % ATPase
  Pflag = 0;
  Kflag = 2;
  % State 1 simulation
  xsim0 = x0;
  xsim0(29) = 5.0e-3 ; % setting Pi, cytoplasm 
  [tsim1,xsim1] = ode15s(@dXdT, [-60 0], xsim0, options1,  T, BX, K_BX, Pflag, x_ATPase ); 
  
  % State 3 simulation
  xsim0 = xsim1(end,:);
  x_ATPase = (i-1)*1.0e-6;       % ATPase
  xsim0(46) = 2.5e-3 ; % setting pyruvate, cytoplasm 
  xsim0(52) = 0.5e-3 ; % setting malate, cytoplasm 
  xsim0(27) = 5.0e-3 ; % setting ATP_c 
  [tsim3,xsim3] = ode15s(@dXdT, [0 150], xsim0, options1,  T, BX, K_BX, Pflag, x_ATPase );  

  tsim_qs{i} = [tsim1; tsim3(2:end); ];
  xsim_qs{i} = [xsim1; xsim3(2:end,:); ];

  NAD  = xsim_qs{i}(end,3);
  NADH = max(1e-12,NADtot - NAD);
  NADPH = xsim_qs{i}(end,67)';
  Nr_qs5(i) = (NADH + NADPH)/(NADtot + NADPtot);
  cytC_qs5(i) = xsim_qs{i}(end,33)/Ctot;
  ADP_qs5(i) = xsim_qs{i}(end,28);
  dPsi_qs5(i) = xsim_qs{i}(end,63);

end

Vmito = 0.001; % [=] l mito (l cuvette)^{-1} 
clear J_ETC4_im_to_matrix
clear J_ROS
% computing fluxes
for i = 1:28
  x_ATPase = (i-1)*1.0e-6; 
  for ii = 1:length(tsim_qs{i})
    [fco,Jco] = dXdT(0,xsim_qs{i}(ii,:)', T,BX,K_BX,Pflag, x_ATPase );
    J_ETC4_im_to_matrix{i}(ii) = Jco(17)/2; % J_O2 in mol / sec / L mito
    J_ROS{i}(ii) = Jco(44);
  end
  % Accounting for electrode response time
  [t,J_electrode] = ode15s(@dXdT_electrode, tsim_qs{i}, 0, options2, tsim_qs{i}, J_ETC4_im_to_matrix{i}+0*J_ROS{i} );  
  J_o2_qs5(i) =  J_electrode(end)*60/674*(1e9)*Vmito;
end

%% Plot simulations of quasi-steady state experiment

JO2 = [0.0105    0.0519    0.0720    0.0869    0.1056;
       0.0145    0.0536    0.0771    0.1086    0.1444];
N   = [0.7284    0.4879    0.4241    0.4063    0.3964;
       0.6642    0.4821    0.3901    0.3688    0.3695];
D   = [0.0000    0.0561    0.1827    0.3589    0.6622;
       0.0000    0.0049    0.0578    0.1648    0.4716];
C   = [0.1862    0.2658    0.3147;
       0.1931    0.2576    0.2959];
Y   = [183.2317  167.5529  162.3817;
       183.4706  167.9643  156.2841];

figure(2); clf;

axes('position',[0.085 0.575 0.40 0.40]); hold on
plot(J_o2_qs1, Nr_qs1,'linewidth',2,'Color','blue'); 
plot(J_o2_qs5, Nr_qs5,'linewidth',2,'Color','red'); 
plot(JO2(1,:)*600,N(1,:),'bo','MarkerFaceColor',[1 1 1],'linewidth',2)
plot(JO2(2,:)*600,N(2,:),'ro','MarkerFaceColor',[1 1 1],'linewidth',2)
ylabel('NAD(P)H (normalized)','interpreter','latex')
axis([0 100 0 1]); box on
set(gca,'xticklabel',[]);

axes('position',[0.575 0.575 0.40 0.40]); hold on
plot(J_o2_qs1, dPsi_qs1,'linewidth',2,'Color','blue'); 
plot(J_o2_qs5, dPsi_qs5,'linewidth',2,'Color','red'); 
plot(JO2(1,[1 2 5])*600,Y(1,:),'bo','MarkerFaceColor',[1 1 1],'linewidth',2)
plot(JO2(2,[1 2 5])*600,Y(2,:),'ro','MarkerFaceColor',[1 1 1],'linewidth',2)
ylabel('$\Delta\Psi$ (mV)','interpreter','latex')
axis([0 100 150 200]); box on
set(gca,'xticklabel',[]);
legend('[Pi]_e = 1 mM', '[Pi]_e = 5 mM','Fontsize',10)

axes('position',[0.085 0.125 0.40 0.40]); hold on
plot(J_o2_qs1, cytC_qs1,'linewidth',2,'Color','blue'); 
plot(J_o2_qs5, cytC_qs5,'linewidth',2,'Color','red'); 
plot(JO2(1,[1 2 5])*600,C(1,:),'bo','MarkerFaceColor',[1 1 1],'linewidth',2)
plot(JO2(2,[1 2 5])*600,C(2,:),'ro','MarkerFaceColor',[1 1 1],'linewidth',2)
xlabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
ylabel('Cyt c$^{2+}$ (normalized)','interpreter','latex')
axis([0 100 0 0.4]); box on

axes('position',[0.575 0.125 0.40 0.40]); hold on
plot(J_o2_qs1, ADP_qs1*1e3,'linewidth',2,'Color','blue'); 
plot(J_o2_qs5, ADP_qs5*1e3,'linewidth',2,'Color','red'); 
plot(JO2(1,:)*600,D(1,:),'bo','MarkerFaceColor',[1 1 1],'linewidth',2)
plot(JO2(2,:)*600,D(2,:),'ro','MarkerFaceColor',[1 1 1],'linewidth',2)
xlabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
ylabel('[ADP]$_e$ (M)','interpreter','latex')
axis([0 100 0 0.7]); box on
