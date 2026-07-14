

%% Setting up simulation

% Set the desired temperature
T = 37;
% Set proton buffer sizes
BX(1) = 0.10;    % matrix
BX(2) = 100;     % buffer
BX(3) = 100;     % intermembrane space DAB
% Proton buffer binding constant(s)
K_BX(1) = 1e-7; % DAB
K_BX(2) = 10^(-7.2);
K_BX(3) = 10^(-7.2);
% Mito volume and water space
VRegion_matrix = 0.001; % [=] l mito (l cuvette)^{-1} 
VWater_matrix = 0.4705; % [=] l matrix water (l mito)^{-1}
VWater_im = 0.2533; % [=] l im water (l mito)^{-1}
% Numerical tolerances
options1 = odeset('NonNegative', 1:62, 'abstol', 1e-7, 'reltol', 1e-8);
options2 = odeset('MaxStep', 1);

% defining total pools
NADtot = 3.0e-3; % NAD + NADH
NADPtot = 0.60e-3; % NADP + NADPH
Qtot = 2.0e-3; % coenzyme Q, corresponds to 10e-3 mol/(L lipid) with lambdaQ = 2; 
Atot = 10e-3; % ADP+ATP
Ctot = 0.632e-3; % cytochrome C 

%% Setting initial conditions
x0 = [];

% MATRIX VARIABLES
x0(1) = 1e-9;           % pyruvate, matrix
x0(2) = 2.0e-3 - 1e-9;  % coenzyme A
x0(3) = NADtot - 1e-9;           % NAD
x0(4) =  2.14e-2;       % CO2 matrix
x0(5) =  1e-9;          % acetyl coA
x0(6) =  1e-9; % NADH
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
x0(20) = Qtot - 1e-9;    % CoQ, oxidised
x0(21) = 1e-9;          % CoQ, reduced
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
x0(56) = 0.050;         % K, matrix

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
x0(27) = 1e-12;         % ATP, buffer
x0(28) = 1e-12;         % ADP, buffer
x0(29) = 1e-12;         % Pi, buffer
x0(46) =  1e-12;        % pyruvate, buffer
x0(47) =  1e-12;        % citrate, buffer
x0(48) =  1e-12;        % oxoglutarate, buffer
x0(49) =  1e-12;        % succinate, buffer
x0(50) =  1e-12;        % glutamate, buffer
x0(51) =  1e-12;        % aspartate, buffer
x0(52) =  1e-12;        % malate, buffer
x0(53) =  0.19e-3;      % O2, buffer
x0(57) = 10^-(7.2);     % H, buffer
x0(58) = 1.0e-9;        % Mg, buffer
x0(59) = 0.130;         % K, buffer

% OTHER VARIABLES
x0(63) = 0; % DPsi_im_to_matrix
x0(64) = 0; % DPsi_buffer_to_im
x0(65) = 0*1.0; % initial PDH activity
x0(66) = NADPtot; % initial NADP_matrix
x0(67) = 0; % initial NADPH_matrix
x0(68) = 0; % initial AMP_buffer
x0(69) = 0; % K+ leak activity

%% Simulations


x_ATPase = 0.40e-6;

% initial simulation to anoxia
xsim0 = x0;
Pflag = 0;
xsim0(46) = 1.0e-3 ; % setting pyruvate, buffer 
xsim0(52) = 0.5e-3 ; % setting malate, buffer 
xsim0(28) = 2.0e-3 ; % setting ADP_c 
xsim0(29) = 5.0e-3 ; % setting Pi, buffer 
[tsim1,xsim1] = ode15s(@dXdT, [-200 900], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase ); 
% reoxygenate
Pflag = 1;
xsim0 = xsim1(end,:);
xsim0(34) =  0.19e-3;      % O2, buffer
xsim0(45) =  0.19e-3;      % O2, buffer
xsim0(53) =  0.19e-3;      % O2, buffer
[tsim2,xsim2] = ode15s(@dXdT, [900 1600], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase ); 

tsim = [tsim1; tsim2(2:end)];
xsim = [xsim1; xsim2(2:end,:)];

Vmito = 0.001; % [=] l mito (l cuvette)^{-1} 
clear J_ETC4_im_to_matrix
clear J_ROS
clear J_o2_el
% computing fluxes
for ii = 1:length(tsim)
  [fco,Jco] = dXdT(0,xsim(ii,:)', T,BX,K_BX,Pflag, x_ATPase );
  J_ETC4_im_to_matrix{1}(ii) = Jco(17)/2; % J_O2 in mol / sec / L mito
  J_ROS{1}(ii) = Jco(44);
end
% Accounting for electrode response time
[t,J_electrode] = ode15s(@dXdT_electrode, tsim, 0, options2, tsim, J_ETC4_im_to_matrix{1}+0*J_ROS{1} );  
J_o2_el =  J_electrode*60/674*(1e9)*Vmito;


%% Plot oxygen data 
    
figure(20); clf; 
    
axes('position',[0.125 0.540 0.40 0.40]); hold on
plot(tsim, J_o2_el, 'b-', 'linewidth', 2);  
ylabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',-200:50:100);
axis([-200 150 0 125]); box on

%%

SUCC =  xsim(:,49)*(1 - Vmito) + Vmito*VWater_matrix*xsim(:,9) + Vmito*VWater_im* xsim(:,44);
figure(22); clf
plot(tsim,SUCC*1e6,'k-','linewidth',2)
axis([-200 1600 0 200]); box on

