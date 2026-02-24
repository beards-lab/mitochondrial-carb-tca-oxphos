clear
Load_Data

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
x0(65) = 1.0; % initial PDH activity
x0(66) = NADPtot; % initial NADP_matrix
x0(67) = 0; % initial NADPH_matrix
x0(68) = 0; % initial AMP_buffer
x0(69) = 0; % K+ leak activity

%% Simulations

% indicies of simulation experiments to run
doi = [1 2 3 4 7 8 17 18 19 20 21 22 25 26 35 36 37 38];

% run once with Pflag = 0 to compute oxygen kinetics (closed system)
Pflag = 0;
for j = 1:length(doi)
  i = doi(j)

  % ATPase 
  if (i == 1) || (i == 2)
    x_ATPase = 0;
  else
    x_ATPase = 0.38e-6;
  end

  % State 1 simulation
  xsim0 = x0;
  xsim0(29) = phosphate(i)*1e-3 ; % setting Pi, buffer 
  [tsim1,xsim1] = ode15s(@dXdT, [0 130], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase ); 
  
  % State 2 simulation
  xsim0 = xsim1(end,:);
  xsim0(46) = pyr(i)*1e-3 ; % setting pyruvate, buffer 
  xsim0(52) = mal(i)*1e-3 ; % setting malate, buffer 
  xsim0(49) = suc(i)*1e-3 ; % setting succinate, buffer 
  xsim0(50) = glu(i)*1e-3 ; % setting glutamate, buffer 
  xsim0(48) = akg(i)*1e-3 ; % setting AKG, buffer 
  xsim0(47) = cit(i)*1e-3 ; % setting citrate, buffer 
  xsim0(51) = asp(i)*1e-3 ; % setting aspartate, buffer 
  [tsim2,xsim2] = ode15s(@dXdT, [130 280], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase );    
   
  % State 3 simulation
  xsim0 = xsim2(end,:);

  Tend = 650;
  xsim0(28) = adp(i)*1e-3 ; % setting ADP_c 
  [tsim3,xsim3] = ode15s(@dXdT, [280 Tend], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase ); 
  
  tsim_ox{i} = [tsim1; tsim2(2:end); tsim3(2:end)];
  xsim_ox{i} = [xsim1; xsim2(2:end,:); xsim3(2:end,:)];

end

Vmito = 0.001; % [=] l mito (l cuvette)^{-1} 
clear J_ETC4_im_to_matrix
clear J_ROS
% computing fluxes
for j = 1:length(doi)
  i = doi(j);
  for ii = 1:length(tsim_ox{i})
    [fco,Jco] = dXdT(0,xsim_ox{i}(ii,:)', T,BX,K_BX,Pflag, x_ATPase );
    J_ETC4_im_to_matrix{i}(ii) = Jco(17)/2; % J_O2 in mol / sec / L mito
    J_ROS{i}(ii) = Jco(44);
  end
  % Accounting for electrode response time
  [t,J_electrode] = ode15s(@dXdT_electrode, tsim_ox{i}, 0, options2, tsim_ox{i}, J_ETC4_im_to_matrix{i}+0*J_ROS{i} );  
  J_o2_el{i} =  J_electrode*60/674*(1e9)*Vmito;

  % NADH calculation
  tsim_N = tsim_ox;
  NAD  = xsim_ox{i}(:,3)';
  NADH = max(1e-12,NADtot - NAD);
  NADPH = xsim_ox{i}(:,67)';
  Nr_sim{i} = (NADH + NADPH)/(NADtot + NADPtot);
  %

end

% % run again with Pflag = 1 to compute NADH kinetics (open system)
% Pflag = 1;
% for j = 1:length(doi)
%   i = doi(j)
% 
%   % State 1 simulation
%   xsim0 = x00;
%   xsim0(29) = phosphate(i)*1e-3 ; % setting Pi, buffer 
%   [tsim1,xsim1] = ode15s(@dXdT, [0 130], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase ); 
%   
%   % State 2 simulation
%   xsim0 = xsim1(end,:);
%   xsim0(46) = pyr(i)*1e-3 ; % setting pyruvate, buffer 
%   xsim0(52) = mal(i)*1e-3 ; % setting malate, buffer 
%   xsim0(49) = suc(i)*1e-3 ; % setting succinate, buffer 
%   xsim0(50) = glu(i)*1e-3 ; % setting glutamate, buffer 
%   xsim0(43) = akg(i)*1e-3 ; % setting AKG, buffer 
%   xsim0(47) = cit(i)*1e-3 ; % setting citrate, buffer 
%   xsim0(51) = asp(i)*1e-3 ; % setting aspartate, buffer 
%   [tsim2,xsim2] = ode15s(@dXdT, [130 280], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase );    
%    
%   % State 3 simulation
%   xsim0 = xsim2(end,:);
% 
%   Tend = 650;
%   xsim0(27) = adp(i)*1e-3 ; % setting ADP_c 
%   [tsim3,xsim3] = ode15s(@dXdT, [280 Tend], xsim0, options1,  T, BX, K_BX,Pflag, x_ATPase ); 
%   
%   tsim_N{i} = [tsim1; tsim2(2:end); tsim3(2:end)];
%   xsim_N{i} = [xsim1; xsim2(2:end,:); xsim3(2:end,:)];
%   
%   NAD  = xsim_N{i}(:,3)';
%   NADH = max(1e-12,NADtot - NAD);
%   NADPH = xsim_N{i}(:,67)';
%   Nr_sim{i} = (NADH + NADPH)/(NADtot + NADPtot);
%   
% end

% toc

%% P+M (low pyr) - low and high Pi (experiments 1, 2) 
    
figure(1); clf; 
    
axes('position',[0.125 0.540 0.40 0.40]); hold on
plot(tN_01,N_01,'g-','linewidth',2);
plot(tsim_N{1}, Nr_sim{1}, 'k-', 'linewidth', 2)
ylabel('NAD(P)H (rel.)','interpreter','latex')
set(gca,'xticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on
text(225,1.075,'Low Pi (0.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.550 0.540 0.40 0.40]); hold on
plot(tN_02,N_02,'g-','linewidth',2);
plot(tsim_N{2}, Nr_sim{2}, 'k-', 'linewidth', 2)
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on   
text(225,1.075,'High Pi (2.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.125 0.090 0.40 0.40]); hold on
fill([tJ_01; flip(tJ_01)],[J_01+eJ_01; flip(J_01-eJ_01)],'g','Linewidth',1,'Edgecolor','g');
plot(tsim_ox{1}, J_o2_el{1}, 'k-', 'linewidth', 2); 
ylabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on

axes('position',[0.550 0.090 0.40 0.40]); hold on
fill([tJ_02; flip(tJ_02)],[J_02+eJ_02; flip(J_02-eJ_02)],'g','Linewidth',1,'Edgecolor','g');
plot(tsim_ox{2}, J_o2_el{2}, 'k-', 'linewidth', 2)
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on

%% P+M (med pyr) - low and high Pi (experiments 3, 4, 7, 8) 
   
figure(3); clf;
    
axes('position',[0.125 0.540 0.40 0.40]); hold on
g = fill([tN_03; flip(tN_03)],[N_03+eN_03; flip(N_03-eN_03)],'g','Linewidth',1,'Edgecolor','g');
plot(tsim_N{3}, Nr_sim{3}, 'k-', 'linewidth', 2)
ylabel('NAD(P)H (rel.)','interpreter','latex')
set(gca,'xticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on
text(225,1.075,'Low Pi (0.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.550 0.540 0.40 0.40]); hold on
fill([tN_04; flip(tN_04)],[N_04+eN_04; flip(N_04-eN_04)],'g','Linewidth',1,'Edgecolor','g');
plot(tsim_N{4}, Nr_sim{4}, 'k-', 'linewidth', 2)
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on   
text(225,1.075,'High Pi (2.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.125 0.090 0.40 0.40]); hold on
fill([tJ_03; flip(tJ_03)],[J_03+eJ_03; flip(J_03-eJ_03)],'g','Linewidth',1,'Edgecolor','g');
fill([tJ_07; flip(tJ_07)],[J_07+eJ_07; flip(J_07-eJ_07)],'y','Linewidth',1,'Edgecolor','y');
plot(tsim_ox{3}, J_o2_el{3}, 'k-', 'linewidth', 2); 
plot(tsim_ox{7}, J_o2_el{7}, 'r-', 'linewidth', 2)
ylabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on
legend('ADP = 0.375','ADP = 1 mM');

axes('position',[0.550 0.090 0.40 0.40]); hold on
fill([tJ_04; flip(tJ_04)],[J_04+eJ_04; flip(J_04-eJ_04)],'g','Linewidth',1,'Edgecolor','g');
fill([tJ_08; flip(tJ_08)],[J_08+eJ_08; flip(J_08-eJ_08)],'y','Linewidth',1,'Edgecolor','y');
plot(tsim_ox{4}, J_o2_el{4}, 'k-', 'linewidth', 2)
plot(tsim_ox{8}(1:length(J_o2_el{8})), J_o2_el{8}, 'r-', 'linewidth', 2)
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on

%% S - low and high Pi (experiments 17 and 18)
%  S + M - low and high Pi (experiments 19 and 20)

figure(17); clf

axes('position',[0.125 0.540 0.40 0.40]); hold on;
plot(tN_17,N_17,'g-','linewidth',2);
% plot(tN_19,N_19,'y-','linewidth',2);
plot(tsim_N{17}, Nr_sim{17}, 'k-', 'linewidth', 2)
% plot(tsim_N{19}, Nr_sim{19}, 'r-', 'linewidth', 2)
ylabel('NAD(P)H (rel.)','interpreter','latex')
set(gca,'xticklabel',[]);
set(gca,'xtick',0:100:600);
axis([100 650 0 1]); box on
text(225,1.075,'Low Pi (0.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.550 0.540 0.40 0.40]);  hold on;
plot(tN_18,N_18,'g-','linewidth',2);
% plot(tN_20,N_20,'y-','linewidth',2);
plot(tsim_N{18}, Nr_sim{18}, 'k-', 'linewidth', 2)
% plot(tsim_N{20}, Nr_sim{20}, 'r-', 'linewidth', 2)
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:600);
axis([100 650 0 1]); box on   
text(225,1.075,'High Pi (2.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.125 0.090 0.40 0.40]); hold on;
pp = plot(tJ_17,J_17,'linewidth',2); pp.Color = 'g'; 
% fill([tJ_17; flip(tJ_17)],[J_17+eJ_17; flip(J_17-eJ_17)],' g','linestyle','none');
fill([tJ_19; flip(tJ_19)],[J_19+eJ_19; flip(J_19-eJ_19)],' y','linestyle','none');
plot(tsim_ox{17}, J_o2_el{17}, 'k-', 'linewidth', 2)
plot(tsim_ox{19}, J_o2_el{19}, 'r-', 'linewidth', 2)
ylabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',0:100:600);
axis([100 650 0 100]); box on

axes('position',[0.550 0.090 0.40 0.40]); hold on;
pp = plot(tJ_18,J_18,'linewidth',2); pp.Color = 'g'; 
% fill([tJ_18; flip(tJ_18)],[J_18+eJ_18; flip(J_18-eJ_18)],'g','linestyle','none');
plot(tJ_20,J_20,'y','linewidth',4); 
plot(tsim_ox{18}, J_o2_el{18}, 'k-', 'linewidth', 2)
plot(tsim_ox{20}, J_o2_el{20}, 'r-', 'linewidth', 2)
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:600);
axis([100 650 0 100]); box on
legend('MAL = 0','MAL = 5 mM');

%% M only - low and high Pi (experiments 21, 22, 23, 24)

figure(21); clf;
    
axes('position',[0.125 0.540 0.40 0.40]); hold on
plot(tN_21,N_21,'g-','linewidth',2);
plot(tsim_N{21}, Nr_sim{21}, 'k-', 'linewidth', 2)
ylabel('NAD(P)H (rel.)','interpreter','latex')
set(gca,'xticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on
text(225,1.075,'Low Pi (0.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.550 0.540 0.40 0.40]); hold on
plot(tN_22,N_22,'g-','linewidth',2);
plot(tsim_N{22}, Nr_sim{22}, 'k-', 'linewidth', 2)
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on   
text(225,1.075,'High Pi (2.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.125 0.090 0.40 0.40]); hold on
fill([tJ_21; flip(tJ_21)],[J_21+eJ_21; flip(J_21-eJ_21)],'g','linewidth',1,'Edgecolor','g');
% plot(tJ_23,J_23,'y-','linewidth',2); 
plot(tsim_ox{21}, J_o2_el{21}, 'k-', 'linewidth', 2)
% plot(tsim_ox{23}, J_o2_el{23}, 'r-', 'linewidth', 2)
ylabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on

axes('position',[0.550 0.090 0.40 0.40]); hold on
fill([tJ_22; flip(tJ_22)],[J_22+eJ_22; flip(J_22-eJ_22)],'g','linewidth',1,'Edgecolor','g');
% plot(tJ_24,J_24,'y-','linewidth',2); 
plot(tsim_ox{22}, J_o2_el{22}, 'k-', 'linewidth', 2)
% plot(tsim_ox{24}, J_o2_el{24}, 'r-', 'linewidth', 2)
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on
% legend('MAL = 5.0','MAL = 1.0 mM','MAL = 0.25 mM');

%% A + M - low and high Pi (experiments 25 and 26)

figure(25); clf;
    
axes('position',[0.125 0.540 0.40 0.40]); hold on;
plot(tN_25,N_25,'linewidth',2,'Color','g'); 
plot(tsim_N{25}, Nr_sim{25}, 'k-', 'linewidth', 2)
ylabel('NAD(P)H (rel.)','interpreter','latex')
set(gca,'xticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on
text(225,1.075,'Low Pi (0.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.550 0.540 0.40 0.40]); hold on;

plot(tN_26,N_26,'linewidth',2,'Color','g'); 
plot(tsim_N{26}, Nr_sim{26}, 'k-', 'linewidth', 2)
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 1]); box on   
text(225,1.075,'High Pi (2.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.125 0.090 0.40 0.40]); hold on;
fill([tJ_25; flip(tJ_25)],[J_25+eJ_25; flip(J_25-eJ_25)],'g','Linewidth',1,'Edgecolor','g');
plot(tsim_ox{25}, J_o2_el{25}, 'k-', 'linewidth', 2)
ylabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on

axes('position',[0.550 0.090 0.40 0.40]); hold on
fill([tJ_26; flip(tJ_26)],[J_26+eJ_26; flip(J_26-eJ_26)],'g','Linewidth',1,'Edgecolor','g');
plot(tsim_ox{26}, J_o2_el{26}, 'k-', 'linewidth', 2)
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:500);
axis([100 600 0 100]); box on

%% S+G - low and high Pi (experiments 35-38)
% also plotting S experiments 17 and 18

figure(35); clf;
   
% axes('position',[0.125 0.575 0.40 0.40]); hold on;
% plot(tN_17,N_17,'g-','linewidth',2);
% plot(tsim_N{17}, Nr_sim{17}, 'k-', 'linewidth', 2)
% ylabel('NAD(P)H (rel.)','interpreter','latex')
% set(gca,'xticklabel',[]);
% set(gca,'xtick',0:100:600);
% axis([100 650 0 1]); box on
% 
% axes('position',[0.550 0.575 0.40 0.40]);  hold on;
% plot(tN_18,N_18,'g-','linewidth',2);
% plot(tsim_N{18}, Nr_sim{18}, 'k-', 'linewidth', 2)
% set(gca,'xticklabel',[]);
% set(gca,'yticklabel',[]);
% set(gca,'xtick',0:100:600);
% axis([100 650 0 1]); box on   

axes('position',[0.125 0.540 0.40 0.40]); hold on;
pp = plot(tJ_17,J_17,'linewidth',2); pp.Color = 'g'; 
pp = plot(tJ_35,J_35,'linewidth',2); pp.Color = 'cyan'; 
pp = plot(tJ_37,J_37,'linewidth',2); pp.Color = 'y'; 
plot(tsim_ox{17}, J_o2_el{17}, 'k-', 'linewidth', 2)
plot(tsim_ox{35}, J_o2_el{35}, 'k-', 'linewidth', 2)
plot(tsim_ox{37}, J_o2_el{37}, 'k-', 'linewidth', 2)
ylabel('$J_{o2}$ (nmol O$_2$ min$^{-1}$ UCS$^{-1}$)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',0:100:600);
axis([100 650 0 100]); box on
legend('GLU = 0','GLU = 1 mM','GLU = 5 mM');
text(225,107.5,'Low Pi (0.5 mM)','fontsize',12,'interpreter','latex')

axes('position',[0.550 0.540 0.40 0.40]); hold on;
pp = plot(tJ_18,J_18,'linewidth',2); pp.Color = 'g'; 
pp = plot(tJ_36,J_36,'linewidth',2); pp.Color = 'cyan'; 
pp = plot(tJ_38,J_38,'linewidth',2); pp.Color = 'y'; 

plot(tsim_ox{18}, J_o2_el{18}, 'k-', 'linewidth', 2)
plot(tsim_ox{36}, J_o2_el{36}, 'k-', 'linewidth', 2)
plot(tsim_ox{38}, J_o2_el{38}, 'k-', 'linewidth', 2)
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:600);
axis([100 650 0 100]); box on
text(225,107.5,'High Pi (2.5 mM)','fontsize',12,'interpreter','latex')

% Plots of metabolite measuremenst for exps. 3 & 4

% total metabolite levels, high and low Pi
ATP_lp =  xsim_ox{3}(:,27)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{3}(:,10) + Vmito*VWater_im* xsim_ox{3}(:,35);
ATP_hp =  xsim_ox{4}(:,27)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{4}(:,10) + Vmito*VWater_im* xsim_ox{4}(:,35);
ADP_lp =  xsim_ox{3}(:,28)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{3}(:,11) + Vmito*VWater_im* xsim_ox{3}(:,37);
ADP_hp =  xsim_ox{4}(:,28)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{4}(:,11) + Vmito*VWater_im* xsim_ox{4}(:,37);
AMP_lp =  xsim_ox{3}(:,68)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{3}(:,13) ;
AMP_hp =  xsim_ox{4}(:,68)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{4}(:,13) ;
PYR_lp =  xsim_ox{3}(:,46)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{3}(:,1) + Vmito*VWater_im* xsim_ox{3}(:,38);
PYR_hp =  xsim_ox{4}(:,46)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{4}(:,1) + Vmito*VWater_im* xsim_ox{4}(:,38);
MAL_lp =  xsim_ox{3}(:,52)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{3}(:,23) + Vmito*VWater_im* xsim_ox{3}(:,42);
MAL_hp =  xsim_ox{4}(:,52)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{4}(:,23) + Vmito*VWater_im* xsim_ox{4}(:,42);
AKG_lp =  xsim_ox{3}(:,48)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{3}(:,16) + Vmito*VWater_im* xsim_ox{3}(:,43);
AKG_hp =  xsim_ox{4}(:,48)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{4}(:,16) + Vmito*VWater_im* xsim_ox{4}(:,43);
SUCC_lp =  xsim_ox{3}(:,49)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{3}(:,9) + Vmito*VWater_im* xsim_ox{3}(:,44);
SUCC_hp =  xsim_ox{4}(:,49)*(1 - Vmito) + Vmito*VWater_matrix* xsim_ox{4}(:,9) + Vmito*VWater_im* xsim_ox{4}(:,44);

figure(40); clf;
   
axes('position',[0.125 0.575 0.40 0.40]); hold on;
plot(tsim_ox{3}, ADP_lp*1e6, 'k', 'linewidth', 2); 
plot(tsim_ox{3}, ATP_lp*1e6, 'b', 'linewidth', 2); 
plot(tsim_ox{3}, AMP_lp*1e6, 'r', 'linewidth', 2);
plot(tm, adp3, 'ko', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]); 
plot(tm, atp3, 'bo', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]);
plot(tm, amp3, 'ro', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]);
ylabel('Conc. ($\mu$M)','interpreter','latex')
set(gca,'xticklabel',[]);
set(gca,'xtick',0:100:600);
text(400,340,'ATP','FontSize',12);
text(400,140,'ADP','FontSize',12);
text(400,30,'AMP','FontSize',12);
axis([100 650 0 400]); box on

axes('position',[0.550 0.575 0.40 0.40]);  hold on;
plot(tsim_ox{4}, ADP_hp*1e6, 'k', 'linewidth', 2); 
plot(tsim_ox{4}, ATP_hp*1e6, 'b', 'linewidth', 2); 
plot(tsim_ox{4}, AMP_hp*1e6, 'r', 'linewidth', 2);
plot(tm, adp4, 'ko', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]); 
plot(tm, atp4, 'bo', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]);
plot(tm, amp4, 'ro', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:600);
% legend('ADP','ATP','AMP')
axis([100 650 0 400]); box on 

axes('position',[0.125 0.125 0.40 0.40]); hold on;
plot(tsim_ox{3}, PYR_lp*1e6 - 35, 'k', 'linewidth', 2); 
plot(tsim_ox{3}, MAL_lp*1e6 - 20, 'b', 'linewidth', 2); 
plot(tsim_ox{3}, AKG_lp*1e6, 'r', 'linewidth', 2); 
plot(tm, pyr3, 'ko', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]); 
plot(tm, mal3, 'bo', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]); 
plot(tm, akg3, 'ro', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]);
ylabel('Conc. ($\mu$M)','interpreter','latex')
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'xtick',0:100:600);
axis([100 650 0 500]); box on

axes('position',[0.550 0.125 0.40 0.40]); hold on;
plot(tsim_ox{4}, PYR_hp*1e6 - 60, 'k', 'linewidth', 2); 
plot(tsim_ox{3}, MAL_lp*1e6 - 20, 'b', 'linewidth', 2); 
plot(tsim_ox{3}, AKG_lp*1e6, 'r', 'linewidth', 2); 
plot(tm, pyr4, 'ko', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]); 
plot(tm, mal4, 'bo', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]); 
plot(tm, akg4, 'ro', 'markersize', 6, 'linewidth',2, 'MarkerFaceColor',[1 1 1]);
xlabel('$t$ (sec)','interpreter','latex')
set(gca,'yticklabel',[]);
set(gca,'xtick',0:100:600);
% legend('pyruvate','malate','ketoglut')
text(440,420,'pyruvate','FontSize',12);
text(440,220,'malate','FontSize',12);
text(440,40,'\alpha-keto.','FontSize',12);

axis([100 650 0 500]); box on


