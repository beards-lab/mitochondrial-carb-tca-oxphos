% function [f,J,] = dXdT(t,x,T,BX,K_BX,Kflag,Pflag,ATPase)
%
%
%
% Output parameters:
%   f     time derivatives of the model
%   J     flux
%
%
% Mandatory input parameters:
%   t     time
%   x     state variables at t=0
%   T     temperature in degreesCelcius
%   BX    buffer sizes
%   K_BX  proton buffer dissociation constants ( matrix buffer im )
%   Pflag = 1 for plate reader, = 0 for oxygraph
%   ATPase buffer ATP hydrolysis activity
%
% State Variables:
% [pyruvate_matrix , COAS_matrix , NAD_matrix , CO2tot_matrix , acetylcoA_matrix , NADH_matrix , succinylcoA_matrix , citrate_matrix , succinate_matrix , ATP_matrix , ADP_matrix , oxaloacetate_matrix , AMP_matrix , Pi_matrix , isocitrate_matrix , ketoglutarate_matrix , C_matrix , GDP_matrix , GTP_matrix , coQ_matrix , coQH2_matrix , fumarate_matrix , malate_matrix , aspartate_matrix , glutamate_matrix , ammonia_matrix , ATP_buffer , ADP_buffer , Pi_buffer , H2O2aq_matrix , SOaq_matrix , cytocox_im , cytocred_im , O2aq_matrix , ATP_im , Pi_im , ADP_im , pyruvate_im , glutamate_im , aspartate_im , citrate_im , malate_im , ketoglutarate_im , succinate_im , O2aq_im , pyruvate_buffer , citrate_buffer , ketoglutarate_buffer , succinate_buffer , glutamate_buffer , aspartate_buffer , malate_buffer , O2aq_buffer ]
%

function [f,J] = dXdT(~,x,T,BX,K_BX,Pflag,x_ATPase)
%% GLOBAL VARIABLES
% temperature 37
% ionic_strength 0.17

% Global experimental values

% minimal parameter value
MinCon = 0;

%% LIST OF STATE VARIABLES
% 1 pyruvate_matrix
% 2 COAS_matrix
% 3 NAD_matrix
% 4 CO2tot_matrix
% 5 acetylcoA_matrix
% 6 NADH_matrix
% 7 succinylcoA_matrix
% 8 citrate_matrix
% 9 succinate_matrix 
% 10 ATP_matrix
% 11 ADP_matrix
% 12 oxaloacetate_matrix
% 13 AMP_matrix
% 14 Pi_matrix
% 15 isocitrate_matrix
% 16 ketoglutarate_matrix
% 17 C_matrix
% 18 GDP_matrix
% 19 GTP_matrix
% 20 coQ_matrix (mol / L mito)
% 21 coQH2_matrix (mol / L mito)
% 22 fumarate_matrix
% 23 malate_matrix
% 24 aspartate_matrix
% 25 glutamate_matrix
% 26 ammonia_matrix
% 27 ATP_buffer
% 28 ADP_buffer
% 29 Pi_buffer
% 30 H2O2aq_matrix
% 31 SOaq_matrix
% 32 cytocox_im
% 33 cytocred_im
% 34 O2aq_matrix
% 35 ATP_im
% 36 Pi_im
% 37 ADP_im
% 38 pyruvate_im
% 39 glutamate_im
% 40 aspartate_im
% 41 citrate_im
% 42 malate_im
% 43 ketoglutarate_im
% 44 succinate_im
% 45 O2aq_im
% 46 pyruvate_buffer
% 47 citrate_buffer
% 48 ketoglutarate_buffer
% 49 succinate_buffer
% 50 glutamate_buffer
% 51 aspartate_buffer
% 52 malate_buffer
% 53 O2aq_buffer
% 54 h_matrix
% 55 m_matrix
% 56 k_matrix
% 57 h_buffer
% 58 m_buffer
% 59 k_buffer
% 60 h_im
% 61 m_im
% 62 k_im
% 63 DPsi_im_to_matrix
% 64 DPsi_buffer_to_im
% 65 PDH activity
% 66 NADP_matrix
% 67 NADPH_matrix
% 68 AMP_buffer
% 69 x_HLE

% PARTIAL VOLUME FRACTIONS
VWater_matrix = 0.4705; % [=] l matrix water (l mito)^{-1}
VRegion_matrix = 0.001; % [=] l mito (l cuvette)^{-1} 
VWater_buffer = 1; % [=] l water (l region)^{-1}
VRegion_buffer = 1 - VRegion_matrix; % [=] l buffer (l cuvette)^{-1} 
VWater_im = 0.2533; % [=] l im water (l mito)^{-1}
VRegion_im = 0.001; % [=] l mito (l cuvette)^{-1}

% % UQ partitioning 
lambdaQ = 5;    % partition coefficient for UQ (L matrix water / L lipid)

%% THERMODYNAMIC STATE and CONSTANTS
RT = 8.314*(T+273.15)/1e3; % kJ  mol^{-1}
F = 0.096484; % kJ mol^{-1} mV^{-1}

%% STATE VARIABLES

% concentration pools
NADtot = 3.0e-3; % NAD + NADH
Qtot = 2.0e-3; % coenzyme Q, corresponds to 20e-3 mol/(L lipid) with lambdaQ = 10;

% Concentrations of Reference Species
pyruvate_matrix = x(1);
COAS_matrix = x(2);
NAD_matrix = x(3);
CO2tot_matrix = x(4);
acetylcoA_matrix = x(5);
NADH_matrix = max(1e-12,NADtot - NAD_matrix);
succinylcoA_matrix = x(7);
citrate_matrix = x(8);
succinate_matrix = x(9);
ATP_matrix = x(10);
ADP_matrix = x(11);
oxaloacetate_matrix = x(12);
AMP_matrix = x(13);
Pi_matrix = x(14);
isocitrate_matrix = x(15);
ketoglutarate_matrix = x(16);
C_matrix = x(17);
GDP_matrix = x(18);
GTP_matrix = x(19);
% coQ_matrix = x(20);
coQH2_matrix = x(21);
coQ_matrix = max(1e-12,Qtot - coQH2_matrix);
fumarate_matrix = x(22);
malate_matrix = x(23);
aspartate_matrix = x(24);
glutamate_matrix = x(25);
ammonia_matrix = x(26);
ATP_buffer = x(27);
ADP_buffer = x(28);
Pi_buffer = x(29);
H2O2aq_matrix = x(30);
SOaq_matrix = x(31);
cytocox_im = x(32);
cytocred_im = x(33);
O2aq_matrix = max(x(34),1e-12);
ATP_im = x(35);
Pi_im = x(36);
ADP_im = x(37);
pyruvate_im = x(38);
glutamate_im = x(39);
aspartate_im = x(40);
citrate_im = x(41);
malate_im = x(42);
ketoglutarate_im = x(43);
succinate_im = x(44);
O2aq_im = x(45);
pyruvate_buffer = x(46);
citrate_buffer = x(47);
ketoglutarate_buffer = x(48);
succinate_buffer = x(49);
glutamate_buffer = x(50);
aspartate_buffer = x(51);
malate_buffer = x(52);
O2aq_buffer = x(53);
% Concentrations of H, Mg, and K
h_matrix = x(54);
m_matrix = x(55);
k_matrix = x(56);
h_buffer = x(57);
m_buffer = x(58);
k_buffer = x(59);
h_im = x(60);
m_im = x(61);
k_im = x(62);

% Membrane potentials
DPsi_im_to_matrix = x(63);
DPsi_buffer_to_im = x(64);

% Kinetic PDH regulation
PDH_activity = x(65);

% NADPH
NADP_matrix  = x(66);
NADPH_matrix = x(67);

% AMP for AK reaction
AMP_buffer = x(68);

% external ammonia
ammonia_im        = 0;

% KLEAK activity
HLE = x(69);

%% DISSOCIATION CONSTANTS
% pyruvate_matrix
Kh(1) = 0.004342603252078945;
Km(1) = 0.09553750450747911;
Kk(1) = 0.11969801449130479;
% COAS_matrix
Kh(2) = 7.5872945883630875e-09;
Km(2) = Inf;
Kk(2) = Inf;
% NAD_matrix
Kh(3) = Inf;
Km(3) = Inf;
Kk(3) = Inf;
% CO2tot_matrix
Kh(4) = 1.514167405186268e-10;
Km(4) = Inf;
Kk(4) = Inf;
% acetylcoA_matrix
Kh(5) = Inf;
Km(5) = Inf;
Kk(5) = Inf;
% NADH_matrix
Kh(6) = Inf;
Km(6) = Inf;
Kk(6) = Inf;
% succinylcoA_matrix
Kh(7) = 0.00010966977350385349;
Km(7) = Inf;
Kk(7) = Inf;
% citrate_matrix
Kh(8) = 2.7546665205956405e-06;
Km(8) = 0.00048047139882426425;
Kk(8) = 0.33133006036859936;
% succinate_matrix
Kh(9) = 6.4527518893172347e-06;
Km(9) = 0.063877481692394447;
Kk(9) = 0.42863973618362095;
% ATP_matrix
Kh(10) = 2.7566016511811682e-07;
Km(10) = 8.4303632255990519e-05;
Kk(10) = 0.09708512798373839;
% ADP_matrix
Kh(11) = 4.1057373286278706e-07;
Km(11) = 0.00071485288102034542;
Kk(11) = 0.1319048728526625;
% oxaloacetate_matrix
Kh(12) = 0.00015141674051862664;
Km(12) = 0.13814901290517506;
Kk(12) = Inf;
% AMP_matrix
Kh(13) = 6.2176978692908864e-07;
Km(13) = 0.019462203648493423;
Kk(13) = 0.23997936127262226;
% Pi_matrix
Kh(14) = 2.3075318323898138e-07;
Km(14) = 0.028149416570991275;
Kk(14) = 0.38034165601214515;
% isocitrate_matrix
Kh(15) = 2.2660048738054451e-06;
Km(15) = 0.0041259283254227218;
Kk(15) = Inf;
% ketoglutarate_matrix
Kh(16) = Inf;
Km(16) = Inf;
Kk(16) = Inf;
% C_matrix
Kh(17) = Inf;
Km(17) = Inf;
Kk(17) = Inf;
% GDP_matrix
Kh(18) = 5.014883398159419e-07;
Km(18) = 0.0028217704376049622;
Kk(18) = 0.038927887780503413;
% GTP_matrix
Kh(19) = 3.7183287636242726e-07;
Km(19) = 0.00031673464679257576;
Kk(19) = 0.019514066748113314;
% coQ_matrix
Kh(20) = Inf;
Km(20) = Inf;
Kk(20) = Inf;
% coQH2_matrix
Kh(21) = Inf;
Km(21) = Inf;
Kk(21) = Inf;
% fumarate_matrix
Kh(22) = 9.7762858862494379e-05;
Km(22) = Inf;
Kk(22) = Inf;
% malate_matrix
Kh(23) = 2.3069280947427119e-05;
Km(23) = 0.025828472176434858;
Kk(23) = 0.79464635030879283;
% aspartate_matrix
Kh(24) = 0.00023023155676210814;
Km(24) = Inf;
Kk(24) = Inf;
% glutamate_matrix
Kh(25) = 8.1530280646559256e-05;
Km(25) = Inf;
Kk(25) = Inf;
% ammonia_matrix
Kh(26) = 5.6234e-10;
Km(26) = Inf;
Kk(26) = Inf;
% ATP_buffer
Kh(27) = 2.7566016511811682e-07;
Km(27) = 8.4303632255990519e-05;
Kk(27) = 0.09708512798373839;
% ADP_buffer
Kh(28) = 4.1057373286278706e-07;
Km(28) = 0.00071485288102034542;
Kk(28) = 0.1319048728526625;
% Pi_buffer
Kh(29) = 2.3075318323898138e-07;
Km(29) = 0.028149416570991275;
Kk(29) = 0.38034165601214515;
% H2O2aq_matrix
Kh(30) = Inf;
Km(30) = Inf;
Kk(30) = Inf;
% SOaq_matrix
Kh(31) = Inf;
Km(31) = Inf;
Kk(31) = Inf;
% cytocox_im
Kh(32) = Inf;
Km(32) = Inf;
Kk(32) = Inf;
% cytocred_im
Kh(33) = Inf;
Km(33) = Inf;
Kk(33) = Inf;
% O2aq_matrix
Kh(34) = Inf;
Km(34) = Inf;
Kk(34) = Inf;
% ATP_im
Kh(35) = 2.7566016511811682e-07;
Km(35) = 8.4303632255990519e-05;
Kk(35) = 0.09708512798373839;
% Pi_im
Kh(36) = 2.3075318323898138e-07;
Km(36) = 0.028149416570991275;
Kk(36) = 0.38034165601214515;
% ADP_im
Kh(37) = 4.1057373286278706e-07;
Km(37) = 0.00071485288102034542;
Kk(37) = 0.1319048728526625;
% pyruvate_im
Kh(38) = 0.004342603252078945;
Km(38) = 0.09553750450747911;
Kk(38) = 0.11969801449130479;
% glutamate_im
Kh(39) = 8.1530280646559256e-05;
Km(39) = Inf;
Kk(39) = Inf;
% aspartate_im
Kh(40) = 0.00023023155676210814;
Km(40) = Inf;
Kk(40) = Inf;
% citrate_im
Kh(41) = 2.7546665205956405e-06;
Km(41) = 0.00048047139882426425;
Kk(41) = 0.33133006036859936;
% malate_im
Kh(42) = 2.3069280947427119e-05;
Km(42) = 0.025828472176434858;
Kk(42) = 0.79464635030879283;
% ketoglutarate_im
Kh(43) = Inf;
Km(43) = Inf;
Kk(43) = Inf;
% succinate_im
Kh(44) = 6.4527518893172347e-06;
Km(44) = 0.063877481692394447;
Kk(44) = 0.42863973618362095;
% O2aq_im
Kh(45) = Inf;
Km(45) = Inf;
Kk(45) = Inf;
% pyruvate_buffer
Kh(46) = 0.004342603252078945;
Km(46) = 0.09553750450747911;
Kk(46) = 0.11969801449130479;
% citrate_buffer
Kh(47) = 2.7546665205956405e-06;
Km(47) = 0.00048047139882426425;
Kk(47) = 0.33133006036859936;
% ketoglutarate_buffer
Kh(48) = Inf;
Km(48) = Inf;
Kk(48) = Inf;
% succinate_buffer
Kh(49) = 6.4527518893172347e-06;
Km(49) = 0.063877481692394447;
Kk(49) = 0.42863973618362095;
% glutamate_buffer
Kh(50) = 8.1530280646559256e-05;
Km(50) = Inf;
Kk(50) = Inf;
% aspartate_buffer
Kh(51) = 0.00023023155676210814;
Km(51) = Inf;
Kk(51) = Inf;
% malate_buffer
Kh(52) = 2.3069280947427119e-05;
Km(52) = 0.025828472176434858;
Kk(52) = 0.79464635030879283;
% O2aq_buffer
Kh(53) = Inf;
Km(53) = Inf;
Kk(53) = Inf;

%% BINDING POLYNOMIALS
P( 1 ) = 1  + h_matrix/Kh(1) + m_matrix/Km(1) + k_matrix/Kk(1);
P( 2 ) = 1  + h_matrix/Kh(2) + m_matrix/Km(2) + k_matrix/Kk(2);
P( 3 ) = 1  + h_matrix/Kh(3) + m_matrix/Km(3) + k_matrix/Kk(3);
P( 4 ) = 1  + h_matrix/Kh(4) + m_matrix/Km(4) + k_matrix/Kk(4);
P( 5 ) = 1  + h_matrix/Kh(5) + m_matrix/Km(5) + k_matrix/Kk(5);
P( 6 ) = 1  + h_matrix/Kh(6) + m_matrix/Km(6) + k_matrix/Kk(6);
P( 7 ) = 1  + h_matrix/Kh(7) + m_matrix/Km(7) + k_matrix/Kk(7);
P( 8 ) = 1  + h_matrix/Kh(8) + m_matrix/Km(8) + k_matrix/Kk(8);
P( 9 ) = 1  + h_matrix/Kh(9) + m_matrix/Km(9) + k_matrix/Kk(9);
P( 10 ) = 1  + h_matrix/Kh(10) + m_matrix/Km(10) + k_matrix/Kk(10);
P( 11 ) = 1  + h_matrix/Kh(11) + m_matrix/Km(11) + k_matrix/Kk(11);
P( 12 ) = 1  + h_matrix/Kh(12) + m_matrix/Km(12) + k_matrix/Kk(12);
P( 13 ) = 1  + h_matrix/Kh(13) + m_matrix/Km(13) + k_matrix/Kk(13);
P( 14 ) = 1  + h_matrix/Kh(14) + m_matrix/Km(14) + k_matrix/Kk(14);
P( 15 ) = 1  + h_matrix/Kh(15) + m_matrix/Km(15) + k_matrix/Kk(15);
P( 16 ) = 1  + h_matrix/Kh(16) + m_matrix/Km(16) + k_matrix/Kk(16);
P( 17 ) = 1  + h_matrix/Kh(17) + m_matrix/Km(17) + k_matrix/Kk(17);
P( 18 ) = 1  + h_matrix/Kh(18) + m_matrix/Km(18) + k_matrix/Kk(18);
P( 19 ) = 1  + h_matrix/Kh(19) + m_matrix/Km(19) + k_matrix/Kk(19);
P( 20 ) = 1  + h_matrix/Kh(20) + m_matrix/Km(20) + k_matrix/Kk(20);
P( 21 ) = 1  + h_matrix/Kh(21) + m_matrix/Km(21) + k_matrix/Kk(21);
P( 22 ) = 1  + h_matrix/Kh(22) + m_matrix/Km(22) + k_matrix/Kk(22);
P( 23 ) = 1  + h_matrix/Kh(23) + m_matrix/Km(23) + k_matrix/Kk(23);
P( 24 ) = 1  + h_matrix/Kh(24) + m_matrix/Km(24) + k_matrix/Kk(24);
P( 25 ) = 1  + h_matrix/Kh(25) + m_matrix/Km(25) + k_matrix/Kk(25);
P( 26 ) = 1  + h_matrix/Kh(26) + m_matrix/Km(26) + k_matrix/Kk(26);
P( 30 ) = 1  + h_matrix/Kh(30) + m_matrix/Km(30) + k_matrix/Kk(30);
P( 31 ) = 1  + h_matrix/Kh(31) + m_matrix/Km(31) + k_matrix/Kk(31);
P( 34 ) = 1  + h_matrix/Kh(34) + m_matrix/Km(34) + k_matrix/Kk(34);
P( 27 ) = 1  + h_buffer/Kh(27) + m_buffer/Km(27) + k_buffer/Kk(27);
P( 28 ) = 1  + h_buffer/Kh(28) + m_buffer/Km(28) + k_buffer/Kk(28);
P( 29 ) = 1  + h_buffer/Kh(29) + m_buffer/Km(29) + k_buffer/Kk(29);
P( 46 ) = 1  + h_buffer/Kh(46) + m_buffer/Km(46) + k_buffer/Kk(46);
P( 47 ) = 1  + h_buffer/Kh(47) + m_buffer/Km(47) + k_buffer/Kk(47);
P( 48 ) = 1  + h_buffer/Kh(48) + m_buffer/Km(48) + k_buffer/Kk(48);
P( 49 ) = 1  + h_buffer/Kh(49) + m_buffer/Km(49) + k_buffer/Kk(49);
P( 50 ) = 1  + h_buffer/Kh(50) + m_buffer/Km(50) + k_buffer/Kk(50);
P( 51 ) = 1  + h_buffer/Kh(51) + m_buffer/Km(51) + k_buffer/Kk(51);
P( 52 ) = 1  + h_buffer/Kh(52) + m_buffer/Km(52) + k_buffer/Kk(52);
P( 53 ) = 1  + h_buffer/Kh(53) + m_buffer/Km(53) + k_buffer/Kk(53);
P( 32 ) = 1  + h_im/Kh(32) + m_im/Km(32) + k_im/Kk(32);
P( 33 ) = 1  + h_im/Kh(33) + m_im/Km(33) + k_im/Kk(33);
P( 35 ) = 1  + h_im/Kh(35) + m_im/Km(35) + k_im/Kk(35);
P( 36 ) = 1  + h_im/Kh(36) + m_im/Km(36) + k_im/Kk(36);
P( 37 ) = 1  + h_im/Kh(37) + m_im/Km(37) + k_im/Kk(37);
P( 38 ) = 1  + h_im/Kh(38) + m_im/Km(38) + k_im/Kk(38);
P( 39 ) = 1  + h_im/Kh(39) + m_im/Km(39) + k_im/Kk(39);
P( 40 ) = 1  + h_im/Kh(40) + m_im/Km(40) + k_im/Kk(40);
P( 41 ) = 1  + h_im/Kh(41) + m_im/Km(41) + k_im/Kk(41);
P( 42 ) = 1  + h_im/Kh(42) + m_im/Km(42) + k_im/Kk(42);
P( 43 ) = 1  + h_im/Kh(43) + m_im/Km(43) + k_im/Kk(43);
P( 44 ) = 1  + h_im/Kh(44) + m_im/Km(44) + k_im/Kk(44);
P( 45 ) = 1  + h_im/Kh(45) + m_im/Km(45) + k_im/Kk(45);

%% THERMODYNAMIC EQUATIONS

% https://pubs.acs.org/doi/pdf/10.1021/jp911381p
DGro_PDH  = +17.50;
DGro_CTS  = +60.32;
DGro_ACON = -5.76;
DGro_IDH  = +93.52;
DGro_AKGDH = +15.28;
DGro_SCS = +56.56;
DGro_SDH = -3.10;
DGro_FUM = -3.38;
DGro_MDH = +71.09;
DGro_NDK = +0.01;
DGro_AAT = +2.8 ; % DAB from Equilibrator 
DGro_ETC4 =-202.524;
DGro_GDH = 134.87; % DAB this yields the same K' as from Equilibrator. 

Keq_PDH_matrix = exp(-DGro_PDH/RT)/P(1)/P(2)/P(3)*P(4)*P(5)*P(6)/h_matrix;
Keq_CTS_matrix = exp(-DGro_CTS/RT)*P(2)/P(5)*P(8)/P(12)/h_matrix^2;
Keq_ACON_matrix = exp(-DGro_ACON/RT)/P(8)*P(15);
Keq_IDH_matrix = exp(-DGro_IDH/RT)/P(3)*P(4)*P(6)/P(15)*P(16)/h_matrix^2;
Keq_AKGDH_matrix = exp(-DGro_AKGDH/RT)/P(2)/P(3)*P(4)*P(6)*P(7)/P(16)/h_matrix;
Keq_SCS_matrix = exp(-DGro_SCS/RT)*P(2)/P(7)*P(9)/P(14)/P(18)*P(19)/h_matrix;
Keq_MDH_matrix = exp(-DGro_MDH/RT)/P(3)*P(6)*P(12)/P(23)/h_matrix;
Keq_NDK_matrix = exp(-DGro_NDK/RT)*P(10)/P(11)*P(18)/P(19);
Keq_AAT_matrix = exp(-DGro_AAT/RT)*P(12)/P(16)/P(24)*P(25);
Keq_GDH_matrix = exp(-DGro_GDH/RT)/P(3)*P(6)*P(16)/P(25)*P(26)/h_matrix^2;

%% FLUX EQUATIONS

% % Equations to compute free NADH
Kn_NADH=0.18120e-3;
Xcp0_NADH=1.0e-3;
NADH_free =(NADH_matrix-Kn_NADH-Xcp0_NADH+sqrt((NADH_matrix-Kn_NADH-Xcp0_NADH)^2+4*Kn_NADH*NADH_matrix))/2;
NAD_free = NAD_matrix;

x_AK = 2;
P_AMPc = 1  + h_buffer/Kh(13) + m_buffer/Km(13) + k_buffer/Kk(13);
Keq_AK = (0.5633)*P_AMPc*P(27)/(P(28)^2);
J_AK_buffer = x_AK*( Keq_AK*ADP_buffer*ADP_buffer - AMP_buffer*ATP_buffer );

% Mitochondrial proton translocating transhydrogenase (MPTT)
% H_i + NADH + NADP = H_x + NAD + NADPH
K = 1.0e-3; % binding constant (DAB); arbitrary
alpha = 0.5;
DGro_MPTTdiv2 = 1.665; % 1/2 of ref free energy for  NADH + NADP = NAD + NADPH
Ke1 = exp(alpha*(F*DPsi_im_to_matrix - DGro_MPTTdiv2)/RT);
Ke2 = exp((alpha-1)*(F*DPsi_im_to_matrix + DGro_MPTTdiv2)/RT);
D = 1 + NADH_free/K + NAD_free/K + NADP_matrix/K + NADPH_matrix/K + ...
        NADH_free*NADP_matrix/(K^2) + NAD_free*NADPH_matrix/(K^2);
J_MPTT = 1.1487*(Ke1*NADH_free*NADP_matrix - NAD_free*NADPH_matrix*Ke2)/D;

% MALIC ENZYME
% NAD(P) + MAL^2- -> NAD(P)H^- + PYR^-
% General rapid by-uni model with p -> 0
Keq_MALIC = 1;
a = malate_matrix;
b = NAD_free;
p = pyruvate_matrix;
q = NADH_free;
KiNADPH = 0.4345e-3;
Ka = 0.16834e-3;  
Kb = (0.5340e-3)*(1 + NADPH_matrix/KiNADPH);  
Kq = (0.1947e-4); % Kmi for NADH
Den = 1 + a/Ka + b/Kb + a*b/(Ka*Kb) + q/Kq; % + o/KdOAA ;
J_MALIC = (5.8400e-05)/(Ka*Kb)*(a*b - p*q/Keq_MALIC)/Den; %/(1 + ATP_matrix/KiATP);

% oxaloacetate decarboxylase Reaction
% OAA^2- + H^+ -> PYR^- + CO2
a = oxaloacetate_matrix;
Ka = 16.535e-6; % apparent OD Km for OAA
D = 1 + a/Ka ;
J_OD = (10.0e-04)*a + (5.2721e-05)*(a/Ka)/D; % 

%PDH_matrix
a=pyruvate_matrix;
b=COAS_matrix;
c=NAD_free;
p=CO2tot_matrix;
q=acetylcoA_matrix;
r = NADH_free;

%KmandKivalues[M_matrix]
KmA=38.3e-6;
KmB=9.9e-6;
KmC=0.82439e-3; 
KiACCOA=40.2e-6;
KiNADH = 0.18608e-3; % 

Vmf_PDH=2.2926e-3;
Vmf_PDH=Vmf_PDH*PDH_activity;

% Kinetic PDH regulation
tau_PDH = 96.4;
Knadh = 0.43542e-3;
Katp  = 1.0718e-3;
alpha_aa = 1/(1 + NADH_matrix/Knadh )/tau_PDH;
beta_aa  = 1.4199/(1 + Knadh/NADH_matrix )/(1 + Katp/ATP_matrix)/tau_PDH;
f(65,:) = alpha_aa*(1-PDH_activity) - beta_aa*PDH_activity;

%Inhibition constants
ai1=1+acetylcoA_matrix/KiACCOA;
ai2=1+NADH_free/KiNADH;
if(a>MinCon)&&(b>MinCon)&&(c>MinCon)
  J_PDH_matrix=Vmf_PDH*(a*b*c-p*q*r/Keq_PDH_matrix)/(KmC*ai2*a*b+KmB*ai1*a*c+KmA*b*c+a*b*c);
else
  J_PDH_matrix=0;
end

%CTS_matrix
KmA = 5.7628e-6; 
KmB = 4.5728e-6; 
KmP = 0.1788e-6; 
KmQ = 2.24e-3;
KeA = 3.0898e-6;
KeB = 0.283;
KeQ = 0.44748e-3; 
KiATP = 39.898e-6;
KiADP = 141.79e-6;
KiAMP = 1.0533e-3;
KiSCOA  = 8.9669e-6;
KiH     = 5.5e-08;
KeP=Keq_CTS_matrix*KeA*KeB/KeQ;
KmP = KmP/(h_matrix/1e-7)^2;

a=oxaloacetate_matrix;
b=acetylcoA_matrix; % 
p=COAS_matrix;
q=citrate_matrix;
uATP=ATP_matrix/P(10);
uADP=ADP_matrix/P(11);
uAMP=AMP_matrix/P(13);

% inhibition factor
I2 = 1 + uATP/KiATP + uADP/KiADP + uAMP/KiAMP + succinylcoA_matrix/KiSCOA;

Vmf_CTS = 0.10354;
n = Vmf_CTS*(a*b - p*q/Keq_CTS_matrix)/(KeA*KmB);
d = (1 + KmA*b/(KeA*KmB) + KmQ*p/(KeQ*KmP)) + ...
    (a/KeA + KmQ*a*p/(KeA*KmP*KeQ) + KmA*p*q/(KeA^2*KmB*Keq_CTS_matrix))*I2 + ...
    ( (1/(KeA*KmB)-KmQ*Keq_CTS_matrix/(KeQ^2*KmP))*a*b + ...
      (1/(KmP*KeQ)-KmA/(KeA^2*Keq_CTS_matrix*KmB))*p*q + ...
      KmQ*a*b*p/(KeA*KmP*KeQ*KeB) + KmA*b*p*q/(KeA*KmB*KeQ*KeP) ) + ...
    (q/KeQ + KmQ*Keq_CTS_matrix*a*b/(KeQ^2*KmP) + KmA*b*q/(KeA*KeQ*KmB));
J_CTS_matrix = (n/d)/(1 + h_matrix/KiH);

%ACON_matrix
a=citrate_matrix;
p=isocitrate_matrix;
KmA=1161e-6;
KmP=434e-6;
Vmf_ACON=0.28387;
Vmr=Vmf_ACON*(KmP/KmA/Keq_ACON_matrix);
J_f=Vmf_ACON*Vmr*a/(KmA*Vmr+Vmr*a+Vmf_ACON/Keq_ACON_matrix*p);
J_ACON_matrix=J_f-Vmf_ACON*Vmr*(p/Keq_ACON_matrix)/(KmA*Vmr+Vmr*a+Vmf_ACON/Keq_ACON_matrix*p);

%IDH_matrix
a = NAD_free/P(3);
b=isocitrate_matrix/P(15);
p=ketoglutarate_matrix/P(16);
r=CO2tot_matrix/P(4);
Ca=C_matrix;
%freeCalciumconcentrationinMolar
q=NADH_free;
KmA=503.3e-6;
KmB=148.9e-6;
Kia=77.6e-6;
KaCaADP=13.3e-6;
KaCaATP=288.6e-6;
KaADP=61.3e-3;
Kiq=4.75;
KiH=0.11e-6;
KiMgEDTA=84.1e-6;
a1=0.0012;
a2=0.0097;
a3=0.0004;
a4=2.21;
n=3;
CaADP=Ca*ADP_matrix/10^(-2.86);
CaATP=Ca*ATP_matrix/10^(-3.86);
a2a=1+(CaADP/KaCaADP)+(CaATP/KaCaATP)+(ADP_matrix/KaADP);
a2b=1+(1/a1)*(CaADP/KaCaADP)+(1/a2)*(CaATP/KaCaATP)+(1/a3)*(ADP_matrix/KaADP);
%N1=1/(1+(1/a4)*(MgEDTA/KiMgEDTA));
%alpha1=a1a/a1b;
alpha2=a2a/a2b;
N1=1;
alpha1=1;
Vmf_IDH=2.2286e-3;
Vmf=Vmf_IDH/(1+h_matrix/KiH);
Jforward=(Vmf*N1)/(1+(KmA/a)*(1+(q/Kiq))*alpha1+alpha2*((KmB/b)^n+(Kia/a)*(KmB/b)^n*(1+q/Kiq)));
J_IDH_matrix=Jforward*(1-(1/Keq_IDH_matrix)*(p*q*r/(a*isocitrate_matrix)));

% DAB code of Qi et al. model
a=ketoglutarate_matrix;
b=COAS_matrix;
c = NAD_free;
p=CO2tot_matrix;
q=succinylcoA_matrix;
Ca=C_matrix;
%freemitoCa
r=NADH_free;

KmA   = 0.273e-3; % 
KmB   = 6.96e-6;
KmC   = 98.6e-6;
Kia   = 75.9e-3;
Kir   = 1.5e-6;  % low Ca value
Kic   = 0.112e-3;
Kiq   = 0.218e-3; 
KaH   = 10^(-6.11);
KaCa  = 0.893e-6;
KiATP = 0.106e-3;
KaADP = 0.305e-3;
KaMg  = 19.49e-6;
aCa   = 0.262;
aATP  = 6.694;
aADP  = 0.173;
aMg   = 1;
bMg   = 4.222;
Vmf_AKGDH=3.5245e-3;

den_alpha_A=(1+Ca/(aCa*KaCa)+(Ca/(aCa*KaCa))^2)*...
  (1+ATP_matrix/(aATP*KiATP))*(1+ADP_matrix/(aADP*KaADP))*...
  (1+m_matrix/(aMg*KaMg)+(m_matrix/(aMg*KaMg))^2);
alpha_A=(1+Ca/KaCa+(Ca/KaCa)^2)*...
  (1+ATP_matrix/KiATP)*(1+ADP_matrix/KaADP)*...
  (1+m_matrix/KaMg+(m_matrix/KaMg)^2)*(KaH/h_matrix)/den_alpha_A;
N= (1 + bMg*m_matrix/(aMg*KaMg) + bMg*(m_matrix/(aMg*KaMg))^2) /...
   (1 + m_matrix/(aMg*KaMg) + (m_matrix/(aMg*KaMg))^2 );
den= a*b*c + KmC*a*b*r/Kir + KmC*a*b + KmB*a*c*q/Kiq + ...
     KmB*a*c + KmB*Kic*a*q*r/(Kiq*Kir) + Kic*KmB*a*q/Kiq + ...
     KmA*b*c*alpha_A + Kia*KmC*b*r/Kir + KmB*Kic*Kia*q*r/(Kiq*Kir) ;
J_for = Vmf_AKGDH*N*a*b*c/den;
if(a>MinCon)&&(b>MinCon)&&(c>MinCon)
  J_AKGDH_matrix = J_for*(1 - (p*q*r)/(a*b*c*Keq_AKGDH_matrix) );
elseif(p>MinCon)&&(q>MinCon)&&(r>MinCon)
  J_AKGDH_matrix = - Vmf_AKGDH*N*(p*q*r*h_matrix)/Keq_AKGDH_matrix/den;
else
  J_AKGDH_matrix = 0;
end


%SCS_matrix--Li et al. model, coded by DAB
I=0.17;
T = 37+273.15;
alphaT=1.10708 - (1.54508e-3)*T + (5.95584e-6)*(T^2);
r1=exp(-alphaT*(1)^2*sqrt(I)/(1+1.6*sqrt(I)));
r2=exp(-alphaT*(2)^2*sqrt(I)/(1+1.6*sqrt(I)));
r3=exp(-alphaT*(3)^2*sqrt(I)/(1+1.6*sqrt(I)));
r4=exp(-alphaT*(4)^2*sqrt(I)/(1+1.6*sqrt(I)));
r5=exp(-alphaT*(5)^2*sqrt(I)/(1+1.6*sqrt(I)));

a = r3*GDP_matrix/P(18);  % GDP3-
b = r1*succinylcoA_matrix/P(7); % SCOA-
c = r2*Pi_matrix/P(14); % PI2-
p = r1*COAS_matrix/P(2); % COAS-
q = r2*succinate_matrix/P(9); % SUC2-
r = r4*GTP_matrix/P(19); % GTP4-

k1  = 8.9553e+00;
km1 = 7.9948e-07;
k2  = 1.0826e-02;
km2 = 4.1850e-07; 
k3  = 1.0332e-03;
km3 = 2.1721e-04;
k4  = 4.2521e-03;
km4 = 1.0107e+02;
k5  = 6.7814e-05;
km5 = 3.8551e-02;
k6  = 3.2078e-06;
km6 = 3.8124e+01;
KIQ = 1.4840e-04;

KmA =(k4*k5*k6*r4^2*r5)/(k1*(k5*k6*r4^2 + k4*k5*r4*r5 + k4*k6*r4*r5));
KmB =(k4*k5*k6*r4^2*r5)/(k2*r2*(k5*k6*r4^2 + k4*k5*r4*r5 + k4*k6*r4*r5));
KmC =(k5*k6*r4^2*(k4*r5 + km3*r5))/(k3*r3*(k5*k6*r4^2 ...
    + k4*k5*r4*r5 + k4*k6*r4*r5));
KmP =(km1*km2*r2*r3*(k4*r5 + km3*r5))/(km4*r4*(km1*km2*r2*r3 ...
    + km1*km3*r2*r5 + km2*km3*r3*r5));
KmQ =(km1*km2*km3*r2*r3*r5)/(km5*r4*(km1*km2*r2*r3 ...
    + km1*km3*r2*r5 + km2*km3*r3*r5));
KmR =(km1*km2*km3*r2*r3*r5)/(km6*(km1*km2*r2*r3 + km1*km3*r2*r5 ...
    + km2*km3*r3*r5));
Kia =(km1*r2)/k1;
Kib =(km2*r3)/(k2*r2);
Kic =(km3*r5)/(k3*r3);
Kip =(k4*r5)/(km4*r4);
Kiq =k5/km5;
Kir =(k6*r4)/km6;
Vmf=15.194e-3;
Vmr=Vmf/Keq_SCS_matrix*KmP*Kiq*Kir/(Kia*Kib*KmC);
J_SCS_matrix = (Vmf*Vmr*a*b*c - Vmf*Vmr*(p*q*r/Keq_SCS_matrix)) ...
  /(Vmr*Kia*Kib*KmC*(1+a*q/Kia/KIQ)+Vmr*Kib*KmC*a+Vmr*Kia*KmB*c*(1+a*q/Kia/KIQ) ...
    + Vmr*KmC*a*b+Vmr*KmB*a*c+Vmr*KmA*b*c*(1+a*q/Kia/KIQ)+Vmr*a*b*c ...
    + Vmf*Kir*KmQ*p*(1+a*q/Kia/KIQ)/Keq_SCS_matrix+Vmf*Kiq*KmP*r/Keq_SCS_matrix+Vmf*KmR*p*q*(1+a*q/Kia/KIQ)/Keq_SCS_matrix+Vmf*KmQ*p*r/Keq_SCS_matrix...
    + Vmf*KmP*q*r/Keq_SCS_matrix+Vmf*p*q*r/Keq_SCS_matrix+Vmf*KmQ*Kir*a*p/Kia/Keq_SCS_matrix+Vmr*Kia*KmB*c*r/Kir...
    + Vmf*KmQ*Kir*a*b*p/Kia/Kib/Keq_SCS_matrix+Vmr*KmA*b*c*r/Kir+Vmf*KmR*a*p*q/Kia/Keq_SCS_matrix...
    + Vmr*Kia*KmB*c*q*r/Kiq/Kir+Vmf*Kir*KmQ*a*b*c*p/Kia/Kib/Kic/Keq_SCS_matrix+Vmf*Kip*KmR*a*b*c*q/Kia/Kib/Kic/Keq_SCS_matrix...
    + Vmf*KmR*a*b*p*q/Kia/Kib/Keq_SCS_matrix+Vmr*KmA*b*c*q*r/Kiq/Kir+Vmr*KmA*Kic*b*p*q*r/Kip/Kiq/Kir...
    + Vmr*Kia*KmB*c*p*q*r/Kip/Kiq/Kir+Vmf*KmR*a*b*c*p*q/Kia/Kib/Kic/Keq_SCS_matrix+Vmr*KmA*b*c*p*q*r/Kip/Kiq/Kir);


%% Bazil SDH_matrix code

% Substrate concencs.
O2 = O2aq_matrix; % oxygen     
SO = 0; % superoxide (zero for now)           
H2O2 = 0; % hydrogen peroxide (zero for now)  
Q = coQ_matrix*lambdaQ; %Q10
QH2 = coQH2_matrix*lambdaQ; % Q10H2
% Q = coQ_matrix; %Q10
% QH2 = coQH2_matrix; % Q10H2
H = h_matrix; % proton
SUC = succinate_matrix; % succinate 
FUM = fumarate_matrix; % fumarate
MALO = 0; % malonate
% MALO = 5e-3; % malonate
OAA = oxaloacetate_matrix; % oxaloacetate

MAL = malate_matrix; % malate
GEA = 0; % general electron acceptor (for TMPD and phenazine)
GED = 0; % no data set has any non-zero general electron donor concentrations

KD_Q    = 2.9047e-10;
KD_QH2  = 1.8800e-10;
KD_SUC  = 3.5507e-04;
KD_FUM  = 1.0000e-03; % 
KD_MAL  = 2.9390e-04;
KD_MALO = 1.4786e-05;
KD_OXA  = 8.2283e-07; %

kf_QH2  = 2.4792e+07;
kf_SUC  = 52.9498;
kf_FADH = 8.1515e+05;
kf_FADH2 = 0;
kf_SO2   = 6.6827e+09;
kf_H2O2  = 2.6096e+03;
kf_SO    = 0;
kf_GEA   = 0;
atpenin      = 0;
beta_atpenin = 14.3598;
KD_atpenin   = 1.6747e-17;
KD_atpenin2  = 6.8349e-08;
KD_H         = 1.6967e-07;
K_FADH  = 1.2882e-08;
K_FADH2 = 1.8621e-08;
Em0_FAD_FADH   = 283;
Em0_FADH_FADH2 = 398;
Em0_FAD_FADH2  = 340.5;
Em0_Q_SQ       = 283.7605;
Em0_Q_QH2      = 464.1;
Em0_Fum_Suc    = 445.0295;
Em_O2_SO       = -160;
Em_O2_H2O2     = 940;
Em0_GEA        = 0;
Em_ISC1        = 0;
Em_ISC2        = -260;
Em_ISC3        = 60;
logHmatrix = log(H);
RTdF = RT/F; % RT/F
FdRT = F/RT; % F/RT

% Binding Polynomials for enzyme, substrates, products and regulators
BP_Q = (1 + Q/KD_Q + QH2/KD_QH2 + atpenin/KD_atpenin);
kf_QH2 = kf_QH2*(1 + atpenin/KD_atpenin2)/(1 + beta_atpenin*atpenin/KD_atpenin2);
kf_SUC = kf_SUC*(1 + atpenin/KD_atpenin2)/(1 + beta_atpenin*atpenin/KD_atpenin2)/(1+H/KD_H);
BP_FAD = (1 + SUC/KD_SUC + FUM/KD_FUM + MAL/KD_MAL + MALO/KD_MALO  + OAA/KD_OXA);

% FMN Occupancy for SO and H2O2 reactions
kf_FADH = kf_FADH/BP_FAD/(1+H/K_FADH);
kf_FADH2 = kf_FADH2/BP_FAD/(1+H/K_FADH2);
kf_H2O2 = kf_H2O2/BP_FAD/(1+H/K_FADH2);

% Adjust redox midpoint potentials
% FAD potential as a pH dependence
Em_FAD_FADH =  Em0_FAD_FADH + (RTdF)*log(H*(1+K_FADH/H)); 
Em_FADH_FADH2 = Em0_FADH_FADH2 + (RTdF) *log(H*(1+H/K_FADH2)/(1+H/K_FADH));
Em_FAD_FADH2 = Em0_FAD_FADH2 + (RTdF/2)*log(H^2*(1+H/K_FADH2)*(1+K_FADH/H)/(1+H/K_FADH));
Em_Q_SQ   =   Em0_Q_SQ; % pH independent
Em_SQ_QH2 =   2*Em0_Q_QH2 + 2*(RTdF)*logHmatrix - Em_Q_SQ; % Em of SQ/QH2 couple at given pH

% Em_Q_QH2 = Em0_Q_QH2 + (RTdF)*logHmatrix;
Em_Fum_Suc = Em0_Fum_Suc + (RTdF)*logHmatrix;

% general electron acceptor (phenazine or TMPD in studies used)
Em_GEA = Em0_GEA  + (RTdF/2) * logHmatrix;

% Midpoint potentials for electron transfer of bound states
SQ_QH2 = Em_SQ_QH2 - (RTdF)*log(KD_QH2); % electron transfer free energy (mV)
Q_SQ = Em_Q_SQ + (RTdF)*log(KD_Q); % electron transfer free energy (mV)
Fum_Suc = Em_Fum_Suc - (RTdF/2)*log(KD_SUC) + (RTdF/2)*log(KD_FUM); % electron transfer free energy (mV)

% Equilibrium constants, Keq
Keq_SO_ISC3            = exp(FdRT*(Em_O2_SO - Em_ISC3));
Keq_SO_SQ              = exp(FdRT*(Em_O2_SO - Q_SQ));
Keq_SO_FADH            = exp(FdRT*(Em_O2_SO - Em_FAD_FADH));
Keq_SO_FADH2           = exp(FdRT*(Em_O2_SO - Em_FADH_FADH2));
Keq_H2O2_FADH2         = exp(2*FdRT*(Em_O2_H2O2 - Em_FAD_FADH2));
Keq_QH2_ISC3           = exp(FdRT*(SQ_QH2 - Em_ISC3)); % should be bound state Em for SQ_QH2
Keq_FUM_FADH2          = exp(2*FdRT*(Em_FAD_FADH2 - Fum_Suc)); % should be Em.FAD - Emb.Suc
Keq_GEA_ISC            = exp(2*FdRT*(Em_GEA - (Em_ISC3 + Em_ISC1))); % since Em of ISC2 so low, using ISC1 as second e- source

% Boltzman redox poise potentials
% compute Eh for Q/SQ reaction
Eh_FAD_FADH = Em_FAD_FADH;
Eh_FADH_FADH2 = Em_FADH_FADH2;
Eh_ISC1 = Em_ISC1;
Eh_ISC2 = Em_ISC2;
Eh_ISC3 = Em_ISC3;    
Eh_Q_SQ = Q_SQ + (RTdF)*log(Q/KD_Q/BP_Q);

% RT normalized DeltaG's of all substrates
% FAD_FADH, FADH_FADH2,FAD_FADH2,ISC1,ISC2,ISC3,Q_SQ;
% Eh.FAD_FADH2
DfG1 = -FdRT* [Eh_FAD_FADH;...
    Eh_FADH_FADH2;...
    Eh_ISC1;...
    Eh_ISC2;...
    Eh_ISC3;...
    Eh_Q_SQ];

% The fractional substrates of redox centers for each elecrtonic states
% E0
s0 = 1;
% E1
substrates_s1 = [1 3 4 5 6];
E1i = DfG1(substrates_s1);
s1 = exp(-E1i)./sum(exp(-E1i));
% E2
substrates_s2 = [1  1  1  1  1  3  3  3  4  4  5;
                 2  3  4  5  6  4  5  6  5  6  6]';
E2i = sum( DfG1(substrates_s2),2 );  
s2 = exp(-E2i)./sum(exp(-E2i));
% E3
substrates_s3 = [1  1  1  1  1  1  1  1  1  1  3  3  3  4;
                 2  2  2  2  3  3  3  4  4  5  4  4  5  5;
                 3  4  5  6  4  5  6  5  6  6  5  6  6  6]';
E3i = sum( DfG1(substrates_s3),2 );  
s3 = exp(-E3i)./sum(exp(-E3i));
% E4
substrates_s4 = [1  1  1  1  1  1  1  1  1  1  3;
                 2  2  2  2  2  2  3  3  3  4  4;
                 3  3  3  4  4  5  4  4  5  5  5;
                 4  5  6  5  6  6  5  6  6  6  6]';
E4i = sum( DfG1(substrates_s4),2 );  
s4 = exp(-E4i)./sum(exp(-E4i));

% State Transition Rates
% commented lines are not used in current model version
s1_FADH = s1(1);
s1_FAD  = sum(s1(2:5));
s1_SQ   = s1(5);
s1_SQempty  = sum(s1(1:4));
s1_SQempty_ISC3ox  = sum(s1(1:3));
% s1_ISC1ox_ISC3ox = sum(s1(p.s1.ISC1ox_ISC3ox,:),1);
s1_ISC3     = s1(4);
% s1_ISC3ox   = sum(s1(p.s1.ISC3ox,:),1);

s2_FAD     = sum(s2(6:11));
s2_FADH    = sum(s2(2:5));
s2_FADH2   = s2(1);
s2_SQ      = sum(s2([5 8 10 11]));
s2_SQempty = sum(s2([1 2 3 4 6 7 9]));
s2_ISC3    = sum(s2([4 7 9]));
% s2_ISC3ox   = sum(s2(p.s2.ISC3ox,:),1);
% s2_ISC1_ISC3 = sum(s2(p.s2.ISC1_ISC3,:),1); 
s2_SQ_ISC3 = s2(11);
s2_SQempty_ISC3ox  = sum(s2([1 2 3 6]));
% s2_ISC1ox_ISC3ox = sum(s2(p.s2.ISC1ox_ISC3ox,:),1);

s3_FAD      = sum(s3([11 12 13 14]));
s3_FADH     = sum(s3(5:10));
s3_FADH2    = sum(s3(1:4));
s3_SQ       = sum(s3([4 7 9 10 12 13 14]));
s3_SQempty  = sum(s3([1 2 3 5 6 8 11]));
s3_ISC3     = sum(s3([3 6 8 11]));
% s3_ISC3ox   = sum(s3(p.s3.ISC3ox,:),1);
% s3_ISC1_ISC3 = sum(s3(p.s3.ISC1_ISC3,:),1);
s3_SQ_ISC3  = sum(s3([10 13 14]));
% s3_ISC1ox_ISC3ox = sum(s3(p.s3.ISC1ox_ISC3ox,:),1);
% s3_SQempty_ISC3ox = sum(s3(p.s2.SQempty_ISC3ox,:),1);

%s4_FAD    = sum(s4(p.s4.FADox,:),1);
s4_FADH     = sum(s4([7 8 9 10]));
s4_FADH2    = sum(s4(1:6));
s4_SQ       = sum(s4([3 5 6 8 9 10 11]));
% s4_SQempty  = sum(s4(p.s4.SQempty,:),1);
s4_ISC3     = sum(s4([2 4 7]));
% s4_ISC3ox   = sum(s4(p.s4.ISC3ox,:),1);
% s4_ISC1_ISC3 = sum(s4(p.s4.ISC1_ISC3,:),1);
s4_SQ_ISC3  = sum(s4([6 9 10 11]));
% s4_ISC1ox_ISC3ox = sum(s4(p.s4.ISC1ox_ISC3ox,:),1);
% s4_SQempty_ISC3ox = sum(s4(p.s4.SQempty_ISC3ox,:),1);

% Unidirectional state transition rates
% 1) E0->E1
% SO + Q -> O2 + SQ
% FAD + H + SO -> FADH + O2

k01_SO_ISC3 = (kf_SO2/Keq_SO_ISC3)* 1/BP_Q * SO * s0;% / (1 + GEA/p.KD_GEA);
k01_SO_Q = (kf_SO/Keq_SO_SQ)* Q/KD_Q/BP_Q * SO * s0;
k01_SO_FAD = (kf_FADH/Keq_SO_FADH)* SO * s0;
k0_1 =   k01_SO_Q + k01_SO_FAD + k01_SO_ISC3;

% 2) E1->E0
% O2 + SQ -> SO + Q
% FADH + O2 -> FAD + H + SO

k10_O2_ISC3 = kf_SO2 * O2 * s1_ISC3 / BP_Q;% / (1 + GEA/p.KD_GEA);
k10_O2_SQ = kf_SO * O2 * s1_SQ;
k10_O2_FADH = kf_FADH * O2 * s1_FADH;
k1_0 = k10_O2_SQ  + k10_O2_FADH + k10_O2_ISC3;

% 3) E0->E2                       
% H2O2 + FAD -> O2 + FADH2
% FAD + Succinate -> FADH2 + fumarate
% QH2 + [4Fe_3S]ox -> SQ + [4Fe_3S]red + 2H
% GED + [4Fe_3S]ox + [3Fe_4S]ox -> GEA  + [4Fe_3S]red + [3Fe_4S]red

k02_H2O2_FAD = (kf_H2O2/Keq_H2O2_FADH2) * H2O2 * s0;
k02_SUC_FAD =  kf_SUC* SUC/KD_SUC/BP_FAD * s0;
k02_QH2_SQempty_ISC3ox =  (kf_QH2/ Keq_QH2_ISC3)* QH2/KD_QH2/BP_Q * s0;
k02_GED = (kf_GEA/Keq_GEA_ISC)*GED;%/(GED+p.KD_GEA);%*s0;
k0_2 = k02_H2O2_FAD + k02_SUC_FAD + k02_QH2_SQempty_ISC3ox +  k02_GED; 

% 4) E2->E0
% O2 + FADH2 ->  H2O2 + FAD
% FADH2 + fumarate -> FAD + Succinate
% SQ-[4Fe_3S]red + 2H -> QH2 + [4Fe_3S]ox
% GEA + [4Fe_3S]red + [3Fe_4S]red = GED + [4Fe_3S]ox + [3Fe_4S]ox

k20_O2_FADH2 = kf_H2O2 * O2 * s2_FADH2;
k20_FUM_FADH2 = (kf_SUC/ Keq_FUM_FADH2)* FUM/KD_FUM/BP_FAD * s2_FADH2; 
k20_SQ_ISC3_QH2 =  kf_QH2 * s2_SQ_ISC3;
k20_GEA = kf_GEA*GEA;%/(GEA+p.KD_GEA);%*s2_ISC1_ISC3;
k2_0 = k20_O2_FADH2 + k20_FUM_FADH2 + k20_SQ_ISC3_QH2 + k20_GEA; 

% 5) E1->E2 
% SO + Q -> O2 + SQ
% SO + FAD + H -> FADH + O2
% SO + FADH + H -> FADH2 + O2

k12_SO_ISC3 = (kf_SO2/Keq_SO_ISC3) * 1/BP_Q * SO * s1_SQempty;% / (1 + GEA/p.KD_GEA);
k12_SO_Q = (kf_SO/Keq_SO_SQ) * Q/KD_Q/BP_Q * SO * s1_SQempty;
k12_SO_FAD = (kf_FADH/Keq_SO_FADH)* SO * s1_FAD;
k12_SO_FADH = (kf_FADH/Keq_SO_FADH2)* SO * s1_FADH;
k1_2 =   k12_SO_Q + k12_SO_FAD + k12_SO_FADH + k12_SO_ISC3;

% 6) E2->E1
% O2 + SQ -> SO + Q
% O2 + FADH -> FAD + H + SO
% O2 + FADH2 -> FADH + H+ SO

k21_O2_ISC3 = kf_SO2 * O2 * s2_ISC3 /  BP_Q;% / (1 + GEA/p.KD_GEA);
k21_O2_SQ = kf_SO * O2 * s2_SQ;
k21_O2_FADH = kf_FADH * O2 * s2_FADH;
k21_O2_FADH2 = kf_FADH2 * O2 * s2_FADH2;
k2_1 = k21_O2_SQ  + k21_O2_FADH + k21_O2_FADH2 + k21_O2_ISC3;

% 7) E1->E3
% H2O2 + FAD -> O2 + FADH2
% FAD + Succinate -> FADH2 + fumarate
% QH2 + [4Fe_3S]ox -> SQ + [4Fe_3S]red + 2H
% GED + [4Fe_3S]ox + [3Fe_4S]ox -> GEA  + [4Fe_3S]red + [3Fe_4S]red

k13_H2O2_FAD = (kf_H2O2/Keq_H2O2_FADH2) * H2O2 * s1_FAD;
k13_SUC_FAD =  kf_SUC* SUC/KD_SUC/BP_FAD * s1_FAD; 
k13_QH2_SQempty_ISC3ox =  (kf_QH2/ Keq_QH2_ISC3)* QH2/KD_QH2/BP_Q * s1_SQempty_ISC3ox;
k13_GED = (kf_GEA/Keq_GEA_ISC)*GED;%/(GED+p.KD_GEA);%*s1_ISC1ox_ISC3ox; 
k1_3 = k13_H2O2_FAD + k13_SUC_FAD + k13_QH2_SQempty_ISC3ox + k13_GED; 

% 8) E3->E1
% O2 + FADH2 ->  H2O2 + FAD
% FADH2 + fumarate -> FAD + Succinate
% SQ-[4Fe_3S]red + 2H -> QH2 + [4Fe_3S]ox
% GEA + [4Fe_3S]red + [3Fe_4S]red = GED + [4Fe_3S]ox + [3Fe_4S]ox

k31_O2_FADH2 = kf_H2O2 * O2 * s3_FADH2;
k31_FUM_FADH2 = (kf_SUC/ Keq_FUM_FADH2)* FUM/KD_FUM/BP_FAD * s3_FADH2;
k31_SQ_ISC3_QH2 =  kf_QH2* s3_SQ_ISC3;
k31_GEA = kf_GEA*GEA;%/(GEA+p.KD_GEA);%*s3_ISC1_ISC3;
k3_1 = k31_O2_FADH2 + k31_FUM_FADH2 + k31_SQ_ISC3_QH2 + k31_GEA;

% 9) E2->E3
% SO + ISC3ox -> O2 + ISC3red
% SO + Q -> O2 + SQ
% SO + FAD + H -> FADH + O2
% SO + FADH + H -> FADH2 + O2

k23_SO_ISC3 = (kf_SO2/Keq_SO_ISC3) * 1/BP_Q * SO * s2_SQempty;% / (1 + GEA/p.KD_GEA);
k23_SO_Q = (kf_SO/Keq_SO_SQ) * Q/KD_Q/BP_Q * SO * s2_SQempty;
k23_SO_FAD = (kf_FADH/Keq_SO_FADH)* SO * s2_FAD;
k23_SO_FADH = (kf_FADH/Keq_SO_FADH2)* SO * s2_FADH;
k2_3 =   k23_SO_Q + k23_SO_FADH + k23_SO_FAD + k23_SO_ISC3;

% 10) E3->E2
% O2 + SQ -> SO + Q
% O2 + FADH -> FAD + H+ + SO
% O2 + FADH2 -> FADH + H + SO

k32_O2_ISC3 = kf_SO2 * O2 * s3_ISC3 / BP_Q;% / (1 + GEA/p.KD_GEA);
k32_O2_SQ = kf_SO * O2 * s3_SQ;
k32_O2_FADH = kf_FADH * O2 * s3_FADH;
k32_O2_FADH2 = kf_FADH2 * O2 * s3_FADH2;
k3_2= k32_O2_SQ + k32_O2_FADH + k32_O2_FADH2 + k32_O2_ISC3;

% 11) E2->E4
% H2O2 + FAD -> O2 + FADH2
% FAD + Succinate -> FADH2 + fumarate
% QH2 + [4Fe_3S]ox -> SQ + [4Fe_3S]red + 2H
% GED + [4Fe_3S]ox + [3Fe_4S]ox -> GEA  + [4Fe_3S]red + [3Fe_4S]red

k24_H2O2_FAD = (kf_H2O2/Keq_H2O2_FADH2) * H2O2 * s2_FAD;
k24_SUC_FAD = kf_SUC* SUC/KD_SUC/BP_FAD * s2_FAD; 
k24_QH2_SQempty_ISC3ox = (kf_QH2/Keq_QH2_ISC3) * QH2/KD_QH2/BP_Q * s2_SQempty_ISC3ox;
k24_GED = (kf_GEA/Keq_GEA_ISC)*GED;%/(GED+p.KD_GEA);%*s2_ISC1ox_ISC3ox; 
k2_4 = k24_H2O2_FAD + k24_SUC_FAD + k24_QH2_SQempty_ISC3ox + k24_GED;

% 12) E4->E2
% O2 + FADH2 ->  H2O2 + FAD
% FADH2 + fumarate -> FAD + Succinate
% SQ-[4Fe_3S]red + 2H -> QH2 + [4Fe_3S]ox
% GEA + [4Fe_3S]red + [3Fe_4S]red = GED + [4Fe_3S]ox + [3Fe_4S]ox

k42_O2_FADH2 = kf_H2O2 * O2 * s4_FADH2;
k42_FUM_FADH2 = (kf_SUC/ Keq_FUM_FADH2)* FUM/KD_FUM/BP_FAD * s4_FADH2; 
k42_SQ_ISC3_QH2 =  kf_QH2* s4_SQ_ISC3;
k42_GEA = kf_GEA*GEA;%/(GEA+p.KD_GEA);%*s4_ISC1_ISC3; 
k4_2 = k42_O2_FADH2 + k42_FUM_FADH2 + k42_SQ_ISC3_QH2 + k42_GEA;

% 13) E3->E4
% SO + Q -> O2 + SQ
% SO + FAD + H -> FADH + O2
% SO + FADH + H -> FADH2 + O2

k34_SO_ISC3 = (kf_SO2/Keq_SO_ISC3) * 1/BP_Q * SO * s3_SQempty;% / (1 + GEA/p.KD_GEA); %%
k34_SO_Q = (kf_SO/Keq_SO_SQ) * Q/KD_Q/BP_Q * SO * s3_SQempty; %%
k34_SO_FAD = (kf_FADH/Keq_SO_FADH)* SO * s3_FAD;
k34_SO_FADH = kf_FADH2 *SO * s3_FADH;
k3_4 =   k34_SO_Q + k34_SO_FAD + k34_SO_FADH + k34_SO_ISC3;

% 14) E4->E3
% O2 + SQ -> SO + Q
% O2 + FADH -> FAD + H + SO
% O2 + FADH2 -> FADH + H+ SO

k43_O2_ISC3 = kf_SO2 * O2 * s4_ISC3 / BP_Q;% / (1 + GEA/p.KD_GEA);
k43_O2_SQ = kf_SO * O2 * s4_SQ;
k43_O2_FADH = kf_FADH * O2 * s4_FADH;
k43_O2_FADH2 =kf_FADH2 * O2 * s4_FADH2;
k4_3 = k43_O2_SQ + k43_O2_FADH + k43_O2_FADH2 + k43_O2_ISC3;

% electron occupancy steady states
D = (k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
E0= (k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2)/D;
E1= (k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3)/D;
E2= (k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2)/D;
E3= (k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3)/D;
E4= 1 -(E0+E1+E2+E3);

% Compute steady-state rates

Vmf_SDH=0.12237e-3 ;

% succinate oxidation rate
J_SUCC = Vmf_SDH*(k02_SUC_FAD*E0 + k13_SUC_FAD*E1 + k24_SUC_FAD*E2 ...
   - k20_FUM_FADH2*E2 - k31_FUM_FADH2*E3 - k42_FUM_FADH2*E4);

% Superoxide formation
J_SDH_SO = Vmf_SDH*((k10_O2_ISC3 + k10_O2_SQ + k10_O2_FADH)*E1 + (k21_O2_ISC3 + k21_O2_SQ + k21_O2_FADH + k21_O2_FADH2)*E2 + (k32_O2_ISC3 + k32_O2_SQ + k32_O2_FADH + k32_O2_FADH2)*E3 + (k43_O2_ISC3 + k43_O2_SQ + k43_O2_FADH + k43_O2_FADH2)*E4...
                  - (k01_SO_ISC3 + k01_SO_Q + k01_SO_FAD)*E0 - (k12_SO_ISC3 + k12_SO_Q + k12_SO_FAD + k12_SO_FADH)*E1 - (k23_SO_ISC3 + k23_SO_Q + k23_SO_FAD + k23_SO_FADH)*E2 - (k34_SO_ISC3 + k34_SO_Q + k34_SO_FAD + k34_SO_FADH)*E3 );
            
% Hydrogen peroxide formation
J_SDH_H2O2 = Vmf_SDH*(k20_O2_FADH2*E2 + k31_O2_FADH2*E3 + k42_O2_FADH2*E4 ...
                    - k02_H2O2_FAD*E0 - k13_H2O2_FAD*E1 - k24_H2O2_FAD*E2);

J_SDH_matrix = J_SUCC;

% J_FUM_matrix: DAB coded 
fum=fumarate_matrix/P(22);
mal=malate_matrix/P(23);
phos=Pi_matrix/P(14);
H = h_matrix;
I=0.17;
alphaT=1.10708-1.54508e-3*(T+273.15)+5.95584e-6*(T+273.15)^2;
gamma1=exp(-sqrt(I)*alphaT*1^2/(1+1.6*sqrt(I)));
gamma2=exp(-sqrt(I)*alphaT*2^2/(1+1.6*sqrt(I)));
gamma3=exp(-sqrt(I)*alphaT*3^2/(1+1.6*sqrt(I)));

k46 = 1.9115e+06;
k45 = 4.4648e+05;
k54 = gamma2*54.1164;
k21f = gamma3*5.8202e+04;
K14m = 0; 
K25 = gamma2/gamma3*1.8077e-10;
K26 = gamma2/gamma3*8.2785e-08;
K14f = (7.6946e-09)/gamma1;
K_aPI =1.4736e-03; % 1/26 E(g2) = 11.07
K_aF = 0.0315;
Kip1m = gamma3/gamma1*0.0460;
Kip4 = gamma2*0.0417;
Kip1f = gamma3/gamma1*0.0143;
Kif4 = gamma2*0.0907;
Kif1f = gamma3/gamma1*0.0040;

% Thermodynamic constraints
Keq = 3.5006;
% k21m = K14m*k12m*k54*Keq/(k45*K25);
k64 = Keq*K26*k46*k54/(k45*K25);
k12f = k45*K25*k21f/(K14f*k54);

% 6 state-model with H binding
E0_FUM=0.22088e-3;
R = (mal*k46  + fum*(k45 + k12f*K14f/H))/(k64*H/K26 + k54*H/K25 + k21f);
E4 = E0_FUM / (K14m/H*(1 + p/Kip1m) + K14f/H*(1 + fum/Kif1f + phos/Kip1f) + ...
    R + (1 + fum/Kif4 + phos/Kip4) + R*H/K25 + R*H/K26 );
E2 = R*E4;
E6 = H/K26*E2;
J_FUM_matrix = ( phos/K_aPI  + fum/K_aF )/(1 + phos/K_aPI  + fum/K_aF )*...
     (k64*E6 - k46*mal*E4 );

% MDH_matrix--DAB version
a = NAD_free;
b = malate_matrix;
p = oxaloacetate_matrix;
q = NADH_free;
AA=1.10708-(1.54508e-3)*(T+273.15)+(5.95584e-6)*(T+273.15).^2;
BB=1.6;
gamma1=exp(-AA*1*sqrt(I)/(1+BB*sqrt(I)));
gamma2=exp(-AA*4*sqrt(I)/(1+BB*sqrt(I)));
Vmf_MDH=0.052998;
k1p=5.80E+04*1e3/60;
%unitsconvertedtoM^-1s^-1
km1p=5.36E+05/60;
%unitsconvertedtos^-1
k2p=6.54E+04*1e3/60;
km2p=2.91E+04/60;
k3p=1.38E+05*gamma1/60;
%km3p=par(8)*gamma1*1e3/60;
k4p=2.24E+05/60;
km4p=4.62E+05*1e3/60;
km4pp=1.81E+05*gamma1*1e3/60;
km4ppp=8.91E+05*gamma2*1e3/60;
kD01=10^-7.60;
kD02=10^-5.88;
kDQ1=10^-9.33;
kDQ2=10^-5.72;
k4pp=(kDQ1/kD01)*(km4pp/km4p)*k4p;
k4ppp=(kDQ2/kD02)*(km4ppp/km4pp)*k4pp;
kD01=h_matrix/kD01;
kD02=h_matrix/kD02;
kDQ1=h_matrix/kDQ1;
kDQ2=h_matrix/kDQ2;
p01=1+kD01+kD01*kD02;
p0Q=1+kDQ1+kDQ1*kDQ2;
k1=k1p/p01;
km1=km1p;
% k2=k2p; % 
k2=k2p; % DAB malate inhibition
km2=km2p;
k3=k3p;
%km3=km3p*(h_matrix/kDQ1)/p0Q;
k4=(k4p+k4pp*kDQ1+k4ppp*kDQ1*kDQ2)/p0Q;
km4=(km4p+km4pp*kD01+km4ppp*kD01*kD02)/p01;
km3=(k1*k2*k3*k4)/(Keq_MDH_matrix*km1*km2*km4);

KmA = k3*k4/(k1*(k3+k4));
KmB = k4*(km2+k3)/(k2*(k3+k4));
KmP = km1*(km2+k3)/(km3*(km1+km2));
KmQ = km1*km2/(km4*(km1+km2));
KeA = km1/k1;
KeB = km2/k2;
KeP = k3/km3;
KeQ = k4/km4;

I1 = 1;
I2 = 1;
I3 = 1;
I4 = 1;

Den = (1 + KmA*b/(KeA*KmB) + KmQ*p/(KeQ*KmP))*I1 + ...
      (a/KeA + KmQ*a*p/(KeA*KmP*KeQ) + KmA*p*q/(KeA^2*KmB*Keq_MDH_matrix))*I2 + ...
      ((1/(KeA*KmB) - KmQ*Keq_MDH_matrix/(KeQ^2*KmP))*a*b + (1/(KmP*KeQ) - KmA/(KeA^2*KmB*Keq_MDH_matrix))*p*q + KmQ*a*b*p/(KeA*KeB*KmP*KeQ) + KmA*b*p*q/(KeA*KmB*KeP*KeQ))*I3 + ...
      (q/KeQ + KmQ*Keq_MDH_matrix*a*b/(KeQ^2*KmP) + KmA*b*q/(KeA*KmB*KeQ))*I4;
J_MDH_matrix = (Vmf_MDH/(KeA*KmB))*(a*b - p*q/Keq_MDH_matrix)/Den; 


% NDK_matrix
a=GTP_matrix;
b=ADP_matrix;
p=GDP_matrix;
q=ATP_matrix;
KmA=111e-6;
KmB=100e-6;
KmP=260e-6;
KmQ=278e-6;
Kia=170e-6;
Kib=143.6e-6;
Kip=146.6e-6;
Kiq=156.5e-6;
KiAMP=650e-6;
fAMP=AMP_matrix*(1+h_matrix/Kh(13))/P(13);
ai=1+fAMP/KiAMP;
Vmf_NDK=0.3857;
Vmr=Vmf_NDK/Keq_NDK_matrix*(KmQ*Kip/Kia/KmB);
if(a>MinCon)&&(b>MinCon)
ab=a*b;
else
ab=0;
end
if(p>MinCon)&&(q>MinCon)
pq=p*q;
else
pq=0;
end
J_NDK_matrix=Vmf_NDK*Vmr*(ab-pq/Keq_NDK_matrix)/ai/(Vmr*KmB*a+Vmr*KmA*b+Vmf_NDK*KmQ/Keq_NDK_matrix*p+Vmf_NDK*KmP/Keq_NDK_matrix*q+Vmr*a*b+Vmf_NDK*KmQ/Kia/Keq_NDK_matrix*a*p+Vmf_NDK/Keq_NDK_matrix*p*q+Vmr*KmA/Kiq*b*q);

%ATPASE_buffer
dG0 = -4.67;
Keq = exp(-dG0)/h_buffer*P(28)*P(29)/P(27);
if ATP_buffer > 1.0e-6
  J_ATPASE_buffer = x_ATPase/( 1 + ADP_buffer/(0.10e-3))*(1 - ADP_buffer*Pi_buffer/(Keq)/ATP_buffer); % from Bazil model
else
  J_ATPASE_buffer = 0;
end

%F1F0ATPASE:im_to_matrix
k1_F1F0 = 0.80089e-3;
dG0 = -4.67;
nH = 8/3;
Keq_F1F0ATPASE = exp(-(dG0 - nH*F*DPsi_im_to_matrix)/RT)*(h_im^nH/h_matrix^(nH-1))*P(10)/P(11)/P(14);
% J_F1F0ATPASE_im_to_matrix = k1_F1F0*(Keq_F1F0ATPASE*ADP_matrix*Pi_matrix-ATP_matrix);

% Random bi-bi model for forward reaction 
a = ADP_matrix;
b = Pi_matrix;
p = ATP_matrix;
Ka = 1.0e-3;
Kb = 2.0e-3;
J_F1F0ATPASE_im_to_matrix = (k1_F1F0/(Ka*Kb))*(Keq_F1F0ATPASE*a*b - p)/(1 + a/Ka + b/Kb + a*b/(Ka*Kb) );

%ETC1:im_to_matrix
ETC1Activity = 0.27496e-3 ;
% NADH_free = NADH_free;
% beta_n=1;
Hp = h_im;
Hn = h_matrix;
NADH1 = NADH_free;
% NADH1 = NADH_matrix;
% NADH_free;
% NAD1 = NAD_matrix;
NAD1 = NAD_free;
Q1 = coQ_matrix*lambdaQ;  % 
QH21 = coQH2_matrix*lambdaQ;  % 
% Q1 = coQ_matrix; 
% QH21 = coQH2_matrix;  
O2 = O2aq_matrix;
H2O2 = H2O2aq_matrix;
SO = SOaq_matrix;
% H2O2 and SO new variables defined
% Update midpoint potentials from thermodynamic data
% dGf_UQH2 = DGfo_coQH2;
% dGf_UQ = DGfo_coQ;
% dGf_NADH = DGfo_NADH;
% dGf_NAD = DGfo_NAD;
% dGf_SO = DGfo_SOaq;
% dGf_H2O2 = DGfo_H2O2aq;
% dGf_O2 = DGfo_O2aq;
% dGf_UQH2 = -19.150559;
% dGf_UQ = 69.2309743;
% dGf_NADH = 40.342360;
% dGf_NAD = 19.230826;
% dGf_SO = 12.4107669;
% dGf_H2O2 = -132.491533;
% dGf_O2 = 16.4;
% CI_Em0_Q_QH2 = (dGf_UQH2 - dGf_UQ)/(-2*F);
% this is 18.67 mV less than value used for CI parameterization
% CI_Em0_NADH = (dGf_NADH - dGf_NAD)/(-2*F);
% this is 2.3 mV higher than value used for CI parameterization
% CI_Em0_SO = (dGf_SO - dGf_O2)/(-1*F);
% this is 8.3 mV higher than value used for CI parameterization
% CI_Em0_H2O2 = (dGf_H2O2 - dGf_O2)/(-2*F);
% this is 8.4 mV less than value used for CI parameterization
CI_Em0_Q_QH2 = 476.6736;
% this is 18.67 mV less than value used for CI parameterization
CI_Em0_NADH = -111.7008;
% this is 2.3 mV higher than value used for CI parameterization
CI_Em0_SO = 33;
% this is 8.3 mV higher than value used for CI parameterization
CI_Em0_H2O2 = 780;
% this is 8.4 mV less than value used for CI parameterization
% Update midpoint potentials for 1e- quinone couples (pH 0)
CI_Kstability = 10;
CI_Em0_Q_SQ = CI_Em0_Q_QH2 + RT/F/2*log(CI_Kstability*1e-14);
% mV
CI_Em0_SQ_QH2 = 2*CI_Em0_Q_QH2 - CI_Em0_Q_SQ;
% mV
% Binding Polynomials for Protonated States
CI_KiH1 = 4.0722e-08;
CI_KiH3 = 3.8958e-07;
PH1 = (1/(1+Hn/CI_KiH1));
PH2 = 1;
PH3 = (Hn/CI_KiH3/(1+Hn/CI_KiH3));
% Binding Polynomials for enzyme, substrates, products and regulators
% NADH binding constants
KdNADHo=(4.6065e-05); % 
KdNADo=7.0545e-04;
KdNADHr=4.9867e-04;
KdNADr=1.1844e-05;
KdNADHrad=1;
KdNADrad=1.5358e-04;
% Quinone binding constants
KdQH2o = 0.1;
KdQo = 0.0175;
KdQH2r = 0.1;
KdQr = 0.0175;
% Binding polynomials
PNo = 1 + NADH1/KdNADHo + NAD1/KdNADo;
% oxidized binding constants
PNr = 1 + NADH1/KdNADHr + NAD1/KdNADr;
% reduced binding constants
PNrad = 1 + NADH1/KdNADHrad + NAD1/KdNADrad;
% reduced binding constants
PQr = 1 + QH21/KdQH2r + Q1/KdQr;
% reduced binding constants
PQo = 1 + QH21/KdQH2o + Q1/KdQo;
% oxidized binding constants
muH = F*DPsi_im_to_matrix + RT*log(Hp/Hn);
% proton chemical potential (kJ/mol)
% Compute pH Corrected Midpoint potentials
% NADH potential
Em_NADH = CI_Em0_NADH + log(10)*RT/F/2*log10(Hn);
% FMN potentials
CI_Em0_FMN = 55.1446;
CI_Em0_FMN2 = 86.7495;
CI_Em0_FMN1 = 23.5369;
CI_FMNrad = 7.9074;
CI_FMNred_pK = 7.0998;
Em_FMNred = CI_Em0_FMN - RT/2/F*log(Hn/10^-CI_FMNred_pK/(1+Hn/10^-CI_FMNred_pK)/Hn.^2);
% FMNred/FMN, pK1 of fully reduced FMN very high (>10)0)
Em_FMNrad = CI_Em0_FMN2 - RT/F*log((Hn/10^-(CI_FMNred_pK)/(1+Hn/10^-CI_FMNred_pK))/((Hn/10^-CI_FMNrad)/(1+Hn/10^-CI_FMNrad))/Hn);
% FMNred/FMNrad
Em_FMNox = CI_Em0_FMN1 - RT/F*log((Hn/10^-CI_FMNrad)/(1+Hn/10^-CI_FMNrad)/Hn);
%FMN/FMNrad
% N2 potential
CI_Em0_N2 = -90;
CI_N2_pKox = 6;
CI_N2_pKred = 8.5;
Em_N2 = CI_Em0_N2 - RT/F*log((Hn+10^-CI_N2_pKox)/(Hn+10^-CI_N2_pKred));
% superoxide and H2O2 potentials
Em_SO = CI_Em0_SO;
% Em0 stays the same for pH > 6, pKa of SO is ~4.8
Em_H2O2 = CI_Em0_H2O2 + log(10)*RT/F*log10(Hn);
% 1st pKa is ~11
% Quinone potentials
Em_Q_QH2 = CI_Em0_Q_QH2 + 2*log(10)*RT/F/2*log10(Hn);
% assuming linked to N-side
Em_Q_SQ = CI_Em0_Q_SQ;
% pH independent
Em_SQ_QH2 = CI_Em0_SQ_QH2 + 2*log(10)*RT/F*log10(Hn);
% assuming linked to N-side
% State Transition Thermodynamics
% NADH-QH2 reductase related
K_NADHFMNred = exp((2*F*(Em_FMNred - Em_NADH))/RT);
% [Fr][NAD]/[F][NADH][H]
K_FMNradN2 = exp((F*(Em_N2 - Em_FMNrad))/RT);
% [Frad][N2r]/[Fr]/[N2]
K_FMNoxN2 = exp((F*(Em_N2 - Em_FMNox))/RT);
% [F][N2r]/[Frad]/[N2]
K_N2SQ = exp((F*(Em_Q_SQ - Em_N2))/RT);
% [N2][SQ]/[N2r][Q]
K_N2QH2 = exp((F*(Em_SQ_QH2 - Em_N2))/RT);
% [N2][QH2]/[N2r][SQ]
% superoxide and hydrogen peroxide
K_FMNredH2O2 = exp((2*F*(Em_H2O2 - Em_FMNred))/RT);
% [F][H2O2]/[Fr][O2]
K_FMNradO2 = exp((F*(Em_SO - Em_FMNrad))/RT);
% [Frad][O2.-]/[Fr][O2]
K_FMNoxO2 = exp((F*(Em_SO - Em_FMNox))/RT);
% [F][O2.-]/[Frad][O2]
K_N2O2 = exp((F*(Em_SO - Em_N2))/RT);
% [N2][O2.-]/[N2r][O2]
K_SQO2 = exp((F*(Em_SO - Em_Q_SQ))/RT);
% [Q][O2.-]/[SQ][O2]
% Substates via Rapid Equilibrium
% substate equilibrium constants
K1 = K_FMNradN2;
%[Frad][N2r]/[Fr][N2]
K2 = K_FMNoxN2;
% [F][N2r]/[Frad][N2]
KSQ = K_N2SQ*KdQr;
% [N2-SQ]/[N2r-Q], KdQr = [N2r][Q]/[N2r-Q]
KQ = (Q1/KdQr/PQr)*KSQ;
% [N2-SQ]/[N2r]
% 0e-
PS0 = 1;
% total
s0_F_N2 = 1/PS0;
% 1e-
PS1 = (1 + K2*(1 + KQ));
% total sum(cell2mat(struct2cell(s1)))
s1_Frad_N2 = 1/PS1;
s1_F_N2r = K2*s1_Frad_N2;
s1_F_N2_SQ = KQ*s1_F_N2r;
% 2e-
PS2 = (1 + K1*(1 + KQ*(1 + K2)));
% total sum(cell2mat(struct2cell(s2)))
s2_Fr_N2 = 1/PS2;
s2_Frad_N2r = K1*s2_Fr_N2;
s2_Frad_N2_SQ = KQ*s2_Frad_N2r;
s2_F_N2r_SQ = K2*s2_Frad_N2_SQ;
% 3e-
PS3 = (1 + KQ*(1 + K1));
% total sum(cell2mat(struct2cell(s3)))
s3_Fr_N2r = 1/PS3;
s3_Fr_N2_SQ = KQ*s3_Fr_N2r;
s3_Frad_N2r_SQ = K1*s3_Fr_N2_SQ;
% 4e-
PS4 = 1;
% total
s4_Fr_N2r_SQ = 1/PS4;
% State Transition Rates
% reverse NADH-Q reductase rates
CI_kfNADH_02 = 1.9642e+03;
CI_kfNADH_24 = 186.3167;
CI_kfNADH_13 = 4.6093;
CI_kfQ_20 = 5.8175e+03;
CI_kfQ_42 = 7.7767e+10;
CI_kfQ_31 = 8.6596;
krN0_02 = CI_kfNADH_02/K_NADHFMNred*KdNADr/KdNADHo*PH1*PH2;
krN0_24 = CI_kfNADH_24/K_NADHFMNred*KdNADr/KdNADHo*PH1*PH2;
krN0_13 = CI_kfNADH_13/K_NADHFMNred*KdNADr/KdNADHo*PH1*PH2;
krQ0_20 = CI_kfQ_20/K_N2QH2*KdQH2o*(PNr*PQo/PNo/PQr)*PH3;
krQ0_42 = CI_kfQ_42/K_N2QH2*KdQH2o*(PNr*PQo/PNo/PQr)*PH3;
krQ0_31 = CI_kfQ_31/K_N2QH2*KdQH2o*(PNr*PQo/PNo/PQr)*PH3;
% reverse ROS rates
% O2 and SQ
CI_kfSQSO_10 = 1.6667e+09;
CI_kfSQSO_21 = 0.0021;
CI_kfSQSO_21b = 0.0021;
CI_kfSQSO_32 = 1.1450e-05;
CI_kfSQSO_32b = 2.3953e-09;
CI_kfSQSO_43 = 2.1602e-07;
krSQSO_10 = CI_kfSQSO_10/K_SQO2;
krSQSO_21 = CI_kfSQSO_21/K_SQO2;
krSQSO_21b = CI_kfSQSO_21b/K_SQO2;
krSQSO_32 = CI_kfSQSO_32/K_SQO2;
krSQSO_32b = CI_kfSQSO_32b/K_SQO2;
krSQSO_43 = CI_kfSQSO_43/K_SQO2;
% O2 and N2
CI_kfN2rSO_10 = 0;
CI_kfN2rSO_21 = 0;
CI_kfN2rSO_21b = 0;
CI_kfN2rSO_32 = 0;
CI_kfN2rSO_32b = 0;
CI_kfN2rSO_43 = 0;
krN2rSO_10 = CI_kfN2rSO_10/K_N2O2;
krN2rSO_21 = CI_kfN2rSO_21/K_N2O2;
krN2rSO_21b = CI_kfN2rSO_21b/K_N2O2;
krN2rSO_32 = CI_kfN2rSO_32/K_N2O2;
krN2rSO_32b = CI_kfN2rSO_32b/K_N2O2;
krN2rSO_43 = CI_kfN2rSO_43/K_N2O2;
% O2 and FMNred
CI_kfFrSO_21 = 451900;
CI_kfFradSO_21b = 1.1943e-07;
CI_kfFrSO_32 = 5.3833e-04;
CI_kfFrSO_32b = 2.8333e-05;
CI_kfFrSO_43 = 1.4005e+05;
krFrSO_21 = CI_kfFrSO_21/K_FMNradO2;
krFradSO_21b = CI_kfFradSO_21b/K_FMNradO2;
krFrSO_32 = CI_kfFrSO_32/K_FMNradO2;
krFrSO_32b = CI_kfFrSO_32b/K_FMNradO2;
krFrSO_43 = CI_kfFrSO_43/K_FMNradO2;
% O2 and FMNrad
CI_kfFradSO_10 = 68055000;
CI_kfFradSO_21 = 4.0417e-11;
CI_kfFradSO_32 = 0.014;
krFradSO_10 = CI_kfFradSO_10/K_FMNoxO2;
krFradSO_21 = CI_kfFradSO_21/K_FMNoxO2;
krFradSO_32 = CI_kfFradSO_32/K_FMNoxO2;
% H2O2 and FMNred
CI_kfH_20 = 7.8243e-07;
CI_kfH_31 = 2.1228e+06;
CI_kfH_31b = 1.7637e+06;
CI_kfH_42 = 62.0583;
krH_20 = CI_kfH_20/K_FMNredH2O2;
krH_31 = CI_kfH_31/K_FMNredH2O2;
krH_31b = CI_kfH_31b/K_FMNredH2O2;
krH_42 = CI_kfH_42/K_FMNredH2O2;
% Nucleotide and pH effects
CI_kfH_20 = CI_kfH_20/PNr;
CI_kfH_42 = CI_kfH_42/PNr;
CI_kfH_31 = CI_kfH_31/PNr;
CI_kfH_31b = CI_kfH_31b/PNr;
CI_kfFrSO_21 = CI_kfFrSO_21/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFrSO_32 = CI_kfFrSO_32/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFrSO_32b = CI_kfFrSO_32b/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFrSO_43 = CI_kfFrSO_43/PNr/(1+Hn/10^-CI_FMNred_pK);
CI_kfFradSO_10 = CI_kfFradSO_10/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfFradSO_21 = CI_kfFradSO_21/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfFradSO_21b = CI_kfFradSO_21b/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfFradSO_32 = CI_kfFradSO_32/PNrad/(1+Hn/10^-CI_FMNrad);
CI_kfN2rSO_10 = CI_kfN2rSO_10/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_21 = CI_kfN2rSO_21/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_21b = CI_kfN2rSO_21b/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_32 = CI_kfN2rSO_32/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_32b = CI_kfN2rSO_32b/(1+Hn/10^-CI_N2_pKred);
CI_kfN2rSO_43 = CI_kfN2rSO_43/(1+Hn/10^-CI_N2_pKred);
% net forward rates % membrane potential dependence
dPsiNf = 1;
dPsiNr = 1;
CI_beta1 = 0.5;
dPsiQf = exp(-4*CI_beta1*muH/RT);
dPsiQr = exp(4*(1-CI_beta1)*muH/RT);
% NADH transitions
kfNADH_02 = CI_kfNADH_02*dPsiNf*NADH1/KdNADHo/PNo*(s0_F_N2)*PH1*PH2;
krNADH_02 = krN0_02*dPsiNr*NAD1/KdNADr/PNr*(s2_Fr_N2);
kfNADH_24 = CI_kfNADH_24*dPsiNf*NADH1/KdNADHo/PNo*(s2_F_N2r_SQ)*PH1*PH2;
krNADH_24 = krN0_24*dPsiNr*NAD1/KdNADr/PNr*(s4_Fr_N2r_SQ);
kfNADH_13 = CI_kfNADH_13*dPsiNf*NADH1/KdNADHo/PNo*(s1_F_N2r + s1_F_N2_SQ)*PH1*PH2;
krNADH_13 = krN0_13*dPsiNr*NAD1/KdNADr/PNr*(s3_Fr_N2r + s3_Fr_N2_SQ);
% QH2 transitions
kfQ_20 = CI_kfQ_20*dPsiQf*(s2_F_N2r_SQ)*PH3;
krQ_20 = krQ0_20*dPsiQr*QH21/KdQH2o/PQo*(s0_F_N2);
kfQ_42 = CI_kfQ_42*dPsiQf*(s4_Fr_N2r_SQ)*PH3;
krQ_42 = krQ0_42*dPsiQr*QH21/KdQH2o/PQo*(s2_Fr_N2);
kfQ_31 = CI_kfQ_31*dPsiQf*(s3_Frad_N2r_SQ)*PH3;
krQ_31 = krQ0_31*dPsiQr*QH21/KdQH2o/PQo*(s1_Frad_N2);
% H2O2 transitions
kfH2O2_20 = CI_kfH_20*O2*s2_Fr_N2;
krH2O2_20 = krH_20*H2O2*s0_F_N2;
kfH2O2_31 = CI_kfH_31*O2*s3_Fr_N2r + CI_kfH_31b*O2*s3_Fr_N2_SQ;
krH2O2_31 = krH_31*H2O2*s1_F_N2r + krH_31b*H2O2*s1_F_N2_SQ;
kfH2O2_42 = CI_kfH_42*O2*s4_Fr_N2r_SQ;
krH2O2_42 = krH_42*H2O2*(s2_F_N2r_SQ);
% SO transitions
kfSO_10 = O2*(CI_kfSQSO_10*s1_F_N2_SQ + CI_kfFradSO_10*s1_Frad_N2 + CI_kfN2rSO_10*s1_F_N2r);
krSO_10 = SO*s0_F_N2*(krSQSO_10*(Q1/KdQo/PQo) + krFradSO_10 + krN2rSO_10);
kfSO_21 = O2*((CI_kfSQSO_21*s2_Frad_N2_SQ + CI_kfSQSO_21b*s2_F_N2r_SQ) + CI_kfFrSO_21*s2_Fr_N2 + (CI_kfFradSO_21*s2_Frad_N2r + CI_kfFradSO_21b*s2_Frad_N2_SQ) + (CI_kfN2rSO_21*s2_Frad_N2r + CI_kfN2rSO_21b*s2_F_N2r_SQ));
krSO_21 = SO*(krFrSO_21*s1_Frad_N2 + krN2rSO_21*s1_Frad_N2 + krN2rSO_21b*s1_F_N2_SQ + krFradSO_21*s1_F_N2r + krFradSO_21b*s1_F_N2_SQ + (krSQSO_21*s1_Frad_N2*(Q1/KdQo/PQo) + krSQSO_21b*s1_F_N2r*(Q1/KdQr/PQr)));
kfSO_32 = O2*((CI_kfSQSO_32*s3_Fr_N2_SQ + CI_kfSQSO_32b*s3_Frad_N2r_SQ) + (CI_kfFrSO_32*s3_Fr_N2r + CI_kfFrSO_32b*s3_Fr_N2_SQ) + CI_kfFradSO_32*s3_Frad_N2r_SQ + (CI_kfN2rSO_32*s3_Fr_N2r + CI_kfN2rSO_32b*s3_Frad_N2r_SQ));
krSO_32 = SO*(krSQSO_32*s2_Fr_N2*(Q1/KdQo/PQo) + krSQSO_32b*s2_Frad_N2r*(Q1/KdQr/PQr) + krFradSO_32*s2_F_N2r_SQ + krFrSO_32*s2_Frad_N2r + krFrSO_32b*s2_Frad_N2_SQ + krN2rSO_32*s2_Fr_N2 + krN2rSO_32b*s2_Frad_N2_SQ);
kfSO_43 = O2*((CI_kfSQSO_43 + CI_kfN2rSO_43)*(s4_Fr_N2r_SQ) + CI_kfFrSO_43*s4_Fr_N2r_SQ);
krSO_43 = SO*(krSQSO_43*s3_Fr_N2r*(Q1/KdQr/PQr) + krFrSO_43*s3_Frad_N2r_SQ + krN2rSO_43*s3_Fr_N2_SQ);
% S0 <-> S2
k0_2 = kfNADH_02 + krQ_20 + krH2O2_20;
k2_0 = kfQ_20 + krNADH_02 + kfH2O2_20;
% S2 <-> S4
k2_4 = kfNADH_24 + krQ_42 + krH2O2_42;
k4_2 = kfQ_42 + krNADH_24 + kfH2O2_42;
% S1 <-> S3
k1_3 = kfNADH_13 + krQ_31 + krH2O2_31;
k3_1 = kfQ_31 + krNADH_13 + kfH2O2_31;
% S0 <-> S1 (superoxide)
k1_0 = kfSO_10;
k0_1 = krSO_10;
% S2 <-> S1 (superoxide)
k2_1 = kfSO_21;
k1_2 = krSO_21;
% S3 <-> S2 (superoxide)
k3_2 = kfSO_32;
k2_3 = krSO_32;
% S4 <-> S3 (superoxide)
k4_3 = kfSO_43;
k3_4 = krSO_43;
% Steady-State Fractional Occupancies (solved analytically)
S0=(k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S1=(k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S2=(k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S3=(k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3)/(k0_2*k1_0*k2_4*k3_1 + k0_1*k1_2*k2_4*k3_1 + k0_1*k1_3*k2_0*k3_4 + k0_2*k1_0*k2_4*k3_2 + k0_1*k1_2*k2_4*k3_2 + k0_1*k1_3*k2_1*k3_4 + k0_2*k1_0*k2_3*k3_4 + k0_2*k1_2*k2_4*k3_1 + k0_1*k1_2*k2_3*k3_4 + k0_1*k1_3*k2_4*k3_2 + k0_2*k1_0*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_2 + k0_2*k1_3*k2_1*k3_4 + k0_1*k1_2*k2_4*k3_4 + k0_1*k1_3*k2_3*k3_4 + k0_2*k1_2*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_2 + k0_1*k1_3*k2_4*k3_4 + k0_2*k1_2*k2_4*k3_4 + k0_2*k1_3*k2_3*k3_4 + k0_2*k1_3*k2_4*k3_4 + k0_1*k1_3*k2_0*k4_2 + k0_1*k1_3*k2_0*k4_3 + k0_1*k1_3*k2_1*k4_2 + k0_2*k1_0*k2_3*k4_2 + k0_1*k1_2*k2_3*k4_2 + k0_1*k1_3*k2_1*k4_3 + k0_2*k1_0*k2_3*k4_3 + k0_2*k1_3*k2_1*k4_2 + k0_1*k1_2*k2_3*k4_3 + k0_1*k1_3*k2_3*k4_2 + k0_2*k1_0*k2_4*k4_3 + k0_2*k1_2*k2_3*k4_2 + k0_2*k1_3*k2_1*k4_3 + k0_1*k1_2*k2_4*k4_3 + k0_1*k1_3*k2_3*k4_3 + k0_2*k1_2*k2_3*k4_3 + k0_2*k1_3*k2_3*k4_2 + k0_1*k1_3*k2_4*k4_3 + k0_2*k1_2*k2_4*k4_3 + k0_2*k1_3*k2_3*k4_3 + k0_2*k1_3*k2_4*k4_3 + k0_2*k1_0*k3_1*k4_2 + k0_1*k1_2*k3_1*k4_2 + k0_2*k1_0*k3_1*k4_3 + k0_2*k1_0*k3_2*k4_2 + k0_1*k1_2*k3_1*k4_3 + k0_1*k1_2*k3_2*k4_2 + k0_2*k1_0*k3_2*k4_3 + k0_2*k1_2*k3_1*k4_2 + k0_1*k1_2*k3_2*k4_3 + k0_1*k1_3*k3_2*k4_2 + k0_2*k1_0*k3_4*k4_2 + k0_2*k1_2*k3_1*k4_3 + k0_2*k1_2*k3_2*k4_2 + k0_1*k1_2*k3_4*k4_2 + k0_1*k1_3*k3_2*k4_3 + k0_2*k1_2*k3_2*k4_3 + k0_2*k1_3*k3_2*k4_2 + k0_1*k1_3*k3_4*k4_2 + k0_2*k1_2*k3_4*k4_2 + k0_2*k1_3*k3_2*k4_3 + k0_2*k1_3*k3_4*k4_2 + k0_1*k2_0*k3_1*k4_2 + k0_1*k2_0*k3_1*k4_3 + k0_1*k2_0*k3_2*k4_2 + k0_1*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_2*k4_3 + k0_1*k2_1*k3_1*k4_3 + k0_1*k2_1*k3_2*k4_2 + k0_2*k2_1*k3_1*k4_2 + k0_1*k2_0*k3_4*k4_2 + k0_1*k2_1*k3_2*k4_3 + k0_1*k2_3*k3_1*k4_2 + k0_2*k2_1*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_2 + k0_1*k2_1*k3_4*k4_2 + k0_1*k2_3*k3_1*k4_3 + k0_2*k2_1*k3_2*k4_3 + k0_2*k2_3*k3_1*k4_2 + k0_1*k2_4*k3_1*k4_3 + k0_2*k2_1*k3_4*k4_2 + k0_2*k2_3*k3_1*k4_3 + k0_2*k2_4*k3_1*k4_3 + k1_0*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_1*k4_3 + k1_0*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_1*k4_2 + k1_0*k2_0*k3_2*k4_3 + k1_0*k2_1*k3_1*k4_3 + k1_0*k2_1*k3_2*k4_2 + k1_2*k2_0*k3_1*k4_2 + k1_0*k2_0*k3_4*k4_2 + k1_0*k2_1*k3_2*k4_3 + k1_0*k2_3*k3_1*k4_2 + k1_2*k2_0*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_2 + k1_0*k2_1*k3_4*k4_2 + k1_0*k2_3*k3_1*k4_3 + k1_2*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_2*k4_2 + k1_0*k2_4*k3_1*k4_3 + k1_2*k2_0*k3_4*k4_2 + k1_3*k2_0*k3_2*k4_3 + k1_3*k2_0*k3_4*k4_2);
S4 = 1 - S0 - S1 - S2 - S3;
% Turnover Rate
JNADH = ETC1Activity*(kfNADH_02*S0 - krNADH_02*S2 ...
+ kfNADH_24*S2 - krNADH_24*S4 ...
+ kfNADH_13*S1 - krNADH_13*S3);
% NADH oxidation releases 1 proton in matrix space
JSO1 = ETC1Activity*(S1*(kfSO_10 - krSO_21) + S2*(kfSO_21 - krSO_32) + S3*(kfSO_32 - krSO_43) + S4*kfSO_43 - krSO_10*S0);
JH2O2 = ETC1Activity*(S2*(kfH2O2_20 - krH2O2_42) + S3*kfH2O2_31 + S4*kfH2O2_42 - krH2O2_31*S1 - krH2O2_20*S0);
J_ETC1_im_to_matrix = JNADH;

%ETC3:im_to_matrix
ETC3_activity = 0.77134e-3;
dPsi = DPsi_im_to_matrix;
Hp = h_im;
Hn = h_matrix;
% Q1 = coQ_matrix;
% QH21 = coQH2_matrix;
Q1 = lambdaQ*coQ_matrix;
QH21 = lambdaQ*coQH2_matrix;
c2 = cytocred_im;
c3 = cytocox_im;
O2 = O2aq_matrix;
SO = SOaq_matrix;
% Set KDQH2 at Qi site
CIII_Kc3 = 1.1193e-06;
CIII_Kc2 = 1.1666e-06;
CIII_KQH2o = 1.0e-3;
CIII_KQo = 0.8e-3;
CIII_KQi = 1.0e-3;
CIII_KQH2i = CIII_KQH2o^2*CIII_KQi*CIII_Kc3^2/CIII_KQo^2/CIII_Kc2^2;
% MR constraint
% Update midpoint potentials from thermodynamic data
CIII_Em0_Q_QH2 = 476.6736;
CIII_Em0_c = 230;
CIII_Em0_SO = 33;
% Update Q Thermodynamics
CIII_Kstabo = 1e-9;
CIII_Kstabi = 0.0078;
CIII_Em0_Q_SQo = CIII_Em0_Q_QH2 + RT/F/2*log(CIII_Kstabo*1e-14);
% mV
CIII_Em0_SQ_QH2o = 2*CIII_Em0_Q_QH2 - CIII_Em0_Q_SQo;
% mV (200 - 300)
CIII_Em0_Q_SQi = CIII_Em0_Q_QH2 + RT/F/2*log(CIII_Kstabi*1e-14);
% mV (assumes Kstabi at pH 7)
CIII_Em0_SQ_QH2i = 2*CIII_Em0_Q_QH2 - CIII_Em0_Q_SQi;
% mV (16-150)
% Binding Polynomials for Protonated States
CIII_pK_ISPox2 = 9.16;
CIII_pK_ISPox1 = 7.63;
CIII_pK_bLox = 5.9;
CIII_pK_bLred = 7.9;
CIII_pK_bHox = 5.7;
CIII_pK_bHred = 7.7;
CIII_pK_QH = 13.2;
CIII_pK_QH2 = 11.3;
CIII_pK_SO = 4.7;
P_ISP = (1 + Hp/10^-CIII_pK_ISPox2 + Hp^2/10^-CIII_pK_ISPox2/10^-CIII_pK_ISPox1);
P_bLox = (1 + Hp/10^-CIII_pK_bLox);
P_bLred = (1 + Hp/10^-CIII_pK_bLred);
P_bHox = (1 + Hn/10^-CIII_pK_bHox);
P_bHred = (1 + Hn/10^-CIII_pK_bHred);
P_QH2o = (1 + Hp/10^-CIII_pK_QH + Hp^2/10^-CIII_pK_QH/10^-CIII_pK_QH2);
P_QH2i = (1 + Hn/10^-CIII_pK_QH + Hn^2/10^-CIII_pK_QH/10^-CIII_pK_QH2);
P_SO = (1 + Hp/10^-CIII_pK_SO);
% Binding Polynomials for enzyme, substrates, products and regulators
% Qo-site
CIII_AA = 0;
CIII_KAA = 1e-10;
P_Qo = (1 + Q1/CIII_KQo + QH21/CIII_KQH2o);
% Cytc c-site
P_c = (1 + c2/CIII_Kc2 + c3/CIII_Kc3);
% Qi-site
P_Qi = (1 + Q1/CIII_KQi + QH21/CIII_KQH2i + CIII_AA/CIII_KAA);
% Midpoint potentials
% Qo-site
CIII_Em0_ISP = 311;
CIII_Em0_bL = 39;
Em_c = CIII_Em0_c;
% pH independent
Em_Q_QH2o = CIII_Em0_Q_QH2 - log(10)*RT/F/2*log10(Hp^2/10^-CIII_pK_QH2/10^-CIII_pK_QH/P_QH2o/Hp^2);
% assuming linked to P-side
Em_Q_SQo = CIII_Em0_Q_SQo;
% pH independent
Em_SQ_QH2o = CIII_Em0_SQ_QH2o - log(10)*RT/F*log10((Hp^2/10^-CIII_pK_QH2/10^-CIII_pK_QH)/P_QH2o/Hp^2);
Em_ISP = CIII_Em0_ISP - log(10)*RT/F*log10(P_ISP/(Hp^2/10^-CIII_pK_ISPox2/10^-CIII_pK_ISPox1));
Em_bL = CIII_Em0_bL - log(10)*RT/F*log10((Hp+10^-CIII_pK_bLox)/(Hp+10^-CIII_pK_bLred));
% assuming linked to P-side (Izrailev et al. 1999, also Crofts)
% Qi-site
CIII_Em0_bH = 160;
Em_bH = CIII_Em0_bH - log(10)*RT/F*log10((Hn+10^-CIII_pK_bHox)/(Hn+10^-CIII_pK_bHred));
% assuming linked to N-side, log10((Hn/10^-CIII_pK_bHred/P_bHred)/(Hn/10^-CIII_pK_bHox/P_bHox))
Em_Q_QH2i = CIII_Em0_Q_QH2 - log(10)*RT/F/2*log10(Hn^2/10^-CIII_pK_QH2/10^-CIII_pK_QH/P_QH2i/Hn^2);
% assuming linked to N-siden
Em_Q_SQi = CIII_Em0_Q_SQi;
% pH independent
Em_SQ_QH2i = CIII_Em0_SQ_QH2i - log(10)*RT/F*log10((Hn^2/10^-CIII_pK_QH2/10^-CIII_pK_QH)/P_QH2i/Hn^2);
% assuming linked to N-side
% Superoxide
Em_SO = CIII_Em0_SO - log(10)*RT/F*log10(1/P_SO/10^-CIII_pK_SO);
% State Transition Thermodynamics
% QH2o<->SQo
dGo = -F*(Em_c - Em_SQ_QH2o);
% dG1 = dGo + RT*log((E2*Q1*c2*Hp^2)/(QH21*c3*E1*1e-14)*(PE1*CIII_KQH2o*CIII_Kc3)/(PE2*CIII_KQo*CIII_Kc2));
% bHred<->bLred
dGi = -F*(Em_bH - Em_bL);
% dG2 = dGi + RT*log((E3)/(E2*Q1)*(PE2*CIII_KQi)/(PE3)) + F*dPsi;
% Substates via Rapid Equilibrium
% ISP protonation state
PISP = Hp/10^-CIII_pK_ISPox2*(1+Hp/10^-CIII_pK_ISPox1)/P_ISP;
% bL reduction by SQo
dGQo = -F*(Em_bL - Em_Q_SQo);
KQo = exp(-dGQo/RT);
% Q-bLred / SQ-bLox
% State 2 fractional species polynomial
% fQoa = 1;
% SQo
rQo = KQo*P_Qo*CIII_KQo/Q1;
% bL_red
PQo = 1 + rQo;
% SQo
% Q reduction by bHred
dGQi = -F*(Em_Q_SQi - Em_bH);
KQi = exp(-dGQi/RT);
% bLox-SQ / bLred-Q
% State 3 fractional species polynomial
% fQia = 1;
% bHred
rQi = KQi*Q1/CIII_KQi/P_Qi;
% bHox-SQi
PQi = 1 + rQi;
% bHred
% SQi reduction by bHred
CIII_beta2 = 0.5;
dGQi2 = -F*(Em_SQ_QH2i - Em_bH) + CIII_beta2*2*F*dPsi;
% protons included in Em_SQ_QH2i
KQi2 = exp(-dGQi2/RT);
% fQi2a = 1;
% bHred-SQi
rQi2 = KQi2*P_Qi*CIII_KQH2i/QH21;
% bHox
PQi2 = 1 + rQi2;
%bHred-SQi
% Superoxide production/consumption rates
% SO thermodynamics
CIII_kSO = 1600;
dG_SO = -F*(Em_SO - Em_Q_SQo);
kSOr = CIII_kSO*exp(dG_SO/RT);
% Reverse Rate Constants
CIII_k120 = 3.0728e+03;
CIII_k230 = 33100000;
CIII_k340 = 3.0728e+03;
CIII_k410 = 14370;
CIII_k450 = 3;
CIII_k520 = 14370;
CIII_k260 = 3.0728e+03;
CIII_k640 = 33100000;
k210 = CIII_k120*exp(dGo/RT);
k320 = CIII_k230*exp(dGi/RT);
k430 = CIII_k340*exp(dGo/RT);
k140 = CIII_k410*exp(dGi/RT);
k540 = CIII_k450*exp(dGo/RT);
k250 = CIII_k520*exp(dGi/RT);
k620 = CIII_k260*exp(dGo/RT);
k460 = CIII_k640*exp(dGi/RT);
% State Transition Rates
% Coupled QH2 oxidation /c3 reduction
kfQH2c3_12 = CIII_k120*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP;
krQH2c3_12 = k210*c2/CIII_Kc2/P_c*PISP/PQo;
kfQH2c3_34 = CIII_k340*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP;
krQH2c3_34 = k430*c2/CIII_Kc2/P_c*PISP/PQo;
kfQH2c3_45 = CIII_k450*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP*rQo/PQo;
krQH2c3_45 = k540*c2/CIII_Kc2/P_c*PISP;
kfQH2c3_26 = CIII_k260*QH21/CIII_KQH2o/P_Qo*c3/CIII_Kc3/P_c*PISP*rQo/PQo;
krQH2c3_26 = k620*c2/CIII_Kc2/P_c*PISP;
% Q uptake at SQ site
CIII_beta1 = 0.5;
kfQ_23 = CIII_k230*rQo/PQo*rQi2/PQi2*exp(-CIII_beta1*F*dPsi/RT/2);
krQ_23 = k320/PQi*exp(CIII_beta1*F*dPsi/RT/2);
kfQ_64 = CIII_k640*rQi2/PQi2*exp(-CIII_beta1*F*dPsi/RT/2);
krQ_64 = k460/PQo/PQi*exp(CIII_beta1*F*dPsi/RT/2);
% QH2 regeneration
kfQH2_41 = CIII_k410*rQo/PQo*rQi/PQi*exp(-CIII_beta1*F*dPsi/RT/2);
krQH2_41 = k140/PQi2*exp(CIII_beta1*F*dPsi/RT/2);
kfQH2_52 = CIII_k520*rQi/PQi*exp(-CIII_beta1*F*dPsi/RT/2);
krQH2_52 = k250/PQo/PQi2*exp(CIII_beta1*F*dPsi/RT/2);
% SO production
kfSO_21 = CIII_kSO*O2/PQo;
krSO_21 = kSOr*SO*Q1/CIII_KQo/P_Qo;
kfSO_43 = CIII_kSO*O2/PQo;
krSO_43 = kSOr*SO*Q1/CIII_KQo/P_Qo;
kfSO_54 = CIII_kSO*O2;
krSO_54 = kSOr*SO*rQo/PQo*Q1/CIII_KQo/P_Qo;
kfSO_62 = CIII_kSO*O2;
krSO_62 = kSOr*SO*rQo/PQo*Q1/CIII_KQo/P_Qo;
% QH2o + c3 -> SQo + c2 + 2H+
k12 = kfQH2c3_12 + krSO_21;
k21 = krQH2c3_12 + kfSO_21;
% bL -> bH
k23 = kfQ_23;
k32 = krQ_23;
% QH2o + c3 -> SQo + c2 + 2H+
k34 = kfQH2c3_34 + krSO_43;
k43 = krQH2c3_34 + kfSO_43;
% bL -> bH
k41 = kfQH2_41;
k14 = krQH2_41;
% QH2o + c3 -> SQo + c2 + 2H+
k45 = kfQH2c3_45 + krSO_54;
k54 = krQH2c3_45 + kfSO_54;
% bL -> bH
k52 = kfQH2_52;
k25 = krQH2_52;
% QH2o + c3 -> SQo + c2 + 2H+
k26 = kfQH2c3_26 + krSO_62;
k62 = krQH2c3_26 + kfSO_62;
% bL -> bH
k64 = kfQ_64;
k46 = krQ_64;
% Steady-State Fractional Occupancies (solved analytically)
E1=(k21*k32*k41*k52*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E2=(k12*k32*k41*k52*k62 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E3=(k12*k23*k41*k52*k62 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k14*k26*k43*k54*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E4=(k14*k21*k32*k52*k62 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k14*k26*k34*k54*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E5=(k12*k25*k32*k41*k62 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k14*k26*k34*k45*k64)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
% E6=(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k14*k26*k34*k46*k54)/(k12*k26*k32*k41*k52 + k12*k26*k32*k41*k54 + k12*k26*k32*k43*k52 + k12*k26*k34*k41*k52 + k14*k21*k32*k46*k52 + k12*k23*k34*k46*k52 + k12*k26*k32*k43*k54 + k12*k26*k32*k45*k52 + k12*k26*k34*k41*k54 + k14*k21*k32*k46*k54 + k14*k21*k34*k46*k52 + k14*k26*k32*k43*k52 + k12*k26*k32*k46*k52 + k12*k23*k34*k46*k54 + k12*k25*k32*k46*k54 + k12*k26*k34*k45*k52 + k14*k21*k34*k46*k54 + k14*k23*k34*k46*k52 + k14*k26*k32*k43*k54 + k14*k26*k32*k45*k52 + k12*k26*k32*k46*k54 + k12*k26*k34*k46*k52 + k14*k26*k32*k46*k52 + k12*k25*k34*k46*k54 + k14*k23*k34*k46*k54 + k14*k25*k32*k46*k54 + k14*k26*k34*k45*k52 + k12*k25*k32*k41*k62 + k12*k26*k34*k46*k54 + k14*k26*k32*k46*k54 + k14*k26*k34*k46*k52 + k14*k25*k34*k46*k54 + k12*k25*k32*k41*k64 + k12*k25*k32*k43*k62 + k12*k25*k34*k41*k62 + k14*k21*k32*k45*k62 + k14*k26*k34*k46*k54 + k12*k23*k34*k45*k62 + k12*k25*k32*k43*k64 + k12*k25*k32*k45*k62 + k12*k25*k34*k41*k64 + k14*k21*k32*k45*k64 + k14*k21*k34*k45*k62 + k14*k25*k32*k43*k62 + k12*k25*k32*k46*k62 + k12*k23*k34*k45*k64 + k12*k25*k32*k45*k64 + k12*k25*k34*k45*k62 + k14*k21*k34*k45*k64 + k14*k23*k34*k45*k62 + k14*k25*k32*k43*k64 + k14*k25*k32*k45*k62 + k12*k25*k34*k46*k62 + k12*k26*k32*k45*k64 + k14*k25*k32*k46*k62 + k12*k25*k34*k45*k64 + k14*k23*k34*k45*k64 + k14*k25*k32*k45*k64 + k14*k25*k34*k45*k62 + k12*k26*k34*k45*k64 + k14*k21*k32*k52*k62 + k14*k25*k34*k46*k62 + k14*k26*k32*k45*k64 + k14*k25*k34*k45*k64 + k12*k23*k34*k52*k62 + k14*k21*k32*k52*k64 + k14*k21*k32*k54*k62 + k14*k21*k34*k52*k62 + k14*k26*k34*k45*k64 + k12*k23*k34*k52*k64 + k12*k23*k34*k54*k62 + k12*k25*k32*k54*k62 + k14*k21*k32*k54*k64 + k14*k21*k34*k52*k64 + k14*k21*k34*k54*k62 + k14*k23*k34*k52*k62 + k12*k26*k32*k52*k64 + k12*k23*k34*k54*k64 + k12*k25*k32*k54*k64 + k12*k25*k34*k54*k62 + k14*k21*k34*k54*k64 + k14*k23*k34*k52*k64 + k14*k23*k34*k54*k62 + k14*k25*k32*k54*k62 + k12*k26*k32*k54*k64 + k12*k26*k34*k52*k64 + k14*k26*k32*k52*k64 + k12*k25*k34*k54*k64 + k14*k23*k34*k54*k64 + k14*k25*k32*k54*k64 + k14*k25*k34*k54*k62 + k12*k23*k41*k52*k62 + k12*k26*k34*k54*k64 + k14*k26*k32*k54*k64 + k14*k26*k34*k52*k64 + k14*k25*k34*k54*k64 + k12*k23*k41*k52*k64 + k12*k23*k41*k54*k62 + k12*k23*k43*k52*k62 + k14*k21*k43*k52*k62 + k14*k26*k34*k54*k64 + k12*k23*k41*k54*k64 + k12*k23*k43*k52*k64 + k12*k23*k43*k54*k62 + k12*k23*k45*k52*k62 + k14*k21*k43*k52*k64 + k14*k21*k43*k54*k62 + k14*k23*k43*k52*k62 + k12*k23*k46*k52*k62 + k12*k23*k43*k54*k64 + k12*k23*k45*k52*k64 + k12*k25*k43*k54*k62 + k14*k21*k43*k54*k64 + k14*k23*k43*k52*k64 + k14*k23*k43*k54*k62 + k14*k23*k45*k52*k62 + k12*k23*k46*k54*k62 + k12*k26*k43*k52*k64 + k14*k23*k46*k52*k62 + k12*k25*k43*k54*k64 + k14*k23*k43*k54*k64 + k14*k23*k45*k52*k64 + k14*k25*k43*k54*k62 + k12*k26*k43*k54*k64 + k12*k32*k41*k52*k62 + k14*k23*k46*k54*k62 + k14*k26*k43*k52*k64 + k14*k25*k43*k54*k64 + k12*k32*k41*k52*k64 + k12*k32*k41*k54*k62 + k12*k32*k43*k52*k62 + k12*k34*k41*k52*k62 + k14*k26*k43*k54*k64 + k12*k32*k41*k54*k64 + k12*k32*k43*k52*k64 + k12*k32*k43*k54*k62 + k12*k32*k45*k52*k62 + k12*k34*k41*k52*k64 + k12*k34*k41*k54*k62 + k14*k32*k43*k52*k62 + k12*k32*k46*k52*k62 + k12*k32*k43*k54*k64 + k12*k32*k45*k52*k64 + k12*k34*k41*k54*k64 + k12*k34*k45*k52*k62 + k14*k32*k43*k52*k64 + k14*k32*k43*k54*k62 + k14*k32*k45*k52*k62 + k12*k32*k46*k54*k62 + k12*k34*k46*k52*k62 + k14*k32*k46*k52*k62 + k12*k34*k45*k52*k64 + k14*k32*k43*k54*k64 + k14*k32*k45*k52*k64 + k14*k34*k45*k52*k62 + k12*k34*k46*k54*k62 + k14*k32*k46*k54*k62 + k14*k34*k46*k52*k62 + k21*k32*k41*k52*k62 + k14*k34*k45*k52*k64 + k14*k34*k46*k54*k62 + k21*k32*k41*k52*k64 + k21*k32*k41*k54*k62 + k21*k32*k43*k52*k62 + k21*k34*k41*k52*k62 + k21*k32*k41*k54*k64 + k21*k32*k43*k52*k64 + k21*k32*k43*k54*k62 + k21*k32*k45*k52*k62 + k21*k34*k41*k52*k64 + k21*k34*k41*k54*k62 + k23*k34*k41*k52*k62 + k21*k32*k46*k52*k62 + k21*k32*k43*k54*k64 + k21*k32*k45*k52*k64 + k21*k34*k41*k54*k64 + k21*k34*k45*k52*k62 + k23*k34*k41*k52*k64 + k23*k34*k41*k54*k62 + k25*k32*k41*k54*k62 + k21*k32*k46*k54*k62 + k21*k34*k46*k52*k62 + k26*k32*k41*k52*k64 + k21*k34*k45*k52*k64 + k23*k34*k41*k54*k64 + k25*k32*k41*k54*k64 + k25*k34*k41*k54*k62 + k21*k34*k46*k54*k62 + k26*k32*k41*k54*k64 + k26*k34*k41*k52*k64 + k25*k34*k41*k54*k64 + k26*k34*k41*k54*k64);
E6 = 1 - E1 - E2 - E3 - E4 - E5;
% Net Turnover Rates
% QH2 consumption
JQH2 = ETC3_activity*((kfQH2c3_12*E1 - krQH2c3_12*E2) + (kfQH2c3_34*E3 - krQH2c3_34*E4) + (kfQH2c3_45*E4 - krQH2c3_45*E5) + (kfQH2c3_26*E2 - krQH2c3_26*E6) - (kfQH2_41*E4 - krQH2_41*E1) - (kfQH2_52*E5 - krQH2_52*E2));
% superoxide production
JSO3 = ETC3_activity*((kfSO_21*E2 - krSO_21*E1) + (kfSO_43*E4 - krSO_43*E3) + (kfSO_54*E5 - krSO_54*E4) + (kfSO_62*E6 - krSO_62*E2)); 
if QH21 < (1e-12)
  JSO3 = 0;
  JQH2 = 0;
end
% for numerical stability, enforce strict stoichiometric coupling
J_ETC3_im_to_matrix = (2*JQH2 - JSO3);

%ETC4:im_to_matrix
Keq_ETC4 = exp(-DGro_ETC4/RT)*exp( -(4*F*DPsi_im_to_matrix)/RT)*P(32)^2/P(33)^2/P(34)^0.5*h_matrix^4/h_im^2;
ETC4_activity = 3.25e-3 ;     
c3 = cytocox_im;
c2 = cytocred_im;
O2 = O2aq_matrix;
dPsi = DPsi_im_to_matrix;
% mV
KO2 = 1e-6;
KM = (1.6115e-04); 
beta = 6.6054e-06; % unitless
nc4 = 2; % Hill coefficient for cytoc red
ter_oxg = (O2/(O2+KO2));
ter1 = (c2^nc4/(c2^nc4 + KM^nc4*(1 + beta*exp(2*dPsi*F/RT))));
ter2 = (1 - (c3/c2)^2/(O2^0.5*Keq_ETC4)); % This is how this term is coded in Bazil code.
J_ETC4_im_to_matrix = ETC4_activity*ter_oxg*ter1*ter2;
if O2 <= 1e-12 || c2 < 1e-9
  J_ETC4_im_to_matrix = 0;
end

% GHK H LEAK 
x_HLE = HLE*5571.8 ;
FRT = DPsi_im_to_matrix*F/RT;
if abs(FRT) > 0.01
  J_HLEAK_im_to_matrix = x_HLE*FRT*(h_im*exp(FRT) - h_matrix)/(exp(FRT)-1);
else
  J_HLEAK_im_to_matrix = x_HLE*(h_im - h_matrix);
end

% GHK K LEAK 
x_KLE = 1.1e-3;
FRT = DPsi_im_to_matrix*F/RT;
if abs(FRT) > 0.01
  J_KLEAK_im_to_matrix = x_KLE*FRT*(k_im*exp(FRT) - k_matrix)/(exp(FRT)-1);
else
  J_KLEAK_im_to_matrix = x_KLE*(k_im - k_matrix);
end

%KH:im_to_matrix
k1_KH = 1.2889e+06;
J_KH_im_to_matrix = k1_KH*(k_im*h_matrix - k_matrix*h_im);

%NH4/H:im_to_matri NH3 transport
J_NH3_im_to_matrix = 83030*(0 - ammonia_matrix/P(26));

% PIH:im_to_matrix
%  based on asymmetric SUC/MAL/Pi model
Kp1 = 1.8e-3; % phosphate (outside) affinity
Kh1 = 0.5e-7; % H (outside) affinity
Kp2 = Kp1; % phosphate (inside) affinity
Kh2 = Kh1; % 
p1 = Pi_im*(h_im/Kh(36))/P(36)/Kp1;
p2 = Pi_matrix*(h_matrix/Kh(14))/P(14)/Kp2;
h1 = h_im/Kh1;
h2 = h_matrix/Kh2;
p1h1 = p1*h1;
p2h2 = p2*h2;
d1 = 1  + p1 + h1 + p1h1;
d2 = 1 +  p2 + h2 + p2h2 ;
x_PIH = 0.25797 ; 
gamma = 1; % 7/25
beta  = 1;
alpha = Kp2*Kh2/(Kp1*Kh1)/gamma;
J_PIH_im_to_matrix = x_PIH*(p1h1 - alpha*beta*p2h2 )/( d1*(gamma + alpha*p2*h2) + d2*(gamma*beta + p1*h1) );

%ANT:im_to_matrix
x_ANT = 0.50092;
ADP_i1 = ADP_im/P(37);
% ADP^3-;
ATP_i1 = ATP_im/P(35);
% ATP^4-;
ADP_x1 = ADP_matrix/P(11);
% ADP^3-;
ATP_x1 = ATP_matrix/P(10);
% ATP^4-;
% F = 0.096485;
del_D = 0.0167;
del_T = 0.0699;
k2_ANT = 9.54/60;
% = 1.59e-1
k3_ANT = 30.05/60;
% = 5.01e-1
K_D_o_ANT = 38.89e-6;
K_T_o_ANT = 56.05e-6;
A = +0.2829;
B = -0.2086;
C = +0.2372;
fi = F*DPsi_im_to_matrix/RT;
k2_ANT_fi = k2_ANT*exp((A*(-3)+B*(-4)+C)*fi);
k3_ANT_fi = k3_ANT*exp((A*(-4)+B*(-3)+C)*fi);
K_D_o_ANT_fi = K_D_o_ANT*exp(3*del_D*fi);
K_T_o_ANT_fi = K_T_o_ANT*exp(4*del_T*fi);
q = k3_ANT_fi*K_D_o_ANT_fi*exp(fi)/(k2_ANT_fi*K_T_o_ANT_fi);
term2 = k2_ANT_fi*ATP_x1*ADP_i1*q/K_D_o_ANT_fi ;
term3 = k3_ANT_fi*ADP_x1*ATP_i1/K_T_o_ANT_fi;
num = term2 - term3;
den = (1 + ATP_i1/K_T_o_ANT_fi + ADP_i1/K_D_o_ANT_fi)*(ADP_x1 + ATP_x1*q);
J_ANT_im_to_matrix = x_ANT*num/den;

%PYRH:im_to_matrix DAB ORDERED Bi-Bi
% based on Vinnakota & Beard DOI: 10.1016/j.bpj.2010.11.079 
x_PYRH = 0.047703;
Kpyr = 7.1631e-05; % 1/26
Khyd   = 10^(-7.0); % 
gamma_pyr = 2.3;
a1 = pyruvate_matrix/P(1)/Kpyr;
b1 = h_matrix/Khyd;
a2 = pyruvate_im/P(38)/Kpyr;
b2 = h_im/Khyd;
den = (a1*b1 + gamma_pyr)*(1 + b2 + a2*b2) + (a2*b2 + gamma_pyr)*(1 + b1 + a1*b1);
J_PYRH_im_to_matrix = x_PYRH*(a2*b2 - a1*b1) / den;

%AAT_matrix--Saito et al. model https://doi.org/10.1113/JP272598
a=aspartate_matrix;
b=ketoglutarate_matrix;
p=oxaloacetate_matrix;
q=glutamate_matrix;
KmA = 1.58*1e-3;
KmB = 0.149*1e-3;
KmP = 0.0399*1e-3;
KmQ = 2.5*1e-3;
Kia = 2.0*1e-3;
Kiq = 1.83*1e-3;
Vmf_GOT=0.5;
% Vmf_GOT=1;
if(a>MinCon)&&(b>MinCon)
  ab=a*b;
else
  ab=0;
end
if(p>MinCon)&&(q>MinCon)
  pq=p*q;
else
  pq=0;
end
% Keq_AAT_matrix = 0.1; % DAB fudge
% J_AAT_matrix = Vmf_GOT*(ab-pq/Keq_AAT_matrix) / ( alpha*KmB*a + alpha*KmA*b + alpha*ab + ...
%     (KmQ*p + KmP*q)/(Keq_AAT_matrix) + (pq + KmQ*a*p/Kia)/(Keq_AAT_matrix) + ...
%     alpha*KmA*b*q/Kiq );
alpha = 0.30/0.58; % KcF/KcR
J_AAT_matrix = Vmf_GOT*(ab-pq/Keq_AAT_matrix) / ( KmB*a + KmA*b + ab + ...
    alpha*(KmQ*p + KmP*q + pq + KmQ*a*p/Kia)/(Keq_AAT_matrix) + KmA*b*q/Kiq );

% %AAT_matrix--Fan Wu model implemented by Santosh
% a=aspartate_matrix;
% b=ketoglutarate_matrix;
% p=oxaloacetate_matrix;
% q=glutamate_matrix;
% KmA=3900e-6;
% KmB=430e-6;
% KmP=88e-6;
% KmQ=8900e-6;
% Kia=3480e-6;
% % Kib=710e-6;
% Kip=50e-6;
% Kiq=8400e-6;
% KiAKG=16.6e-3;
% ai=1+b/KiAKG;
% Vmf_GOT=0.12935;
% Vmr=Vmf_GOT/Keq_AAT_matrix*(KmQ*Kip/Kia/KmB);
% 
% if(a>MinCon)&&(b>MinCon)
%   ab=a*b;
% else
%   ab=0;
% end
% if(p>MinCon)&&(q>MinCon)
%   pq=p*q;
% else
%   pq=0;
% end
% J_AAT_matrix=Vmf_GOT*Vmr*(ab-pq/Keq_AAT_matrix)/ (Vmr*KmB*a+Vmr*KmA*ai*b+Vmf_GOT*KmQ/Keq_AAT_matrix*ai*p+Vmf_GOT*KmP/Keq_AAT_matrix*q+Vmr*a*b+Vmf_GOT*KmQ/Kia/Keq_AAT_matrix*a*p+Vmf_GOT/Keq_AAT_matrix*p*q+Vmr*KmA/Kiq*b*q);

% GDH1_matrix
a = glutamate_matrix;
b = NAD_matrix;
r = NADH_free;
p = ketoglutarate_matrix;
q = ammonia_matrix; 
KmA=3.5e-3;
KmB=80e-6;
KmP=1.1e-3;
KmR=40e-6;
KiA=3.5e-3;
KiB=1e-3;
KiP=0.25e-3;
KmQ=6e-3;
KiQ=6e-3;
KiR=4e-6;
% Keq_GDH_matrix=3e-6;
Vmf_GDH1 = 0.84903 ;
Den=KiB*KmA+KmA*b+KmB*a+a*b+a*b*q/KiP+KiB*KmA*r/KiR+...
KmA*b*q/KiQ+KiB*KmA*KmR*p*q/KmQ/KiP/KiR+KmB*a*r/KiR+a*b*p*q/KiP/KiQ+...
KiB*KmA*KmP*q*r/KmQ/KiP/KiR+KiB*KmA*KmP*p*q*r/KmQ/KiP/KiR+a*b*p/KiP+...
KiB*KmA*KmP*q/KmQ/KiP+KiB*KmA*a*p*r/KiA/KiP/KiR+KiB*KmA*p*r/KiP/KiR+...
KmR*KmA*p*q*b/KmQ/KiP/KiR+KiB*KmA*a*p*q*r/KmQ/KiA/KiP/KiR;
J_GDH1_matrix = Vmf_GDH1*(a*b-p*q*r/Keq_GDH_matrix)/Den;
if a <= 0
  J_GDH1_matrix = min(J_GDH1_matrix,0);
end
if q <= 0
  J_GDH1_matrix = max(J_GDH1_matrix,0);
end

% ASPGLU
Keq_ASPGLU = exp( +(1*F*DPsi_im_to_matrix)/RT)/P(24)*P(25)/P(39)*P(40);
Kg1 = 1.0718e-3; % glutamate (outside) affinity
Ka1 = 0.69225e-3; % aspartate (outside) affinity
Kg2 = 0.94083e-3; % glutamate (inside) affinity
Ka2 = Kg2*Ka1/Kg1; % thermo constraint
g1 = glutamate_im*h_im/Kh(39)/P(39)/Kg1;
g2 = glutamate_matrix*h_matrix/Kh(25)/P(25)/Kg2;
a1 = aspartate_im/P(40)/Ka1;
a2 = aspartate_matrix/P(24)/Ka2;
d1 = 1 + g1 + a1;
d2 = 1 + g2 + a2;
x_ASPGLU = 0.10;
if (g1 + a1 >1e-6)||(g2 + a2 >1e-6)
  J_ASPGLU_im_to_matrix = x_ASPGLU*(g1*a2*Keq_ASPGLU - g2*a1 )/( d1*(g2+a2) + d2*(g1 + a1) );
else
  J_ASPGLU_im_to_matrix = 0;
end

%GLUH:im_to_matrix
x_GLUH = 7.4769e+05;
glutamate1 = glutamate_matrix/P(25);
glutamate2 = glutamate_im/P(39);
J_GLUH_im_to_matrix = x_GLUH * (glutamate2*h_im - glutamate1*h_matrix);

%CITMAL:im_to_matrix
J_CITMAL_im_to_matrix = 0;

% AKGMAL:im_to_matrix
alpha = 1.0;
Km1 = 1.0e-3; % malate (outside) affinity
Ka1 = 1.0e-3; % akg (outside) affinity
Km2 = alpha*Km1; % thermo constraint
Ka2 = alpha*Ka1; % thermo constraint
m1 = malate_im/P(42)/Km1;
m2 = malate_matrix/P(23)/Km2;
a1 = ketoglutarate_im/P(43)/Ka1;
a2 = ketoglutarate_matrix/P(16)/Ka2;
d1 = 1 + m1 + a1;
d2 = 1 + m2 + a2;
x_AKGMAL = 0.0070174;
J_AKGMAL_im_to_matrix = -x_AKGMAL*(m1*a2 - m2*a1 )/( d1*(m2+a2) + d2*(m1 + a1) );
  
% Dicarboxylate transporter DAB August 2025
alpha = 1;
Km1 = 3.4183e-3; % im malate affinity (adjustable)
Ks1 = 3.1379e-3; % im succinate affinity (adjustable)
Kp1 = 6.9112e-3; % im Pi affinity (adjustable)
Kp2 = alpha*Kp1;  % thermo constraint
Ks2 = alpha*Ks1;  % thermo constraint
Km2 = alpha*Km1; % thermo constraint
Ko2 = 1.8716e-6; % g1
Kf2 = 1.4062e-3; % 

m1 = malate_im/P(42)/Km1;
m2 = malate_matrix/P(23)/Km2;
p1 = Pi_im/P(36)/Kp1;
p2 = Pi_matrix/P(14)/Kp2;
s1 = succinate_im/P(44)/Ks1;
s2 = succinate_matrix/P(9)/Ks2;
o2 = oxaloacetate_matrix/P(12)/Ko2;
f2 = fumarate_matrix/P(22)/Kf2;
d1 = 1 + s1 + m1 + p1;  
d2 = 1 + s2 + m2 + p2 + o2 + f2;
D = alpha*d1*(s2 + m2 + p2) + d2*(s1 + m1 + p1);

x_DICARBOX = 8.2904e-3;
J_DICARBPI_im_to_matrix  = x_DICARBOX*alpha*(p1*m2 + p1*s2 - p2*m1 - p2*s1)/D;
J_DICARBMAL_im_to_matrix = x_DICARBOX*alpha*(m1*s2 + m1*p2 - m2*s1 - m2*p1)/D;
J_DICARBSUC_im_to_matrix = x_DICARBOX*alpha*(s1*m2 + s1*p2 - s2*m1 - s2*p1)/D;
J_DICARBOAA_im_to_matrix = 0;

%O2PERM:im_to_matrix
x_O2perm = 1000;
J_O2PERM_im_to_matrix = x_O2perm*(O2aq_im - O2aq_matrix);
%ADPPERM:buffer_to_im
x_ADPPERM = 400;
J_ADPPERM_buffer_to_im = x_ADPPERM * (ADP_buffer - ADP_im);
%ATPPERM:buffer_to_im
x_ATPPERM = 400;
J_ATPPERM_buffer_to_im = x_ATPPERM * (ATP_buffer - ATP_im);
%PIPERM:buffer_to_im
x_PIPERM = 400;
J_PIPERM_buffer_to_im = x_PIPERM * (Pi_buffer - Pi_im);
%HPERM:buffer_to_im
x_HPERM = 1e8;
J_HPERM_buffer_to_im = x_HPERM * (h_buffer - h_im);
%KPERM:buffer_to_im
x_KPERM = 400;
J_KPERM_buffer_to_im = x_KPERM * (k_buffer - k_im);
%MPERM:buffer_to_im
x_MPERM = 400;
J_MPERM_buffer_to_im = x_MPERM * (m_buffer - m_im);
%PYRPERM:buffer_to_im
x_PYRPERM = 400;
J_PYRPERM_buffer_to_im = x_PYRPERM * (pyruvate_buffer - pyruvate_im);
%CITPERM:buffer_to_im
x_CITPERM = 400;
J_CITPERM_buffer_to_im = x_CITPERM * (citrate_buffer - citrate_im);
%AKGPERM:buffer_to_im
x_AKGPERM = 400;
J_AKGPERM_buffer_to_im = x_AKGPERM * (ketoglutarate_buffer - ketoglutarate_im);
%SUCPERM:buffer_to_im
x_SUCPERM = 400;
J_SUCPERM_buffer_to_im = x_SUCPERM * (succinate_buffer - succinate_im);
%GLUPERM:buffer_to_im
x_GLUPERM = 400;
J_GLUPERM_buffer_to_im = x_GLUPERM * (glutamate_buffer - glutamate_im);
%ASPPERM:buffer_to_im
x_ASPPERM = 400;
J_ASPPERM_buffer_to_im = x_ASPPERM * (aspartate_buffer - aspartate_im);
%MALPERM:buffer_to_im
x_MALPERM = 400;
J_MALPERM_buffer_to_im = x_MALPERM * (malate_buffer - malate_im);
%O2PERM2:buffer_to_im
x_O2perm2 = 1000;
J_O2PERM2_buffer_to_im = x_O2perm2*(O2aq_buffer - O2aq_im);

%% REACTANT TIME DERIVATIVES
f(1,:) = ( 0  - 1*J_PDH_matrix + 1*J_PYRH_im_to_matrix + J_OD + J_MALIC ) / VWater_matrix; % pyruvate_matrix
f(2,:) = ( 0  - 1*J_PDH_matrix + 1*J_CTS_matrix - 1*J_AKGDH_matrix + 1*J_SCS_matrix ) / VWater_matrix; % COAS_matrix
f(3,:) = ( 0  - 1*J_PDH_matrix - 1*J_IDH_matrix - 1*J_AKGDH_matrix - 1*J_MDH_matrix - J_GDH1_matrix + 1*J_ETC1_im_to_matrix + 1*J_MPTT - J_MALIC) / VWater_matrix; % NAD_matrix 
f(4,:) = ( 0  ); % [clamped] % CO2tot_matrix
f(5,:) = ( 0  + 1*J_PDH_matrix - 1*J_CTS_matrix ) / VWater_matrix; % acetylcoA_matrix
f(6,:) = ( 0  + 1*J_PDH_matrix + 1*J_IDH_matrix + 1*J_AKGDH_matrix + 1*J_MDH_matrix + J_GDH1_matrix - 1*J_ETC1_im_to_matrix - 1*J_MPTT + J_MALIC ) / VWater_matrix; % NADH_matrix 
f(7,:) = ( 0  + 1*J_AKGDH_matrix - 1*J_SCS_matrix ) / VWater_matrix; % succinylcoA_matrix
f(8,:) = ( 0  + 1*J_CTS_matrix - 1*J_ACON_matrix + 1*J_CITMAL_im_to_matrix ) / VWater_matrix; % citrate_matrix
f(9,:) = ( 0  + 1*J_SCS_matrix - 1*J_SDH_matrix + 1*J_DICARBSUC_im_to_matrix) / VWater_matrix; % succinate_matrix 
f(10,:) = ( 0  + 1*J_NDK_matrix + 1*J_F1F0ATPASE_im_to_matrix - 1*J_ANT_im_to_matrix ) / VWater_matrix; % ATP_matrix
f(11,:) = ( 0  - 1*J_NDK_matrix - 1*J_F1F0ATPASE_im_to_matrix + 1*J_ANT_im_to_matrix ) / VWater_matrix; % ADP_matrix
B = 0.25e-3;
K = 1.0e-6;
beta = 1 + B*K/(K + oxaloacetate_matrix)^2;
f(12,:) = ( 0  - 1*J_CTS_matrix + 1*J_MDH_matrix + 1*J_AAT_matrix - J_OD + J_DICARBOAA_im_to_matrix) / VWater_matrix / beta ; % oxaloacetate_matrix
f(13,:) = ( 0  ) / VWater_matrix; % AMP_matrix
f(14,:) = ( 0  - 1*J_SCS_matrix - 1*J_F1F0ATPASE_im_to_matrix + 1*J_PIH_im_to_matrix + 1*J_DICARBPI_im_to_matrix ) / VWater_matrix; % Pi_matrix; 
f(15,:) = ( 0  + 1*J_ACON_matrix - 1*J_IDH_matrix ) / VWater_matrix; % isocitrate_matrix
f(16,:) = ( 0  + 1*J_IDH_matrix - 1*J_AKGDH_matrix - 1*J_AAT_matrix + 1*J_GDH1_matrix + 1*J_AKGMAL_im_to_matrix ) / VWater_matrix; % ketoglutarate_matrix
f(17,:) = ( 0  ) / VWater_matrix; % C_matrix
f(18,:) = ( 0  - 1*J_SCS_matrix + 1*J_NDK_matrix ) / VWater_matrix; % GDP_matrix
f(19,:) = ( 0  + 1*J_SCS_matrix - 1*J_NDK_matrix ) / VWater_matrix; % GTP_matrix
f(20,:) = ( 0  - 1*J_SDH_matrix - 1*J_ETC1_im_to_matrix + 1*J_ETC3_im_to_matrix  ); % coQ_matrix (mole / L mito)
f(21,:) = ( 0  + 1*J_SDH_matrix + 1*J_ETC1_im_to_matrix - 1*J_ETC3_im_to_matrix  ); % coQH2_matrix (mole / L mito)
f(22,:) = ( 0  + 1*J_SDH_matrix - 1*J_FUM_matrix ) / VWater_matrix; % fumarate_matrix
f(23,:) = ( 0  + 1*J_FUM_matrix - 1*J_MDH_matrix - J_MALIC - 1*J_CITMAL_im_to_matrix - 1*J_AKGMAL_im_to_matrix + 1*J_DICARBMAL_im_to_matrix ) / VWater_matrix; % malate_matrix
f(24,:) = ( 0  - 1*J_AAT_matrix - 1*J_ASPGLU_im_to_matrix ) / VWater_matrix; % aspartate_matrix
f(25,:) = ( 0  + 1*J_AAT_matrix - 1*J_GDH1_matrix + 1*J_ASPGLU_im_to_matrix + 1*J_GLUH_im_to_matrix ) / VWater_matrix; % glutamate_matrix
f(26,:) = ( + 1*J_GDH1_matrix + 1*J_NH3_im_to_matrix ) / VWater_matrix; % ammonia_matrix
f(30,:) = ( 0  ); % [clamped] % H2O2aq_matrix
f(31,:) = ( 0  ); % [clamped] % SOaq_matrix
f(27,:) = ( 0  - 1*J_ATPASE_buffer - 1*J_ATPPERM_buffer_to_im/VRegion_buffer*VRegion_im + J_AK_buffer) / VWater_buffer; % ATP_buffer
f(28,:) = ( 0  + 1*J_ATPASE_buffer - 1*J_ADPPERM_buffer_to_im/VRegion_buffer*VRegion_im - 2*J_AK_buffer) / VWater_buffer; % ADP_buffer
f(29,:) = ( 0  + 1*J_ATPASE_buffer - 1*J_PIPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % Pi_buffer
f(46,:) = ( 0  - 1*J_PYRPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % pyruvate_buffer
f(47,:) = ( 0  - 1*J_CITPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % citrate_buffer
f(48,:) = ( 0  - 1*J_AKGPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % ketoglutarate_buffer
f(49,:) = ( 0  - 1*J_SUCPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % succinate_buffer
f(50,:) = ( 0  - 1*J_GLUPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % glutamate_buffer
f(51,:) = ( 0  - 1*J_ASPPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % aspartate_buffer
f(52,:) = ( 0  - 1*J_MALPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % malate_buffer
f(32,:) = ( 0  - 2*J_ETC3_im_to_matrix/VRegion_im*VRegion_matrix + 2*J_ETC4_im_to_matrix/VRegion_im*VRegion_matrix ) / VWater_im; % cytocox_im 
f(33,:) = ( 0  + 2*J_ETC3_im_to_matrix/VRegion_im*VRegion_matrix - 2*J_ETC4_im_to_matrix/VRegion_im*VRegion_matrix ) / VWater_im; % cytocred_im 
f(35,:) = ( 0  + 1*J_ANT_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_ATPPERM_buffer_to_im ) / VWater_im; % ATP_im
f(36,:) = ( 0  - 1*J_PIH_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_DICARBPI_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_PIPERM_buffer_to_im ) / VWater_im; % Pi_im 
f(37,:) = ( 0  - 1*J_ANT_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_ADPPERM_buffer_to_im ) / VWater_im; % ADP_im
f(38,:) = ( 0  - 1*J_PYRH_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_PYRPERM_buffer_to_im ) / VWater_im; % pyruvate_im
f(39,:) = ( 0  - 1*J_ASPGLU_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_GLUH_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_GLUPERM_buffer_to_im ) / VWater_im; % glutamate_im
f(40,:) = ( 0  + 1*J_ASPGLU_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_ASPPERM_buffer_to_im ) / VWater_im; % aspartate_im
f(41,:) = ( 0  - 1*J_CITMAL_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_CITPERM_buffer_to_im ) / VWater_im; % citrate_im
f(42,:) = ( 0  + 1*J_CITMAL_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_AKGMAL_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_DICARBMAL_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_MALPERM_buffer_to_im ) / VWater_im; % malate_im
f(43,:) = ( 0  - 1*J_AKGMAL_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_AKGPERM_buffer_to_im ) / VWater_im; % ketoglutarate_im
f(44,:) = ( 0  + 1*J_SUCPERM_buffer_to_im - 1*J_DICARBSUC_im_to_matrix/VRegion_im*VRegion_matrix ) / VWater_im; % succinate_im  

% oxygen kinetics
if Pflag
  f(34,:) = 0; % O2aq_matrix
  f(53,:) = 0; % O2aq_buffer
  f(45,:) = 0; % O2aq_im
else
  f(34,:) = ( 0  - 0.5*J_ETC4_im_to_matrix + 1*J_O2PERM_im_to_matrix - 0*JSO1 - 0*JSO3 ) / VWater_matrix; % O2aq_matrix
  f(53,:) = ( 0  - 1*J_O2PERM2_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer; % O2aq_buffer
  f(45,:) = ( 0  - 1*J_O2PERM_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_O2PERM2_buffer_to_im ) / VWater_im; % O2aq_im
end

% added (DAB) for NADP(H)
f(66,:) = ( 0  - 1*J_MPTT ) / VWater_matrix; % NADP_matrix
f(67,:) = ( 0  + 1*J_MPTT ) / VWater_matrix; % NADPH_matrix

% added (NLC) for AK reaction
f(68,:) = ( + J_AK_buffer) / VWater_buffer; % AMP_buffer

% HLeak kinetics
J_ROS = JSO1 + JH2O2 + JSO3 + J_SDH_SO + J_SDH_H2O2;
J0 = 0.5e-4;
J1 = 0.5e-4;
tau_HLE1 = 5;
tau_HLE2 = 469.18;
HLE_inf = 0.5*(1 + tanh((J_ROS - J0)/J1) );
if HLE_inf > HLE 
  f(69,:) = ( HLE_inf - HLE ) / tau_HLE1 ;
else
  f(69,:) = ( HLE_inf - HLE ) / tau_HLE2 ;
end


%% ION EQUATIONS
% COMPARTMENT matrix:
ii = [1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  30  31  34]; % Indices of SVs in compartment matrix
% PARTIAL DERIVATIVES
pHBpK = -sum( (h_matrix*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpM = -sum( (h_matrix*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_matrix./Km(ii)+k_matrix./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_matrix*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_matrix*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_matrix./Kh(ii)+k_matrix./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_matrix*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_matrix*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_matrix./Kh(ii)+m_matrix./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 - 1*J_PDH_matrix + 2*J_CTS_matrix + 0*J_ACON_matrix + 0*J_IDH_matrix - 1*J_AKGDH_matrix + 1*J_SCS_matrix + 0*J_SDH_matrix + 0*J_FUM_matrix + 1*J_MDH_matrix + 0*J_NDK_matrix + 0*J_AAT_matrix + 2*J_GDH1_matrix + 1.66667*J_F1F0ATPASE_im_to_matrix + 0*J_F1F0ATPASE_im_to_matrix - 5*J_ETC1_im_to_matrix - 2*J_ETC3_im_to_matrix  - 4*J_ETC4_im_to_matrix + 1*J_HLEAK_im_to_matrix + 2*J_PIH_im_to_matrix - 1*J_KH_im_to_matrix + 1*J_PYRH_im_to_matrix + 1*J_GLUH_im_to_matrix + 1*J_ASPGLU_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_CITMAL_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_OD + 1*J_MPTT ) / VWater_matrix; 

J_M = (0) / VWater_matrix;
J_K = (0 + 1*J_KH_im_to_matrix + 1*J_KLEAK_im_to_matrix) / VWater_matrix;
Phi_H = J_H - sum( h_matrix*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = J_M - sum( m_matrix*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K - sum( k_matrix*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(1)/K_BX(1)/(1+h_matrix/K_BX(1))^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(54,:) =  ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
f(55,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(56,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;
% COMPARTMENT buffer:
ii = [27  28  29  46  47  48  49  50  51  52  53]; % Indices of SVs in compartment buffer
% PARTIAL DERIVATIVES
pHBpK = -sum( (h_buffer*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpM = -sum( (h_buffer*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_buffer./Km(ii)+k_buffer./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_buffer*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_buffer*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_buffer./Kh(ii)+k_buffer./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_buffer*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_buffer*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_buffer./Kh(ii)+m_buffer./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 + 1*J_ATPASE_buffer - 1*J_HPERM_buffer_to_im/VRegion_buffer*VRegion_im) / VWater_buffer;
J_M = (0 - 1*J_MPERM_buffer_to_im/VRegion_buffer*VRegion_im) / VWater_buffer;
J_K = (0 - 1*J_KPERM_buffer_to_im/VRegion_buffer*VRegion_im ) / VWater_buffer;
Phi_H = J_H - sum( h_buffer*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = J_M - sum( m_buffer*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K - sum( k_buffer*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(2)/K_BX(2)/(1+h_buffer/K_BX(2))^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(57,:) =  ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
f(58,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(59,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;
% COMPARTMENT im:
ii = [32  33  35  36  37  38  39  40  41  42  43  44  45]; % Indices of SVs in compartment im
% PARTIAL DERIVATIVES
pHBpK = -sum( (h_im*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpM = -sum( (h_im*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_im./Km(ii)+k_im./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_im*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_im*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_im./Kh(ii)+k_im./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_im*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_im*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_im./Kh(ii)+m_im./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 - 2.66667*J_F1F0ATPASE_im_to_matrix/VRegion_im*VRegion_matrix + 4*J_ETC1_im_to_matrix/VRegion_im*VRegion_matrix + (2 + 2)*J_ETC3_im_to_matrix/VRegion_im*VRegion_matrix + (4-2)*J_ETC4_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_HLEAK_im_to_matrix/VRegion_im*VRegion_matrix - 2*J_PIH_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_KH_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_PYRH_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_GLUH_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_HPERM_buffer_to_im - 1*J_ASPGLU_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_CITMAL_im_to_matrix/VRegion_im*VRegion_matrix  - 1*J_MPTT/VRegion_im*VRegion_matrix  ) / VWater_im;
J_M = (0 + 1*J_MPERM_buffer_to_im) / VWater_im;
J_K = (0 - 1*J_KH_im_to_matrix/VRegion_im*VRegion_matrix - 1*J_KLEAK_im_to_matrix/VRegion_im*VRegion_matrix + 1*J_KPERM_buffer_to_im ) / VWater_im;
Phi_H = J_H - sum( h_im*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = J_M - sum( m_im*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K - sum( k_im*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(3)/K_BX(3)/(1+h_im/K_BX(3))^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(60,:) =  ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
f(61,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(62,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;

%% ELECTROPHYS EQUATIONS
C_im_to_matrix = 3.0e-6;%*(1+abs(DPsi_im_to_matrix/10));
f(63) = ( 0 - 2.66667*J_F1F0ATPASE_im_to_matrix  + 4*J_ETC1_im_to_matrix  + 2*J_ETC3_im_to_matrix  + 4*J_ETC4_im_to_matrix  - 1*J_HLEAK_im_to_matrix - 1*J_KLEAK_im_to_matrix - 1*J_ANT_im_to_matrix  + 1*J_ASPGLU_im_to_matrix - 1*J_MPTT ) / C_im_to_matrix ;
% C_buffer_to_im = Inf;
% f(64) = ( 0 + 1*J_HPERM_buffer_to_im  + 1*J_KPERM_buffer_to_im  + 2*J_MPERM_buffer_to_im ) / C_buffer_to_im;
f(64) = 0;
%% FLUX VECTOR:
J(1) =  J_PDH_matrix;
J(2) =  J_CTS_matrix;
J(3) =  J_ACON_matrix;
J(4) =  J_IDH_matrix;
J(5) =  J_AKGDH_matrix;
J(6) =  J_SCS_matrix;
J(7) =  J_SDH_matrix;
J(8) =  J_FUM_matrix;
J(9) =  J_MDH_matrix;
J(10) =  J_NDK_matrix;
J(11) =  J_AAT_matrix;
J(12) =  J_GDH1_matrix;
J(13) =  J_ATPASE_buffer;
J(14) =  J_F1F0ATPASE_im_to_matrix;
J(15) =  J_ETC1_im_to_matrix;
J(16) =  J_ETC3_im_to_matrix;
J(17) =  J_ETC4_im_to_matrix;
J(18) =  J_HLEAK_im_to_matrix;
J(19) =  J_PIH_im_to_matrix;
J(20) =  J_ANT_im_to_matrix;
J(21) =  J_KH_im_to_matrix;
J(22) =  J_PYRH_im_to_matrix;
J(23) =  J_ASPGLU_im_to_matrix;
J(24) =  J_GLUH_im_to_matrix;
J(25) =  J_CITMAL_im_to_matrix;
J(26) =  J_AKGMAL_im_to_matrix;
J(27) =  0;
J(28) =  J_DICARBSUC_im_to_matrix;
J(29) =  J_O2PERM_im_to_matrix;
J(30) =  J_ADPPERM_buffer_to_im;
J(31) =  J_ATPPERM_buffer_to_im;
J(32) =  J_PIPERM_buffer_to_im;
J(33) =  J_HPERM_buffer_to_im;
J(34) =  J_KPERM_buffer_to_im;
J(35) =  J_MPERM_buffer_to_im;
J(36) =  J_PYRPERM_buffer_to_im;
J(37) =  J_CITPERM_buffer_to_im;
J(38) =  J_AKGPERM_buffer_to_im;
J(39) =  J_SUCPERM_buffer_to_im;
J(40) =  J_GLUPERM_buffer_to_im;
J(41) =  J_ASPPERM_buffer_to_im;
J(42) =  J_MALPERM_buffer_to_im;
J(43) =  J_O2PERM2_buffer_to_im;
J(44) =  JSO1 ; 
J(45) =  JH2O2;
J(46) =  JSO3;
J(47) =  J_SDH_SO;
J(48) =  J_SDH_H2O2;
