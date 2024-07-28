function [P, T, Rho, c, mu, k, g] = USAtmo(Alt)

%**********************************************************
%|Input|
% 1)Alt：Geometric Altitude of desired altitude(m) 
%
%|Output|
% 1)P：Pressure(Pa)
% 2)T：Temperature(K)
% 3)Rho：Density(kg/m^3)
% 4)c：Speed of Sound(m/s)
% 5)mu：Dynamic Viscosity(N*s/m^2)
% 6)k：Thermal Conductivity Coefficient(W/(m*K))
% 7)g：Gravity(m/s^2)
%**********************************************************

Alt = Alt / 1000;


%% 数据校验
if length(Alt) ~= 1
    error('*** Error：Alt Input ***')
elseif Alt < 0 || Alt > 86
    error('*** Error：Program only valid for 0 <= Alt <= 86km ***')
end


%% 初始化
%Conversion Factor Used in 80<alt<86 km
Z_M = 80:0.5:86;
M_M_0 = [1, 0.999996, 0.999989, 0.999971, 0.999941, 0.999909, 0.999870,...
            0.999829, 0.999786, 0.999741, 0.999694, 0.999641, 0.999579];

%Constants
M_0 = 28.9644;
Beta = 1.458e-6;
Gamma = 1.4;
g_0 = 9.80665;
R = 8.31432e3;
r_E = 6.356766e3;
S = 110.4;


%% 计算
%Temperature
T = USAtmo_T(Alt);
T_M = T;
if Alt > 80 && Alt < 86
    T = T * interp1(Z_M,M_M_0, Alt);
end

%Pressure
P = USAtmo_P(Alt);

%Density
Rho = M_0 * P / (R*T_M);

%Speed of Sound
c = sqrt(Gamma * R * T_M / M_0);

%Dynamic Viscosity
mu = Beta * T.^1.5 / (T + S);

%Thermal Conductivity Coefficient
k = 2.64638e-3 * T.^1.5 / (T + 245 * 10.^(-12/T));

%Gravity
g = g_0 * (r_E / (r_E+Alt)).^2;


end


%% 
function P = USAtmo_P(Alt)

%Constants
N_A = 6.022169e26;
g_0 = 9.80665;
M_0 = 28.9644;
R = 8.31432e3;
r_E = 6.356766e3;

%Geopotential/Geometric Altitudes used for Geometric Altitudes < 86 km
H = [0, 11, 20, 32, 47, 51, 71, 84.852];
Z = r_E * H ./ (r_E-H);
Z(8) = 86;

%Defined temperatures/lapse rates/pressures/density at each layer
T_M_B = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65];
L = [-6.5, 0, 1, 2.8, 0, -2.8, -2]/1e3;
P_ref = [1.01325e5, 2.2632e4, 5.4748e3, 8.6801e2, 1.1090e2, 6.6938e1, 3.9564];

P = zeros(size(Alt));
for i = 1 : length(Alt)
    Z_i = Alt(i);
    index = find(Z>=Z_i) - 1 + double(Z_i==0);
    index = index(1);
    Z_H = r_E * Z_i / (r_E+Z_i);
    if L(index) == 0
        P(i) = P_ref(index) * exp(-g_0*M_0*(Z_H-H(index))*1e3/(R*T_M_B(index)));
    else
        P(i) = P_ref(index) *...
              (T_M_B(index)/(T_M_B(index)+L(index)*(Z_H-H(index))*1e3))^...
              (g_0*M_0/(R*L(index)));
    end
end

end


%%
function T = USAtmo_T(Alt)

%Constants
r_E = 6.356766e3;

%Defined temperatures at each layer
T0 = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65,...
     214.65, 186.95, 186.8673, 240, 360, 1000];
L = [-6.5, 0, 1, 2.8, 0, -2.8, -2, 0, 0, 12, 0];

%Geopotential/Geometric Altitudes used for Geometric Altitudes < 86 km
H = [0, 11, 20, 32, 47, 51, 71];
Z = r_E * H ./ (r_E-H);

%Geometric Altitudes used for Altitudes >86 km
Z(8:12) = [86, 91, 110, 120, 1000];

%Temperature Calculation with Molecular Temperature below 86 km and Kinetic Temperature above
if Alt >= Z(1) && Alt <= Z(8)
    T = interp1(Z,T0,Alt);
    
elseif Alt > Z(8) && Alt <= Z(9)
    T = T0(9);
    
elseif Alt > Z(9) && Alt <= Z(10)
    a = 19.9429;
    A = -76.3232;
    T_c = 263.1905;
    T = T_c + A * sqrt(1-((Alt-Z(9))/a)^2);
    
elseif  Alt > Z(10) && Alt <= Z(11)
    T = interp1(Z,T0,Alt);
    
elseif Alt > Z(11)
    Lambda = L(10) / (T0(12)-T0(11));
    xi = (Alt-Z(11)) * (r_E+Z(11)) / (r_E+Alt);
    T = T0(12) - (T0(12)-T0(11)) * exp(-Lambda*xi);
end

end



















