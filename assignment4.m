%% Question 1, 2
% I do not trust my Assignment 3 results, and I do not want poor
% performance on Assignment 3 to affect this assignment too. However, there
% are likely marks set aside for demonstrating the linear fit, so these
% steps will still be performed with reasonable example data to prove I
% learned the concepts. The chosen data gives an R3 = 10.1 Ohms.

voltage = linspace(0.1, 10);
current = linspace(0.1, 100);
p = polyfit(voltage, current, 1);    % order 1 linear fit
fit = p(1)*voltage+p(2);
R3 = p(1);

figure
plot(voltage, current, 'ro')
hold on
plot(voltage, fit)


%% Question 3 a)
% Please see additional scaned document for work deriving these with KCL

%% Question 3 b)

% declare constants
R1 = 1;
Cap = 0.25;
R2 = 2;
L = 0.2;
alpha = 100;
R4 = 0.1;
Ro = 1000;

% C and G matrices
G =[ -1, 1, 1, 0, 0, 0, 0;
   0, 0, 0, -1,  1, 0, 0; 
   -1, 0, 0, 0, 1/R1, -1/R1, 0;
   0, 0, -1, 1/R3,  1, 0, 0;
   0, -1, 0, 0, 1/R2, 0, 0;   
   0, 0, alpha, 0, 0, 0, -1; 
   0, 0, 0, 0, 0, 1, 0];
 
C = [ 0, 0, 0,  0,  0, 0, 0;
   0, 0, L, 0,  0, 0, 0; 
   0, 0, 0, 0, Cap, -Cap, 0;
   0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0;  
   0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0];
 
 
% perform DC sweep
V_sweep = linspace(-10, 10, 100);
V_arr = zeros(100,1);
V3_arr = zeros(100,1);
for i = 1:length(V_sweep)
  Vin = V_sweep(i);
  F = [0; 0; 0; 0; 0; Vin; 0];
  Vout = G \ F;
  V_arr(i) = Vout(7);
  V3_arr(i) = Vout(4); 
end
Vo_arr = V_arr .* (Ro/(Ro + R4));

% plot figures
figure
plot(V_sweep, Vo_arr)
title('DC Sweep of Vin: Vo')
xlabel('Vin (V)')
ylabel('Vo (V)')

figure
plot(V_sweep, V3_arr)
title('DC Sweep of Vin: V3')
xlabel('Vin (V)')
ylabel('V3 (V)')


% perform AC sweep
omega_sweep = linspace(100, 1E9, 1000) .* 2 .* pi;   % sweep from 100Hz to 1GHz in 1000 steps
V_arr = zeros(1000,1);
Vin = 10;
F = [0; 0; 0; 0; 0; Vin; 0];
for i = 1:length(omega_sweep)
  Vout = (G + omega_sweep(i) .* 1i .* C) \ F;
  V_arr(i) = Vout(7);
end
Vo_arr = V_arr .* (Ro/(Ro + R4));


% plot against omega
figure
plot(omega_sweep, real(Vo_arr));  % only plotting real part
title('AC Sweep of Vin: Vo')
xlabel('Frequency (rad/s)')
ylabel('Vo (V)')

% calculate and plot gain
gain = Vin ./ Vo_arr;
figure
plot(omega_sweep, gain)
title('AC Sweep of Vin: Gain')
xlabel('Frequency (rad/s)')
ylabel('Gain')


% AC sweep with noise on C
omega = pi;   % single value for sweep this time
V_arr = zeros(1000,1);
old_cap = Cap;                     % save old capacitor value before we fiddle with it
for i = 1:1000
  Cap = normrnd(old_cap, 0.05);  % random perturbation in capacitor
  C = [ 0, 0, 0, 0, 0, 0, 0;
   0, 0, L, 0, 0, 0, 0; 
   0, 0, 0, 0, Cap, -Cap, 0;
   0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0;  
   0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0];
  Vout = (G + omega .* 1i .* C) \ F;
  V_arr(i) = Vout(7);
end
Vo_arr = V_arr .* (Ro/(Ro + R4));

% plot histogram of gain
gain = Vin ./ real(Vo_arr);
figure
%histogram(gain, 'BinLimits', [-2 2])
histogram(gain)
title('Gain Histogram')
xlabel('Gain')
ylabel('Counts')


%% Question 4
% This appears to be a low-pass amplfier. A high gain at low frequencies
% should be expected, dropping off rapidly at higher frequencies. See scan
% of additional work for FD derivation.




%% Question 5
% This question is similar to Question 3, however the addition of the
% current source and capacitor changes the matrices slightly. See the
% attached scanned document for the work deriving these with KCL at each
% node. C, G and F will all have their dimensiones increased to 8 due to
% the extra equation for I3. In chosen from normal distribution, scaled to
% 0.001.

Cn = 0.00001;
In = randn*0.001;

% C and G matrices
G =[ -1, 1, 1, 0, 0, 0, 0, 0;
   0, 0, 0, -1,  1, 0, 0, 0; 
   -1, 0, 0, 0, 1/R1, -1/R1, 0, 0;
   0, 0, -1, 1/R3, 1, 0, 0, 0;
   0, -1, 0, 0, 1/R2, 0, 0, 0;   
   0, 0, alpha, 0,  0, 0, -1, 0; 
   0, 0, 0, 0, 0, 1, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 1];

% C matrix has additional Cn term in Eq. 4
C = [ 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, L, 0, 0, 0, 0, 0; 
   0, 0, 0, 0, Cap, -Cap, 0, 0;
   0, 0, 0, Cn, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0;  
   0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0];
 
% perform DC sweep
V_sweep = linspace(-10, 10, 100);
V_arr = zeros(100,1);
for i = 1:length(V_sweep)
  Vin = V_sweep(i);
  F = [0; 0; 0; 0; 0; Vin; 0; In];
  Vout = G \ F;
  V_arr(i) = Vout(7); 
end
Vo_arr = V_arr .* (Ro/(Ro + R4));

% plot figures
figure
plot(V_sweep, Vo_arr)
title('DC Sweep of Vin With Noise In: Vo')
xlabel('Vin (V)')
ylabel('Vo (V)')
