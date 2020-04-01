%% Question 1, 2
% Voltage sweep from 0.1V to 10V gives an R3 of 14 Ohms.

R3 = 14;

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

% C, G and F matrices
G =[ -1,  1,  1,    0,    0,     0,  0;
      0,  0,  0,   -1,    1,     0,  0; 
     -1,  0,  0,    0, 1/R1, -1/R1,  0;
      0,  0, -1, 1/R3,    1,     0,  0;
      0, -1,  0,    0, 1/R2,     0,  0;     
      0,  0,  alpha,    0,    0,     0, -1; 
      0,  0,  0,    0,    0,     1,  0];
  
C = [ 0,  0,  0,    0,    0,     0,  0;
      0,  0,  L,    0,    0,     0,  0; 
      0,  0,  0,    0,  Cap,  -Cap,  0;
      0,  0,  0,    0,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0;    
      0,  0,  0,    0,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0];
  
 
% perform DC sweep
V_sweep = linspace(-10, 10, 100);
Vo_arr = zeros(100,1);
V3_arr = zeros(100,1);
for i = 1:length(V_sweep)
    Vin = V_sweep(i);
    F = [0; 0; 0; 0; 0; 0; Vin];
    Vout = G \ F;
    Vo_arr(i) = Vout(7);       % not sure if it's actually 7, but this is a good placeholder for now
    V3_arr(i) = Vout(3);       % not sure if actually V3, come back to this   
end

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
omega_sweep = linspace(100, 1E9, 1000) .* 2 .* pi;     % sweep from 100Hz to 1GHz in 1000 steps
Vo_arr = zeros(1000,1);
Vin = 10;
F = [0; 0; 0; 0; 0; 0; Vin];
for i = 1:length(omega_sweep)
    Vout = (G + omega_sweep(i) .* 1i .* C) \ F;
    Vo_arr(i) = Vout(7);            % again, not sure if 7 is right. Placeholder
end


% plot against omega
figure
plot(omega_sweep, real(Vo_arr));    % only plotting real part
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
omega = pi;      % single value for sweep this time
Vo_arr = zeros(1000,1);
old_cap = Cap;                                          % save old capacitor value before we fiddle with it
for i = 1:1000
    Cap = normrnd(old_cap, 0.05);   % random perturbation in capacitor
    C = [ 0,  0,  0,    0,    0,     0,  0;
      0,  0,  L,    0,    0,     0,  0; 
      0,  0,  0,    0,  Cap,  -Cap,  0;
      0,  0,  0,    0,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0;    
      0,  0,  0,    0,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0];
    Vout = (G + omega .* 1i .* C) \ F;
    Vo_arr(i) = Vout(7);            % again, not sure if 7 is right. Placeholder
end

% plot histogram of gain
gain = Vin ./ real(Vo_arr);
figure
histogram(gain, 'BinLimits', [-2 2])
title('Gain Histogram')
xlabel('Gain')
ylabel('Counts')


%% Question 4
% This appears to be a low-pass amplfier. A high gain at low frequencies
% should be expected, dropping off rapidly at higher frequencies. See scan
% of additional work for FD derivation.


