% =========================================================================
% Power System Analysis - Gauss-Seidel Version
% Part 1: Gauss-Seidel Load Flow
% Part 2: Symmetrical Component Fault Analysis (LG and LLG)
%
% This script uses the Gauss-Seidel method for load flow and is
% designed to work with the standardized IEEE 5-bus and 9-bus system data.
% =========================================================================

clc;
clear all;
close all;

%% 1. Data Input
% -------------------------------------------------------------------------
% Select the system to analyze by uncommenting one of the lines below.
fileName = 'IEEE_5.xlsx';
% fileName = 'IEEE_9.xlsx';

% --- Check if file exists in the current directory ---
if ~isfile(fileName)
    fprintf('Error: The file "%s" was not found.\n', fileName);
    fprintf('Please make sure the Excel file is saved in the same directory as this script.\n');
    fprintf('Your current MATLAB folder is: %s\n', pwd);
    return;
end

% --- Read data using readmatrix ---
try
    busData = readmatrix(fileName, 'Sheet', 'BusData');
    lineData = readmatrix(fileName, 'Sheet', 'LineData');
    fprintf('Successfully read data from %s.\n', fileName);
catch ME
    fprintf('\nError reading data from "%s".\n', fileName);
    fprintf('Check that the file contains sheets named "BusData" and "LineData" (case-sensitive).\n');
    fprintf('MATLAB Error Message: %s\n', ME.message);
    return;
end

% --- Base MVA ---
BaseMVA = 100;

%% 2. Process Line Data and Build Positive Sequence Ybus
% -------------------------------------------------------------------------
fprintf('Building positive sequence Ybus matrix...\n');
fromBus = lineData(:, 1);
toBus = lineData(:, 2);
R1 = lineData(:, 3);
X1 = lineData(:, 4);
B1_half = 1i * lineData(:, 5);
taps = lineData(:, 6);
numBuses = max(max(fromBus), max(toBus));
numLines = length(fromBus);
Z1_series = R1 + 1i * X1;
y1_series = 1 ./ Z1_series;
Ybus1 = zeros(numBuses, numBuses);
for k = 1:numLines
    p = fromBus(k);
    q = toBus(k);
    Ybus1(p, q) = Ybus1(p, q) - y1_series(k) / taps(k);
    Ybus1(q, p) = Ybus1(p, q);
    Ybus1(p, p) = Ybus1(p, p) + y1_series(k) / (taps(k)^2) + B1_half(k);
    Ybus1(q, q) = Ybus1(q, q) + y1_series(k) + B1_half(k);
end

%% 3. Process Bus Data and Initialize for Load Flow
% -------------------------------------------------------------------------
busNo = busData(:, 1);
busType = busData(:, 2);
Vmag = busData(:, 3);
Vang_deg = busData(:, 4);
PL = busData(:, 5) / BaseMVA;
QL = busData(:, 6) / BaseMVA;
PG = busData(:, 7) / BaseMVA;
QG = busData(:, 8) / BaseMVA;
Psp = PG - PL;
Qsp = QG - QL;
slackBus = find(busType == 1);
pvBuses = find(busType == 2);
pqBuses = find(busType == 3);
Vang_rad = deg2rad(Vang_deg);
Vmag(pqBuses) = 1.0; % Flat start for PQ buses
V = Vmag .* exp(1i * Vang_rad);

%% 4. Gauss-Seidel (G-S) Solver
% -------------------------------------------------------------------------
fprintf('Starting Gauss-Seidel iterations...\n');
iter = 0;
tolerance = 1;
max_iter = 500; % G-S may require more iterations

startTime = tic;
V_prev = V;
while (tolerance > 1e-5 && iter < max_iter)
    iter = iter + 1;
    
    % Iterate through all non-slack buses
    for i = 2:numBuses
        
        sum_YV = 0;
        for j = 1:numBuses
            if i ~= j
                sum_YV = sum_YV + Ybus1(i, j) * V(j);
            end
        end
        
        % For PQ buses, P and Q are specified
        if busType(i) == 3 % PQ Bus
            S_sp = Psp(i) + 1i * Qsp(i);
            V(i) = (1/Ybus1(i,i)) * (conj(S_sp)/conj(V_prev(i)) - sum_YV);
        
        % For PV buses, P and |V| are specified
        elseif busType(i) == 2 % PV Bus
            % First, calculate the new reactive power Q
            Q_calc = -imag(conj(V_prev(i)) * (sum_YV + Ybus1(i,i)*V_prev(i)));
            
            % For this simplified model, we assume Q is within limits.
            % In a full model, you would check Q against Qmin and Qmax.
            
            S_sp = Psp(i) + 1i * Q_calc;
            V_new_complex = (1/Ybus1(i,i)) * (conj(S_sp)/conj(V_prev(i)) - sum_YV);
            
            % Keep the specified voltage magnitude, but update the angle
            V(i) = Vmag(i) * (V_new_complex / abs(V_new_complex));
        end
    end
    
    % Check for convergence
    tolerance = max(abs(abs(V) - abs(V_prev)));
    V_prev = V; % Update for next iteration
end
solveTime = toc(startTime);

% Store pre-fault voltage for fault calculations
V_prefault = V;

%% 5. Display Load Flow Results
% -------------------------------------------------------------------------
Vmag = abs(V);
Vang_rad = angle(V);
S_final_pu = V .* conj(Ybus1 * V);
P_final_pu = real(S_final_pu);
PG_final_MW = (P_final_pu + PL) * BaseMVA;
QG_final_MVAr = (imag(S_final_pu) + QL) * BaseMVA;
Vang_deg_final = rad2deg(Vang_rad);
total_loss_MW = sum(P_final_pu) * BaseMVA;
fprintf('\n----------------- LOAD FLOW RESULTS (Gauss-Seidel) -----------------\n');
if iter >= max_iter
    fprintf('\nWARNING: Load flow did not converge.\n');
else
    fprintf('\nConverged in %d iterations (%.4f seconds).\n', iter, solveTime);
end
fprintf('\n==================================================================\n');
fprintf('                         Final Bus Data\n');
fprintf('==================================================================\n');
fprintf(' Bus | Type | V (p.u.) | Angle (deg) | P Gen (MW) | Q Gen (MVAr)\n');
fprintf('-----|------|----------|-------------|------------|--------------\n');
for i = 1:numBuses
    fprintf(' %3d | %4d |  %7.4f |   %8.4f  |  %9.2f |   %9.2f\n', ...
        busNo(i), busType(i), Vmag(i), Vang_deg_final(i), PG_final_MW(i), QG_final_MVAr(i));
end
fprintf('------------------------------------------------------------------\n');
fprintf('Total Real Power Loss: %.4f MW\n', total_loss_MW);
fprintf('------------------------------------------------------------------\n');


%% 6. Symmetrical Component and Fault Analysis Setup
% -------------------------------------------------------------------------
fprintf('\nBuilding sequence impedance matrices for fault analysis...\n');
Zbus1 = inv(Ybus1);
Zbus2 = Zbus1;
R0 = R1;
X0 = 3 * X1;
Z0_series = R0 + 1i * X0;
y0_series = 1 ./ Z0_series;
Ybus0 = zeros(numBuses, numBuses);
for k = 1:numLines
    if taps(k) == 1
        p = fromBus(k);
        q = toBus(k);
        Ybus0(p, q) = Ybus0(p, q) - y0_series(k);
        Ybus0(q, p) = Ybus0(p, q);
        Ybus0(p, p) = Ybus0(p, p) + y0_series(k);
        Ybus0(q, q) = Ybus0(q, q) + y0_series(k);
    end
end
genBuses = find(busType == 2 | busType == 1);
Xg0 = 0.05;
for i = 1:length(genBuses)
    bus_idx = genBuses(i);
    Ybus0(bus_idx, bus_idx) = Ybus0(bus_idx, bus_idx) + 1/(1i*Xg0);
end
Zbus0 = pinv(Ybus0);
a = exp(1i*120*pi/180);
A = [1 1 1; 1 a^2 a; 1 a a^2];

%% 7. Interactive Fault Calculation Loop
% -------------------------------------------------------------------------
while true
    fprintf('\n----------------- FAULT ANALYSIS -----------------\n');
    faultBus = input('Enter the bus number to apply fault (or 0 to exit): ');
    if faultBus == 0, break; end
    if faultBus < 1 || faultBus > numBuses
        fprintf('Invalid bus number. Please try again.\n');
        continue;
    end
    faultType = input('Select fault type (1 for LG, 2 for LLG): ');
    Zf = input('Enter fault impedance Zf (p.u.) [e.g., 0]: ');
    Z1 = Zbus1(faultBus, faultBus);
    Z2 = Zbus2(faultBus, faultBus);
    Z0 = Zbus0(faultBus, faultBus);
    Vf_prefault = V_prefault(faultBus);
    if faultType == 1
        Ia1 = Vf_prefault / (Z1 + Z2 + Z0 + 3*Zf);
        Ia2 = Ia1;
        Ia0 = Ia1;
    elseif faultType == 2
        Z_parallel = (Z2 * (Z0 + 3*Zf)) / (Z2 + Z0 + 3*Zf);
        Ia1 = Vf_prefault / (Z1 + Z_parallel);
        Ia2 = -Ia1 * (Z0 + 3*Zf) / (Z2 + Z0 + 3*Zf);
        Ia0 = -Ia1 * Z2 / (Z2 + Z0 + 3*Zf);
    else
        fprintf('Invalid fault type. Please try again.\n');
        continue;
    end
    I_seq = [Ia0; Ia1; Ia2];
    I_phase = A * I_seq;
    if faultType == 1
        I_fault_mag_pu = abs(I_phase(1));
    else
        I_fault_mag_pu = abs(I_phase(2));
    end
    I_base = BaseMVA / (sqrt(3) * 13.8);
    I_fault_mag_kA = I_fault_mag_pu * I_base;
    MVA_symmetrical = abs(Vf_prefault) * I_fault_mag_pu * BaseMVA;
    fprintf('\n--- Fault Results for Bus %d ---\n', faultBus);
    fprintf('Sequence Currents (p.u.):\n');
    fprintf('  I0: %.4f /_ %.2f deg\n', abs(Ia0), rad2deg(angle(Ia0)));
    fprintf('  I1: %.4f /_ %.2f deg\n', abs(Ia1), rad2deg(angle(Ia1)));
    fprintf('  I2: %.4f /_ %.2f deg\n', abs(Ia2), rad2deg(angle(Ia2)));
    fprintf('\nPhase Currents (p.u.):\n');
    fprintf('  Ia: %.4f /_ %.2f deg\n', abs(I_phase(1)), rad2deg(angle(I_phase(1))));
    fprintf('  Ib: %.4f /_ %.2f deg\n', abs(I_phase(2)), rad2deg(angle(I_phase(2))));
    fprintf('  Ic: %.4f /_ %.2f deg\n', abs(I_phase(3)), rad2deg(angle(I_phase(3))));
    fprintf('\n--- Breaker Rating Calculation ---\n');
    fprintf('Total Fault Current: %.4f kA\n', I_fault_mag_kA);
    fprintf('Symmetrical Interrupting MVA: %.2f MVA\n', MVA_symmetrical);
end
fprintf('\nFault analysis complete.\n');
