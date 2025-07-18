classdef PowerSystemAnalyzerGUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        FileSelectionLabel              matlab.ui.control.Label
        SelectDataFileButton            matlab.ui.control.Button
        FileNameField                   matlab.ui.control.EditField
        RunAnalysisButton               matlab.ui.control.Button
        StatusLabel                     matlab.ui.control.Label
        StatusField                     matlab.ui.control.EditField
        RightPanel                      matlab.ui.container.Panel
        LoadFlowResultsTabGroup         matlab.ui.container.TabGroup
        LoadFlowResultsTab              matlab.ui.container.Tab
        LoadFlowUITable                 matlab.ui.control.Table
        LoadFlowSummaryLabel            matlab.ui.control.Label
        FaultAnalysisTab                matlab.ui.container.Tab
        FaultAnalysisGridLayout         matlab.ui.container.GridLayout
        FaultInputsPanel                matlab.ui.container.Panel
        FaultBusDropDownLabel           matlab.ui.control.Label
        FaultBusDropDown                matlab.ui.control.DropDown
        FaultTypeDropDownLabel          matlab.ui.control.Label
        FaultTypeDropDown               matlab.ui.control.DropDown
        FaultImpedancepuEditFieldLabel  matlab.ui.control.Label
        FaultImpedancepuEditField       matlab.ui.control.NumericEditField
        CalculateFaultButton            matlab.ui.control.Button
        FaultResultsPanel               matlab.ui.container.Panel
        SequenceCurrentsLabel           matlab.ui.control.Label
        SequenceCurrentsTextArea        matlab.ui.control.TextArea
        PhaseCurrentsLabel              matlab.ui.control.Label
        PhaseCurrentsTextArea           matlab.ui.control.TextArea
        BreakerRatingLabel              matlab.ui.control.Label
        BreakerRatingTextArea           matlab.ui.control.TextArea
    end

    
    properties (Access = private)
        % --- Properties to store data and results ---
        busData
        lineData
        V_prefault
        Zbus0
        Zbus1
        Zbus2
        numBuses
        BaseMVA = 100;
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SelectDataFileButton
        function SelectDataFileButtonPushed(app, event)
            [file, path] = uigetfile('*.xlsx', 'Select IEEE Data File');
            if isequal(file, 0)
                app.StatusField.Value = 'File selection cancelled.';
            else
                fullPath = fullfile(path, file);
                app.FileNameField.Value = fullPath;
                app.StatusField.Value = ['Selected file: ' file];
                app.RunAnalysisButton.Enable = 'on';
            end
        end

        % Button pushed function: RunAnalysisButton
        function RunAnalysisButtonPushed(app, event)
            app.StatusField.Value = 'Running analysis... Please wait.';
            drawnow; % Update the UI to show the message
            
            try
                fileName = app.FileNameField.Value;
                
                % --- Part 1: Load Flow Analysis ---
                app.StatusField.Value = 'Reading data from Excel file...';
                app.busData = readmatrix(fileName, 'Sheet', 'BusData');
                app.lineData = readmatrix(fileName, 'Sheet', 'LineData');
                
                app.StatusField.Value = 'Building Ybus matrix...';
                [Ybus1, R1, X1, taps, fromBus, toBus] = app.buildYbus(app.lineData);
                
                app.StatusField.Value = 'Running FDLF solver...';
                [V, Vang_deg_final, PG_final_MW, QG_final_MVAr, total_loss_MW, iter, solveTime] = app.runFDLF(Ybus1, app.busData);
                app.V_prefault = V; % Store pre-fault voltage
                
                app.StatusField.Value = 'Populating results...';
                app.displayLoadFlowResults(V, Vang_deg_final, PG_final_MW, QG_final_MVAr, total_loss_MW, iter, solveTime);
                
                % --- Part 2: Prepare for Fault Analysis ---
                app.StatusField.Value = 'Building sequence impedance matrices...';
                app.buildZbusMatrices(Ybus1, R1, X1, taps, fromBus, toBus);
                
                % Enable fault analysis components
                app.FaultBusDropDown.Items = cellstr(num2str((1:app.numBuses)'));
                app.FaultBusDropDown.Value = '1'; % Default to bus 1
                app.FaultAnalysisTab.Parent.SelectedTab = app.FaultAnalysisTab;
                app.CalculateFaultButton.Enable = 'on';
                app.FaultBusDropDown.Enable = 'on';
                app.FaultTypeDropDown.Enable = 'on';
                app.FaultImpedancepuEditField.Enable = 'on';
                
                app.StatusField.Value = 'Analysis complete. Ready for fault calculation.';
                
            catch ME
                app.StatusField.Value = ['ERROR: ' ME.message];
                uialert(app.UIFigure, ['An error occurred: ' ME.message], 'Analysis Error');
            end
        end

        % Button pushed function: CalculateFaultButton
        function CalculateFaultButtonPushed(app, event)
            try
                % Get user inputs from the GUI
                faultBus = str2double(app.FaultBusDropDown.Value);
                faultTypeStr = app.FaultTypeDropDown.Value;
                Zf = app.FaultImpedancepuEditField.Value;

                if strcmp(faultTypeStr, 'Line-to-Ground (LG)')
                    faultType = 1;
                else
                    faultType = 2;
                end

                % Get pre-calculated data
                Z1 = app.Zbus1(faultBus, faultBus);
                Z2 = app.Zbus2(faultBus, faultBus);
                Z0 = app.Zbus0(faultBus, faultBus);
                Vf_prefault = app.V_prefault(faultBus);
                
                % Calculate sequence currents
                if faultType == 1 % LG
                    Ia1 = Vf_prefault / (Z1 + Z2 + Z0 + 3*Zf);
                    Ia2 = Ia1;
                    Ia0 = Ia1;
                else % LLG
                    Z_parallel = (Z2 * (Z0 + 3*Zf)) / (Z2 + Z0 + 3*Zf);
                    Ia1 = Vf_prefault / (Z1 + Z_parallel);
                    Ia2 = -Ia1 * (Z0 + 3*Zf) / (Z2 + Z0 + 3*Zf);
                    Ia0 = -Ia1 * Z2 / (Z2 + Z0 + 3*Zf);
                end
                
                % Calculate phase currents
                a = exp(1i*120*pi/180);
                A = [1 1 1; 1 a^2 a; 1 a a^2];
                I_seq = [Ia0; Ia1; Ia2];
                I_phase = A * I_seq;
                
                % Calculate breaker ratings
                if faultType == 1
                    I_fault_mag_pu = abs(I_phase(1));
                else
                    I_fault_mag_pu = abs(I_phase(2));
                end
                I_base = app.BaseMVA / (sqrt(3) * 13.8); % Assuming 13.8kV base
                I_fault_mag_kA = I_fault_mag_pu * I_base;
                MVA_symmetrical = abs(Vf_prefault) * I_fault_mag_pu * app.BaseMVA;
                
                % Display results in text areas
                app.displayFaultResults(Ia0, Ia1, Ia2, I_phase, I_fault_mag_kA, MVA_symmetrical);
                app.StatusField.Value = sprintf('Fault calculation for Bus %d complete.', faultBus);

            catch ME
                app.StatusField.Value = ['ERROR: ' ME.message];
                uialert(app.UIFigure, ['An error occurred during fault calculation: ' ME.message], 'Calculation Error');
            end
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 950 650];
            app.UIFigure.Name = 'Power System Load Flow and Fault Analyzer';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {300, '1x'};
            app.GridLayout.RowHeight = {'1x'};

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Title = 'Setup & Control';
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create FileSelectionLabel
            app.FileSelectionLabel = uilabel(app.LeftPanel);
            app.FileSelectionLabel.HorizontalAlignment = 'center';
            app.FileSelectionLabel.FontSize = 14;
            app.FileSelectionLabel.FontWeight = 'bold';
            app.FileSelectionLabel.Position = [28 580 245 22];
            app.FileSelectionLabel.Text = '1. Select System Data File';

            % Create SelectDataFileButton
            app.SelectDataFileButton = uibutton(app.LeftPanel, 'push');
            app.SelectDataFileButton.ButtonPushedFcn = createCallbackFcn(app, @SelectDataFileButtonPushed, true);
            app.SelectDataFileButton.FontSize = 14;
            app.SelectDataFileButton.Position = [28 540 120 30];
            app.SelectDataFileButton.Text = 'Browse...';

            % Create FileNameField
            app.FileNameField = uieditfield(app.LeftPanel, 'text');
            app.FileNameField.Editable = 'off';
            app.FileNameField.Position = [28 505 245 25];
            app.FileNameField.Value = 'No file selected.';

            % Create RunAnalysisButton
            app.RunAnalysisButton = uibutton(app.LeftPanel, 'push');
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);
            app.RunAnalysisButton.BackgroundColor = [0.4706 0.6706 0.1882];
            app.RunAnalysisButton.FontSize = 14;
            app.RunAnalysisButton.FontColor = [1 1 1];
            app.RunAnalysisButton.Enable = 'off';
            app.RunAnalysisButton.Position = [28 450 245 35];
            app.RunAnalysisButton.Text = '2. Run Load Flow & Prepare';

            % Create StatusLabel
            app.StatusLabel = uilabel(app.LeftPanel);
            app.StatusLabel.HorizontalAlignment = 'center';
            app.StatusLabel.FontSize = 14;
            app.StatusLabel.FontWeight = 'bold';
            app.StatusLabel.Position = [28 100 245 22];
            app.StatusLabel.Text = 'Status';

            % Create StatusField
            app.StatusField = uieditfield(app.LeftPanel, 'text');
            app.StatusField.Editable = 'off';
            app.StatusField.BackgroundColor = [0.902 0.902 0.902];
            app.StatusField.Position = [28 40 245 50];
            app.StatusField.Value = 'Ready. Please select a data file.';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Title = 'Analysis Results';
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create LoadFlowResultsTabGroup
            app.LoadFlowResultsTabGroup = uitabgroup(app.RightPanel);
            app.LoadFlowResultsTabGroup.Position = [20 20 600 500];

            % Create LoadFlowResultsTab
            app.LoadFlowResultsTab = uitab(app.LoadFlowResultsTabGroup);
            app.LoadFlowResultsTab.Title = 'Load Flow Results';

            % Create LoadFlowUITable
            app.LoadFlowUITable = uitable(app.LoadFlowResultsTab);
            app.LoadFlowUITable.ColumnName = {'Bus'; 'Type'; 'V (p.u.)'; 'Angle (deg)'; 'P Gen (MW)'; 'Q Gen (MVAr)'};
            app.LoadFlowUITable.RowName = {};
            app.LoadFlowUITable.Position = [15 50 570 400];

            % Create LoadFlowSummaryLabel
            app.LoadFlowSummaryLabel = uilabel(app.LoadFlowResultsTab);
            app.LoadFlowSummaryLabel.FontSize = 14;
            app.LoadFlowSummaryLabel.Position = [15 510 570 50];
            app.LoadFlowSummaryLabel.Text = 'Load Flow Summary';

            % Create FaultAnalysisTab
            app.FaultAnalysisTab = uitab(app.LoadFlowResultsTabGroup);
            app.FaultAnalysisTab.Title = 'Fault Analysis';

            % Create FaultAnalysisGridLayout
            app.FaultAnalysisGridLayout = uigridlayout(app.FaultAnalysisTab);
            app.FaultAnalysisGridLayout.RowHeight = {200, '1x'};
            app.FaultAnalysisGridLayout.ColumnWidth = {'1x'};

            % Create FaultInputsPanel
            app.FaultInputsPanel = uipanel(app.FaultAnalysisGridLayout);
            app.FaultInputsPanel.Title = 'Fault Inputs';
            app.FaultInputsPanel.Layout.Row = 1;
            app.FaultInputsPanel.Layout.Column = 1;

            % Create FaultBusDropDownLabel
            app.FaultBusDropDownLabel = uilabel(app.FaultInputsPanel);
            app.FaultBusDropDownLabel.HorizontalAlignment = 'right';
            app.FaultBusDropDownLabel.Position = [30 130 65 22];
            app.FaultBusDropDownLabel.Text = 'Fault Bus';

            % Create FaultBusDropDown
            app.FaultBusDropDown = uidropdown(app.FaultInputsPanel);
            app.FaultBusDropDown.Items = {'N/A'};
            app.FaultBusDropDown.Enable = 'off';
            app.FaultBusDropDown.Position = [110 130 100 22];
            app.FaultBusDropDown.Value = 'N/A';

            % Create FaultTypeDropDownLabel
            app.FaultTypeDropDownLabel = uilabel(app.FaultInputsPanel);
            app.FaultTypeDropDownLabel.HorizontalAlignment = 'right';
            app.FaultTypeDropDownLabel.Position = [30 90 65 22];
            app.FaultTypeDropDownLabel.Text = 'Fault Type';

            % Create FaultTypeDropDown
            app.FaultTypeDropDown = uidropdown(app.FaultInputsPanel);
            app.FaultTypeDropDown.Items = {'Line-to-Ground (LG)', 'Double Line-to-Ground (LLG)'};
            app.FaultTypeDropDown.Enable = 'off';
            app.FaultTypeDropDown.Position = [110 90 200 22];

            % Create FaultImpedancepuEditFieldLabel
            app.FaultImpedancepuEditFieldLabel = uilabel(app.FaultInputsPanel);
            app.FaultImpedancepuEditFieldLabel.HorizontalAlignment = 'right';
            app.FaultImpedancepuEditFieldLabel.Position = [10 50 85 22];
            app.FaultImpedancepuEditFieldLabel.Text = 'Fault Imp (p.u.)';

            % Create FaultImpedancepuEditField
            app.FaultImpedancepuEditField = uieditfield(app.FaultInputsPanel, 'numeric');
            app.FaultImpedancepuEditField.Enable = 'off';
            app.FaultImpedancepuEditField.Position = [110 50 100 22];

            % Create CalculateFaultButton
            app.CalculateFaultButton = uibutton(app.FaultInputsPanel, 'push');
            app.CalculateFaultButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateFaultButtonPushed, true);
            app.CalculateFaultButton.BackgroundColor = [0.851 0.3255 0.098];
            app.CalculateFaultButton.FontSize = 14;
            app.CalculateFaultButton.FontColor = [1 1 1];
            app.CalculateFaultButton.Enable = 'off';
            app.CalculateFaultButton.Position = [350 70 150 50];
            app.CalculateFaultButton.Text = 'Calculate Fault';

            % Create FaultResultsPanel
            app.FaultResultsPanel = uipanel(app.FaultAnalysisGridLayout);
            app.FaultResultsPanel.Title = 'Fault Results';
            app.FaultResultsPanel.Layout.Row = 2;
            app.FaultResultsPanel.Layout.Column = 1;

            % Create SequenceCurrentsLabel
            app.SequenceCurrentsLabel = uilabel(app.FaultResultsPanel);
            app.SequenceCurrentsLabel.FontSize = 14;
            app.SequenceCurrentsLabel.FontWeight = 'bold';
            app.SequenceCurrentsLabel.Position = [20 280 150 22];
            app.SequenceCurrentsLabel.Text = 'Sequence Currents';

            % Create SequenceCurrentsTextArea
            app.SequenceCurrentsTextArea = uitextarea(app.FaultResultsPanel);
            app.SequenceCurrentsTextArea.Editable = 'off';
            app.SequenceCurrentsTextArea.FontName = 'Courier New';
            app.SequenceCurrentsTextArea.Position = [20 190 250 80];

            % Create PhaseCurrentsLabel
            app.PhaseCurrentsLabel = uilabel(app.FaultResultsPanel);
            app.PhaseCurrentsLabel.FontSize = 14;
            app.PhaseCurrentsLabel.FontWeight = 'bold';
            app.PhaseCurrentsLabel.Position = [300 280 150 22];
            app.PhaseCurrentsLabel.Text = 'Phase Currents';

            % Create PhaseCurrentsTextArea
            app.PhaseCurrentsTextArea = uitextarea(app.FaultResultsPanel);
            app.PhaseCurrentsTextArea.Editable = 'off';
            app.PhaseCurrentsTextArea.FontName = 'Courier New';
            app.PhaseCurrentsTextArea.Position = [300 190 250 80];

            % Create BreakerRatingLabel
            app.BreakerRatingLabel = uilabel(app.FaultResultsPanel);
            app.BreakerRatingLabel.FontSize = 14;
            app.BreakerRatingLabel.FontWeight = 'bold';
            app.BreakerRatingLabel.Position = [20 120 200 22];
            app.BreakerRatingLabel.Text = 'Breaker Rating Calculation';

            % Create BreakerRatingTextArea
            app.BreakerRatingTextArea = uitextarea(app.FaultResultsPanel);
            app.BreakerRatingTextArea.Editable = 'off';
            app.BreakerRatingTextArea.FontName = 'Courier New';
            app.BreakerRatingTextArea.Position = [20 30 530 80];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end
    
    % App specific helper functions
    methods (Access = private)
        
        % Function to build Ybus
        function [Ybus1, R1, X1, taps, fromBus, toBus] = buildYbus(app, lineData)
            fromBus = lineData(:, 1);
            toBus = lineData(:, 2);
            R1 = lineData(:, 3);
            X1 = lineData(:, 4);
            B1_half = 1i * lineData(:, 5);
            taps = lineData(:, 6);
            app.numBuses = max(max(fromBus), max(toBus));
            numLines = length(fromBus);
            Z1_series = R1 + 1i * X1;
            y1_series = 1 ./ Z1_series;
            Ybus1 = zeros(app.numBuses, app.numBuses);
            for k = 1:numLines
                p = fromBus(k);
                q = toBus(k);
                Ybus1(p, q) = Ybus1(p, q) - y1_series(k) / taps(k);
                Ybus1(q, p) = Ybus1(p, q);
                Ybus1(p, p) = Ybus1(p, p) + y1_series(k) / (taps(k)^2) + B1_half(k);
                Ybus1(q, q) = Ybus1(q, q) + y1_series(k) + B1_half(k);
            end
        end
        
        % Function to run FDLF
        function [V, Vang_deg_final, PG_final_MW, QG_final_MVAr, total_loss_MW, iter, solveTime] = runFDLF(app, Ybus1, busData)
            busNo = busData(:, 1);
            busType = busData(:, 2);
            Vmag = busData(:, 3);
            Vang_deg = busData(:, 4);
            PL = busData(:, 5) / app.BaseMVA;
            QL = busData(:, 6) / app.BaseMVA;
            PG = busData(:, 7) / app.BaseMVA;
            QG = busData(:, 8) / app.BaseMVA;
            Psp = PG - PL;
            Qsp = QG - QL;
            pvBuses = find(busType == 2);
            pqBuses = find(busType == 3);
            Vang_rad = deg2rad(Vang_deg);
            Vmag(pqBuses) = 1.0;
            V = Vmag .* exp(1i * Vang_rad);
            p_indices = [pvBuses; pqBuses];
            q_indices = pqBuses;
            iter = 0;
            tolerance = 1;
            max_iter = 100;
            B_prime = -imag(Ybus1(p_indices, p_indices));
            B_double_prime = -imag(Ybus1(q_indices, q_indices));
            inv_B_prime = inv(B_prime);
            inv_B_double_prime = inv(B_double_prime);
            startTime = tic;
            while (tolerance > 1e-5 & iter < max_iter)
                iter = iter + 1;
                S_calc = V .* conj(Ybus1 * V);
                P_calc = real(S_calc);
                dQ = Qsp(q_indices) - imag(S_calc(q_indices));
                dP = Psp(p_indices) - P_calc(p_indices);
                mismatch = [dP; dQ];
                d_ang = inv_B_prime * (dP ./ abs(V(p_indices)));
                d_mag = inv_B_double_prime * (dQ ./ abs(V(q_indices)));
                Vang_rad(p_indices) = Vang_rad(p_indices) + d_ang;
                Vmag(q_indices) = Vmag(q_indices) + d_mag;
                V = Vmag .* exp(1i * Vang_rad);
                tolerance = max(abs(mismatch));
            end
            solveTime = toc(startTime);
            S_final_pu = V .* conj(Ybus1 * V);
            P_final_pu = real(S_final_pu);
            PG_final_MW = (P_final_pu + PL) * app.BaseMVA;
            QG_final_MVAr = (imag(S_final_pu) + QL) * app.BaseMVA;
            Vang_deg_final = rad2deg(Vang_rad);
            total_loss_MW = sum(P_final_pu) * app.BaseMVA;
        end
        
        % Function to display load flow results
        function displayLoadFlowResults(app, V, Vang_deg_final, PG_final_MW, QG_final_MVAr, total_loss_MW, iter, solveTime)
            busNo = app.busData(:,1);
            busType = app.busData(:,2);
            Vmag = abs(V);
            
            T = table(busNo, busType, Vmag, Vang_deg_final, PG_final_MW, QG_final_MVAr);
            T.Properties.VariableNames = {'Bus', 'Type', 'V_pu', 'Angle_deg', 'P_Gen_MW', 'Q_Gen_MVAr'};
            app.LoadFlowUITable.Data = T;
            
            summary_text = sprintf('Converged in %d iterations (%.4f seconds).\nTotal Real Power Loss: %.4f MW', iter, solveTime, total_loss_MW);
            app.LoadFlowSummaryLabel.Text = summary_text;
        end

        % Function to build Zbus matrices
        function buildZbusMatrices(app, Ybus1, R1, X1, taps, fromBus, toBus)
            app.Zbus1 = inv(Ybus1);
            app.Zbus2 = app.Zbus1;
            
            R0 = R1;
            X0 = 3 * X1;
            Z0_series = R0 + 1i * X0;
            y0_series = 1 ./ Z0_series;
            Ybus0 = zeros(app.numBuses, app.numBuses);
            numLines = length(fromBus);
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
            busType = app.busData(:,2);
            genBuses = find(busType == 2 | busType == 1);
            Xg0 = 0.05;
            for i = 1:length(genBuses)
                bus_idx = genBuses(i);
                Ybus0(bus_idx, bus_idx) = Ybus0(bus_idx, bus_idx) + 1/(1i*Xg0);
            end
            app.Zbus0 = pinv(Ybus0);
        end
        
        % Function to display fault results
        function displayFaultResults(app, Ia0, Ia1, Ia2, I_phase, I_fault_mag_kA, MVA_symmetrical)
            seq_text = {
                sprintf('I0: %.4f /_ %.2f deg', abs(Ia0), rad2deg(angle(Ia0)));
                sprintf('I1: %.4f /_ %.2f deg', abs(Ia1), rad2deg(angle(Ia1)));
                sprintf('I2: %.4f /_ %.2f deg', abs(Ia2), rad2deg(angle(Ia2)))
            };
            app.SequenceCurrentsTextArea.Value = seq_text;
            
            phase_text = {
                sprintf('Ia: %.4f /_ %.2f deg', abs(I_phase(1)), rad2deg(angle(I_phase(1))));
                sprintf('Ib: %.4f /_ %.2f deg', abs(I_phase(2)), rad2deg(angle(I_phase(2))));
                sprintf('Ic: %.4f /_ %.2f deg', abs(I_phase(3)), rad2deg(angle(I_phase(3))))
            };
            app.PhaseCurrentsTextArea.Value = phase_text;
            
            breaker_text = {
                sprintf('Total Fault Current: %.4f kA', I_fault_mag_kA);
                sprintf('Symmetrical Interrupting MVA: %.2f MVA', MVA_symmetrical)
            };
            app.BreakerRatingTextArea.Value = breaker_text;
        end
        
    end
    
    % App initialization and construction
    methods (Access = public)

        % Construct app
        function app = PowerSystemAnalyzerGUI
            % Create components
            createComponents(app)

            % Register the app object
            registerApp(app, app.UIFigure)
            
            % No startup function needed for this app
            % runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
