function fmcw_radar_sim
    clc; close all;  % Clear screen and close figures, set up basic timing
    %% ---------------- Constants & Config ----------------
   
    c = physconst('LightSpeed');   
    dt = 0.05; isRunning = false; t = 0;
    
    % Basic FMCW settings (not super strict, just for demo)
    fc = 10e9; Pt = 100; Gt = 30; BW = 150e6; Tchirp = 1e-2;
    lambda = c / fc; scanPeriod = 2.0;
    
    % Chirp slope (how fast freq changes)
    slope_S = BW / Tchirp;
    
    %% ---------------- Target Initial ----------------
    % Start target a bit far, with some heading and speed
    R = 70e3; theta = deg2rad(30); v = 200/3.6; psi = deg2rad(150);
    targetRadius = 2; 
    
    %% ---------------- Figure Setup ----------------
    % Make a figure with two panels: map and freq-time
    fig = figure('Name', 'FMCW Radar: Frequency-Time Analysis', ...
        'NumberTitle', 'off', 'Position', [100 100 1200 650]);
    
    % Left panel: simple map view
    axMap = subplot(1, 2, 1);
    limitVal = 120e3;
    axis(axMap, [-limitVal limitVal -limitVal limitVal]);
    axis(axMap, 'equal', 'manual'); grid(axMap, 'on'); hold(axMap, 'on');
    title(axMap, 'Radar Map (Fixed Center)');
    
    % Beam line, target marker, and trail (initially empty)
    hBeam = plot(axMap, [0 0], [0 0], 'g-', 'LineWidth', 1.2);
    hTarget = plot(axMap, NaN, NaN, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hTrail = plot(axMap, NaN, NaN, 'r--');
    plot(axMap, 0, 0, 'ks', 'MarkerFaceColor', 'k'); % radar center
    
    % Right Panel: show TX and RX frequency vs time
    axFreq = subplot(1, 2, 2);
    hold(axFreq, 'on'); grid(axFreq, 'on');
    xlabel(axFreq, 'Time [ms]'); ylabel(axFreq, 'Frequency [GHz]');
    title(axFreq, 'FMCW Triangular Chirp (TX vs RX)');
    
    hTxFreq = plot(axFreq, NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName', 'TX Signal');
    hRxFreq = plot(axFreq, NaN, NaN, 'r--', 'LineWidth', 1.5, 'DisplayName', 'RX Signal');
    legend(axFreq, 'Location', 'northeast');
    
    %% ---------------- UI Controls ----------------
    % Small control panel to change target and scan settings
    uip = uipanel('Title','Target & Radar Settings','Position',[.05 .01 .9 .15]);
    
    editR = createEdit(uip, '70', 10, 20);      createLabel(uip, 'Dist [km]', 10, 45);
    editTheta = createEdit(uip, '30', 90, 20);  createLabel(uip, 'Azim [deg]', 90, 45);
    editV = createEdit(uip, '200', 170, 20);    createLabel(uip, 'Vel [km/h]', 170, 45);
    editPsi = createEdit(uip, '150', 250, 20);  createLabel(uip, 'Course [deg]', 250, 45);
    editScan = createEdit(uip, '2.0', 330, 20); createLabel(uip, 'Scan [s]', 330, 45);
    
    btnStart = uicontrol(uip, 'Style','togglebutton','String','Start',...
        'Position', [420 20 100 40], 'Callback', @startStopCallback);
        
    %% ---------------- Main Loop ----------------
    % Trail arrays to show where target was
    trailX = nan(1, 5000); trailY = nan(1, 5000); trailIdx = 0;
    
    while ishandle(fig)
        if isRunning
            t = t + dt;
            
            % 1. Rotate the beam and move the target a bit
            beamAngle = mod(2 * pi * (t / scanPeriod), 2 * pi);
            set(hBeam, 'XData', [0 limitVal*cos(beamAngle)], 'YData', [0 limitVal*sin(beamAngle)]);
            
            % Simple kinematics: move target by small step 
            vx = v * cos(psi); vy = v * sin(psi);
            x_real = R * cos(theta) + vx * dt;
            y_real = R * sin(theta) + vy * dt;
            
            % Update polar coords from new position
            R = sqrt(x_real^2 + y_real^2);
            theta = atan2(y_real, x_real);
            
            % 2. Check if beam is pointing at the target (rough check)
            targetAngleNorm = mod(theta, 2*pi);
            angleDiff = abs(beamAngle - targetAngleNorm);
            if angleDiff > pi, angleDiff = 2*pi - angleDiff; end
            
            if angleDiff < (2 * pi * dt / scanPeriod)
                % Beam "hits" target, do some FMCW math (simple model)
                tau = 2 * R / c; % round-trip delay
                vr = v * cos(psi - theta); % radial speed (approx)
                fd = -2 * vr / lambda;     % Doppler shift (signs are messy sometimes)
                
                %% -------- BEAT FREQUENCY CALCULATION --------
                % Sampling for beat signal 
                fs_beat = 20e6;                 
                
                % time vector for one chirp (baseband model)
                t_beat = 0 : 1/fs_beat : Tchirp - 1/fs_beat;
                
                % TX baseband phase (quadratic because of chirp)
                phase_tx = pi * slope_S * t_beat.^2;
                tx_sig = exp(1j * phase_tx);
                
                % RX: delayed and Doppler-shifted (simple approx)
                phase_rx = pi * slope_S * (t_beat - tau).^2 + 2*pi*fd*t_beat;
                rx_sig = exp(1j * phase_rx);
                
                % Mix TX and RX to get beat signal (dechirp)
                beat_sig = tx_sig .* conj(rx_sig);
                
                % FFT for beat spectrum (use next power of 2 for speed)
                N = length(beat_sig);
                N_FFT = 2^nextpow2(N);
                
                % Apply a Hann window to make peak nicer, then FFT
                win = hann(N)'; 
                BEAT_FFT = abs(fft(beat_sig .* win, N_FFT));
                
                % Frequency axis for the FFT
                f_axis = (0:N_FFT-1)*(fs_beat/N_FFT);
                
                % Find the main peak in the first half (positive freqs)
                [~, idx] = max(BEAT_FFT(1:floor(N_FFT/2)));
                f_beat = f_axis(idx);

                %% --- Estimate range and radial velocity from peaks ---
                % Range from beat freq (simple linear relation)
                R_est = (c * f_beat) / (2 * slope_S);
                
                % Radial velocity from Doppler (signs may be flipped)
                v_r_est = -fd * lambda / 2;
                
                %% --- Rough received power estimate ---
                % Radar equation, just to show numbers
                sigma = pi * targetRadius^2;
                Gt_lin = 10^(Gt/10);
                Pr = (Pt * Gt_lin^2 * lambda^2 * sigma) / ((4*pi)^3 * R_est^4);   

                % ---------------------------------------------------------
                
                % Update map marker and trail
                set(hTarget, 'XData', x_real, 'YData', y_real);
                trailIdx = trailIdx + 1;
                if trailIdx <= 5000
                    trailX(trailIdx) = x_real; trailY(trailIdx) = y_real;
                end
                set(hTrail, 'XData', trailX(1:trailIdx), 'YData', trailY(1:trailIdx));
                
                % --- Update chirp visualization (low res so it's smooth) ---
                t_plot = linspace(0, 2*Tchirp, 1000);
                freq_tx = zeros(size(t_plot));
                
                for i = 1:length(t_plot)
                    t_mod = mod(t_plot(i), 2*Tchirp);
                    if t_mod <= Tchirp
                        freq_tx(i) = fc + slope_S * t_mod;
                    else
                        freq_tx(i) = (fc + BW) - slope_S * (t_mod - Tchirp);
                    end
                end
                
                % RX freq is just TX shifted by Doppler for the plot
                freq_rx = freq_tx + fd; 
                shift_samples = round(tau / (t_plot(2)-t_plot(1)));
                % Put NaNs in front to show delay visually
                if shift_samples > 0 && shift_samples < length(freq_rx)
                     freq_rx = [nan(1, shift_samples), freq_rx(1:end-shift_samples)];
                end

                set(hTxFreq, 'XData', t_plot*1000, 'YData', freq_tx/1e9);
                set(hRxFreq, 'XData', t_plot*1000, 'YData', freq_rx/1e9);
                
                % Keep y-limits reasonable so plot doesn't jump around
                ylim(axFreq, [(fc - 0.1*BW)/1e9, (fc + 1.1*BW)/1e9]);
                
                % Show some numbers in the map title (quick debug info)
                title(axMap, sprintf( ...
                    'Detected | f_{beat} = %.1f kHz| R = %.2f km | v_r = %.1f m/s | P_r=%.2e W', ...
                    f_beat/1e3, R_est/1e3, v_r_est, Pr ));
            end
        end
        drawnow limitrate;
        pause(dt);
    end
    
    %% ---------------- Callbacks & Helpers ----------------
    function startStopCallback(src, ~)
        isRunning = src.Value;
        if isRunning
            % Read all the UI fields when starting (user inputs)
            R = str2double(editR.String)*1e3;
            theta = deg2rad(str2double(editTheta.String));
            v = str2double(editV.String)/3.6;
            psi = deg2rad(str2double(editPsi.String));
            scanPeriod = str2double(editScan.String);
            
            % Reset time and trail
            t = 0; trailIdx = 0; trailX(:)=NaN; trailY(:)=NaN;
            set(hTarget, 'XData', NaN);
            src.String = 'Stop';
        else
            src.String = 'Start';
        end
    end
    function h = createEdit(parent, def, x, y)
        % Tiny helper to make edit boxes
        h = uicontrol(parent, 'Style','edit','String',def,'Position',[x y 70 25]);
    end
    function createLabel(parent, txt, x, y)
        % Tiny helper to make labels
        uicontrol(parent, 'Style','text','String',txt,'Position',[x y 70 15],'FontSize',8);
    end
end
