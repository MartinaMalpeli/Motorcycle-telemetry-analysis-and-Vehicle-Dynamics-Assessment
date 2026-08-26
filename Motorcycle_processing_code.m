%% MOTORCYCLE DYNAMIC ANALYSIS: TIME-LAT/LON SEGMENTATION & WEAVE STABILITY
%  Course: Vehicle Dynamics / Experimental Mechanics
%  Description: 
%     This script performs a comprehensive analysis of motorcycle dynamics data.
%     The workflow includes:
%     1. Pre-processing and High-Pass filtering to remove signal drift.
%     2. GPS-based segmentation to separate Outbound and Return runs.
%     3. Interactive peak selection for identifying specific oscillation events.
%     4. Advanced Modal Analysis: FFT for frequency identification and 
%        Hilbert Transform for damping ratio extraction.
%     5. Statistical clustering by speed (Mean +/- Std Dev) and Gauss Diagram visualization.
%     6. Spectral comparison between three speed values.

clear; clc; close all;

%% 1. DATA LOADING AND PRE-PROCESSING
filename = 'experimental_data.csv';

% --- Configuration Parameters ---
MIN_RUN_DURATION = 30;      % Minimum seconds to consider a valid run
THRESHOLD_MULT = 3.5;       % Multiplier for trigger threshold
MIN_EVENT_DIST_SEC = 4;     % Minimum time gap between trigger events

disp('1. Loading Data...');
if ~exist(filename, 'file')
    error('File not found. Please check the filename or path.');
end

opts = detectImportOptions(filename);
opts.Delimiter = ';'; 
opts.VariableNamesLine = 1; 
opts.DataLines = [3 Inf];
opts.PreserveVariableNames = true;
T = readtable(filename, opts);

% Cleaning and Conversion (String to Double if necessary)
varNames = T.Properties.VariableNames;
for i = 1:length(varNames)
    if iscell(T.(varNames{i})) || isstring(T.(varNames{i}))
        T.(varNames{i}) = str2double(strrep(string(T.(varNames{i})), ',', '.'));
    end
end

% Variable Assignment
Time = T.Time; 
V_Sat = T.V_Sat; V_Sat(isnan(V_Sat))=0;
Gyro_Z = T.GYRO_Z; 
Lon = T.Longitude; Lon(isnan(Lon))=0; 
Lat = T.Latitude; Lat(isnan(Lat))=0;
Steer = (T.GYRO_STEER - 2.5) * 40; 

% Physical Quantities Calculation
Fs = 1 / mean(diff(Time));
V_ms = V_Sat / 3.6; 
Acc_Raw = gradient(V_ms) ./ gradient(Time); 

%% 2. SIGNAL FILTERING (HIGH-PASS > 0.5 Hz)
disp('2. Signal Filtering (High-Pass > 1 Hz)...');
% A high-pass filter is applied to remove DC offset and slow drift,
% facilitating peak detection.
Fc = 2;             % Cut-off frequency [Hz]
[b_hp, a_hp] = butter(2, Fc/(Fs/2), 'high');

Steer_HP = filtfilt(b_hp, a_hp, Steer);
Gyro_Z_HP = filtfilt(b_hp, a_hp, Gyro_Z);

%% 3. PRELIMINARY SEGMENTATION (GRADIENT METHOD)
disp('3. Performing Preliminary Segmentation...');
V_Proc = smoothdata(V_Sat, 'movmean', 20);

% Calculate radial distance gradient from start point
% Positive gradient = Moving away (Outbound), Negative = Coming back (Return)
Dist_Rate = smoothdata(gradient(sqrt(((Lat-Lat(1))*111320).^2 + ((Lon-Lon(1))*40075e3*cosd(Lat(1))/360).^2))./gradient(Time), 'movmean', 300);

Is_Moving = fillmissing(V_Proc > 5, 'previous');
diff_mov = diff([0; Is_Moving; 0]);
idx_s = find(diff_mov == 1); 
idx_e = find(diff_mov == -1)-1;

Run_ID_Vector = zeros(size(Time));
glob_cnt = 0;

for k = 1:length(idx_s)
    if (Time(idx_e(k)) - Time(idx_s(k))) > MIN_RUN_DURATION
        blk_rate = Dist_Rate(idx_s(k):idx_e(k));
        blk_ids = zeros(size(blk_rate));
        curr_st = 0; 
        for t=1:length(blk_rate)
            if blk_rate(t)>0.5, n_st=1; elseif blk_rate(t)<-0.5, n_st=2; else, n_st=curr_st; end
            if n_st~=0 && n_st~=curr_st
                glob_cnt = glob_cnt+1;
                % Logic to ensure odd IDs are Outbound, even IDs are Return
                if (n_st==1 && mod(glob_cnt,2)==0) || (n_st==2 && mod(glob_cnt,2)~=0), glob_cnt=glob_cnt+1; end
            end
            curr_st=n_st; if curr_st~=0, blk_ids(t)=glob_cnt; end
        end
        Run_ID_Vector(idx_s(k):idx_e(k)) = blk_ids;
    end
end

%% 4. IMPULSE DETECTION (TRIGGER)
disp('4. Detecting Impulses...');
Steer_Diff = [0; abs(diff(smoothdata(Steer,'movmean',10)))];
[Trigger_Pks, locs] = findpeaks(Steer_Diff, 'MinPeakHeight', mean(Steer_Diff)+THRESHOLD_MULT*std(Steer_Diff), 'MinPeakDistance', round(MIN_EVENT_DIST_SEC/0.04));
Events_Speed = V_Sat(locs);

%% 5. ROBUST RADIAL SEGMENTATION & GEOPLOT
disp('5. Radial Segmentation & Geoplotting...');
deg_to_m_lat = 111320; 
deg_to_m_lon = 40075000 * cosd(Lat(1)) / 360;
dist_radial = sqrt(((Lat - Lat(1)) * deg_to_m_lat).^2 + ((Lon - Lon(1)) * deg_to_m_lon).^2); 
dist_smooth = smoothdata(dist_radial, 'gaussian', 500);

min_prom = max(dist_smooth) * 0.3; 
[~, idx_ta] = findpeaks(dist_smooth, 'MinPeakProminence', min_prom);
[~, idx_base] = findpeaks(-dist_smooth, 'MinPeakProminence', min_prom);
cut_pts = sort(unique([1; idx_ta; idx_base; length(Time)]));

segments = [];
for k = 1:length(cut_pts)-1
    is = cut_pts(k); ie = cut_pts(k+1);
    delta = dist_smooth(ie) - dist_smooth(is);
    if abs(delta) > (max(dist_smooth)*0.2)
        % 1 = OUTBOUND (Go), 2 = RETURN
        segments = [segments; is, ie, (delta>0)*1 + (delta<=0)*2]; 
    end
end
num_plots = min(6, size(segments, 1));

% --- Longitudinal Trajectory ---
figure('Name','Longitude-Time Trajectory','Color','w'); 
plot(Time, Lon, 'r'); 
xlabel('Time [s]'); ylabel('Longitude [deg]'); 
grid on;

% ---  Detection Outbound & Return ---
figure('Name', 'Detection Outbound & Return', 'Color', 'w');
for k = 1:num_plots
    idx = segments(k,1):segments(k,2);
    subplot(3, 2, k);
    
    geoscatter(Lat(idx), Lon(idx), 15, Time(idx), 'filled');
    try geobasemap('satellite'); catch, end
    
    colormap('jet'); 
    caxis([0 max(Time)]);
    
    if segments(k,3) == 1
        lbl = ": GO";
    else
        lbl = ": RETURN";
    end
    title(['Segment ' num2str(k) char(lbl)]);
end
c = colorbar('Position', [0.93 0.1 0.02 0.8]); 
c.Label.String = 'Global Time [s]';

%% 6. INTERACTIVE ZONE SELECTION
fprintf('\n--- STARTING INTERACTIVE SELECTION (Fig 6) ---\n');
fprintf('Instructions:\n1. Click 2 diagonal points to draw a box around the oscillation.\n2. Confirm via dialog.\n');

Segment_Analysis = struct('Tratto', [], 'Type', [], 'PeakTimes', [], 'PeakAmps', []);

for k = 1:num_plots
    % Setup Data
    idx = segments(k,1):segments(k,2);
    t_seg = Time(idx); v_seg = V_Sat(idx); st_seg = Steer(idx); gz_seg = Gyro_Z(idx);
    
    if segments(k,3) == 1
        type_str = "OUTBOUND";
    else
        type_str = "RETURN";
    end

    hFig = figure('Name', ['Zone Analysis Segment ' num2str(k)], 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    
    subplot(2,1,1); plot(t_seg, v_seg, 'k'); grid on; title(['Speed Profile - Segment ' num2str(k)]); xlim([min(t_seg) max(t_seg)]);
    ax_sig = subplot(2,1,2); hold on; grid on;
    plot(t_seg, st_seg, 'r'); plot(t_seg, gz_seg, 'b'); legend('Steer','Gyro Z');
    title({['Segment ' num2str(k) ' (' char(type_str) ')']; '\color{blue}Click 2 points (Diagonal) to select area. Use dialog box to confirm.'});
    xlim([min(t_seg) max(t_seg)]);

    all_peaks_t = []; all_peaks_a = [];
    selecting_active = true;

    while selecting_active
        try [tx, ty] = ginput(2); catch, selecting_active=false; break; end
        
        % Check for valid input
        if length(tx) < 2
            continue; 
        end
        
        t_start = min(tx); t_stop = max(tx); thresh = min(abs(ty));
        
        % Temporary Graphics
        hRect = rectangle('Position', [t_start, -thresh, t_stop-t_start, 2*thresh], 'EdgeColor', 'm', 'LineStyle', '--');
        hLine = yline([thresh, -thresh], 'm:');
        
        % Find Peaks within the box
        mask = t_seg >= t_start & t_seg <= t_stop;
        
        % Safety check
        if sum(mask) == 0
            delete(hRect); delete(hLine);
            continue;
        end
        
        warning('off', 'signal:findpeaks:largeMinPeakHeight');
        [~, l_st] = findpeaks(abs(st_seg(mask)), 'MinPeakHeight', thresh);
        [~, l_gz] = findpeaks(abs(gz_seg(mask)), 'MinPeakHeight', thresh);
        warning('on', 'signal:findpeaks:largeMinPeakHeight');
        
        sub_t = t_seg(mask); sub_st = st_seg(mask); sub_gz = gz_seg(mask);
        cand_t = [sub_t(l_st); sub_t(l_gz)]; cand_a = [sub_st(l_st); sub_gz(l_gz)];
        
        if isempty(cand_t)
            title('\color{red}No peaks found above threshold in this area!');
            pause(1);
            delete(hRect); delete(hLine);
            title({['Segment ' num2str(k) ' (' char(type_str) ')']; '\color{blue}Click 2 points (Diagonal) to select area.'});
            continue;
        end
        
        hCand = plot(cand_t, cand_a, 'co', 'LineWidth', 2);
        
        % Dialog Box
        choice = questdlg(sprintf('Found %d peaks. Accept?', length(cand_t)), 'Confirm Selection', 'Save', 'Retry', 'Finish Segment', 'Save');
        
        switch choice
            case 'Save'
                delete(hRect); delete(hLine); delete(hCand);
                
                all_peaks_t = [all_peaks_t; cand_t]; all_peaks_a = [all_peaks_a; cand_a];
                plot(sub_t(l_st), sub_st(l_st), 'ro', 'MarkerFaceColor','r');
                plot(sub_t(l_gz), sub_gz(l_gz), 'bo', 'MarkerFaceColor','b');
                
            case 'Retry'
                delete(hRect); delete(hLine); delete(hCand);
                
            case 'Finish Segment'
                if ~isempty(cand_t)
                    if strcmp(questdlg('Save the last selection?','Exit','Yes','No','Yes'),'Yes')
                        all_peaks_t=[all_peaks_t; cand_t]; all_peaks_a=[all_peaks_a; cand_a];
                    end
                end
                selecting_active = false;
                
            otherwise
                selecting_active = false;
        end
    end
    Segment_Analysis(k) = struct('Tratto', k, 'Type', type_str, 'PeakTimes', all_peaks_t, 'PeakAmps', all_peaks_a);
    close(hFig);
end

%% 7. ADVANCED ANALYSIS: MULTI-EVENT HILBERT (CLUSTERING) - CORRECTED
fprintf('\n-------------------------------------------------\n');
fprintf('ADVANCED ANALYSIS: Separazione Eventi Multipli & Hilbert Ibrida\n');

% --- PARAMETRI CHIAVE ---
MIN_GAP_SEC = 1.5;      % Se due picchi distano più di questo, sono eventi diversi
PAD_DUR = 2.0;          % Secondi di "padding" per stabilizzare la Hilbert ai bordi
Fs = 1/mean(diff(Time));

% Filtro Weave Standard
[b_wea, a_wea] = butter(2, [0.4 5.0]/(Fs/2), 'bandpass'); 

Global_Speed = []; 
Global_Damp  = []; 
Global_Freq  = []; 
Global_Time  = [];

% Prepariamo la figura per i risultati
figure('Name', 'Multi-Event Analysis Results', 'Color', 'w');
plot_idx = 1;
max_plots = 9; % Numero massimo di grafici da mostrare

for k = 1:num_plots
    % Recuperiamo i picchi che hai selezionato manualmente nel blocco precedente
    if isempty(Segment_Analysis(k).PeakTimes), continue; end
    
    % --- 1. CLUSTERING DEGLI EVENTI ---
    % Ordiniamo i picchi temporalmente
    [p_times_sorted, sort_idx] = sort(Segment_Analysis(k).PeakTimes);
    
    % Calcoliamo la distanza temporale tra picchi consecutivi
    dt_peaks = diff(p_times_sorted);
    
    % Se la distanza > MIN_GAP_SEC, è un nuovo evento.
    % Creiamo un vettore ID (es: [1 1 1 2 2 2 3 3])
    cluster_ids = [1; cumsum(dt_peaks > MIN_GAP_SEC) + 1];
    unique_clusters = unique(cluster_ids);
    
    % Dati completi del tratto corrente
    idx_tratto = segments(k,1):segments(k,2);
    t_full_tratto = Time(idx_tratto); 
    gz_full_tratto = Gyro_Z_HP(idx_tratto); 
    v_full_tratto = V_Sat(idx_tratto);
    
    % --- 2. CICLO SU OGNI SOTTO-EVENTO (BURST) ---
    for i_cls = 1:length(unique_clusters)
        cls_id = unique_clusters(i_cls);
        
        % Tempi dei picchi appartenenti a questo singolo burst
        current_peaks_t = p_times_sorted(cluster_ids == cls_id);
        
        % Se ho meno di 2 picchi selezionati in questo cluster, è rumore o errore -> Salta
        if length(current_peaks_t) < 2, continue; end
        
        % Definiamo una finestra temporale "larga" attorno a questo evento
        % per dare spazio alla Hilbert di assestarsi
        t_start_cut = min(current_peaks_t) - PAD_DUR;
        t_stop_cut  = max(current_peaks_t) + PAD_DUR;
        
        mask_ev = t_full_tratto >= t_start_cut & t_full_tratto <= t_stop_cut;
        if sum(mask_ev) < 50, continue; end
        
        % Estrazione dati evento singolo
        t_ev = t_full_tratto(mask_ev);
        gz_ev = gz_full_tratto(mask_ev);
        v_ev_mean = mean(v_full_tratto(mask_ev));
        
        % Zero temporale locale per i grafici
        t_ev_zero = t_ev - t_ev(1);
        
        % --- 3. ELABORAZIONE SEGNALE (FILTRO + HILBERT) ---
        gz_clean = filtfilt(b_wea, a_wea, gz_ev);
        
        % Hilbert Transform
        analytic = hilbert(gz_clean);
        env_raw = abs(analytic);
        
        % Smoothing dell'inviluppo (cruciale per togliere il ripple)
        env_smooth = smoothdata(env_raw, 'gaussian', 15);
        
        % --- 4. IDENTIFICAZIONE PUNTO DI START (Hybrid Method) ---
        % Cerchiamo il picco MASSIMO del segnale reale (non dell'inviluppo)
        % all'interno della zona centrale (escludiamo il padding ai bordi)
        padding_samples = round(PAD_DUR * Fs * 0.5); % Margine di sicurezza in campioni
        valid_range = (padding_samples) : (length(gz_clean) - padding_samples);
        
        % Fix per array vuoti o indici non validi
        if isempty(valid_range) || valid_range(end) > length(gz_clean)
             valid_range = 1:length(gz_clean); 
        end
        
        [max_val, idx_max_rel] = max(abs(gz_clean(valid_range)));
        idx_max = idx_max_rel + valid_range(1) - 1;
        
        % --- 5. FITTING SULLA CODA DI DECADIMENTO ---
        % Prendiamo l'inviluppo dal picco massimo in poi
        y_decay = env_smooth(idx_max:end);
        t_decay = t_ev_zero(idx_max:end);
        
        % Tagliamo quando il segnale diventa rumore (sotto il 10% del picco)
        noise_floor = max_val * 0.10;
        mask_valid = y_decay > noise_floor;
        
        y_fit = y_decay(mask_valid);
        t_fit = t_decay(mask_valid);
        
        % CORREZIONE QUI: Sostituito range() con max-min
        if (max(t_fit) - min(t_fit)) < 0.3, continue; end
        
        % Stima Frequenza (FFT sulla parte attiva)
        L = length(gz_clean); nfft = 2^nextpow2(L);
        f_ax = (Fs/2)*linspace(0,1,nfft/2+1);
        Y = fft(gz_clean, nfft);
        [~, imax] = max(abs(Y(1:nfft/2+1)));
        freq_dom = f_ax(imax);
        
        % Se la frequenza è fuori range moto (troppo bassa o alta), scarta
        if freq_dom < 0.5 || freq_dom > 5, continue; end
        
        % REGRESSIONE LOGARITMICA: ln(y) = -sigma*t + C
        p = polyfit(t_fit, log(y_fit), 1);
        sigma = -p(1);
        
        wn = sqrt(sigma^2 + (2*pi*freq_dom)^2);
        damp_perc = (sigma / wn) * 100;
        
        % Filtro risultati assurdi (Damping negativo o > 50%)
        if damp_perc > 0 && damp_perc < 50
            % Salvataggio Risultati Globali
            Global_Speed = [Global_Speed; v_ev_mean];
            Global_Damp  = [Global_Damp; damp_perc];
            Global_Freq  = [Global_Freq; freq_dom];
            Global_Time  = [Global_Time; t_ev(idx_max)]; % Tempo assoluto evento
            
            % --- PLOTTING ---
            if plot_idx <= max_plots
                subplot(3, 3, plot_idx); 
                hold on; grid on;
                
                % Segnale Reale
                plot(t_ev_zero, gz_clean, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
                
                % Inviluppo Hilbert
                plot(t_ev_zero, env_smooth, 'b', 'LineWidth', 1.5);
                
                % Ricostruzione del Fit (Rosso tratteggiato)
                % Disegnata solo sulla parte usata per il calcolo
                y_recon = exp(polyval(p, t_fit));
                plot(t_fit, y_recon, 'r--', 'LineWidth', 2);
                
                % Marker inizio decadimento
                plot(t_ev_zero(idx_max), env_smooth(idx_max), 'ro', 'MarkerFaceColor', 'r');
                
                title(sprintf('V=%.0f | D=%.1f%%', v_ev_mean, damp_perc), 'FontSize', 9);
                xlim([min(t_fit)-0.5, max(t_fit)+0.5]);
                
                plot_idx = plot_idx + 1;
            end
        end
    end
end
%% 8. SPEED vs DAMPING SCATTER (Weave Only)
mask_weave = ~isnan(Global_Damp);

if sum(mask_weave) > 0
    figure('Name', 'Weave Stability Map', 'Color', 'w');
    
    scatter(Global_Speed(mask_weave), Global_Damp(mask_weave), 20, Global_Time(mask_weave), 'filled', 'MarkerEdgeColor', 'k');
    
    colormap('jet');
    c = colorbar;
    c.Label.String = 'Absolute Time [s]';
    caxis([0 max(Time)]);
    
    grid on;
    xlabel('Riding Speed [km/h]');
    ylabel('Weave Damping Ratio [%]');
    title('Weave Stability Map (Decay Fit)');
    
    xlim([0 max(V_Sat)+5]);
    ylim([0 50]); 
else
    disp('No Weave events analyzable for damping.');
end

%% 9. SPEED vs FREQUENCY SCATTER 
if ~isempty(Global_Freq)
    figure('Name', 'Speed vs Frequency Map', 'Color', 'w');
    
    scatter(Global_Speed, Global_Freq, 20, Global_Time, 'filled', 'MarkerEdgeColor', 'k');
    
    colormap('jet');
    c = colorbar;
    c.Label.String = 'Absolute Time [s]';
    caxis([0 max(Time)]);
    
    grid on;
    xlabel('Riding Speed [km/h]');
    ylabel('Detected Frequency [Hz]');
    title('Frequency Map ');
    
    % Highlight Weave Band
    yline(2, 'k--'); yline(4.5, 'k--');
    %text(5, 3.25, 'Standard Weave Band', 'Color', 'k','Location', 'Best');
    
    xlim([0 max(V_Sat)+5]);
    ylim([0 12]); 
else
    disp('No frequencies detected.');
end

%% 10. CLUSTERED STATISTICAL ANALYSIS 
if ~isempty(Global_Speed)
    fprintf('\n======================================================\n');
    fprintf('   STATISTICAL TRENDS (MEAN +/- STD DEV)    \n');
    fprintf('======================================================\n');
    
    % Filter valid data (Weave only)
    mask_stat = ~isnan(Global_Damp);
    stat_speed = Global_Speed(mask_stat);
    stat_damp  = Global_Damp(mask_stat);
    stat_freq  = Global_Freq(mask_stat);
    
    % Define Clusters: Fixed 10 km/h steps starting from 0
    max_v = ceil(max(stat_speed)/10)*10;
    if max_v < 10, max_v = 10; end
    edges = 0:10:max_v;
    
    [~, ~, bins] = histcounts(stat_speed, edges);
    bin_centers = (edges(1:end-1) + edges(2:end))' / 2;
    
    % Pre-allocation
    num_bins = length(bin_centers);
    stats_mean_damp = NaN(num_bins, 1); stats_std_damp  = NaN(num_bins, 1);
    stats_mean_freq = NaN(num_bins, 1); stats_std_freq  = NaN(num_bins, 1);
    stats_count     = zeros(num_bins, 1);
    
    % Computation loop
    for b = 1:max(bins)
        idx_bin = (bins == b);
        count = sum(idx_bin);
        
        if count > 0
            stats_count(b)      = count;
            stats_mean_damp(b)  = mean(stat_damp(idx_bin));
            stats_mean_freq(b)  = mean(stat_freq(idx_bin));
            
            if count > 1
                stats_std_damp(b) = std(stat_damp(idx_bin));
                stats_std_freq(b) = std(stat_freq(idx_bin));
            else
                stats_std_damp(b) = 0;
                stats_std_freq(b) = 0;
            end
        end
    end
    
    % --- 1. PRINT TABLE ---
    StatsTable = table(bin_centers, stats_count, stats_mean_freq, stats_std_freq, stats_mean_damp, stats_std_damp, ...
        'VariableNames', {'Speed_Center', 'Count', 'Freq_Mean', 'Freq_Std', 'Damp_Mean', 'Damp_Std'});
    disp(StatsTable);
    
    % --- 2. PLOT ERRORBARS ---
    figure('Name', 'Statistical Trends vs Speed', 'Color', 'w');
    valid_bins = stats_count > 0;
    x_plot = bin_centers(valid_bins);
    
    % Subplot 1: Frequency Trend
    subplot(2,1,1);
    errorbar(x_plot, stats_mean_freq(valid_bins), stats_std_freq(valid_bins), '-s', 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', [0.2 0.6 0.2], 'Color', [0.2 0.6 0.2], 'CapSize', 8);
    grid on; 
    ylabel('Frequency [Hz]', 'FontWeight', 'bold');
    title('Frequency vs Riding Speed (Mean \pm Std Dev)');
    xlim([0 max_v]); ylim([0 6]); 
    
    % Subplot 2: Damping Trend
    subplot(2,1,2);
    errorbar(x_plot, stats_mean_damp(valid_bins), stats_std_damp(valid_bins), '-s', 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'Color', 'b', 'CapSize', 8);
    grid on; 
    xlabel('Riding Speed [km/h]', 'FontWeight', 'bold'); 
    ylabel('Damping Ratio [%]', 'FontWeight', 'bold'); 
    title('Modal Damping vs Riding Speed (Mean \pm Std Dev)');
    xlim([0 max_v]); 
    ylim([0 50]); 
end

%% 11. GAUSS DIAGRAM & DECAY RATE

try
    vel_vec = Global_Speed;   
    freq_vec = Global_Freq;   % Già in Hz
    damp_vec = Global_Damp;   
    
    mask_valid = ~isnan(freq_vec) & ~isnan(damp_vec) & ~isnan(vel_vec);
    vel_vec = vel_vec(mask_valid);
    freq_vec = freq_vec(mask_valid);
    damp_vec = damp_vec(mask_valid);
catch
    warning('Variabili Globali non trovate. Esegui prima la Sezione 7.');
    return;
end

if isempty(vel_vec)
    disp('Nessun dato valido per il Root Locus.');
else
    % 1. Conversione in Hz
    zeta = damp_vec / 100;
    % Calcoliamo Sigma in Hz: sigma_hz = (zeta * freq_hz) / sqrt(1 - zeta^2)
    sigma_hz = zeros(size(zeta));
    for k = 1:length(zeta)
        z = zeta(k);
        if z < 0.99
            % freq_vec(k) è già in Hz, quindi sigma_hz sarà in Hz
            sigma_hz(k) = (z * freq_vec(k)) / sqrt(1 - z^2);
        else
            sigma_hz(k) = NaN;
        end
    end
    
    % 2. Creazione Grafico
    figure('Name', 'Weave Root Locus - Hz Scale', 'Color', 'w');
    hold on; grid on; box on;
    
    real_part_hz = -sigma_hz; 
    imag_part_hz = freq_vec; % Hz
    
    % Scatter Plot
    scatter(real_part_hz, imag_part_hz, 20, vel_vec, 'filled', 'MarkerEdgeColor', 'k');
    
    % % 3. Linee di Riferimento (Iso-Damping) - Adattate per Hz
    % max_mod_hz = max(sqrt(real_part_hz.^2 + imag_part_hz.^2), [], 'omitnan') * 1.1;
    % if isempty(max_mod_hz) || max_mod_hz == 0, max_mod_hz = 5; end
    % 
    % zetas_ref = [0.05, 0.1, 0.15, 0.2]; 
    % for z = zetas_ref
    %     x_L = -z * max_mod_hz;
    %     y_L = sqrt(1 - z^2) * max_mod_hz;
    %     plot([0, x_L], [0, y_L], 'k:', 'LineWidth', 1);
    %     text(x_L, y_L, [num2str(z*100) '%'], 'FontSize', 8, 'HorizontalAlignment', 'right');
    % end
    
    % 4. Formattazione
    xline(0, 'r--', 'Stability Limit', 'LineWidth', 2);
    xlabel('Real Part (Damping) [Hz]', 'FontSize', 12);
    ylabel('Frequency (Imaginary Part) [Hz]', 'FontSize', 12);
    title('Weave Root Locus (Gauss Diagram in Hz)', 'FontSize', 14);
    
    c = colorbar;
    c.Label.String = 'Speed [km/h]';
    colormap(jet);
    
    % Limiti assi
    xlim([min(real_part_hz)*1.2, 0.5]); 
    ylim([0, max(imag_part_hz)*1.2]);
    hold off;
end
%% PLOT: Correct Experimental Root Locus (Gauss Diagram)
% ADATTATO: Usa le variabili Globali calcolate nella Sezione 7

% 1. Recupero Variabili Corrette (Global_...)
try
    % QUI C'ERA L'ERRORE: Ora punta alle variabili generate dalla Sez. 7
    vel_vec = Global_Speed;   % Velocità [km/h]
    freq_vec = Global_Freq;   % Frequenza [Hz]
    damp_vec = Global_Damp;   % Damping [%]
    
    % Filtriamo eventuali NaN per sicurezza
    mask_valid = ~isnan(freq_vec) & ~isnan(damp_vec) & ~isnan(vel_vec);
    vel_vec = vel_vec(mask_valid);
    freq_vec = freq_vec(mask_valid);
    damp_vec = damp_vec(mask_valid);
    
catch
    warning('Variabili Globali non trovate. Esegui prima la Sezione 7.');
    return;
end

if isempty(vel_vec)
    disp('Nessun dato valido per il Root Locus.');
else
    % 2. Ricostruzione di Sigma (Decay Rate) e Wd (Rad/s)
    % Formula inversa: sigma = (zeta * wd) / sqrt(1 - zeta^2)
    zeta = damp_vec / 100;           % Da % a decimale (0.1, 0.2...)
    wd_rad = freq_vec * 2 * pi;      % Da Hz a rad/s (Parte Immaginaria)
    
    sigma_rec = zeros(size(zeta));
    for k = 1:length(zeta)
        z = zeta(k);
        if z < 0.99 % Evita divisioni per zero o numeri complessi se zeta >= 1
            sigma_rec(k) = (z * wd_rad(k)) / sqrt(1 - z^2);
        else
            sigma_rec(k) = NaN;
        end
    end
    
    % 3. Creazione Grafico (Assi omogenei rad/s vs rad/s)
    figure('Name', 'Weave Root Locus - True Scale', 'Color', 'w');
    hold on; grid on; box on;
    
    % Asse X: Parte Reale (-Sigma)
    % Asse Y: Parte Immaginaria (Wd)
    real_part = -sigma_rec; 
    
    % Scatter Plot Colorato in base alla velocità
    scatter(real_part, wd_rad, 60, vel_vec, 'filled', 'MarkerEdgeColor', 'k');
    
    % 4. Linee di Riferimento (Iso-Damping Lines)
    % Calcoliamo il raggio massimo per disegnare le linee
    max_mod = max(sqrt(real_part.^2 + wd_rad.^2), [], 'omitnan');
    if isempty(max_mod) || max_mod == 0, max_mod = 20; end
    max_mod = max_mod * 1.1; 
    
    zetas_ref = [0.05, 0.1, 0.15, 0.2]; % Linee per 5%, 10%, 15%, 20%
    
    for z = zetas_ref
        % Coordinate linea: x = -z*R, y = sqrt(1-z^2)*R
        x_L = -z * max_mod;
        y_L = sqrt(1 - z^2) * max_mod;
        plot([0, x_L], [0, y_L], 'k:', 'LineWidth', 1);
        
        % Etichetta percentuale
        text(x_L, y_L, [num2str(z*100) '%'], ...
             'BackgroundColor', 'w', 'FontSize', 8, 'HorizontalAlignment', 'right');
    end
    
    % 5. Formattazione Finale
    xline(0, 'r--', 'Unstable Limit', 'LineWidth', 2);
    
    xlabel('Real Part ($\sigma$) [rad/s]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Imaginary Part ($\omega_d$) [rad/s]', 'Interpreter', 'latex', 'FontSize', 12);
    title('Weave Root Locus ', 'FontSize', 14);
    
    c = colorbar;
    c.Label.String = 'Speed [km/h]';
    colormap(jet);
    
    % FONDAMENTALE: axis equal forza X e Y ad avere la stessa scala
    axis normal;
    
    % Limiti intelligenti
    min_x = min(real_part);
    if isempty(min_x), min_x = -5; end
    xlim([min_x*1.2, 1]); 
    ylim([0, max(wd_rad)*1.1]);
    
    hold off;
end
%% 12. SPECTRAL COMPARISON AT USER-SELECTED SPEEDS (WEAVE MAGNITUDE)
if ~isempty(Global_Speed)
    fprintf('\n======================================================\n');
    fprintf('   SPECTRAL COMPARISON: 3 SPEED TARGETS \n');
    fprintf('======================================================\n');

    % 1. Input Dialog for 3 Speeds
    prompt = {'Target Speed 1 [km/h]:', 'Target Speed 2 [km/h]:', 'Target Speed 3 [km/h]:'};
    dlgtitle = 'Select Speeds for Frequency Analysis';
    dims = [1 40];
    % Default values based on data range
    definput = {'80', '120', '160'}; 
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    if isempty(answer)
        disp('Selection cancelled.');
    else
        target_speeds = str2double(answer);
        
        figure('Name', 'Weave Intensity Comparison (FFT)', 'Color', 'w');
        
        for i = 1:3
            t_spd = target_speeds(i);
            
            % 2. Find the closest analyzed event
            [diff_val, idx_best] = min(abs(Global_Speed - t_spd));
            
            % Check if the match is reasonable (e.g., within 5 km/h)
            if diff_val > 10
                subplot(3,1,i);
                text(0.5, 0.5, sprintf('No event found near %d km/h', t_spd), ...
                     'HorizontalAlignment', 'center');
                continue;
            end
            
            % 3. Extract Raw Signal for that specific event
            % We use Global_Time to locate the event in the full time history
            event_t_center = Global_Time(idx_best);
            real_v = Global_Speed(idx_best);
            
            % Window: +/- 1.5 seconds around the peak
            t_win_start = event_t_center - 1.5;
            t_win_end   = event_t_center + 1.5;
            
            mask_win = Time >= t_win_start & Time <= t_win_end;
            sig_chunk = Gyro_Z_HP(mask_win); % Use HP signal to remove DC
            
            % 4. Compute FFT (Amplitude Spectrum)
            L = length(sig_chunk);
            if L > 0
                nfft = 2^nextpow2(L);
                Y = fft(sig_chunk, nfft);
                f = (Fs/2) * linspace(0,1,nfft/2+1);
                P1 = 2*abs(Y(1:nfft/2+1)/L); % Single-sided spectrum
                
                % 5. Plot
                subplot(3,1,i);
                plot(f, P1, 'b', 'LineWidth', 1.5);
                grid on;
                
                % Highlight Weave Band (0.5 - 4.5 Hz)
                xline(0.5, 'g--'); xline(4.5, 'g--');
                
                % Formatting
                ylabel('Amplitude [deg/s]');
                title(sprintf('Target: %d km/h | Actual: %.1f km/h', t_spd, real_v));
                xlim([0 12]); % Show up to 12 Hz
                
                % Find peak within Weave Band for annotation
                mask_band = f >= 0.5 & f <= 4.5;
                [max_amp, idx_peak] = max(P1(mask_band));
                f_band = f(mask_band);
                freq_peak = f_band(idx_peak);
                
                text(freq_peak, max_amp, sprintf('  %.1f Hz', freq_peak), 'Color', 'r', 'FontWeight', 'bold');
            end
        end
        xlabel('Frequency [Hz]');
        sgtitle('Frequency Content Evolution at Different Speeds');
    end
end
%% 13. WOBBLE ANALYSIS (RMS + FREQ vs SPEED CORRELATION)
%  Description: 
%     Calculates Wobble intensity and correlates Frequency with Speed.
%     Restored 'Critical' threshold and added direct Scatter Plot.
fprintf('\n======================================================\n');
fprintf('   WOBBLE ANALYSIS (RMS & FREQ vs SPEED)    \n');
fprintf('======================================================\n');

% --- Configuration ---
F_WOBBLE_MIN = 6;       % [Hz] Lower bound
F_WOBBLE_MAX = 12;      % [Hz] Upper bound
FREQ_THRESH  = 8;       % [Hz] Threshold of interest
RMS_WINDOW_SEC = 0.5;   % [s] Window size
Fs = 1 / mean(diff(Time)); 

% 1. Band-Pass Filtering
[b_wob, a_wob] = butter(2, [F_WOBBLE_MIN F_WOBBLE_MAX]/(Fs/2), 'bandpass');
Steer_Wobble_Filt = filtfilt(b_wob, a_wob, Steer);

% 2. RMS Calculation
win_samples = round(RMS_WINDOW_SEC * Fs);
Steer_Wobble_Sq = Steer_Wobble_Filt .^ 2;
Mean_Sq = smoothdata(Steer_Wobble_Sq, 'movmean', win_samples);
Wobble_RMS = sqrt(Mean_Sq);

% 3. Frequency Extraction (Hilbert)
analytic_sig = hilbert(Steer_Wobble_Filt);
inst_phase = unwrap(angle(analytic_sig));
inst_freq = [0; diff(inst_phase)] / (2*pi) * Fs; 
inst_freq_smooth = smoothdata(inst_freq, 'movmean', Fs/5);

% 4. Visualization
figure('Name', 'Wobble Detailed Analysis', 'Color', 'w');

% --- A. Time History ---
subplot(3,1,1);
plot(Time, Steer_Wobble_Filt, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;
plot(Time, Wobble_RMS, 'r', 'LineWidth', 1.5);
grid on;
ylabel('Steer [deg]');
title(sprintf('A. Time History (Band-Pass %.0f-%.0f Hz)', F_WOBBLE_MIN, F_WOBBLE_MAX));
legend('Raw Signal', 'RMS Intensity', 'Location', 'best');
xlim([min(Time) max(Time)]);

% Filters for Scatter Plots (Straight line & Significant Wobble)
Steer_Mean = smoothdata(Steer, 'movmean', Fs*2); 
Is_Straight = abs(Steer_Mean) < 1.5; 
Is_Significant = Wobble_RMS > 0.4; % Show only if vibration > 0.4 deg
mask_clean = (V_Sat > 10) & Is_Straight & Is_Significant;

% --- B. Intensity vs Speed ---
subplot(3,1,2);
if sum(mask_clean) > 0
    scatter(V_Sat(mask_clean), Wobble_RMS(mask_clean), 20, Time(mask_clean), 'filled');
    colormap('jet');
    c = colorbar; c.Label.String = 'Time [s]';
    grid on;
    xlabel('Speed [km/h]'); ylabel('RMS Intensity [deg]');
    title('B. Wobble Intensity vs Speed');
    xlim([0 max(V_Sat)+5]); ylim([0 max(Wobble_RMS)*1.2]);
    
    % Reference Lines
    yline(1.0, 'k--', 'Noticeable'); 
    yline(2.5, 'r--', 'Critical', 'LineWidth', 1.5); % <--- REINSERITA
else
    text(0.5, 0.5, 'No significant wobble detected.', 'Units', 'normalized', 'Color', 'r');
end

% --- C. Frequency vs Speed (NEW CORRELATION PLOT) ---


fprintf('Wobble Analysis Update: Critical line restored. Freq vs Speed plot generated.\n');
%% 14. SECTOR-BASED SPECTRAL ANALYSIS (WEAVE IDENTIFICATION)
%  Description: 
%     Calculates the Auto-Spectrum (PSD) for each identified segment individually.
%     Uses RAW signal + Detrending to preserve low-frequency Weave content (1-4 Hz)
%     which might have been attenuated by the global High-Pass filter.

fprintf('\n======================================================\n');
fprintf('   SECTOR SPECTRAL ANALYSIS (WEAVE)    \n');
fprintf('======================================================\n');

% --- Configuration ---
% Weave is typically between 1.5 Hz and 4.5 Hz for heavy motorcycles
F_MAX_PLOT = 6;  % [Hz] Limit the plot to see low freq clearly
Fs = 1 / mean(diff(Time)); 

figure('Name', 'Weave Auto-Spectrum per Sector', 'Color', 'w');

num_plots = min(6, size(segments, 1)); % Handle cases with fewer than 6 segments

for k = 1:num_plots
    % 1. Extract Data for the specific sector
    idx = segments(k,1):segments(k,2);
    
    % DECISION: Use RAW Gyro_Z to avoid filter attenuation at 2Hz.
    % We apply 'detrend' to remove DC offset/Drift locally.
    sig_raw = Gyro_Z(idx); 
    sig_proc = detrend(sig_raw, 'linear'); 
    
    % Get average speed for context
    avg_speed = mean(V_Sat(idx));
    
    % 2. Compute FFT (Auto-Spectrum)
    L = length(sig_proc);
    if L > 0
        nfft = 2^nextpow2(L);
        Y = fft(sig_proc, nfft);
        
        % Compute Single-Sided Spectrum (P1)
        P2 = abs(Y/L);
        P1 = P2(1:nfft/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        f_axis = Fs * (0:(nfft/2)) / nfft;
        
        % 3. Plotting
        subplot(3, 2, k);
        plot(f_axis, P1, 'b', 'LineWidth', 1.5); hold on;
        grid on;
        
        % Highlight Weave Band (Typical 2-4 Hz)
        xregion_h = xline(2, 'g:', 'LineWidth', 1); 
        xline(4, 'g:', 'LineWidth', 1);
        
        % Find Peak in Weave Band for automatic labeling
        mask_weave = f_axis >= 1.5 & f_axis <= 4.5;
        [max_amp, idx_max] = max(P1(mask_weave));
        f_weave_band = f_axis(mask_weave);
        
        if ~isempty(f_weave_band)
            peak_freq = f_weave_band(idx_max);
            
            % If peak is significant enough (> 0.5 deg/s for Gyro)
            if max_amp > 0.5
                plot(peak_freq, max_amp, 'rv', 'MarkerFaceColor','r');
                text(peak_freq, max_amp*1.1, sprintf('%.1f Hz', peak_freq), ...
                     'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
            end
        end
        
        % Formatting
        title(sprintf('Seg %d (Avg V: %.0f km/h)', k, avg_speed));
        if k > 4, xlabel('Frequency [Hz]'); end
        if mod(k,2) ~= 0, ylabel('|Yaw Rate| [deg/s]'); end
        
        xlim([0 F_MAX_PLOT]); % Focus on Low Frequencies
        
        % Auto-scale Y but keep a minimum floor to avoid empty looking plots
        ylim([0 max(2, max(P1)*1.2)]); 
    end
end
sgtitle('Weave Identification by Sector (Auto-Spectrum)');


%% 16. SECTION-BASED CLUSTERING (POINTS ONLY, +/- 3 KM/H)
fprintf('\n======================================================\n');
fprintf('   SECTION CLUSTERING: 10 KM/H STEPS (+/- 3 TOLERANCE)    \n');
fprintf('======================================================\n');

% 1. CONTROLLO ESISTENZA DATI
if ~exist('Global_Speed', 'var') || isempty(Global_Speed)
    error('Errore: Esegui prima la Sezione 7 per generare i dati "Global_...".');
end

% 2. MAPPATURA EVENTI -> SEZIONI
% Associa ogni singolo evento (picco analizzato) alla sezione GPS di appartenenza
Global_SectionID = NaN(size(Global_Time));

for i = 1:length(Global_Time)
    t_ev = Global_Time(i);
    for k = 1:size(segments, 1)
        % Tempo inizio/fine della sezione k
        t_start_seg = Time(segments(k,1));
        t_end_seg   = Time(segments(k,2));
        
        if t_ev >= t_start_seg && t_ev <= t_end_seg
            Global_SectionID(i) = k;
            break; 
        end
    end
end

% 3. IMPOSTAZIONE CLUSTER
speed_centers = 0:10:180; % Target: 0, 10, 20 ... 180 km/h
tolerance = 3;            % Tolleranza: +/- 3 km/h
num_bins = length(speed_centers);
num_secs = size(segments, 1);

% Matrici per salvare le medie: [Velocità x Sezione]
Mean_Freq = NaN(num_bins, num_secs);
Mean_Damp = NaN(num_bins, num_secs);
Mean_Real = NaN(num_bins, num_secs); % Parte Reale (Decay Rate)

% 4. CALCOLO DELLE MEDIE
for s = 1:num_secs
    % Indici dei dati che appartengono alla sezione 's'
    idx_sec = find(Global_SectionID == s);
    
    if isempty(idx_sec), continue; end
    
    % Dati grezzi della sezione
    sec_v = Global_Speed(idx_sec);
    sec_f = Global_Freq(idx_sec);
    sec_d = Global_Damp(idx_sec); 
    
    % Calcolo sigma (parte reale) per ogni singolo evento
    % sigma = (zeta * wn) / sqrt(1-zeta^2) -> approssimato a wd per zeta piccoli
    sec_zeta = sec_d / 100;
    sec_wd   = sec_f * 2 * pi; % rad/s
    sec_sigma = zeros(size(sec_zeta));
    
    for j = 1:length(sec_zeta)
        z = sec_zeta(j);
        if z < 0.99
            sec_sigma(j) = (z * sec_wd(j)) / sqrt(1 - z^2);
        else
            sec_sigma(j) = NaN;
        end
    end
    
    % Loop sui bin di velocità (10, 20, 30...)
    for b = 1:num_bins
        v_target = speed_centers(b);
        
        % Maschera: prendi solo velocità nel range [target-3, target+3]
        mask_bin = sec_v >= (v_target - tolerance) & sec_v <= (v_target + tolerance);
        
        if sum(mask_bin) > 0
            Mean_Freq(b, s) = mean(sec_f(mask_bin));
            Mean_Damp(b, s) = mean(sec_d(mask_bin));
            Mean_Real(b, s) = mean(sec_sigma(mask_bin));
        end
    end
end

% 5. PLOTTING
% Definizione colori fissi per le 6 sezioni per coerenza tra i grafici
custom_colors = [
    0 0.4470 0.7410; % Blu
    0.8500 0.3250 0.0980; % Arancione
    0.9290 0.6940 0.1250; % Giallo ocra
    0.4940 0.1840 0.5560; % Viola
    0.4660 0.6740 0.1880; % Verde
    0.6350 0.0780 0.1840  % Rosso scuro
];
% Se ci sono più di 6 sezioni, estende la colormap
if num_secs > 6
    custom_colors = lines(num_secs); 
end

legend_list = cell(1, num_secs);
for k=1:num_secs, legend_list{k} = ['Section ' num2str(k)]; end

% --- GRAFICO 1: FREQUENZA vs VELOCITÀ ---
figure('Name', 'Mean Frequency vs Speed (Points)', 'Color', 'w');
hold on; grid on;
for s = 1:num_secs
    % 'LineStyle', 'none' toglie la linea di collegamento
    % 'Marker', 'o' mette il pallino
    % 'MarkerFaceColor' riempie il pallino col colore della sezione
    plot(speed_centers, Mean_Freq(:, s), 'o', ...
        'LineStyle', 'none', ...
        'Color', custom_colors(s,:), ...
        'MarkerFaceColor', custom_colors(s,:), ...
        'MarkerSize', 6); 
end
xlabel('Speed [km/h]'); ylabel('Frequency [Hz]');
title({'Frequency vs Speed', ['(Mean Values @ \pm' num2str(tolerance) ' km/h)']});
legend(legend_list, 'Location', 'best');
xlim([0 180]); ylim([0 6]);

% --- GRAFICO 2: SMORZAMENTO vs VELOCITÀ ---
figure('Name', 'Mean Damping vs Speed (Points)', 'Color', 'w');
hold on; grid on;
for s = 1:num_secs
    plot(speed_centers, Mean_Damp(:, s), 's', ... % Quadrato per variare
        'LineStyle', 'none', ...
        'Color', custom_colors(s,:), ...
        'MarkerFaceColor', custom_colors(s,:), ...
        'MarkerSize', 6);
end
yline(0, 'k--');
xlabel('Speed [km/h]'); ylabel('Damping Ratio [%]');
title({'Damping Ratio vs Speed', ['(Mean Values @ \pm' num2str(tolerance) ' km/h)']});
legend(legend_list, 'Location', 'best');
xlim([0 180]); ylim([0 40]);

% --- GRAFICO 3: ROOT LOCUS (CLUSTERIZZATO) ---
% Qui plottiamo la "Parte Reale Media" vs "Frequenza Media" per ogni step di velocità
figure('Name', 'Root Locus Clustered by Speed', 'Color', 'w');
hold on; grid on;
for s = 1:num_secs
    % Asse X: -Mean_Real (Sigma negativo è stabile)
    % Asse Y: Mean_Freq (in Hz, per semplicità di lettura, o rad/s)
    
    % Nota: Uso Frequenza in Hz sull'asse Y come richiesto spesso nei report sperimentali
    % Se vuoi rad/s, moltiplica Y per 2*pi
    x_val = -Mean_Real(:, s);
    y_val = Mean_Freq(:, s); % Hz
    
    plot(x_val, y_val, 'd', ... % Rombo
        'LineStyle', 'none', ...
        'Color', custom_colors(s,:), ...
        'MarkerFaceColor', custom_colors(s,:), ...
        'MarkerSize', 6);
end
xline(0, 'r--', 'Stability Limit');
xlabel('Real Part (Decay Rate) [rad/s]'); 
ylabel('Frequency [Hz]');
title('Root Locus (Mean Operating Points)');
legend(legend_list, 'Location', 'best');
axis normal; 
% Setto i limiti per vedere bene la zona stabile
xlim([-8 1]); ylim([0 6]);

disp('Grafici a punti (Cluster 10 km/h) generati.');