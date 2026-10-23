% co2model_batch.m — Headless batch OCIM2-48L CO2 perturbation model
%
% Runs 4 scenarios: baseline (no dredging), 1x, 10x, 100x dredging forcing.
% No figures, no interactive output — designed for SLURM / matlab -batch.
%
% Dredging-forcing window, protocol of Atwood et al. (2024, §2.2, p. 3).
%   DREDGE_YEAR_START, DREDGE_YEAR_END: emission years (inclusive) where the source is applied; zero elsewhere.
%   DREDGE_SCALE: factor applied to the forcing read (0.1: step7 writes the decadal total of C_ri, imposed every
%                 year; 0.1 gives the annual mass M_ijk of Eq. 10 of the manuscript).
%   With a window, the dredging scenarios restart from the control state at the first step of the window, as
%   co2model_dredge.m did ("starts from baseline state at 2012"): before the window, the source is zero
%   and the state is that of the control.
%   Without these variables: previous behaviour, identical (source at every step, 1780-2100, scale 1).
%
% Usage:
%   matlab -batch "run('co2model_batch.m')"
%   matlab -batch "forcing_file = '/path/to/jdredge.mat'; run('co2model_batch.m')"
%
% Dependencies (all in working directory or on path):
%   OCIM2_48L_CTL.mat, woa09si.mat, woa09po4.mat, netemission.txt
%   schmidt.m, inpaint_nans.m, sw_pres.m, eqco2.m, eqdic.m,
%   nsnew.m, nsgmres.m, ns_step_eb.m, d0.m, mfactor.m, CO2SYS.m

fprintf('=== CO2MODEL_BATCH: OCIM2-48L DREDGING PERTURBATION ===\n');
fprintf('Start: %s\n\n', datestr(now));
t_wall_start = tic;

% ── Section 0: Configuration ──────────────────────────────────────────────

% OCIM_OUT_DIR: replay without overwriting the published outputs. Default unchanged.
out_dir = getenv('OCIM_OUT_DIR');
if isempty(out_dir)
    out_dir = fullfile(getenv('HOME'), 'scratch', 'output_V6');
end
disp(['  out_dir = ' out_dir]);
if ~exist(out_dir, 'dir')
    out_dir = '.';
end
run_timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% ── Section 1: Load OCIM ─────────────────────────────────────────────────

fprintf('--- Section 1: Load OCIM2-48L ---\n');
load OCIM2_48L_CTL.mat output sol
grid = output.grid;
M3d = output.M3d;
Area3 = grid.DXT3d .* grid.DYT3d;
VOL = Area3 .* grid.DZT3d;

spyr = 365.25 * 24 * 3600;
TR = output.TR;

iocn = find(M3d(:) == 1);
isurf = find(M3d(:,:,1) == 1);
ideep = setdiff(iocn, isurf);
ns = length(isurf);
m = size(TR, 1);

fprintf('  Ocean points: m = %d, surface = %d\n', m, ns);

% ── Section 2: Gas exchange ──────────────────────────────────────────────

fprintf('--- Section 2: Gas exchange ---\n');
u10 = mean(output.data.u10cfc11, 2);
fice = mean(output.data.ficecfc11, 2);

temp = output.data.dens.tstar(1:ns);
salt = output.data.dens.sstar(1:ns);
Sc_co2 = schmidt(temp, 'co2');

a = sol(length(sol));

Kw = a * (1 - fice) .* (u10.^2) .* (Sc_co2 ./ 660).^(-.5);
clear output
fprintf('  Kw range: [%.4e, %.4e] m/yr\n', min(Kw), max(Kw));

% ── Section 3: Nutrients ─────────────────────────────────────────────────

fprintf('--- Section 3: Nutrients ---\n');
load woa09si.mat si_an LON LAT DEPTH
load woa09po4.mat po4_an
Si = interp3(LON, LAT, DEPTH, si_an, grid.XT3d, grid.YT3d, grid.ZT3d);
Po4 = interp3(LON, LAT, DEPTH, po4_an, grid.XT3d, grid.YT3d, grid.ZT3d);
for k = 1:size(M3d, 3)
    si(:,:,k) = inpaint_nans(Si(:,:,k));
    po4(:,:,k) = inpaint_nans(Po4(:,:,k));
end
fprintf('  Nutrients interpolated and gap-filled.\n');

% ── Section 4: CO2 system parameters ────────────────────────────────────

fprintf('--- Section 4: CO2 system ---\n');
co2syspar.pres = sw_pres(grid.ZT3d(isurf), grid.YT3d(isurf));
co2syspar.si = si(isurf);
co2syspar.po4 = po4(isurf);
co2syspar.temp = temp;
co2syspar.salt = salt;

rhobar = 1025;
sbar = VOL(isurf)' * salt(1:ns) ./ sum(VOL(isurf));

alk = 2310 * salt ./ sbar;
co2syspar.alk = alk(1:ns);

pco2atm = 280;
fprintf('  Pre-industrial pCO2: %d uatm\n', pco2atm);

% ── Section 5: Equilibrium DIC ───────────────────────────────────────────

fprintf('--- Section 5: Equilibrium DIC ---\n');
fprintf('Solving for equilibrium DIC distribution...\n');
tic
dic0 = zeros(m, 1) + 2100;
options.atol = 5e-6; options.rtol = 1e-16; options.iprint = 0;
[dic, ierr] = nsnew(dic0, @(x)eqdic(x, grid, M3d, TR, Kw, pco2atm, co2syspar), options);
t_eq = toc;
fprintf('  Equilibrium DIC solved in %.1f s (ierr=%d)\n', t_eq, ierr);

[~, ~, Jco2pre] = eqdic(dic, grid, M3d, TR, Kw, pco2atm, co2syspar);
[co2swpre, Rpre, k0pre, pHpre] = eqco2(dic(1:ns), co2syspar);
fprintf('  Equilibrium pH range: [%.3f, %.3f]\n', min(pHpre), max(pHpre));

% ── Section 6: Anthropogenic emissions ───────────────────────────────────

fprintf('--- Section 6: Emissions ---\n');
f = load('netemission.txt');
year_em = f(:,1);
em_gcb = f(:,2);  % umol C/yr

% Match original co2model.m: 1780-2100 starting from pre-industrial equilibrium.
% IC = equilibrium at 280 uatm, so we must spin up through historical emissions
% to get a physically consistent baseline at present day.
dt = 1;
year = 1780:2100;
Jem = zeros(length(year) + 1, 1);
% netemission.txt maps year -> emission; fill Jem(2:end) for each year
Jem(2:length(em_gcb) + 1) = em_gcb;
% Extend last emission value for years beyond emission file
Jem(length(em_gcb) + 2:length(year) + 1) = em_gcb(end);

t_sim = 0:dt:length(year);
nt = length(t_sim);
fprintf('  Simulation: %d-%d (%d timesteps, dt=%d yr)\n', year(1), year(end), nt, dt);
em_pos = Jem(Jem > 0);
if ~isempty(em_pos)
    fprintf('  Emissions range: [%.4e, %.4e] umol C/yr\n', min(em_pos), max(em_pos));
else
    fprintf('  Emissions: all zero (pre-industrial only)\n');
end

% ── Section 6b: Dredging-forcing window ─────────────────────────────────
% Step i applies Jem(i), the emissions of year year(1) + t_sim(i) - dt (Jem(2) = first line of
% netemission.txt). annee_pas(1) is the initial state, without a step.
annee_pas = year(1) + t_sim - dt;
dredge_y0 = str2double(getenv('DREDGE_YEAR_START'));
dredge_y1 = str2double(getenv('DREDGE_YEAR_END'));
dredge_scale = str2double(getenv('DREDGE_SCALE'));
window = ~isnan(dredge_y0) || ~isnan(dredge_y1);
if isnan(dredge_y0), dredge_y0 = -Inf; end
if isnan(dredge_y1), dredge_y1 = Inf; end
if isnan(dredge_scale), dredge_scale = 1; end
actif = (annee_pas >= dredge_y0) & (annee_pas <= dredge_y1);
actif(1) = false;
i0 = find(actif, 1);
if window
    assert(year_em(1) == year(1), 'netemission.txt starts in %d, not in %d', year_em(1), year(1));
    assert(~isempty(i0) && i0 > 2, 'Empty window, or open at the first step');
    fprintf('  Dredging window: %d-%d, steps %d to %d (%d active steps) ; scale %g\n', dredge_y0, dredge_y1, ...
            i0, find(actif, 1, 'last'), nnz(actif), dredge_scale);
end

% ── Section 7: Gas exchange operator + preconditioner ────────────────────

fprintf('--- Section 7: Preconditioner ---\n');

Q = 0 * M3d;
Q(isurf) = Kw ./ grid.dzt(1);
Qgas = d0(Q(iocn));

V = grid.DXT3d(iocn) .* grid.DYT3d(iocn) .* grid.DZT3d(iocn);

% Preconditioner: factored Jacobian at equilibrium (reused across scenarios)
Jtrawl_zero = zeros(m, 1);
xi = [dic; pco2atm];
[F, Jac] = ns_step_eb(xi, TR, M3d, V, Qgas, pco2atm, co2syspar, Jem(1), Jtrawl_zero, dt, dic);
fprintf('  Factoring preconditioner...\n');
tic
FM = mfactor(Jac);
t_pre = toc;
precond = @(x)mfactor(FM, x);
fprintf('  Preconditioner factored in %.1f s\n', t_pre);

% ── Section 8: Load Jdredge forcing ──────────────────────────────────────

fprintf('\n--- Section 8: Load Jdredge forcing ---\n');

if ~exist('forcing_file', 'var')
    files = dir(fullfile(out_dir, 'jdredge_ocim2_48l_*.mat'));
    if isempty(files)
        error('No jdredge_ocim2_48l_*.mat found in %s. Set forcing_file manually.', out_dir);
    end
    [~, newest] = max([files.datenum]);
    forcing_file = fullfile(out_dir, files(newest).name);
end

fprintf('  Loading forcing: %s\n', forcing_file);
jdata = load(forcing_file);

if isfield(jdata, 'Jdredge')
    Jdredge_1x = jdata.Jdredge(:);
elseif isfield(jdata, 'Jtrawl')
    Jdredge_1x = jdata.Jtrawl(:);
else
    error('Neither Jdredge nor Jtrawl found in forcing file.');
end

% Use pre-computed amplified scenarios if available, else scale
if isfield(jdata, 'Jdredge_10x')
    Jdredge_10x = jdata.Jdredge_10x(:);
else
    Jdredge_10x = Jdredge_1x * 10;
end
if isfield(jdata, 'Jdredge_100x')
    Jdredge_100x = jdata.Jdredge_100x(:);
else
    Jdredge_100x = Jdredge_1x * 100;
end
if dredge_scale ~= 1   % Section 6b
    Jdredge_1x = Jdredge_1x * dredge_scale;
    Jdredge_10x = Jdredge_10x * dredge_scale;
    Jdredge_100x = Jdredge_100x * dredge_scale;
end

assert(length(Jdredge_1x) == m, 'Jdredge length (%d) != m (%d)', length(Jdredge_1x), m);

% G1: Units contract
total_forcing_gC = sum(Jdredge_1x .* V .* rhobar .* 12e-6);
fprintf('  Forcing total: %.4f TgC/yr\n', total_forcing_gC / 1e12);
fprintf('  Jdredge unit: umol C kg-1 yr-1\n');
fprintf('  Active cells: %d / %d\n', nnz(Jdredge_1x), m);
assert(total_forcing_gC > 0, 'Total Jdredge forcing must be positive (source to ocean).');

% Load metadata if available
if isfield(jdata, 'conservation_error_pct')
    fprintf('  Conservation error (from step7): %.4f%%\n', jdata.conservation_error_pct);
end

% ── Section 9: Transient simulation — 4 scenarios ───────────────────────

fprintf('\n--- Section 9: Transient simulation ---\n');

scenario_names = {'baseline', 'dredge_1x', 'dredge_10x', 'dredge_100x'};
scenario_forcings = {Jtrawl_zero, Jdredge_1x, Jdredge_10x, Jdredge_100x};
n_scenarios = length(scenario_names);

% Storage for final-step summary
pco2a_final = zeros(1, n_scenarios);
flux_final = zeros(1, n_scenarios);

for s = 1:n_scenarios
    fprintf('\n  === Scenario %d/%d: %s ===\n', s, n_scenarios, scenario_names{s});
    t_scen_start = tic;

    Jforcing = scenario_forcings{s};

    % Initialize state arrays
    DIC = zeros(m, nt);
    pco2a = zeros(1, nt);
    pH_out = zeros(ns, nt);
    Revelle = zeros(ns, nt);
    co2_sw = zeros(ns, nt);
    Jco2 = zeros(ns, nt);
    Jgas = zeros(m, nt);

    % Initial conditions (equilibrium)
    DIC(:,1) = dic;
    pco2a(1) = pco2atm;
    Jco2(:,1) = Jco2pre(1:ns);
    pH_out(:,1) = pHpre;
    co2_sw(:,1) = co2swpre;
    Revelle(:,1) = Rpre;

    % Time-stepping
    options.atol = 2e-6; options.rtol = 1e-16; options.iprint = 0;

    % Section 6b: with a window, a dredging scenario takes the control state (scenario 1) up to the step
    % preceding the window, where the source is zero, and integrates only from i0.
    i_start = 2;
    if window && s > 1
        DIC(:,1:i0-1) = DIC_b(:,1:i0-1);
        pco2a(1:i0-1) = pco2a_b(1:i0-1);
        Jco2(:,1:i0-1) = Jco2_b(:,1:i0-1);
        pH_out(:,1:i0-1) = pH_b(:,1:i0-1);
        co2_sw(:,1:i0-1) = co2sw_b(:,1:i0-1);
        Revelle(:,1:i0-1) = Revelle_b(:,1:i0-1);
        Jgas(:,1:i0-1) = Jgas_b(:,1:i0-1);
        i_start = i0;
        fprintf('    Control state reused up to step %d ; integration of steps %d to %d\n', i0 - 1, i0, nt);
    end

    for i = i_start:nt
        if actif(i)
            Jf = Jforcing;
        else
            Jf = Jtrawl_zero;
        end
        xi = [DIC(:,i-1); pco2a(i-1)];
        [x2, ierr] = nsgmres(xi, ...
            @(x)ns_step_eb(x, TR, M3d, V, Qgas, pco2a(i-1), co2syspar, ...
                           Jem(i), Jf, dt, DIC(:,i-1)), ...
            precond, options);

        DIC(:,i) = x2(1:m);
        pco2a(i) = x2(m+1);

        [~, ~, Jgas(:,i)] = ns_step_eb(x2, TR, M3d, V, Qgas, ...
            pco2a(i-1), co2syspar, Jem(i), Jf, dt, DIC(:,i-1));

        [co2, R, ~, ph] = eqco2(DIC(1:ns,i), co2syspar);
        Revelle(:,i) = R;
        pH_out(:,i) = ph;
        co2_sw(:,i) = co2;
        Jco2(:,i) = Jgas(1:ns,i);

        if mod(i, 50) == 0 || i == nt
            flux_GtC = -V' * Jgas(:,i) * 1e-6 * rhobar * 12 / 1e15;
            fprintf('    step %3d/%d  year=%d  pCO2=%.2f uatm  flux=%.4f GtC/yr\n', ...
                    i, nt, year(1) + t_sim(i), pco2a(i), flux_GtC);
        end
    end

    t_scen = toc(t_scen_start);

    if window && s == 1   % Section 6b: control state, reused by the dredging scenarios
        DIC_b = DIC; pco2a_b = pco2a; Jco2_b = Jco2; pH_b = pH_out;
        co2sw_b = co2_sw; Revelle_b = Revelle; Jgas_b = Jgas;
    end

    % Store final values for comparison
    pco2a_final(s) = pco2a(end);
    flux_final(s) = -V' * Jgas(:,end) * 1e-6 * rhobar * 12 / 1e15;

    % Save per-scenario output
    out_file = fullfile(out_dir, sprintf('ocim_%s_%s.mat', scenario_names{s}, run_timestamp));
    scenario_name = scenario_names{s};
    Jdredge_forcing = Jforcing;

    save(out_file, 't_sim', 'year', 'DIC', 'pco2a', 'Jco2', 'co2_sw', ...
         'Revelle', 'pH_out', 'Jem', 'Jdredge_forcing', 'grid', 'M3d', ...
         'scenario_name', 'annee_pas', 'actif', 'dredge_y0', 'dredge_y1', 'dredge_scale', '-v7.3');

    fprintf('  Saved: %s (%.1f min elapsed)\n', out_file, t_scen / 60);
end

% ── Section 10: Summary ─────────────────────────────────────────────────

fprintf('\n=== SCENARIO COMPARISON ===\n');
fprintf('%-20s  %12s  %12s  %16s\n', 'Scenario', 'pCO2_final', 'delta_pCO2', 'flux_2100_GtC_yr');
delta_pco2 = pco2a_final - pco2a_final(1);
for s = 1:n_scenarios
    fprintf('%-20s  %12.2f  %+12.4f  %+16.4e\n', ...
        scenario_names{s}, pco2a_final(s), delta_pco2(s), flux_final(s));
end

% Linearity check
if delta_pco2(2) ~= 0
    ratio_10 = delta_pco2(3) / delta_pco2(2);
    ratio_100 = delta_pco2(4) / delta_pco2(2);
    fprintf('\nLinearity check:\n');
    fprintf('  delta_10x / delta_1x  = %.2f (expect ~10)\n', ratio_10);
    fprintf('  delta_100x / delta_1x = %.2f (expect ~100)\n', ratio_100);
else
    fprintf('\nWARNING: delta_pCO2 for 1x scenario is zero — cannot check linearity.\n');
end

t_total = toc(t_wall_start);
fprintf('\nTotal wall time: %.1f min\n', t_total / 60);
fprintf('=== CO2MODEL_BATCH COMPLETE: %s ===\n', datestr(now));
