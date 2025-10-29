% RESFOURIONET GROUND TRUTH DATA GENERATOR (Monte Carlo Simulation)
%
% This script runs a Monte Carlo simulation of light propagation in tissue
% to generate absorbed energy distributions (source functions) for
% subsequent thermal modeling. This generated data serves as the ground
% truth for testing the ResFouriONet architecture.
%
% The script initializes simulation parameters, iterates over focal depths,
% and saves the resulting absorbed energy matrices.
%
%
% Author: Aditya Roy
% Date: January 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;
tic

% --- 1. SIMULATION INPUT PARAMETERS ---
% These parameters are based on the experimental setup used in the paper.
wavelength = [920]; % Wavelength in nm (Fixed for this implementation)
fov = [0.1];        % Field of View (Linear scan range for photon launch) in mm
w0 = [10.86];       % 1/e^2 beam radius at objective back aperture in mm


% --- 2. MAIN EXECUTION LOOP (Iterates over all parameter combinations) ---
fprintf('Starting Monte Carlo Ground Truth Generation...\n');

for i = 1:length(fov)
    for j = 1:length(w0)            
        fprintf('\n--- Running Simulation for: Wavelength=%d nm, FOV=%.2f mm, w0=%.2f mm, T_surf=%d C ---\n', ...
                wavelength(1), fov(i), w0(j));
            
        % Call the main simulation function
        RunMonteCarloDepthScan(fov(i), w0(j), wavelength);
          
    end
end

toc
fprintf('\nSimulation finished. Data saved in the generated directory.\n');


% =========================================================================
% --- 3. CORE SIMULATION FUNCTIONS (Local Functions) ---
% =========================================================================


function [absPortion, lostPortions] = RunMonteCarloDepthScan(FOV, w0, wavelength)
% RunMonteCarloDepthScan: Initializes geometry/thermal parameters and
% iterates the Monte Carlo simulation across a range of tissue depths.
%
% Inputs:
% FOV (mm): Field of View for point scan (used for spatial averaging)
% w0 (mm): 1/e^2 beam radius at the objective back aperture
% wavelength (nm): Light source wavelength
% T_surf (C): Surface temperature boundary condition

global params;

% --- A. PARAMETER INITIALIZATION (Based on pointScan920.m) ---

% Optical properties of media (Model 3: Custom/Fixed parameters for 920nm)
params.opt.wavelength = wavelength; % nm
params.opt.absTissue = 0.039;       % 1/mm. Absorption coefficient (mu_a)
params.opt.scatterTissue = 6.7;     % 1/mm. Scattering coefficient (mu_s)
params.opt.gTissue = 0.9;           % Anisotropy factor
params.opt.nTissue = 1.36;          % Refractive index of brain
params.opt.nWater = 1.328;          % Refractive index of water


% Geometry and Focusing Parameters
params.geo.NA = 1.05;           % Numerical Aperture
params.geo.f = 7.2;             % Focal length of the objective: mm
params.geo.w0 = w0;             % 1/e^2 beam radius at objective back aperture
params.geo.FOV = FOV;           % Linear FOV of focal scan: mm
params.geo.d_glass = 0.17;      % Thinkness of cover glass: mm
params.geo.r_glass = 2;         % Radius of cover glass: mm

% Parameters updated during MC run (initialized to zero/empty)
params.geo.focalDepthTissue = 0;
params.geo.dstep = []; 
params.geo.zrange = []; 
params.geo.rmax = [];


% --- B. DEPTH SCAN SETUP ---

utot = params.opt.scatterTissue + params.opt.absTissue; % Total attenuation coefficient (mu_t)
% Depth array: [0, 1, 2, 3, 4, 5, 6] scattering lengths
depthArray = [0:6] ./ utot; 
depthArray(1) = 0.01; % Replace zero depth with a small non-zero value (0.01 mm)

% Initialize storage for simulation results
absPortion = zeros(size(depthArray));
lostPortions = zeros(4, length(depthArray));

% Create output directory
outputDir = ['./920_mcData_fov',num2str(FOV),'_w0',num2str(w0)];
if ~exist(outputDir, 'dir')
    mkdir (outputDir);
end
fprintf('Saving results to: %s\n', outputDir);


% --- C. DEPTH SCAN LOOP ---

for i = 1:length(depthArray)    
    params.geo.focalDepthTissue = depthArray(i); % Set the current focal depth
    
    % Run the Monte Carlo simulation
    [frac_abs, frac_trans, r1, d1, catcher, nlaunched, lostphotons] = MonteCarlo(500, [0 0 2], 1, 0.02); 
    
    % Store and display results
    focal_depth_um = utot * depthArray(i);
    absorbed_fraction = sum(sum(catcher)) / nlaunched;
    disp(['  Attenuation depth ',  num2str(depthArray(i), '%.2f'), ' l_t: Absorbed frac=', num2str(absorbed_fraction, '%.6f')]);
    
    absPortion(i) = absorbed_fraction;
    lostPortions(:,i) = lostphotons' / nlaunched;
    
    % Save the raw output data
    save([outputDir, '/light_920nm_', num2str(i-1),'attleng.mat'], 'frac_abs', 'frac_trans', 'r1', 'd1', 'catcher');
end

end % End of RunMonteCarloDepthScan function


% =========================================================================
% --- 4. MONTE CARLO ENGINE (from original MonteCarloLight_V3.m) ---
% =========================================================================

function [frac_abs, frac_trans, r, depth, catcher, nlaunched, lostphotons] = MonteCarlo(nphotonpackets, zrange, rmax, dr)
% MonteCarloLight_V3: Main Monte Carlo simulation engine.
% Simulates photon packet propagation, absorption, and scattering in a
% semi-infinite medium using the Henyey-Greenstein phase function.
%
% Inputs (Optional):
% nphotonpackets: Multiplier for total launched photons (100000 * n)
% zrange (mm): [z_min, z_center, z_max] for simulation volume
% rmax (mm): Maximum radial extent of simulation volume
% dstep (mm): Discretization step (dr=dz)

global params;

tic % Initialize timing

% --- A. INPUT & PARAMETER SETUP ---

dz = dr; 

% Parameters from the global structure
fll = params.geo.focalDepthTissue; % Focal depth inside the tissue
zmin = zrange(1); zcenter = zrange(2); zmax = zrange(3); 

% Tissue optical properties
absorption = params.opt.absTissue;
scattering = params.opt.scatterTissue;
g = params.opt.gTissue;


% Initialize MC variables
nphotonstotal_touse = 100000 * nphotonpackets; 
nphotons = 100000; % Number of photons "in play" at any particular time
nlaunched = nphotons;
lostphotons = [0 0 0 0]; % [Radial, Too Deep, Back to Glass, Back to Skull]

catcher = zeros(length(0:dr:rmax), length(zmin:dz:zmax)); % Energy deposition matrix

% --- B. MONTE CARLO LOOP INITIALIZATION ---

fprintf('  Monte Carlo running with N_total = %d packets... \n', nphotonstotal_touse);

% Launch first batch of photons using Gaussian beam profile
[cors, dirs] = photoninitGauss(nphotons, params.geo.w0, zcenter, fll);

w = ones(1, nphotons); % Initial weight matrix
newplotthreshold = nphotons; % Progress update threshold
stringout = ['Progress: 0.00%'];
fprintf(stringout);

% --- C. THE MAIN MONTE CARLO LOOP ---

while 1
    % 1. Choose step size (s)
    s = -log(rand(1, nphotons)) / (absorption + scattering);
    
    % 2. Move photon packets
    cors = cors + repmat(s, 3, 1) .* dirs;
    
    % 3. Check for out-of-bounds
    rcors = sqrt(sum(cors(1:2,:).^2));
    
    % Count photons lost to different boundaries (before setting weights to 0)
    lostphotons(1) = lostphotons(1) + sum(w(rcors >= rmax + dr));      
    lostphotons(2) = lostphotons(2) + sum(w(cors(3,:) >= zmax + dz));  
    lostphotons(3) = lostphotons(3) + sum(w(cors(3,:) < zmin & rcors < params.geo.r_glass)); 
    lostphotons(4) = lostphotons(4) + sum(w(cors(3,:) < zmin & rcors >= params.geo.r_glass)); 
    
    % Kill out-of-bounds packets
    outofbounds = rcors >= rmax + dr | cors(3,:) >= zmax + dz | cors(3,:) < zmin;    
    w(outofbounds) = 0;
    
    % Return killed packets to center (for indexing safety, weight remains 0)
    cors(:, w == 0) = [zeros(2, sum(w == 0)); ones(1, sum(w == 0)) * zcenter];
    rcors = sqrt(sum(cors(1:2,:).^2)); 
    
    % 4. Deposit absorbed energy into 'catcher' matrix
    zs = floor((cors(3,:) - zmin) / dz) + 1; % Z-index
    rs = floor(rcors / dr) + 1; % R-index
    
    % Filter for valid indices (inside the catcher bounds)
    valid_indices = rs >= 1 & rs <= size(catcher, 1) & zs >= 1 & zs <= size(catcher, 2);
    
    ins = sub2ind(size(catcher), rs(valid_indices), zs(valid_indices));
    wacumm = accumarray(ins', w(valid_indices));
    ns = (wacumm ~= 0);
    catcher(ns) = catcher(ns) + wacumm(ns) * (absorption / (absorption + scattering));
    
    % 5. Reduce photon packet weight by absorbed energy
    w = w * (scattering / (absorption + scattering));
    
    % 6. Kill packets near boundaries
    outofbounds2 = rcors >= rmax | cors(3,:) >= zmax;
    w(outofbounds2) = 0;
    
    % 7. Change direction (Scattering via Henyey-Greenstein)
    if g > 0 && g <= 1
        costh = ((1 + g.^2 - ((1 - g.^2) ./ (1 - g + 2 * g .* rand(1, nphotons))).^2) ./ (2 .* g));
    elseif g == 0
        costh = 2 * rand(1, nphotons) - 1;
    else
        error('g must be between 0 and 1');
    end
    phi = 2 * pi * rand(1, nphotons);
    sinth = sqrt(1 - costh.^2);
    
    % Direction update logic (using vector rotation)
    temp = sqrt(1 - dirs(3,:).^2);
    tofix = abs(dirs(1,:)) < 0.0001 & abs(dirs(2,:)) < 0.0001; % Vertically directed photons
    
    % Standard rotation (Henyey-Greenstein)
    uxyz = repmat(sinth, 3, 1) .* ...
           [(dirs(1:2,:).*repmat(dirs(3,:).*cos(phi),2,1)+ ...
             [-dirs(2,:);dirs(1,:)].*repmat(sin(phi),2,1)) ./ repmat(temp,2,1); ...
            -cos(phi).*temp] + ...
            dirs .* repmat(costh, 3, 1);
            
    % Special case for near-vertical photons
    if any(tofix)
        uxyz(:, tofix) = [repmat(sinth(tofix), 2, 1) .* [cos(phi(tofix)); sin(phi(tofix))]; ...
                         costh(tofix) .* sign(dirs(3, tofix))];
    end
    dirs = uxyz;
    
    % Re-normalize direction vector
    mag = sqrt(sum(dirs.^2));
    dirs = dirs ./ repmat(mag, 3, 1);
    
    % 8. Roulette
    chance = rand(1, nphotons);
    
    % Preserve energy: 10% of low-weight packets are kept, gaining 10x weight
    w(w < 1e-4 & chance <= 0.1 & w > 0) = w(w < 1e-4 & chance <= 0.1 & w > 0) / 0.1; 
    
    % Determine packets to be destroyed
    todestroy = (w < 1e-4 & chance > 0.1 & w > 0) | outofbounds | outofbounds2;
    ntodestroy = sum(todestroy);
    
    % 9. Replace destroyed photon packets
    if ntodestroy > 0
        if ntodestroy + nlaunched <= nphotonstotal_touse 
            [cors(:, todestroy), dirs(:, todestroy)] = photoninitGauss(ntodestroy, params.geo.w0, zcenter, fll);
            w(todestroy) = 1;
            nlaunched = nlaunched + ntodestroy;
        elseif nlaunched < nphotonstotal_touse 
            % Launch only the remaining needed photons
            which = find(todestroy);
            replaceins = (which(1:nphotonstotal_touse - nlaunched));
            [cors(:, replaceins), dirs(:, replaceins)] = photoninitGauss(length(replaceins), params.geo.w0, zcenter, fll);
            w(replaceins) = 1;
            nlaunched = nlaunched + length(replaceins);
            w(which(nphotonstotal_touse - nlaunched + 1:end)) = 0; 
        else
            w(w < 1e-4 & chance >= 0.1 & w > 0) = 0; % Kill remaining
        end
    end
    
    % 10. Termination condition
    if all(w == 0)
        break
    end
    
    % 11. Progress output
    if nlaunched >= max(newplotthreshold, nphotons + 1)
        tout = toc;
        timeremaining = tout * ((nphotonstotal_touse / (nlaunched - nphotons)) - 1); 
        
        newplotthreshold = newplotthreshold + nphotons;
        removeout = repmat('\b', 1, length(stringout));
        stringout = ['  Progress: ', num2str(nlaunched / nphotonstotal_touse * 100, '%.2f'), ...
            '%% of photons launched. Est. Time Remaining: ', num2str(timeremaining, '%.1f'), ' s.'];
        fprintf([removeout stringout])
    end
end

% --- D. CLEAN UP & OUTPUT GENERATION ---
totaltime = toc;
fprintf(['\n  Simulation Time: ', num2str(totaltime, '%.2f'), ' s\n'])

% Calculate absorbed energy density (fraction of launched energy per volume)
r = 0:dr:rmax;
depth = zmin:dz:zmax;

% Correct for cylindrical geometry (normalization by radial area)
frac_abs = (catcher ./ repmat(2 * pi * ((1:size(catcher, 1)) - 0.5)', 1, size(catcher, 2)))';
frac_abs = frac_abs ./ (dr.^2 * dz) ./ (nphotonstotal_touse); 

% The fraction of total power converted to temperature rise
frac_trans = frac_abs / absorption;

% Save final discretization parameters back to global for reference
params.geo.dstep = dr;
params.geo.zrange = zrange;
params.geo.rmax = rmax;


end % End of MonteCarlo function


% =========================================================================
% --- 5. PHOTON INITIALIZATION ---
% =========================================================================

function [cors, dirs] = photoninitGauss(num, r, zcenter, fl)
% photoninitGauss: Initializes the position and direction vectors for
% 'num' photons packets, simulating a focused Gaussian beam
% that is clipped by the objective's back aperture.
%
% Inputs:
% num: Number of photons to initialize
% r (mm): 1/e^2 radius of the Gaussian beam at the objective back aperture (params.geo.w0)
% zcenter (mm): Z-position where the light ray simulation starts (e.g., z=0)
% fl (mm): Focal depth inside the tissue (params.geo.focalDepthTissue)

global params;

% --- A. RADIAL POSITION/ANGLE CALCULATION ---

% Radial location for 2D Gaussian with 1/e^2 radius of r. 
% Generate redundant photons as some will be clipped.
rinit_redundant = sqrt(-0.5 * log(rand(1, 3 * num))) * r; 

% Clipping by the back aperture (NA restriction)
back_aperture_radius = params.geo.f * params.geo.NA;
rinit = rinit_redundant(rinit_redundant <= back_aperture_radius); 
rinit = rinit(1:num); 

% Sine condition for infinitely conjugated system (to calculate angle with the z-axis)
thinit = asind(rinit / params.geo.f / params.opt.nWater); 
thinitpos = -180 + 2 * rand(1, num) * 180; % Initial angular position (phi)

% --- B. DIRECTION VECTORS (DIRS) ---
% Directions are encoded in spherical coordinates [xdir, ydir, zdir] (unit vector)
zdir = cosd(thinit);
temp = sqrt(1 - zdir.^2);
xdir = -temp .* cosd(thinitpos); 
ydir = -temp .* sind(thinitpos);

dirs = [xdir; ydir; zdir];

% --- C. POSITION VECTORS (CORS) ---

% Simulating point scanning by giving each photon a random translation in xy.
FOV = params.geo.FOV;
dx = (FOV * rand(1, num) - 0.5 * FOV); % Uniform random translation in [-FOV/2, FOV/2]
dy = (FOV * rand(1, num) - 0.5 * FOV);

% Initializing x, y, and z coordinates at zcenter (based on focal depth 'fl')
% Beam size at z=zcenter (focusing geometry)
xcor = fl * tand(thinit) .* cosd(thinitpos);
ycor = fl * tand(thinit) .* sind(thinitpos);
zcor = ones(1, num) * zcenter;

% Apply the scan translation
cors = [xcor + dx; ycor + dy; zcor];

end % End of photoninitGauss function
