% MAIN  Gravitational microlensing lab - analysis driver.
%
% Put OGLE light curves ("phot.dat" files: time, magnitude, error) in the
% data/ folder, list them in the table below and run this script.
% If none of the listed files exist, a synthetic event is analysed instead,
% so the script always runs.
%
% Columns of the table:
%   name     label used in the printout and figure file names
%   file     path of the photometry file
%   blended  true -> also fit the blending fraction fbl (4 parameters)
%   ref      published parameters to compare with (may be struct())

addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));

events = {
    "SMC001", "data/SMC001.dat", false, struct('T0', 2454476.072, 'tau', 51.032, 'uMin', 0.121)
    "BLG002", "data/BLG002.dat", true,  struct('T0', 2454517.540, 'tau', 16.779, 'uMin', 0.411, 'fbl', 0.512)
    "BLG004", "data/BLG004.dat", false, struct('T0', 2454521.168, 'tau', 18.108, 'uMin', 0.145)
    "BLG005", "data/BLG005.dat", false, struct('T0', 2454534.171, 'tau', 17.577, 'uMin', 0.505)
    "BLG007", "data/BLG007.dat", false, struct('T0', 2454534.013, 'tau', 28.010, 'uMin', 1.397)
    "BLG012", "data/BLG012.dat", false, struct('T0', 2454501.089, 'tau',  9.322, 'uMin', 0.576)
    "BLG015", "data/BLG015.dat", false, struct('T0', 2454510.512, 'tau',  6.520, 'uMin', 1.105)
    };

saveDir = "results";
results = struct();
found = false;

for i = 1:size(events, 1)
    [name, file, blended, ref] = events{i, :};
    if ~isfile(file)
        fprintf('Skipping %s: %s not found.\n', name, file);
        continue
    end
    found = true;
    lc = loadLightCurve(file);
    results.(name) = analyzeEvent(lc, 'Name', name, 'Blended', blended, ...
        'Reference', ref, 'SaveDir', saveDir);
end

if ~found
    fprintf('\nNo data files found - running the synthetic demo instead.\n');
    truth = struct('T0', 2454517.5, 'tau', 17, 'uMin', 0.4, 'fbl', 0.5);
    lc = makeSyntheticLightCurve(truth, 'NumPoints', 800);
    results.demo = analyzeEvent(lc, 'Name', "demo", 'Blended', true, ...
        'Reference', truth, 'SaveDir', saveDir);
end
