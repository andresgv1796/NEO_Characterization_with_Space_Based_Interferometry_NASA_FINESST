function [IC, tF] = load_JPL_orbit_data(filename)
%LOAD_JPL_ORBIT_DATA Load vertical orbit initial conditions and periods.
%   [IC, tF] = LOAD_JPL_ORBIT_DATA(filename) reads the CSV file 
%   specified by filename, extracts:
%     • IC – the initial conditions (columns 2–7) 
%     • tF – the full orbit period (column 9)
%
%   Inputs:
%     filename  – string, name of CSV file (with or without path)
%
%   Outputs (you choose names when calling):
%     IC        – N×6 matrix of [x0 y0 z0 vx0 vy0 vz0]
%     tF        – N×1 vector of orbit periods
%
%   Example:
%     [IC_JPL_L3_Vertical, t_F_JPL_L3_Vertical] = ...
%         load_JPL_orbit_data('L3_vertical_orbits_JPL_IC.csv');

    % Input validation
    if nargin < 1 || (~ischar(filename) && ~isstring(filename))
        error('You must supply a valid filename as a string.');
    end

    % Verify file exists
    if ~isfile(filename)
        error('File "%s" not found in the current folder or MATLAB path.', filename);
    end

    % Read numeric data from CSV
    data = readmatrix(filename);

    % Extract outputs
    IC = data(:, 2:7);
    tF = data(:,    9);
end
