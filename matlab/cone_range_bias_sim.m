function results = cone_range_bias_sim(options)
%CONE_RANGE_BIAS_SIM Validate cone-range nonlinear least-squares bias.
%
% RESULTS = CONE_RANGE_BIAS_SIM(OPTIONS) runs a Monte Carlo experiment for a
% simple cone-range positioning problem (e.g., GNSS-A ship track over a
% seafloor transponder). The script estimates the empirical bias of the
% nonlinear least-squares (NLS) solution and compares it with a
% second-order analytical approximation based on Hessian terms.
%
% OPTIONS fields (all optional):
%   .sigma       Scalar range noise standard deviation in meters (default 1.0)
%   .nTrials     Number of Monte Carlo trials           (default 2000)
%   .symmetric   Logical flag: true for 360° coverage, false for a 160°
%                unilateral arc that exaggerates vertical bias (default false)
%   .depth       True transponder depth in meters (positive number, default 2000)
%   .radius      Horizontal radius of ship track in meters (default 1500)
%   .seed        RNG seed for repeatability (default 1)
%
% OUTPUT struct fields:
%   biasMC       Empirical bias [dx; dy; dz] from Monte Carlo (meters)
%   biasTheory   Second-order bias approximation (meters)
%   covMC        Empirical covariance of the estimator (m^2)
%   geometry     Positions of ship samples (m)
%   trueState    Ground-truth [x; y; z] (m)
%   sigma        Noise sigma used in simulation (m)
%
% Example:
%   results = cone_range_bias_sim(struct('sigma',0.5,'nTrials',5000,'symmetric',false));
%   disp(results.biasMC), disp(results.biasTheory);
%
% The Gauss-Newton solver is implemented without external toolboxes so that the
% script can run in base MATLAB or GNU Octave.

arguments
    options.sigma (1,1) double {mustBePositive} = 1.0
    options.nTrials (1,1) double {mustBeInteger, mustBePositive} = 2000
    options.symmetric (1,1) logical = false
    options.depth (1,1) double {mustBePositive} = 2000
    options.radius (1,1) double {mustBePositive} = 1500
    options.seed (1,1) double = 1
end

rng(options.seed);
trueState = [0; 0; -abs(options.depth)];

shipXYZ = build_ship_track(options.radius, options.depth, options.symmetric);
trueRanges = predict_ranges(trueState, shipXYZ);

% Analytical bias approximation (second-order Taylor around truth)
[J, Hessians] = range_jacobians(trueState, shipXYZ);
F = J.' * J;
Sigma_delta = options.sigma^2 * inv(F); % Covariance of first-order solution
traceTerms = arrayfun(@(k) trace(Hessians(:,:,k) * Sigma_delta), 1:size(J,1));
theoryBias = -0.5 * (F \ (J.' * traceTerms.'));

% Monte Carlo estimation
estimates = zeros(3, options.nTrials);
for k = 1:options.nTrials
    noisyRanges = trueRanges + options.sigma * randn(size(trueRanges));
    estimates(:, k) = solve_gn(noisyRanges, shipXYZ, trueState);
end

errors = estimates - trueState;
results.biasMC     = mean(errors, 2);
results.covMC      = cov(errors.');
results.biasTheory = theoryBias;
results.geometry   = shipXYZ;
results.trueState  = trueState;
results.sigma      = options.sigma;

if nargout == 0
    fprintf('--- Cone-range NLS bias validation ---\n');
    fprintf('Trials: %d | Noise sigma: %.3f m | Geometry: %s\n', ...
        options.nTrials, options.sigma, ternary(options.symmetric,'symmetric','unilateral'));
    fprintf('True state     : [%.3f, %.3f, %.3f] m\n', trueState);
    fprintf('Empirical bias : [%.4f, %.4f, %.4f] m\n', results.biasMC);
    fprintf('Theory bias    : [%.4f, %.4f, %.4f] m\n', results.biasTheory);
    fprintf('RMS (x,y,z)    : [%.4f, %.4f, %.4f] m\n', sqrt(diag(results.covMC)));
end
end

% -------------------------------------------------------------------------
function shipXYZ = build_ship_track(radius, depth, symmetric)
%BUILD_SHIP_TRACK Generate surface sampling positions.
%   When symmetric=false, use an arc (unilateral geometry) to emphasize bias.
z = zeros(1, 12);
if symmetric
    thetas = linspace(0, 2*pi, 12+1);
    thetas(end) = [];
else
    thetas = linspace(-20, 180, 12) * pi/180;
end
shipXYZ = [radius * cos(thetas(:)), radius * sin(thetas(:)), z(:)];
end

function ranges = predict_ranges(state, shipXYZ)
diffs = shipXYZ - state.';
ranges = sqrt(sum(diffs.^2, 2));
end

function [J, Hessians] = range_jacobians(state, shipXYZ)
%RANGE_JACOBIANS Compute Jacobian and Hessians of range observations.
nObs = size(shipXYZ,1);
J = zeros(nObs, 3);
Hessians = zeros(3,3,nObs);
for i = 1:nObs
    d = state.' - shipXYZ(i,:);
    r = norm(d);
    J(i,:) = d / r;
    I = eye(3);
    Hessians(:,:,i) = (I / r) - (d.' * d) / r^3;
end
end

function est = solve_gn(ranges, shipXYZ, initial)
%SOLVE_GN Plain Gauss-Newton solver for the range-based positioning problem.
maxIter = 30;
tol = 1e-8;
x = initial;
for iter = 1:maxIter
    pred = predict_ranges(x, shipXYZ);
    residual = ranges - pred;
    [J, ~] = range_jacobians(x, shipXYZ);
    step = (J.' * J) \ (J.' * residual);
    x = x + step;
    if norm(step) < tol
        break;
    end
end
est = x;
end

function out = ternary(condition, a, b)
if condition, out = a; else, out = b; end
end
