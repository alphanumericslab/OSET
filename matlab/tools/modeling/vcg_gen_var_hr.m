function [vcg, phi] = vcg_gen_var_hr(HR, fs, alpha, b, theta, theta0, varargin)
% 
% vcg_gen_var_hr - Synthetic cardiac dipole/vectorcardiogram (VCG) generator with variable heart rate
%
% Usage:
%   [vcg, phi] = vcg_gen_var_hr(HR, fs, alpha, b, theta, theta0, method)
%
% This function generates a synthetic cardiac dipole/VCG signal with a variable
%   heart rate. The dipole/VCG signal is modeled using the 'differential form' of
%   the sum of Gaussian dipole/VCG model.
%
% Inputs:
%   HR: Heart rate vector in beats per minute (BPM) over time
%   fs: Sampling rate
%   alpha: Structure containing the amplitudes of Gaussian functions for
%       modeling x, y, and z coordinates of the cardiac dipole/VCG
%   b: Structure containing the widths of Gaussian functions for modeling
%       x, y, and z coordinates of the cardiac dipole/VCG
%   theta: Structure containing the phase of Gaussian functions for modeling
%       x, y, and z coordinates of the cardiac dipole/VCG
%   theta0: Initial phase of the synthetic dipole/VCG
%   method: Implementation method: 'LOOP' (uses a while loop over the HR
%       vector), or 'MATRIX' (uses a matrix implementation). Default is the
%       'MATRIX' method
%
% Outputs:
%   vcg: Structure containing the x, y, and z coordinates of the cardiac dipole/VCG
%   phi: Vector containing the vcg phase
%
% References:
%   - Clifford, G. D., Nemati, S., & Sameni, R. (2010). An artificial vector
%     model for generating abnormal electrocardiographic rhythms. Physiological
%     Measurement, 31(5), 595–609.
%   - Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%     Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%     ECG Signals. EURASIP Journal on Advances in Signal Processing, 2007(1).
%   - McSharry, P. E., Clifford, G. D., Tarassenko, L., & Smith, L. A.
%     (2003). A dynamical model for generating synthetic electrocardiogram signals.
%     IEEE Transactions on Biomedical Engineering, 50(3), 289–294.
%
% Revision History:
%   2012: First release
%   2023: Renamed from the deprecated version DipoleGenerator3

if nargin > 6 && ~isempty(varargin{1})
    method = varargin{1};
else
    method = 'MATRIX';
end

switch method
    case 'LOOP'
        % Initialize parameters
        dt = 1 / fs;
        X = [];
        Y = [];
        Z = [];
        phi = theta0;
        i = 1;

        % Generate dipole/VCG signal for variable heart rate
        while true
            w = 2 * pi * HR(i) / 60; % the beat-wise angular velocity

            % Compute phase differences
            dtetaix = mod(phi(end) - theta.x + pi, 2 * pi) - pi;
            dtetaiy = mod(phi(end) - theta.y + pi, 2 * pi) - pi;
            dtetaiz = mod(phi(end) - theta.z + pi, 2 * pi) - pi;

            % Compute initial values for the first iteration
            if i == 1
                X = sum(alpha.x .* exp(-dtetaix .^ 2 ./ (2 * b.x .^ 2)));
                Y = sum(alpha.y .* exp(-dtetaiy .^ 2 ./ (2 * b.y .^ 2)));
                Z = sum(alpha.z .* exp(-dtetaiz .^ 2 ./ (2 * b.z .^ 2)));
            end

            % Update x, y, and z state variables
            X = cat(2, X, X(end) - dt * sum(w * alpha.x ./ (b.x .^ 2) .* dtetaix .* exp(-dtetaix .^ 2 ./ (2 * b.x .^ 2))));
            Y = cat(2, Y, Y(end) - dt * sum(w * alpha.y ./ (b.y .^ 2) .* dtetaiy .* exp(-dtetaiy .^ 2 ./ (2 * b.y .^ 2))));
            Z = cat(2, Z, Z(end) - dt * sum(w * alpha.z ./ (b.z .^ 2) .* dtetaiz .* exp(-dtetaiz .^ 2 ./ (2 * b.z .^ 2))));

            % Update dipole/VCG phase
            phi = cat(2, phi, mod(phi(end) + w * dt + pi, 2 * pi) - pi);

            % Check for a phase transition
            if abs(phi(end) - phi(end - 1)) > pi
                i = i + 1;
            end

            % Check for completion of the HR vector
            if i > length(HR)
                break;
            end
        end

        % Store the dipole/VCG coordinates in the output structure
        vcg.x = X;
        vcg.y = Y;
        vcg.z = Z;

    case 'MATRIX'
        dt = 1/fs; % time step
        w = 2 * pi * HR / 60; % the beat-wise angular velocity
        samples_per_beat = round(fs./(HR/60)); % number of samples per beat
        samples_per_beat(1) =  round((1 - (theta0+pi)/(2*pi))*samples_per_beat(1));
        cumlen = cumsum(samples_per_beat);

        N = sum(samples_per_beat);
        ww = zeros(1,N);
        ww(1:cumlen(1)) = repmat(w(1),1,samples_per_beat(1));
        for i = 2:length(w)
            %     ww(cumlen(i-1)+1:cumlen(i)) = repmat(w(i),1,samples_per_beat(i));
            ww(cumlen(i-1)+1:cumlen(i)) = w(i);
        end
        ww(cumlen(length(w)):end) = w(end);

        phi = cumsum([theta0+pi,ww*dt]);
        phi = mod(phi,2*pi)-pi;
        phi = phi(2:end);

        dtetaix = mod(phi(ones(length(theta.x),1),:)' - theta.x(ones(1,N),:) + pi , 2*pi) - pi;
        dtetaiy = mod(phi(ones(length(theta.y),1),:)' - theta.y(ones(1,N),:) + pi , 2*pi) - pi;
        dtetaiz = mod(phi(ones(length(theta.z),1),:)' - theta.z(ones(1,N),:) + pi , 2*pi) - pi;

        X = sum(alpha.x(ones(1,N),:) .* exp(-dtetaix .^2 ./ (2*b.x(ones(1,N),:) .^ 2)),2);
        Y = sum(alpha.y(ones(1,N),:) .* exp(-dtetaiy .^2 ./ (2*b.y(ones(1,N),:) .^ 2)),2);
        Z = sum(alpha.z(ones(1,N),:) .* exp(-dtetaiz .^2 ./ (2*b.z(ones(1,N),:) .^ 2)),2);

        vcg.x = X';
        vcg.y = Y';
        vcg.z = Z';

    otherwise
        error('Undefined implementation method');
end
