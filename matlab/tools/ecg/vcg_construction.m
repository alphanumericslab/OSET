function vcg = vcg_construction(ecg, lead_names, varargin)
% VCG_CONSTRUCTION - Constructs the Vectorcardiogram (VCG) from ECG leads.
%
% This function constructs the Vectorcardiogram (VCG) from a set of ECG
% leads. It uses a weighted sum of lead vectors corresponding to the
% available ECG leads to create the VCG. The lead vectors can be provided
% as inputs, or the function will use default approximate lead vector
% directions
%
% Usage:
%   vcg = vcg_construction(ecg, lead_names)
%   vcg = vcg_construction(ecg, lead_names, lead_vectors)
%   vcg = vcg_construction(ecg, lead_names, lead_vectors, plot_flag)
%
% Inputs:
%   - ecg: The ECG signal matrix, where each row represents an ECG lead and
%     each column represents a time sample.
%   - lead_names: A cell array containing the names of the ECG leads. These
%     names are used to identify which lead corresponds to which vector in
%     the VCG construction. The first dimension should match ecg channels
%   - lead_vectors (optional): A structure containing pre-defined lead
%     vectors for the ECG leads. If not provided, default vectors will be
%     used based on common lead positions.
%   - plot_flag (optional): A logical value (true or false) indicating
%     whether to plot lead vectors and the reconstructed VCG. Default is
%     false.
%
% Output:
%   - vcg: The reconstructed Vectorcardiogram (VCG), 3 x number of samples
%
% Supported ECG Lead Names (case-sensitive for lead_vectors, and case-insensitive with potential _ in the names for lead_names):
%   - I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6, Vx, Vy, Vz
%
% Revision History:
%   2023: First release.
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for additional inputs
if nargin > 2 && ~isempty(varargin{1})
    lead_vectors = varargin{1};
else
    lead_vectors = [];
end

if nargin > 3 && ~isempty(varargin{2})
    plot_flag = varargin{2};
else
    plot_flag = false;
end

% ECG lead vectors (in cm). See: https://www.3decgleads.com/ecg-leads,
% from: Eldrige , et al., 2014: https://scst.org.uk/wp-content/uploads/2020/02/SCST_ECG_Recording_Guidelines_2017am.pdf
% and: https://ecgwaves.com/topic/ekg-ecg-leads-electrodes-systems-limb-chest-precordial/
if isfield(lead_vectors, 'V1') && isempty(lead_vectors.V1)
    V1 = lead_vectors.V1;
else
    % V1 = [6, -4, 3]'; % 4th Intercostal Space on the Right Border of the Sternum
    V1 = [10, -10, 60]'; % 4th Intercostal Space on the Right Border of the Sternum
end
V1 = V1/norm(V1);

if isfield(lead_vectors, 'V2')  && isempty(lead_vectors.V2)
    V2 = lead_vectors.V2;
else
    % V2 = [6, 2, 3]'; % 4th Intercostal Space on the Left Border of the Sternum
    V2 = [10, 10, 60]'; % 4th Intercostal Space on the Left Border of the Sternum
end
V2 = V2/norm(V2);

if isfield(lead_vectors, 'V3')  && isempty(lead_vectors.V3)
    V3 = lead_vectors.V3;
else
    % V3 = [7, 5.5, 2.5]';  % midway between ECG chest lead V2(C2) and ECG chest lead V4(C4)
    V3 = [15, 15, 60]';  % midway between ECG chest lead V2(C2) and ECG chest lead V4(C4)
end
V3 = V3/norm(V3);

if isfield(lead_vectors, 'V4')  && isempty(lead_vectors.V4)
    V4 = lead_vectors.V4;
else
    % V4 = [8, 9, 2]';  % 5th Intercostal Space on the Left Mid-Clavicular Line, the vertical line which passes through the center of the Clavicle
    V4 = [20, 25, 60]';  % 5th Intercostal Space on the Left Mid-Clavicular Line, the vertical line which passes through the center of the Clavicle
end
V4 = V4/norm(V4);

if isfield(lead_vectors, 'V5')  && isempty(lead_vectors.V5)
    V5 = lead_vectors.V5;
else
    % V5 = [8, 13, 2]';  % positioned at the same level as the ECG chest lead V4(C4) on the Left Anterior Axillary Line, the vertical line that goes along the Anterior Axillary Fold
    V5 = [20, 35, 60]';  % positioned at the same level as the ECG chest lead V4(C4) on the Left Anterior Axillary Line, the vertical line that goes along the Anterior Axillary Fold
end
V5 = V5/norm(V5);

if isfield(lead_vectors, 'V6')  && isempty(lead_vectors.V6)
    V6 = lead_vectors.V6;
else
    % V6 = [8, 16, 2]'; % positioned at the same level as ECG chest lead V4(C4) and ECG chest lead V5(C5) on the Left Mid-Axillary Line, the vertical line that goes midway between the Anterior Axillary Fold and the Posterior Axillary Fold
    V6 = [20, 45, 60]'; % positioned at the same level as ECG chest lead V4(C4) and ECG chest lead V5(C5) on the Left Mid-Axillary Line, the vertical line that goes midway between the Anterior Axillary Fold and the Posterior Axillary Fold
end
V6 = V6/norm(V6);

if isfield(lead_vectors, 'I')  && isempty(lead_vectors.I)
    I = lead_vectors.I;
else
    I = [0, 1, 0]'; % normalized electrical potential difference between the left arm (LA) and right arm (RA)
end
I = I/norm(I);

if isfield(lead_vectors, 'II')  && isempty(lead_vectors.II)
    II = lead_vectors.II;
else
    II = [0, cos(pi/3), -sin(pi/3)]'; % normalized electrical potential difference between the right arm (RA) and the left leg (LL)
end
II = II/norm(II);

if isfield(lead_vectors, 'III')  && isempty(lead_vectors.III)
    III = lead_vectors.III;
else
    III = II - I; % Einthoven s law
end
III = III/norm(III);

if isfield(lead_vectors, 'aVL')  && isempty(lead_vectors.aVL)
    aVL = lead_vectors.aVL;
else
    aVL = (I - III) / 2; % Goldberger's equation
end
aVL = aVL/norm(aVL);

if isfield(lead_vectors, 'aVR')  && isempty(lead_vectors.aVR)
    aVR = lead_vectors.aVR;
else
    aVR = -(I + II) / 2; % Goldberger's equation
end
aVR = aVR/norm(aVR);

if isfield(lead_vectors, 'aVF')  && isempty(lead_vectors.aVF)
    aVF = lead_vectors.aVF;
else
    aVF = (II + III) / 2; % Goldberger's equation
end
aVF = aVF/norm(aVF);

if isfield(lead_vectors, 'Vx')  && isempty(lead_vectors.Vx)
    Vx = lead_vectors.Vx;
else
    Vx = [1, 0, 0]'; % anteroposterior direction, Frank lead
end
Vx = Vx/norm(Vx);

if isfield(lead_vectors, 'Vy')  && isempty(lead_vectors.Vy)
    Vy = lead_vectors.Vy;
else
    Vy = [0, 1, 0]'; % horizontal direction, Frank lead
end
Vy = Vy/norm(Vy);

if isfield(lead_vectors, 'Vz')  && isempty(lead_vectors.Vz)
    Vz = lead_vectors.Vz;
else
    Vz = [0, 0, 1]'; % frontal or vertical direction direction, Frank lead
end
Vz = Vz/norm(Vz);

vcg = zeros(3, size(ecg, 2));
mean_lead = zeros(3, 1);
% construct weighted sum of all leads
for ch = 1 : size(ecg, 1)
    leadname = lower(erase(lead_names{ch}, '_'));
    switch leadname
        case 'i'
            vcg = vcg + I * ecg(ch, :);
            mean_lead = mean_lead + I;
        case 'ii'
            vcg = vcg + II * ecg(ch, :);
            mean_lead = mean_lead + II;
        case 'iii'
            vcg = vcg + III * ecg(ch, :);
            mean_lead = mean_lead + III;
        case 'avr'
            vcg = vcg + aVR * ecg(ch, :);
            mean_lead = mean_lead + aVR;
        case 'avl'
            vcg = vcg + aVL * ecg(ch, :);
            mean_lead = mean_lead + aVL;
        case 'avf'
            vcg = vcg + aVF * ecg(ch, :);
            mean_lead = mean_lead + aVF;
        case 'v1'
            vcg = vcg + V1 * ecg(ch, :);
            mean_lead = mean_lead + V1;
        case 'v2'
            vcg = vcg + V2 * ecg(ch, :);
            mean_lead = mean_lead + V2;
        case 'v3'
            vcg = vcg + V3 * ecg(ch, :);
            mean_lead = mean_lead + V3;
        case 'v4'
            vcg = vcg + V4 * ecg(ch, :);
            mean_lead = mean_lead + V4;
        case 'v5'
            vcg = vcg + V5 * ecg(ch, :);
            mean_lead = mean_lead + V5;
        case 'v6'
            vcg = vcg + V6 * ecg(ch, :);
            mean_lead = mean_lead + V6;
        case 'vx'
            vcg = vcg + Vx * ecg(ch, :);
            mean_lead = mean_lead + Vx;
        case 'vy'
            vcg = vcg + Vy * ecg(ch, :);
            mean_lead = mean_lead + Vy;
        case 'vz'
            vcg = vcg + Vz * ecg(ch, :);
            mean_lead = mean_lead + Vz;
        otherwise
            error(['undefined lead name: ' lead_names{ch}])
    end
end
vcg = vcg ./ mean_lead; % weighted mean of the available leads

if plot_flag
    figure
    hold on
    h_v1 = quiver3(0, 0, 0, V1(1), V1(2), V1(3), 'linewidth', 2);
    h_v2 = quiver3(0, 0, 0, V2(1), V2(2), V2(3), 'linewidth', 2);
    h_v3 = quiver3(0, 0, 0, V3(1), V3(2), V3(3), 'linewidth', 2);
    h_v4 = quiver3(0, 0, 0, V4(1), V4(2), V4(3), 'linewidth', 2);
    h_v5 = quiver3(0, 0, 0, V5(1), V5(2), V5(3), 'linewidth', 2);
    h_v6 = quiver3(0, 0, 0, V6(1), V6(2), V6(3), 'linewidth', 2);
    h_i = quiver3(0, 0, 0, I(1), I(2), I(3));
    h_ii = quiver3(0, 0, 0, II(1), II(2), II(3));
    h_iii = quiver3(0, 0, 0, III(1), III(2), III(3));
    h_avl = quiver3(0, 0, 0, aVL(1), aVL(2), aVL(3), '--');
    h_avr = quiver3(0, 0, 0, aVR(1), aVR(2), aVR(3), '--');
    h_avf = quiver3(0, 0, 0, aVF(1), aVF(2), aVF(3), '--');
    h_vx = quiver3(0, 0, 0, Vx(1), Vx(2), Vx(3), ':', 'linewidth', 2);
    h_vy = quiver3(0, 0, 0, Vy(1), Vy(2), Vy(3), ':', 'linewidth', 2);
    h_vz = quiver3(0, 0, 0, Vz(1), Vz(2), Vz(3), ':', 'linewidth', 2);
    legend([h_v1, h_v2, h_v3, h_v4, h_v5, h_v6, h_i, h_ii, h_iii, h_avl, h_avr, h_avf, h_vx, h_vy, h_vz], {'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'I', 'II', 'III', 'aVL', 'aVR', 'aVF', 'Vx', 'Vy', 'Vz'});
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid
    set(gca, 'fontsize', 16)
    title('ECG lead vectors');
    xlim([0, 2])
    ylim([-1, 1])
    zlim([-1, 1])

    figure
    plot3(vcg(1, :), vcg(2, :), vcg(3, :))
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid
    set(gca, 'fontsize', 16)
    title('Reconstructed VCG');

end
