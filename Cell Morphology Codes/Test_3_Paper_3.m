% Load and convert image
filename = '1_20x.jpg';
img = imread(filename);

if size(img, 3) == 3
    gray = rgb2gray(img);
else
    gray = img;
end

% Enhance dark elongated features (cells) via inversion and local contrast
gray_inv = imcomplement(gray);                     % dark -> light
gray_eq = adapthisteq(gray_inv);                   % local contrast adjustment

% Threshold and clean binary image
bw = imbinarize(gray_eq);                          % threshold image
bw = imopen(bw, strel('disk', 2));                 % remove small connections
bw = imfill(bw, 'holes');                          % fill small holes
bw = bwareaopen(bw, 100);                          % remove small regions

% Extract region properties
stats = regionprops(bw, 'Area', 'Eccentricity', ...
    'MajorAxisLength', 'MinorAxisLength', ...
    'PixelList', 'Centroid');

% Apply filtering based on morphology
ecc = [stats.Eccentricity];
area = [stats.Area];
aspect_ratio = [stats.MajorAxisLength] ./ [stats.MinorAxisLength];
valid_idx = (ecc > 0.75) & (area > 150) & (aspect_ratio < 8);

stats = stats(valid_idx);  % Filtered cell candidates
fprintf('Found %d valid regions\n', numel(stats));

% Fit ellipses
E = cell(numel(stats), 1);
for k = 1:numel(stats)
    xy = double(stats(k).PixelList);
    E{k} = fit_ellipse(xy(:,1), xy(:,2));
end

% Plot results
imshow(gray); hold on;
title('Final Fitted Ellipses on Elongated Dark Cells');

for k = 1:numel(E)
    if isempty(E{k}) || ~isfield(E{k}, 'Xdata'), continue; end
    plot(E{k}.Xdata, E{k}.Ydata, 'r-', 'LineWidth', 1.5);
    text(E{k}.X0, E{k}.Y0, sprintf('%d', k), ...
         'Color', 'y', 'FontSize', 8, 'HorizontalAlignment', 'center');
end

% --- Fit ellipse function ---
function ellipse_t = fit_ellipse(x, y)
    x = x(:);
    y = y(:);
    D = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
    S = D'*D;
    C = zeros(6); C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

    [gevec, geval] = eig(S,C);
    I = find(real(diag(geval)) > 0 & ~isinf(diag(geval)));
    if isempty(I), ellipse_t = []; return; end
    a = real(gevec(:, I(1)));

    b = 0.5 * atan2(a(2), a(1) - a(3));
    cos_phi = cos(b); sin_phi = sin(b);
    Ao = a(6);
    Au = a(4)*cos_phi + a(5)*sin_phi;
    Av = -a(4)*sin_phi + a(5)*cos_phi;
    Auu = a(1)*cos_phi^2 + a(3)*sin_phi^2 + a(2)*cos_phi*sin_phi;
    Avv = a(1)*sin_phi^2 + a(3)*cos_phi^2 - a(2)*cos_phi*sin_phi;

    tuCentre = -0.5 * [Au/Auu; Av/Avv];
    xy = [cos_phi, sin_phi; -sin_phi, cos_phi]' * tuCentre;

    ellipse_t.X0 = xy(1);
    ellipse_t.Y0 = xy(2);
    ellipse_t.phi = b;
    ellipse_t.a = sqrt(abs(-Ao/Auu));
    ellipse_t.b = sqrt(abs(-Ao/Avv));

    theta = linspace(0, 2*pi, 100);
    ellipse_t.Xdata = ellipse_t.X0 + ellipse_t.a * cos(theta) * cos(b) - ellipse_t.b * sin(theta) * sin(b);
    ellipse_t.Ydata = ellipse_t.Y0 + ellipse_t.a * cos(theta) * sin(b) + ellipse_t.b * sin(theta) * cos(b);
end