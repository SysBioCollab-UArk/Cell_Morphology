% --- Improved Ellipse Fitting for Dark Elongated Epithelial Cells ---

% Load and preprocess
img = imread('1_20x.jpg');
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Enhance dark cells (bottom-hat filter emphasizes dark features)
se = strel('disk', 12);
tophatFiltered = imbothat(img, se);

% Contrast adjustment
adjusted = imadjust(tophatFiltered);

% Thresholding and cleanup
bw = imbinarize(adjusted);
bw = bwareaopen(bw, 100);           % remove tiny noise
bw = imclose(bw, strel('disk', 3)); % close small gaps
bw = imfill(bw, 'holes');

% Label and extract region stats
labeled = bwlabel(bw);
stats = regionprops(labeled, 'Area', 'Eccentricity', ...
    'MajorAxisLength', 'MinorAxisLength', 'PixelList', 'Centroid');

% Filter good candidates
minArea = 500;
minEcc = 0.8;
minMinor = 5;

validIdx = find([stats.Area] > minArea & ...
                [stats.Eccentricity] > minEcc & ...
                [stats.MinorAxisLength] > minMinor);

% Visual debug overlay of regions
figure; imshow(img); hold on; title('Filtered Regions');
for i = validIdx
    rectangle('Position', regionprops(labeled==i, 'BoundingBox').BoundingBox, ...
              'EdgeColor', 'g', 'LineWidth', 1);
end
hold off;

% Fit ellipses only to valid regions
E = {};
validStats = stats(validIdx);

for k = 1:numel(validStats)
    pixels = double(validStats(k).PixelList);
    ell = fit_ellipse(pixels(:,1), pixels(:,2));

    % Sanity check on result
    if isfield(ell, 'a') && isfield(ell, 'b') && all(isfinite([ell.a ell.b])) && ...
       ell.a > 5 && ell.a < 300 && ell.b > 3 && ell.b < 300
        E{end+1} = ell;
    end
end

% Display ellipses
figure; imshow(img); hold on;
title('Final Fitted Ellipses on Dark Elongated Cells');
for k = 1:length(E)
    if isfield(E{k}, 'Xdata')
        plot(E{k}.Xdata, E{k}.Ydata, 'r-', 'LineWidth', 1.5);
        plot(E{k}.X0, E{k}.Y0, 'b+');
        text(E{k}.X0, E{k}.Y0, num2str(k), 'Color', 'y', 'FontSize', 9);
    end
end
hold off;

% -------------------------- Helper Function --------------------------
function ellipse_t = fit_ellipse(x, y)
    x = x(:); y = y(:);
    D = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
    S = D'*D;
    C = zeros(6); C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

    [gevec, geval] = eig(S,C);
    I = find(real(diag(geval)) > 0 & ~isinf(diag(geval)));
    if isempty(I)
        ellipse_t = struct(); return;
    end
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
    ellipse_t.X0_in = ellipse_t.X0;
    ellipse_t.Y0_in = ellipse_t.Y0;

    theta = linspace(0, 2*pi, 100);
    ellipse_t.Xdata = ellipse_t.X0 + ellipse_t.a * cos(theta) * cos(b) - ellipse_t.b * sin(theta) * sin(b);
    ellipse_t.Ydata = ellipse_t.Y0 + ellipse_t.a * cos(theta) * sin(b) + ellipse_t.b * sin(theta) * cos(b);
end
