grayImage = imread('1_20x.jpg');

if size(grayImage, 3) == 3
    grayImage = rgb2gray(grayImage);
end

% Step 1: Preprocess
bw = imbinarize(grayImage);
bw = bwareaopen(bw, 50);
labeled = bwlabel(bw);
stats = regionprops(labeled, 'PixelList');

% Step 2: Extract segments and fit ellipses
L = {stats.PixelList};
G = length(L);
E = cell(G, 1);
for k = 1:G
    E{k} = fit_ellipse(double(L{k}(:,1)), double(L{k}(:,2)));
end

% Step 3: Merge ellipses based on rules
% [L, E] = mergeEllipses(L, E);

% Step 4: Refine the ellipses based on fitting error
% E = refineEllipses(L, E, 0.2);

% Optional: Visualize results
% imshow(grayImage); hold on;
% for k = 1:length(E)
%     plot(E{k}.Xdata, E{k}.Ydata, 'LineWidth', 1.5);
% end
% hold off;

% --- Visualize Ellipses on Top of the Grayscale Image ---
figure;
imshow(grayImage); hold on;
title('Fitted Ellipses');

for k = 1:length(E)
    if isempty(E{k}) || ~isfield(E{k}, 'Xdata') || isempty(E{k}.Xdata)
        continue;
    end
    plot(E{k}.Xdata, E{k}.Ydata, 'r-', 'LineWidth', 1.5);  % ellipse outline
    plot(E{k}.X0, E{k}.Y0, 'b+', 'MarkerSize', 6, 'LineWidth', 1);  % ellipse center
    text(E{k}.X0, E{k}.Y0, num2str(k), 'Color', 'y', 'FontSize', 10);  % index label
end

hold off;



% --- Function Definitions Below ---

function ellipse_t = fit_ellipse(x, y)
x = x(:);
y = y(:);
D = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
S = D'*D;
C = zeros(6);
C(1,3) = -2;
C(2,2) = 1;
C(3,1) = -2;

[gevec, geval] = eig(S,C);
I = find(real(diag(geval)) > 0 & ~isinf(diag(geval)));
a = real(gevec(:, I(1)));

b = 0.5 * atan2(a(2), a(1) - a(3));
cos_phi = cos(b);
sin_phi = sin(b);
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
ellipse_t.Xdata = ellipse_t.X0 + ellipse_t.a * cos(theta) * cos(ellipse_t.phi) - ellipse_t.b * sin(theta) * sin(ellipse_t.phi);
ellipse_t.Ydata = ellipse_t.Y0 + ellipse_t.a * cos(theta) * sin(ellipse_t.phi) + ellipse_t.b * sin(theta) * cos(ellipse_t.phi);
end

function [L, E] = mergeEllipses(L, E)
G = length(L);
i = 1;

while i <= G
    merged = false;
    for j = i+1:G
        Li = L{i};
        Lj = L{j};
        Ei = E{i};
        Ej = E{j};
        allPoints = [Li; Lj];
        Eij = fit_ellipse(double(allPoints(:,1)), double(allPoints(:,2)));

        if satisfies_case1(Ei, Ej, Eij) && ...
           ~satisfies_case2(Ei, Ej, Eij) && ...
           ~satisfies_case3(Ei, Ej, Eij)
            continue;
        elseif satisfies_case2(Ei, Ej, Eij) || ...
               satisfies_case3(Ei, Ej, Eij)
            L{i} = allPoints;
            E{i} = Eij;
            L(j) = [];
            E(j) = [];
            G = G - 1;
            i = 1;
            merged = true;
            break;
        end
    end
    if ~merged
        i = i + 1;
    end
end
end

function E = refineEllipses(L, E, disThRe)
M = length(L);
for i = 1:M
    Li = L{i};
    bestError = Inf;
    bestE = [];

    for j = 1:length(E)
        Ej = E{j};
        if isempty(Ej)
            continue;
        end
        PLi = [Ej.Xdata', Ej.Ydata'; Li];
        newE = fit_ellipse(double(PLi(:,1)), double(PLi(:,2)));
        disj = ellipse_fit_error(newE, PLi);
        
        if disj < bestError
            bestError = disj;
            bestE = newE;
        end
    end

    if bestError < disThRe
        E{i} = bestE;
    end
end
end

function isTrue = satisfies_case1(Ei, Ej, Eij)
isTrue = abs(Ei.a - Eij.a) < 5 && abs(Ej.b - Eij.b) < 5;
end

function isTrue = satisfies_case2(Ei, Ej, Eij)
isTrue = Eij.a > 100;
end

function isTrue = satisfies_case3(Ei, Ej, Eij)
isTrue = Eij.b < 5;
end

function error = ellipse_fit_error(E, points)
xc = E.X0_in;
yc = E.Y0_in;
a = E.a;
b = E.b;
theta = E.phi;
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
ptsRot = (R' * ([points(:,1) - xc, points(:,2) - yc])')';

ellipseEq = (ptsRot(:,1)/a).^2 + (ptsRot(:,2)/b).^2;
error = mean(abs(ellipseEq - 1));
end