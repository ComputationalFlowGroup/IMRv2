% Compare multiple alternative graded functions
% Define parameters
G_0 = 1;          % Example value for G_0
G_1 = 2;          % Example value for G_1
a = 2;            % Example value for a
n = 0.3;            % Example value for n
f = linspace(0, 1, 100);  % Create a range of f from 0 to 1

% Function 1: G = G0 + (G1 - G0)
G = G_0 + (G_1 - G_0);

% Function 2: y = G * (1 / (1 + exp(-a(2f - 1))))
y1 = G * (1 ./ (1 + exp(-a * (2 * f - 1))));

% Function 3: y = G * (1/2 * (1 + tanh(a(2f - 1))))
y2 = G * (1/2 * (1 + tanh(a * (2 * f - 1))));

% Function 4: y = G * (3f^2 - 2f^3)
y3 = G * (3 * f.^2 - 2 * f.^3);

% Function 5: y = G * (6f^5 - 15f^4 + 10f^3)
y4 = G * (6 * f.^5 - 15 * f.^4 + 10 * f.^3);

% Function 6: y = G * (f^a / (f^n + (1 - f)^a))
y5 = G * (f.^a ./ (f.^n + (1 - f).^a));

% Function 7: y = G * (1/pi * atan(a(2f - 1)) + 1/2)
y6 = G * (1/pi * atan(a * (2 * f - 1)) + 1/2);

% Function 8: y = G * (1/2 * (1 + erf(a(2f - 1))))
y7 = G * (1/2 * (1 + erf(a * (2 * f - 1))));

% Function 9: y = G * f^n
y8 = G * f.^n;

% Function 10: y = G * f^a
y9 = G * f.^a;

% Function 11: y = G * exp(-(f/n)^a)
y10 = G * exp(-((f / n).^a));

% Function 12: y = G * ((ln(1 + f^a)) / ln(1 + 1^a))^n
y11 = G * ((log(1 + f.^a) ./ log(1 + 1.^a)).^n);

% Function 13: y = G * ((ln(1 + f^n)) / ln(1 + 1^n))^a
y12 = G * ((log(1 + f.^n) ./ log(1 + 1.^n)).^a);

% Function 14: y = G * f^a * (3 - 2f)
y13 = G * f.^a .* (3 - 2 * f);

% Function 15: y = (1 + f^a)^((n-1)/a)
y14 = G * (1 + f.^a).^( (n-1)/a );


% so don't have to use default colors
RBG = get(groot,"FactoryAxesColorOrder");
H = compose("#%02X%02X%02X",round(RBG*255));

% Plotting all functions for comparison
figure;

plot(f, y1, 'r', 'LineWidth', 2); hold on;
plot(f, y2, 'g', 'LineWidth', 2);
plot(f, y3, 'b', 'LineWidth', 2);
plot(f, y4, 'm', 'LineWidth', 2);
plot(f, y5, 'c', 'LineWidth', 2);
plot(f, y6, 'k', 'LineWidth', 2);
plot(f, y7, 'y', 'LineWidth', 2);
plot(f, y8, 'b--', 'LineWidth', 2);
plot(f, y9, 'g--', 'LineWidth', 2);
plot(f, y10, 'r--', 'LineWidth', 2);
plot(f, y11, 'm--', 'LineWidth', 2);
plot(f, y12, 'c--', 'LineWidth', 2);
plot(f, y13, 'k--', 'LineWidth', 2);
plot(f, y14, 'y--', 'LineWidth', 2);

% Labeling the plot
title('Comparison of Functions');
xlabel('f');
ylabel('y');
legend('y = G * (1/(1 + exp(-a(2f-1))))', 'y = G * (1/2 * (1 + tanh(a(2f-1))))', ...
    'y = G * (3f^2 - 2f^3)', 'y = G * (6f^5 - 15f^4 + 10f^3)', ...
    'y = G * (f^a / (f^n + (1 - f)^a))', 'y = G * (1/pi * atan(a(2f - 1)) + 1/2)', ...
    'y = G * (1/2 * (1 + erf(a(2f - 1))))', 'y = G * f^n', ...
    'y = G * f^a', 'y = G * exp(-(f/n)^a)', 'y = G * ((ln(1 + f^a)) / ln(1 + 1^a))^n', ...
    'y = G * ((ln(1 + f^n)) / ln(1 + 1^n))^a', 'y = G * f^a * (3 - 2f)', ...
    'y = G * (1 + f^a)^((n-1)/a)','Location', 'northeastoutside');

grid on;
hold off;

%%
RGB = get(groot, 'FactoryAxesColorOrder');
H = compose("#%02X%02X%02X", round(RGB*255));

% Display the colors in HEX format
disp(H);

% Now using these colors in the plot
figure;

% Data and color settings (using custom colors)
y_data = {y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14};
line_styles = repmat({'-'}, 1, 14);  % Default style '-'
line_styles(8:14) = {'--'};  % Change to dashed lines from y8 to y14

% Plot all in a loop with custom colors
for i = 1:length(y_data)
    plot(f, y_data{i}, 'Color', H{i}, 'LineWidth', 2, 'LineStyle', line_styles{i});
    hold on;
end
%%
% Define custom color set (14 distinct colors in RGB format)
custom_colors = [
    [0.5, 0.0, 0.0];  % Dark Red
    [0.0, 0.5, 0.0];  % Dark Green
    [0.0, 0.0, 0.5];  % Dark Blue
    [0.5, 0.0, 0.5];  % Dark Magenta
    [0.0, 0.5, 0.5];  % Dark Cyan
    [0.5, 0.5, 0.0];  % Olive
    [0.5, 0.5, 0.5];  % Gray
    [1.0, 0.0, 0.0];  % Bright Red
    [0.0, 1.0, 0.0];  % Bright Green
    [0.0, 0.0, 1.0];  % Bright Blue
    [1.0, 0.0, 1.0];  % Bright Magenta
    [0.0, 1.0, 1.0];  % Bright Cyan
    [1.0, 0.5, 0.0];  % Orange
    [0.5, 0.0, 1.0];  % Purple
    [0.8, 0.8, 0.2];  % Light Olive
    [0.2, 0.8, 0.8];  % Light Cyan
    [0.8, 0.2, 0.8];  % Light Magenta
    [0.8, 0.5, 0.5];  % Light Red
];

% Plot using the custom colors
figure;
y_data = {y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14};
line_styles = repmat({'-'}, 1, 14);  % Default style '-'
line_styles(8:14) = {'--'};  % Change to dashed lines from y8 to y14

for i = 1:length(y_data)
    plot(f, y_data{i}, 'Color', custom_colors(i,:), 'LineWidth', 2, 'LineStyle', line_styles{i});
    hold on;
end

legend('y = G * (1/(1 + exp(-a(2f-1))))', 'y = G * (1/2 * (1 + tanh(a(2f-1))))', ...
    'y = G * (3f^2 - 2f^3)', 'y = G * (6f^5 - 15f^4 + 10f^3)', ...
    'y = G * (f^a / (f^n + (1 - f)^a))', 'y = G * (1/pi * atan(a(2f - 1)) + 1/2)', ...
    'y = G * (1/2 * (1 + erf(a(2f - 1))))', 'y = G * f^n', ...
    'y = G * f^a', 'y = G * exp(-(f/n)^a)', 'y = G * ((ln(1 + f^a)) / ln(1 + 1^a))^n', ...
    'y = G * ((ln(1 + f^n)) / ln(1 + 1^n))^a', 'y = G * f^a * (3 - 2f)', ...
    'y = G * (1 + f^a)^((n-1)/a)','Location', 'northeastoutside');
%ylim([G_0-0.5 G_1+0.5]);

%%
% Parameters
a = 2;
n = 0.3;
f = linspace(0, 10, 1000);  % f from 0 to ∞ (approximated by 0 to 10)

% Function 1: Sigmoid
y1 = 1 ./ (1 + exp(-a * (2 * f - 1)));  % Already in [0,1]

% Function 2: Tanh-based sigmoid
y2 = 0.5 * (1 + tanh(a * (2 * f - 1)));  % Already in [0,1]

% Function 3: Smoothstep (classic)
s = f ./ (1 + f);  % Rescale to [0,1] for input to smoothstep
y3 = 3 * s.^2 - 2 * s.^3;

% Function 4: Smootherstep
y4 = 6 * s.^5 - 15 * s.^4 + 10 * s.^3;

% Function 5: Normalized Rational Sigmoid
y5 = (f.^a) ./ (f.^a + (1 + f).^n);  % Adjusted denominator to guarantee max < ∞

% Function 6: Arctangent-based
y6 = (1/pi) * atan(a * f);  % Range (0 to 0.5), need scaling
y6 = y6 / max(y6);  % Normalize to [0,1]

% Function 7: Erf-based sigmoid
y7 = 0.5 * (1 + erf(a * (2 * f - 1)));  % Already in [0,1]

% Function 8: Power law (normalized)
y8 = f.^n ./ (1 + f.^n);  % Always in [0,1]

% Function 9: Power law (stronger)
y9 = f.^a ./ (1 + f.^a);  % Always in [0,1]

% Function 10: Exponential decay (flipped)
y10 = 1 - exp(- (f / n).^a);  % Range [0,1)

% Function 11: Normalized Logarithmic Growth (n on outside)
y11 = (log(1 + f.^a) ./ log(1 + Inf.^a)).^n;  % log(1 + Inf) → 1 at ∞
y11 = (log(1 + f.^a) ./ log(1 + 1^a)).^n;
y11 = y11 / max(y11);  % Normalize

% Function 12: Same but a outside
y12 = (log(1 + f.^n) ./ log(1 + 1^n)).^a;
y12 = y12 / max(y12);  % Normalize

% Function 13: f^a * (3 - 2f), normalized
numerator = f.^a .* (3 - 2 * s);
y13 = numerator / max(numerator);

% Function 14: y = (1 + f^a)^((n-1)/a), normalized
y14 = (1 + f.^a).^((n - 1) / a);
y14 = (y14 - min(y14)) / (max(y14) - min(y14));  % Normalize to [0,1]

% Plotting all functions for comparison
figure;

plot(f, y1, 'r', 'LineWidth', 2); hold on;
plot(f, y2, 'g', 'LineWidth', 2);
plot(f, y3, 'b', 'LineWidth', 2);
plot(f, y4, 'm', 'LineWidth', 2);
plot(f, y5, 'c', 'LineWidth', 2);
plot(f, y6, 'k', 'LineWidth', 2);
plot(f, y7, 'y', 'LineWidth', 2);
plot(f, y8, 'b--', 'LineWidth', 2);
plot(f, y9, 'g--', 'LineWidth', 2);
plot(f, y10, 'r--', 'LineWidth', 2);
plot(f, y11, 'm--', 'LineWidth', 2);
plot(f, y12, 'c--', 'LineWidth', 2);
plot(f, y13, 'k--', 'LineWidth', 2);
plot(f, y14, 'y--', 'LineWidth', 2);

% Labeling the plot
title('Comparison of Functions');
xlabel('f');
ylabel('y');
legend('y = G * (1/(1 + exp(-a(2f-1))))', 'y = G * (1/2 * (1 + tanh(a(2f-1))))', ...
    'y = G * (3f^2 - 2f^3)', 'y = G * (6f^5 - 15f^4 + 10f^3)', ...
    'y = G * (f^a / (f^n + (1 - f)^a))', 'y = G * (1/pi * atan(a(2f - 1)) + 1/2)', ...
    'y = G * (1/2 * (1 + erf(a(2f - 1))))', 'y = G * f^n', ...
    'y = G * f^a', 'y = G * exp(-(f/n)^a)', 'y = G * ((ln(1 + f^a)) / ln(1 + 1^a))^n', ...
    'y = G * ((ln(1 + f^n)) / ln(1 + 1^n))^a', 'y = G * f^a * (3 - 2f)', ...
    'y = G * (1 + f^a)^((n-1)/a)','Location', 'northeastoutside');

grid on;
hold off;