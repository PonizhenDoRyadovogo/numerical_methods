function [result_x, result_y, count_iteration] = systemSimpleIterations(psi1, psi2, x0, y0, epsilon)
    current_x = psi1(x0, y0);
    current_y = psi2(x0, y0);
    next_x = psi1(current_x, current_y);
    next_y = psi2(current_x, current_y);
    count_iteration = 1;
    while(max(abs(next_x - current_x), abs(next_y - current_y)) >= epsilon)
        current_x = next_x;
        current_y = next_y;
        next_x = psi1(current_x, current_y);
        next_y = psi2(current_x, current_y);
        count_iteration = count_iteration + 1;
    end
    result_x = next_x;
    result_y = next_y;
end

function [result_x, result_y, count_iteration] = systemNewton

end

example_psi1 = @(x, y) (0.7 - cos(y + 1)) / 3;
example_psi2 = @(x,y) (2 - sin(x)) / 2;
epsilon = 0.5 * 10 ^ (-4);
[result_x, result_y, count_itertion] = systemSimpleIterations(example_psi1, example_psi2, 0.2, 0.7, epsilon);
fprintf('Пример из учебника: |f1(x,y) = sin(x) + 2y - 2 = 0\n                    |f2(x,y) = cos(y+1) + 3x - 0.7 = 0')
fprintf('\nМетодом простых итераций при точке начального приближения [0.2;0.7]:\n');
fprintf('x = %.5f', result_x);
fprintf('\ny = %.5f', result_y);
fprintf('\nколичество итераций: %d', count_itertion);
equ1 = @(x, y) sin(x) + 2 * y - 2;
res = equ1(result_x, result_y);
fprintf('\nпроверка подстановкой в одно из уравнений: %.5f\n', res);

psi1 = @(x, y) 0.5 + 0.5 * sin(y - 0.5);
psi2 = @(x, y) 0.375 - 0.25 * cos(x);
[result_x, result_y, count_itertion] = systemSimpleIterations(psi1, psi2, 0.3229, 0.138, epsilon);
fprintf('\nПример: |f1(x,y) = cos(x) + 4y - 1.5 = 0\n        |f2(x,y) = 2x - sin(y-0.5) - 1 = 0')
fprintf('\nМетодом простых итераций при точке начального приближения [0.3229;0.138]:\n');
fprintf('x = %.5f', result_x);
fprintf('\ny = %.5f', result_y);
fprintf('\nколичество итераций: %d', count_itertion);

[result_x, result_y, count_itertion] = systemSimpleIterations(psi1, psi2, -0.1771, 0.638, epsilon);
fprintf('\nМетодом простых итераций при точке начального приближения [0.3229;0.138]:\n');
fprintf('x = %.5f', result_x);
fprintf('\ny = %.5f', result_y);
fprintf('\nколичество итераций: %d', count_itertion);
equ2 = @(x, y) cos(x) + 4 * y - 1.5;
res = equ2(result_x, result_y);
fprintf('\nпроверка подстановкой в одно из уравнений: %.5f', res);

x = linspace(0, 1, 10000);
plot(x, 0.375 - (cos(x)/4), 'r', 'LineWidth', 2);
hold on;
plot(x, asin(2 * x - 1) + 0.5, 'b', 'LineWidth', 2);
hold off;
grid on;
legend('0.375 - cos(x)/4', 'arcsin(2x-1) + 0.5');
xlabel('x');
ylabel('y');
title('Графики функций y = 0.375 - cos(x)/4 и y = arcsin(2x-1) + 0.5');
