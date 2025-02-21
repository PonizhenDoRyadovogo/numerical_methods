function sufficientCondition(A)
    norm = max(sum(abs(A), 2));
    fprintf('Максимальная сумма модулей элементов в строках: %g ', norm);
    if norm > 1
        fprintf('> 1 => Необходимое условие не выполняется!\n');
    else
        fprintf('< 1 => Необходимое условие выполняется!\n');
    end
end

function [result_X, count_iteration] = simpleIteration(B, b, X0, epsilon)
    sufficientCondition(B);
    current_X = X0;
    next_X = B * X0 + b;
    count_iteration = 1;
    while max(abs(next_X - current_X)) > epsilon
        current_X = next_X;
        next_X = B * current_X + b;
        count_iteration = count_iteration + 1;
    end
    result_X = next_X;
end

function [result_X, count_iteration] = zeidelMethod(B, b, X0, epsilon)
    sufficientCondition(B);
    H = tril(B, -1);
    F = B - H;
    E = eye(length(b));
    current_X = X0;
    next_X = inv(E - H) * F * current_X + inv(E - H) * b;
    count_iteration = 1;
    while max(abs(next_X - current_X)) > epsilon
        current_X = next_X;
        next_X = inv(E - H) * F * current_X + inv(E - H) * b;
        count_iteration = count_iteration + 1;
    end
    result_X = next_X;
end

% x1 = 0.42*x1 - 0.52*x2 + 0.03*x3 + 0.44
% x2 = 0.31*x1 - 0.26*x2 - 0.36*x3 + 1.42
% x3 = 0.12*x1 + 0.08*x2 -0.14*x3 -0.24*x4 - 0.83
% x4 = 0.15*x1 -0.35*x2 - 0.18*x3 - 1.42

B = [0.42, -0.52, 0.03, 0; 0.31, -0.26, -0.36, 0; 0.12, 0.08, -0.14, -0.24; 0.15, -0.35, -0.18, 0];
b = [0.44; 1.42; -0.83; -1.42];
initial_X = [0; 0; 0; 0];

[X, k] = simpleIteration(B, b, initial_X, 1e-6);
fprintf('При точке начального приближения X0 = \n');
disp(initial_X);
fprintf('решение с помощью метода простых итераций:\n');
disp(X);
fprintf('количество итераций: %d\n', k);

otherB = [-0.58, -0.52, 0.03, 0; 0.31, -1.26, -0.36, 0; 0.12, 0.08, -1.14, -0.24; 0.15, -0.35, -0.18, -1];
check = otherB * X + b;
fprintf('Проверка:\n');
for i = 1:length(b)
    fprintf("%.4f\n", check(i));
end

[X, k] = simpleIteration(B, b, b, 1e-6);
fprintf('При точке начального приближения X0 = \n');
disp(b);
fprintf('решение с помощью метода простых итераций:\n');
disp(X);
fprintf('количество итераций: %d\n', k);

otherB = [-0.58, -0.52, 0.03, 0; 0.31, -1.26, -0.36, 0; 0.12, 0.08, -1.14, -0.24; 0.15, -0.35, -0.18, -1];
check = otherB * X + b;
fprintf('Проверка:\n');
for i = 1:length(b)
    fprintf("%.4f\n", check(i));
end
disp("==============================================================");


[X, k] = zeidelMethod(B, b, initial_X, 1e-6);
fprintf('При точке начального приближения X0 = \n');
disp(initial_X);
fprintf('решение с помощью метода Зейделя:\n');
disp(X);
fprintf('количество итераций: %d\n', k);
check = otherB * X + b;
fprintf('Проверка:\n');
for i = 1:length(b)
    fprintf("%.4f\n", check(i));
end

[X, k] = zeidelMethod(B, b, b, 1e-6);
fprintf('При точке начального приближения X0 = \n');
disp(b);
fprintf('решение с помощью метода Зейделя:\n');
disp(X);
fprintf('количество итераций: %d\n', k);
check = otherB * X + b;
fprintf('Проверка:\n');
for i = 1:length(b)
    fprintf("%.4f\n", check(i));
end
disp("==============================================================");