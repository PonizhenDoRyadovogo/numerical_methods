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
    lambda = eig(B);
    H = tril(B, -1);
    F = B - H;
    E = eye(length(b));
    null_space = null(F + lambda(2) * H - lambda(2) * E);
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

function result_X = sqrtMethod(A, f)
    n = length(f);
    S = zeros(n,n);
    D = zeros(n,n);
    D(1,1) = sign(A(1,1));
    S(1,1) = sqrt(A(1,1));
    for j = 2:n
        S(1, j) = A(1, j) / S(1,1);
    end
    for i = 2:n
        tmp_sum = 0;
        tmp_sum2 = 0;
        for p = 1:i-1
            tmp_sum = tmp_sum + S(p, i) ^ 2;
            tmp_sum2 = tmp_sum2 + (S(p,i) ^ 2) * D(p,p);
        end
        S(i, i) = sqrt(A(i,i) - tmp_sum);
        D(i,i) = sign(A(i,i) - tmp_sum2);
        for j = i+1:n
            tmp_sum = 0;
            for p = 1:i-1
                tmp_sum = tmp_sum + S(p, i) * S(p, j);
            end
            S(i, j) = (A(i, j) - tmp_sum) / S(i,i);
        end
    end

    B = S' * D;
    y = zeros(n, 1);
    y(1) = f(1) / B(1,1);
    for k = 2:n
        tmp_sum = 0;
        for p = 1:k-1
            tmp_sum = tmp_sum + B(k,p) * y(p);
        end
        y(k) = (f(k) - tmp_sum) / B(k,k);
    end

    result_X = zeros(n,1);
    result_X(n) = y(n) / S(n,n);
    for k = n-1:-1:1
        tmp_sum = 0;
        for p = k+1:n
            tmp_sum = tmp_sum + S(k,p) * result_X(p);
        end
        result_X(k) = (y(k) - tmp_sum) / S(k,k);
    end
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
disp("====================Метод квадратного корня===========================");

% 3.23*x1 + 1.62*x2 - 0.65*x3 = 1.28
% 1.62*x1 - 2.33*x2 - 1.43*x3 = 0.87
% -0.65*x1 - 1.43*x2 + 2.18*x3 = -2.87

A = [3.23, 1.62, -0.65; 1.62, -2.33, -1.43; -0.65, -1.43, 2.18];
f = [1.28; 0.87; -2.87];

X = sqrtMethod(A, f);
fprintf('Система:\n');
fprintf('3.23*x1 + 1.62*x2 - 0.65*x3 = 1.28\n');
fprintf('1.62*x1 - 2.33*x2 - 1.43*x3 = 0.87\n');
fprintf('-0.65*x1 - 1.43*x2 + 2.18*x3 = -2.87\n');
fprintf('решение X = \n');
disp(X);
fprintf('Проверка:\n');
check = A * X - f;
for i = 1:length(f)
    fprintf("%.4f\n", check(i));
end