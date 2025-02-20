function [lambda, X] = danilevski(A)
    n = size(A,1);
    if size(A,2) ~= n
        error('Матрица A должна быть квадратной.');
    end
    M = A; 
    % S = S_1 * S_2 * ... * S_{n-1}.
    S = eye(n);    
    for i = 1:(n-1)
        pivot = M(i,n);
        
        if pivot == 0
            % Поищем ненулевой элемент ниже в том же столбце n
            rowSwap = i+1;
            while (rowSwap <= n) && (M(rowSwap,n) == 0)
                rowSwap = rowSwap + 1;
            end
            
            if rowSwap > n
                warning('Шаг %d: в последнем столбце все элементы = 0 от строки %d вниз. Нужны более сложные перестановки!', i, i);
                break; 
            else
                % Меняем строки i и rowSwap
                M([i rowSwap], :) = M([rowSwap i], :);
                % То же самое надо отразить в общем преобразовании S
                S([i rowSwap], :) = S([rowSwap i], :);
                pivot = M(i,n);
            end
        end
        
        if pivot == 0
            warning('Шаг %d: опорный элемент остался 0, пропускаем преобразование.', i);
            continue;
        end
 
        Si = eye(n);
        Si(i,:) = M(i,:) / pivot;  % замена i-й строки
        
        Si_inv = inv(Si);
        M = Si_inv * M;
        M = M * Si;
        
        S = S * Si; 
    end
    
    coeffs = poly(M);               
    lambda = roots(coeffs);
    
    X = zeros(n,n);
    
    for i = 1:n
        lam = lambda(i);
        
        % Найдём y, ненулевое вектор-решение (M - lam*I)*y = 0.
        y_basis = null(M - lam*eye(n));
        if isempty(y_basis)
            warning('Не удалось найти собственный вектор для λ = %g (численные особенности).', lam);
            continue;
        end
       
        y = y_basis(:,1);
        y = y / norm(y);
        
        % Теперь x = S*y
        x = S * y;
        x = x / norm(x);
        
        X(:,i) = x;
    end
end

A = [  1.21   0.45   -0.17   -0.12;
       0.45   1.37   -0.11   0.38;
       -0.17   -0.11   1.44   -0.17;
       -0.12   0.38   -0.17   1.79];

[lambda, V] = danilevski(A);
res = eig(A);
disp('Собственные значения, найденные с помощью метода Данилевского:');
disp(lambda);
disp('Собственные векторы (по столбцам), найденные с помощью метода Данилевского:');
disp(V);
disp('Собственные значения, найденные с помощью встроенных функций:')
disp(res);
% проверка собственных векторов
residual = [];
for i = 1:length(V)
    residual = A * V(:, i) - lambda(i) * V(:, i);
end
disp('Проверка собственных векторов Ax - lambda*x = 0:');
for i = 1:length(residual)
    fprintf('Вектор номер %d: %.10f\n', i, residual(i));
end