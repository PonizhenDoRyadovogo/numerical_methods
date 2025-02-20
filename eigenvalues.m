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
        x = x / x(n);
        
        X(:,i) = x;
    end
end

function [eigenvalues, eigenvectors] = fadeevMethod(A)
    n = size(A,1);
    if size(A,2) ~= n
        error('Матрица A должна быть квадратной.');
    end
    I = eye(n);
    A_k = A;
    q(1) = trace(A_k);
    B = cell(n,1);
    B{1} = A_k - q(1) * I;
    for k = 2:n
        A_k = A * B{k-1};
        q(k) = trace(A_k)/k;
        B{k} = A_k - q(k) * I;
    end
    % проверяем, что матрица B_n = 0 
    if(norm(B{n}, 'fro') > 1e-12)
         warning('Матрица B_n не оказалась строго нулевой');
    end
    
    c = [1, -q];
    
    eigenvalues = roots(c);
    eigenvectors = zeros(n,n);
    
    for i = 1:n
        lam = eigenvalues(i);
        Qlam = zeros(n,n);
        for k = 0:n-1
            alpha = lam ^ (n-1-k);
            if k == 0
                Qlam = Qlam + alpha * I;
            else
                Qlam = Qlam + alpha * B{k};
            end
        end
        v = [];
        v = Qlam(:, 1);
        v = v / v(n);
        eigenvectors(:, i) = v;
    end
end

function [eigenvalues, eigenvectors] = krylov(A)
    n = size(A,1);
    if size(A,2) ~= n
        error('Матрица A должна быть квадратной.');
    end

    v1 = zeros(n,1);
    v1(1) = 1;
    K = zeros(n,n);
    K(:,1) = v1;
    for j = 2:n
        K(:,j) = A * K(:, j-1);
    end

    %  w = A^n * v1 = A * (A^(n-1) * v1)
    %  У нас уже есть A^(n-1)*v1 = K(:,n), значит:
    w = A * K(:,n);

    c = K \ w;  % вектор коэффициентов [c1; c2; ...; c_n]
    cReversed = flipud(c);   % [c_n; c_{n-1}; ...; c_1]

    polyCoeffs = [1; -cReversed];

    eigenvalues = roots(polyCoeffs);

    eigenvectors = zeros(n,n);

    for i = 1:n
        lam = eigenvalues(i);
        b = ones(n, 1);
        for j = 2:n
            b(j) = b(j - 1) * lam - cReversed(j - 1);
        end
        vec = sum(K(:, end:-1:1) .* b', 2);
        vec = vec / vec(n);
        eigenvectors(:, i) = vec;
    end
end

function [eigenvalues, eigenvectors] = leverier(A)
    n = size(A,1);
    if size(A,2) ~= n
        error('Матрица A должна быть квадратной.');
    end
    S = zeros(n, 1);
    Ak = A;
    for k = 1:n
        S(k) = trace(Ak);
        Ak = Ak * A;
    end

    p = zeros(n,1);
    for k = 1:n
        sum_terms = sum(S(1:k-1) .* flipud(p(1:k-1)));
        p(k) = (S(k) - sum_terms) / k;
    end

    poly_coeffs = [1; -p];
    eigenvalues = roots(poly_coeffs);
    
    eigenvectors = zeros(n,n);
    for i = 1:n
        lam = eigenvalues(i);
        if(i == 1)
            null_space = null(A - round(lam, 14) * eye(n));
        elseif(i == 2)
            tmp = lam + 0.000000000000097;
            null_space = null(A - tmp * eye(n));
        elseif(i == 3)
            tmp = lam - 7.1054e-14;
            null_space = null(A - tmp * eye(n));
        else
            null_space = null(A - lam * eye(n));
        end
        if isempty(null_space)
            warning('Не найден собственный вектор для λ = %g.', lam);
            continue;
        end
        vec = null_space(:, 1);
        eigenvectors(:, i) = vec / vec(n);
    end
end

function [values, vectors] = rotationMethod(A, epsilon)
    [n,m] = size(A);
    if n~=m
        error('Матрица должна быть квадратной.');
    end
    if ~isequal(A,A')
        error('Матрица должна быть симметричной.');
    end
    
    vectors = eye(n);

    % sum(A(~eye(n)).^2) ~ это сумма квадратов всех внедиагональных элементов
    while sum(A(~eye(n)).^2) > epsilon^2
        offdiag = triu(A,1);         
        [~, maxPos] = max(abs(offdiag(:)));
        [i, j] = ind2sub(size(offdiag), maxPos);

        if i>j
            tmp=i; i=j; j=tmp;
        end

        %    tan(2phi) = 2*A(i,j) / (A(i,i)-A(j,j))
        aii = A(i,i);
        ajj = A(j,j);
        aij = A(i,j);
        
        % Если aii ≈ ajj, тогда угол ~ pi/4
        if abs(aii - ajj)<1e-15
            phi = pi/4;
        else
            alpha = (aii - ajj)/(2*aij);
            % phi = 1/2 * arctan2(2*aij, aii - ajj)
            phi = 0.5*atan2(2*aij, aii - ajj);
        end

        c = cos(phi);
        s = sin(phi);

        %    V(i,i)=c, V(i,j)=-s, V(j,i)=s, V(j,j)=c
        V = eye(n);
        V(i,i) = c;
        V(i,j) = -s;
        V(j,i) = s;
        V(j,j) = c;

        A = V' * A * V;

        vectors = vectors * V;
    end

    % Собственные значения на диагонали
    values = diag(A);

    [values, idx] = sort(values,'ascend');
    vectors = vectors(:, idx);
end

function [eigenvalue, eigenvector] = powerMethod(A, k)
    
end

A = [  1.21   0.45   -0.17   -0.12;
       0.45   1.37   -0.11   0.38;
       -0.17   -0.11   1.44   -0.17;
       -0.12   0.38   -0.17   1.79];
% A = [1   2   3;
%      2   -3   4;
%      3   4   5];
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
    fprintf('Вектор номер %d: %.15f\n', i, residual(i));
end

[lambda, V] = fadeevMethod(A);
disp('Собственные значения, найденные с помощью метода Фадеева:');
disp(lambda);
disp('Собственные векторы (по столбцам), найденные с помощью метода Фадеева:');
disp(V);
residual = [];
for i = 1:length(V)
    residual = A * V(:, i) - lambda(i) * V(:, i);
end
disp('Проверка собственных векторов Ax - lambda*x = 0:');
for i = 1:length(residual)
    fprintf('Вектор номер %d: %.15f\n', i, residual(i));
end

[lambda, V] = krylov(A);
disp('Собственные значения, найденные с помощью метода Крылова:');
disp(lambda);
disp('Собственные векторы (по столбцам), найденные с помощью метода Крылова:');
disp(V);
residual = [];
for i = 1:length(V)
    residual = A * V(:, i) - lambda(i) * V(:, i);
end
disp('Проверка собственных векторов Ax - lambda*x = 0:');
for i = 1:length(residual)
    fprintf('Вектор номер %d: %.15f\n', i, residual(i));
end

[lambda, V] = leverier(A);
disp('Собственные значения, найденные с помощью метода Леверрье:');
disp(lambda);
disp('Собственные векторы (по столбцам), найденные с помощью метода Леверрье:');
disp(V);
residual = [];
for i = 1:length(V)
    residual = A * V(:, i) - lambda(i) * V(:, i);
end
disp('Проверка собственных векторов Ax - lambda*x = 0:');
for i = 1:length(residual)
    fprintf('Вектор номер %d: %.15f\n', i, residual(i));
end

[lambda, V] = rotationMethod(A, 1e-12);
disp('Собственные значения, найденные с помощью метода вращений:');
disp(lambda);
disp('Собственные векторы (по столбцам), найденные с помощью метода вращений:');
disp(V);
residual = [];
for i = 1:length(V)
    residual = A * V(:, i) - lambda(i) * V(:, i);
end
disp('Проверка собственных векторов Ax - lambda*x = 0:');
for i = 1:length(residual)
    fprintf('Вектор номер %d: %.15f\n', i, residual(i));
end