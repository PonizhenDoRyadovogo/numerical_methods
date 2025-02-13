function roots = lobachevckyMethod(poly, degree_accuracy)
    coeffs_vector = coeffs(poly, 'ALL');
    degree_poly = length(coeffs_vector) - 1;
    quadr_matrix = zeros(degree_accuracy, length(coeffs_vector));
    roots = zeros(1, degree_poly);
    %заполнили первую строку в матрице квадрирования
    for column = 1:length(coeffs_vector)
        quadr_matrix(1,column) = coeffs_vector(column);
    end
    %заполняем матрцу квадрирования
    for row = 2:degree_accuracy % отдельно обрабатываем крайние коэффициенты
        quadr_matrix(row, 1) = quadr_matrix(row - 1, 1) * quadr_matrix(row - 1, 1);
        quadr_matrix(row, degree_poly + 1) = quadr_matrix(row - 1, degree_poly + 1) * quadr_matrix(row - 1, degree_poly + 1);
    end
    for row = 2:degree_accuracy
        for column = 2:degree_poly
            a = 0;
            k = 1;
            while((column - k >= 0) && (column + k)<= degree_poly)
                k = k + 1;
            end
            k = k - 1;
            for p = 1:k
                a = a + ((-1) ^ k) * 2 * quadr_matrix(row - 1, column - p) * quadr_matrix(row - 1, column + p);
            end
            a = (quadr_matrix(row - 1, column) * quadr_matrix(row - 1, column)) + a;
            quadr_matrix(row, column) = a;
        end
    end
    %просматриваем матрицу квадрирования, чтобы понять какие у корни у
    %нашего уравнения
    numbers_columns = [];
    is_complex = 0; 
    for column = 2:degree_poly
        for row = 2:degree_accuracy
            if(is_complex == 1)
                is_complex = 0;
                break;
            end
            if(quadr_matrix(row, column) < 0)
                numbers_columns(end+1) = column - 1;
                is_complex = 1;
                break;
            end
        end
    end
    % вычисляем корни
    m = 2 ^ (degree_accuracy);
    for k = 1:degree_poly
        if(ismember(k, numbers_columns)) % комплексно-сопряженные корни
            r = (quadr_matrix(degree_accuracy, 4)/quadr_matrix(degree_accuracy, 1)) ^ (1/(2 * m));
            %roots(1)
        else % только вещественные корни
            roots(k) = (abs((quadr_matrix(degree_accuracy, k))/(quadr_matrix(degree_accuracy, k - 1)))) ^ (1/m);
            if(subs(poly, roots(k)) > subs(poly, -roots(k)))
                roots(k) = -roots(k);
            end
        end
    end
end

flag = 1;
if(flag == 0)
    syms x;
    poly_str = input('Введите многочлен (например, 2*x^3 - 3*x^2 - x - 1.5): ', 's');
    poly_expr = str2sym(poly_str);
    disp('Ваш многочлен:');
    disp(poly_expr);
else
    syms x;
    % poly = 2*x^3 - 3*x^2 - x - 1.5;
    poly = x^3 + 3*x^2 - x - 1;
    lobachevckyMethod(poly, 5);
end
