function [x1, x2, x3] = lobachevckyMethod(poly, degree_accuracy)
    coeffs_vector = coeffs(poly, 'ALL');
    quadr_matrix = zeros(degree_accuracy + 1, length(coeffs_vector));
    
    %заполнили первую строку в матрице квадрирования
    for column = 1:length(coeffs_vector)
        quadr_matrix(1,column) = coeffs_vector(column);
    end
    %заполняем матрцу квадрирования
    for row = 2:(degree_accuracy + 1) % отдельно обрабатываем крайние коэффициенты
        quadr_matrix(row, 1) = quadr_matrix(row - 1, 1) * quadr_matrix(row - 1, 1);
        quadr_matrix(row, 4) = quadr_matrix(row - 1, 4) * quadr_matrix(row - 1, 4);
    end
    for row = 2:(degree_accuracy + 1)
        for column = 2:3
            quadr_matrix(row, column) = quadr_matrix(row - 1, column) * quadr_matrix(row - 1, column) - 2 * quadr_matrix(row - 1, column - 1) * quadr_matrix(row - 1, column + 1);
        end
    end
    %просматриваем матрицу квадрирования, чтобы понять какие у корни у
    %нашего уравнения
    numbers_columns = [];
    is_complex = 0; 
    for column = 2:3
        for row = 2:(degree_accuracy + 1)
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
    if(isempty(numbers_columns))
        x1 = (abs((quadr_matrix(degree_accuracy + 1, 2))/(quadr_matrix(degree_accuracy + 1, 1)))) ^ (1/m);
        if(abs(subs(poly, x1)) > abs(subs(poly, -x1)))
            x1 = -x1;
        end
        x2 = (abs((quadr_matrix(degree_accuracy + 1, 3))/(quadr_matrix(degree_accuracy + 1, 2)))) ^ (1/m);
        if(abs(subs(poly, x2)) > abs(subs(poly, -x2)))
            x2 = -x2;
        end
        x3 = (abs((quadr_matrix(degree_accuracy + 1, 4))/(quadr_matrix(degree_accuracy + 1, 3)))) ^ (1/m);
        if(abs(subs(poly, x3)) > abs(subs(poly, -x3)))
            x3 = -x3;
        end
    else
        
    end
    % for k = 1:3
    %     if(ismember(k, numbers_columns)) % комплексно-сопряженные корни
    %         r = (quadr_matrix(degree_accuracy, 4)/quadr_matrix(degree_accuracy, 1)) ^ (1/(2 * m));
    %         %roots(1)
    %     else % только вещественные корни
    %         x1 = (abs((quadr_matrix(degree_accuracy, k))/(quadr_matrix(degree_accuracy, k - 1)))) ^ (1/m);
    %         if(subs(poly, roots(k)) > subs(poly, -roots(k)))
    %             roots(k) = -roots(k);
    %         end
    %     end
    % end
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
