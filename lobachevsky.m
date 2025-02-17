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
    numbers_columns = 0;
    is_complex = 0; 
    for column = 2:3
        for row = 2:(degree_accuracy + 1)
            if(quadr_matrix(row, column) < 0)
                numbers_columns = column;
                is_complex = 1;
                break;
            end
        end
        if(is_complex == 1)
            break;
        end
    end
    % вычисляем корни
    m = 2 ^ (degree_accuracy);
    if(numbers_columns == 0)%все корни вещественные
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
        if(numbers_columns == 2)
            x3 = (abs((quadr_matrix(degree_accuracy + 1, 4))/(quadr_matrix(degree_accuracy + 1, 3)))) ^ (1/m);
            if(abs(subs(poly, x3)) > abs(subs(poly, -x3)))
                x3 = -x3;
            end
            r = (abs(quadr_matrix(degree_accuracy + 1, 4) / quadr_matrix(degree_accuracy + 1, 2))) ^ (1 / (2 * m));
            cosinus = (-(coeffs_vector(2)/coeffs_vector(1)) + (-x3)) / 2 * r;
            sinus = sqrt(1 - (cosinus * cosinus));
            x1 = r * cosinus + (r * sinus)*1i;
            x2 = r * cosinus - (r * sinus)*1i;
        else
            x1 = (abs((quadr_matrix(degree_accuracy + 1, 2))/(quadr_matrix(degree_accuracy + 1, 1)))) ^ (1/m);
            if(abs(subs(poly, x1)) > abs(subs(poly, -x1)))
                x1 = -x1;
            end
            r = (abs(quadr_matrix(degree_accuracy + 1, 4) / quadr_matrix(degree_accuracy + 1, 2))) ^ (1 / (2 * m));
            tmp = -(quadr_matrix(1, 2)/quadr_matrix(1, 1));
            tmp2 = tmp - x1;
            tmp3 = 2 * r;
            cosinus = tmp2 / tmp3;
            sinus = sqrt(1 - (cosinus * cosinus));
            x2 = r * cosinus + (r * sinus)*1i;
            x3 = r * cosinus - (r * sinus)*1i;
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
    poly = 2*x^3 - 3*x^2 - x - 1.5;
    coefs = [2 -3 -1 -1.5];

    % poly = x^3 + 3*x^2 - x - 1;
    % coefs = [1 3 -1 -1];
    roots_values = roots(coefs);
    fprintf('Точное решение:\n');
    disp(roots_values);

    % [x1, x2, x3] = lobachevckyMethod(poly, 5);
    % fprintf('Решение методом Лобачевского c числом шагов k=%d:\n', 5);
    % fprintf('x1 = %.7f\n', x1);
    % if(imag(x2) >= 0)
    %     fprintf('x2 = %.7f + %.7fi\n', real(x2), imag(x2));
    % else
    %     fprintf('x2 = %.7f - %.7fi\n', real(x2), abs(imag(x2)));
    % end
    % if(imag(x3) >= 0)
    %     fprintf('x3 = %.7f + %.7fi\n', real(x3), imag(x3));
    % else
    %     fprintf('x3 = %.7f - %.7fi\n', real(x3), abs(imag(x3)));
    % end

    [x1, x2, x3] = lobachevckyMethod(poly, 9);
    fprintf('Решение методом Лобачевского c числом шагов k=%d:\n', 9);
    fprintf('x1 = %.7f\n', x1);
    if(imag(x2) >= 0)
        fprintf('x2 = %.7f + %.7fi\n', real(x2), imag(x2));
    else
        fprintf('x2 = %.7f - %.7fi\n', real(x2), abs(imag(x2)));
    end
    if(imag(x3) >= 0)
        fprintf('x3 = %.7f + %.7fi\n', real(x3), imag(x3));
    else
        fprintf('x3 = %.7f - %.7fi\n', real(x3), abs(imag(x3)));
    end

    polinom = @(x) 2*x^3 - 3*x^2 - x - 1.5;
    result_poly = polinom(x3);
    fprintf('\nresult = %.10f', result_poly);
end
