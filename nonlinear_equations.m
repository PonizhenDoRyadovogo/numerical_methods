function [result, count_iteration, array_x] = simpleIteration(func, initial_x, epsilon)
    current_x = func(initial_x);
    array_x = [];
    array_x(end+1) = current_x;
    next_x = func(current_x);
    array_x(end+1) = next_x;
    count_iteration = 1;
    while(abs(next_x - current_x) >= epsilon)
        current_x = next_x;
        next_x = func(current_x);
        count_iteration = count_iteration + 1;
        array_x(end+1) = next_x;
    end
    result = next_x;
end

function [result, count_iteration] = secantMethod(func, initial_x, epsilon)
    current_x = func(initial_x);
    next_x = (initial_x * func(current_x) - current_x * func(initial_x))/(func(current_x) - current_x - func(initial_x) + initial_x);
    count_iteration = 1;
    while(abs(next_x - current_x) >= epsilon)
        tmp_x = current_x;
        current_x = next_x;
        next_x = (tmp_x * func(current_x) - current_x * func(tmp_x))/(func(current_x) - current_x - func(tmp_x) + tmp_x);
        count_iteration = count_iteration + 1;
    end
    result = next_x;
end

function [result, count_iteration] = steffensen(func, initial_x, epsilon)
    current_x = initial_x;
    next_x = (current_x * func(func(current_x)) - (func(current_x) * func(current_x)))/(func(func(current_x)) - 2 * func(current_x) + current_x);
    count_iteration = 1;
    while(abs(next_x - current_x) >= epsilon)
        current_x = next_x;
        next_x = (current_x * func(func(current_x)) - (func(current_x) * func(current_x)))/(func(func(current_x)) - 2 * func(current_x) + current_x);
        count_iteration = count_iteration + 1;
    end
    result = next_x;
end

function [result, count_iteration] = newton(func, deriative, initial_x, epsilon)
    current_x = initial_x - (func(initial_x))/(deriative(initial_x));
    next_x = current_x - (func(current_x))/(deriative(current_x));
    count_iteration = 1;
    while(abs(next_x - current_x) >= epsilon)
        current_x = next_x;
        next_x = current_x - (func(current_x))/(deriative(current_x));
        count_iteration = count_iteration + 1;
    end
    result = next_x;
end

function [result, count_iteration] = secantMethodNewton(func, initial_x, epsilon)
    current_x = func(initial_x);
    next_x = (initial_x * func(current_x) - current_x * func(initial_x))/(func(current_x) - func(initial_x));
    count_iteration = 1;
    while(abs(next_x - current_x) >= epsilon)
        tmp_x = current_x;
        current_x = next_x;
        next_x = (tmp_x * func(current_x) - current_x * func(tmp_x))/(func(current_x) - func(tmp_x));
        count_iteration = count_iteration + 1;
    end
    result = next_x;
end

function [result, count_iteration] = newtonConstDeriative(func, der, initial_x, epsilon)
    current_x = initial_x - (func(initial_x))/(der);
    next_x = current_x - (func(current_x))/(der);
    count_iteration = 1;
    while(abs(next_x - current_x) >= epsilon)
        current_x = next_x;
        next_x = current_x - (func(current_x))/(der);
        count_iteration = count_iteration + 1;
    end
    result = next_x;
end

phi = @(x) exp((3 * cos(x))/10);
initial_x = 1.136;
epsilon = 0.5 * 10 ^ (-4);
func = @(x) 10 * log(x) - 3 * cos(x);
der_func = @(x) (10/x) + 3 * sin(x);

der = der_func(initial_x);

[origin_x, count_iteration, array_x] = simpleIteration(phi, initial_x, epsilon);
fprintf('Метод простых итераций:\n');
fprintf('при точке начального приближения x0 = %.3f значение x = %.5f, количество итераций: %d', initial_x, origin_x, count_iteration);
[origin_x, count_iteration, array_x2] = simpleIteration(phi, 0, epsilon);
fprintf('\nпри точке начального приближения x0 = 0 значение x = %.5f, количество итераций: %d', origin_x, count_iteration);

[origin_x, count_iteration] = secantMethod(phi, initial_x, epsilon);
fprintf('\nМетод секущих:\n');
fprintf('при точке начального приближения x0 = %.3f значение x = %.5f, количество итераций: %d', initial_x, origin_x, count_iteration);
[origin_x, count_iteration] = secantMethod(phi, 0, epsilon);
fprintf('\nпри точке начального приближения x0 = 0 значение x = %.5f, количество итераций: %d', origin_x, count_iteration);

[origin_x, count_iteration] = steffensen(phi, initial_x, epsilon);
fprintf('\nМетод Стеффенсена:\n');
fprintf('при точке начального приближения x0 = %.3f значение x = %.5f, количество итераций: %d', initial_x, origin_x, count_iteration);
[origin_x, count_iteration] = steffensen(phi, 0, epsilon);
fprintf('\nпри точке начального приближения x0 = 0 значение x = %.5f, количество итераций: %d', origin_x, count_iteration);

[origin_x, count_iteration] = newton(func, der_func, initial_x, epsilon);
fprintf('\nМетод Ньютона:\n');
fprintf('при точке начального приближения x0 = %.3f значение x = %.5f, количество итераций: %d', initial_x, origin_x, count_iteration);
[origin_x, count_iteration] = newton(func, der_func, 0.5, epsilon);
fprintf('\nпри точке начального приближения x0 = 0.5 значение x = %.5f, количество итераций: %d', origin_x, count_iteration);

[origin_x, count_iteration] = secantMethodNewton(func, initial_x, epsilon);
fprintf('\nМетод секущих Ньютона:\n');
fprintf('при точке начального приближения x0 = %.3f значение x = %.5f, количество итераций: %d', initial_x, origin_x, count_iteration);
[origin_x, count_iteration] = secantMethodNewton(func, 0.75, epsilon);
fprintf('\nпри точке начального приближения x0 = 0.75 значение x = %.5f, количество итераций: %d', origin_x, count_iteration);

[origin_x, count_iteration] = newtonConstDeriative(func, der, initial_x, epsilon);
fprintf('\nМетод Ньютона с постоянной производной:\n');
fprintf('при точке начального приближения x0 = %.3f значение x = %.5f, количество итераций: %d', initial_x, origin_x, count_iteration);
[origin_x, count_iteration] = newtonConstDeriative(func, der_func(0.5), initial_x, epsilon);
fprintf('\nпри точке начального приближения x0 = 0.75 значение x = %.5f, количество итераций: %d', origin_x, count_iteration);

%для показа графического нахождения точки начального приближения x0
x = linspace(0.1, 10, 10000); 
psi1 = 10 * log(x);  
psi2 = 3 * cos(x); 
origin_func = 10 * log(x) - 3 * cos(x);

subplot(2,1,1);
plot(x, psi1, 'r', 'LineWidth', 2);
hold on; 
plot(x, psi2, 'b', 'LineWidth', 2); 
hold off;
grid on;
legend('10ln(x)', '3cos(x)');
xlabel('x');
ylabel('y');
title('Графики функций psi1(x)=10ln(x) и psi2(x)=3cos(x)');

x = linspace(0.1, 10, 10000);
subplot(2,1,2);
plot(x, origin_func, 'g', 'LineWidth',2);
hold on;

%рисуем "след"
array_y = 10 * log(array_x) - 3 * cos(array_x);
plot(array_x, array_y, 'mo-', 'MarkerFaceColor', 'm', 'MarkerSize', 5, 'LineWidth', 3);
array_y2 = 10 * log(array_x2) - 3 * cos(array_x2);
plot(array_x2, array_y2, 'bo-', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'LineWidth', 1.5);
xline(0, 'k','LineWidth',1.5);
yline(0, 'k', 'LineWidth',1.5);
legend('10ln(x) - 3cos(x)','След метода простых итераций при точке начального приближения x0=1.136','след метода простых итераций при точке начального приближения x0=0','ось X', 'ось Y');
hold off;
grid on;
xlabel('x');
ylabel('y');
title('График функции f(x)=10ln(x)-3cos(x)');




