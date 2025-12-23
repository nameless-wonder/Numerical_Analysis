
clc;
clear;
close all;
lab2_exercise1(1, 0.1, 10); % Call the function for y0 = 1
lab2_exercise1(0.99, 0.1, 10); % Call the function for y0 = 0.99
lab2_exercise1(1.01, 0.1, 10); % Call the function for y0 = 1.01
function[] = lab2_exercise1(y0,h,N)
    yE = zeros(1,N+1);
    yBE = zeros(1,N+1);

    t = 0:h:N*h;
    exact_solution = @(t) t + exp(-100*t) + 1;
    yE(1) = y0;
    yBE(1) = y0;
    fun = @(t,y) -100*y + 100*t + 101;
    for n = 1:N
        yE(n+1) = yE(n) + h*fun(t(n),yE(n));
    end

    for n=1:N
        yBE(n+1) = (yBE(n) + h * (100*t(n+1) + 101)) / (1 + 100*h);
    end


    figure;
    plot(t, exact_solution(t), 'k-', 'LineWidth', 1.5); hold on;
    plot(t, yE, 'b--o', 'LineWidth', 1.2);
    title('Forward Euler Method');
    xlabel('Time t');
    ylabel('Solution y');
    legend('Exact Solution', 'Forward Euler Approximation');
    grid on;

    % Plot Backward Euler solution and exact solution
    figure;
    plot(t, exact_solution(t), 'k-', 'LineWidth', 1.5); hold on;
    plot(t, yBE, 'r--o', 'LineWidth', 1.2);
    title('Backward Euler Method');
    xlabel('Time t');
    ylabel('Solution y');
    legend('Exact Solution', 'Backward Euler Approximation');
    grid on;

    % Display results
    disp('Time steps t:');
    disp(t);
    disp('Forward Euler solution yE:');
    disp(yE);
    disp('Backward Euler solution yBE:');
    disp(yBE);



end
