clc;
clear;
close all;

lab1_exercise2(1000,25,1,0.5,0.1,0.02,100,10);

function [] = lab1_exercise2(N,T,a1,a2,b1,b2,y10,y20)
    y = zeros(2,N+1);
    y(:,1) = [y10;y20];

    h = T/N;

    for n = 1:N
        y(:,n+1) = y(:,n) + [y(1,n)*(a1 - b1*y(2,n));y(2,n)*(-a2+b2*y(1,n))]*h;
    end
    timer = 0:h:T;
    plot(timer,y);
    % Plot populations y1 and y2 as functions of time
    figure;
    plot(timer, y(1, :), 'r', 'LineWidth', 1.5); hold on;
    plot(timer, y(2, :), 'b', 'LineWidth', 1.5);
    title('Population Dynamics Over Time');
    xlabel('Time');
    ylabel('Population');
    legend('Prey (y_1)', 'Predator (y_2)');
    grid on;

    % Plot the trajectory of (y1(t), y2(t)) in the phase plane
    figure;
    plot(y(1, :), y(2, :), 'k', 'LineWidth', 1.5);
    title('Phase Plane Trajectory');
    xlabel('Prey Population (y_1)');
    ylabel('Predator Population (y_2)');
    grid on;


end