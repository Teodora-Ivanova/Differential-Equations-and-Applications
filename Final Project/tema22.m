function task22
% The case where a = 0.6
clf 
clc

% Define the variables given in the task
L = 10;
a = 0.6;
x = 0:L / 100:L;
t(1) = 0;
t(2) = 2;
t(3) = 5;

% Plot the movement of the string in the given moments t(1) = 0, t(2) = 2, t(3) = 5
for i = 1 : length(t)
    % Plot the graphics in the moment t(0) = 0
    subplot(3, 1, 1)
    plot(x, fourier(x, t(1)), 'r', 'LineWidth', 2)
    title('t = 0', 'Color', 'r')
    % Constraints for y: -11*pi and 11*pi are chosen 
    % in order to get better view of the graphics
    axis([0 L -11*pi 11*pi]) 
    grid on
    grid minor
    
    % Plot the graphics in the moment t(1) = 2
    subplot(3, 1, 2)
    plot(x, fourier(x, t(2)), 'b', 'LineWidth', 2)
    title('t = 2', 'Color', 'b')
    axis([0 L -11*pi 11*pi])
    grid on
    grid minor
   
    % Plot the graphics in the moment t(2) = 5
    subplot(3, 1, 3)
    plot(x, fourier(x, t(3)), 'm', 'LineWidth', 2)
    title('t = 5', 'Color', 'm')
    axis([0 L -11*pi 11*pi])
    grid on
    grid minor
    
    M(i) = getframe;
end
movie(M, 3);

% Define Fourier function
function y = fourier(x, t)
    y = 0;
    for k = 1:35 % 35th partial sum
        % Define Xk
        Xk = sin(((2*k+1)*pi/(2*L)).*x);
        % Define Ak
        Ak = (2/L)*trapz(x, phi(x).*Xk);
        % Define Bk
        Bk = (4/(2*k+1)*pi*a)*trapz(x, psi(x).*Xk);
        % Define Tk
        Tk = Ak*cos((2*k+1)*pi*a*t/(2*L)) + Bk*sin((2*k+1)*pi*a*t/(2*L));
        y = y + Xk.*Tk;
    end
end

% Define function phi(x)
function y = phi(x)
    for j = 1:length(x)
        if pi <= x(j) && x(j) <= 2*pi
            y(j) = -5*((x(j)-pi)^2)*((x(j)-2*pi)^2);
        elseif (0 <= x(j) && x(j) < pi) || (2*pi < x(j) && x(j) <= L)
            y(j) = 0;
        end
    end
end

% Define function psi(x)
function y = psi(x)
    y = 0;
end

end