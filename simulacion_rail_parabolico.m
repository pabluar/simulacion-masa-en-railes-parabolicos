clear;clc;
close all;
%SIMULACIÓN MOI. MASA EN RAILES PARABOLICOS CON Y SIN GIRO
%PABLO ARIAS RODRIGUEZ Y SARA MOHAMED ACEVEDO
% Datos del problema

m = input('Masa [kg] (recomendado para visualizar 1kg): ');
a = input('Parametro a de la parabola [1/m] (recomendado 1m^-1): ');

s0 = 1.4;       % Posicion horizontal inicial [m]
v0 = 0;         % Velocidad inicial plano xy ds/dt [m/s]
tf = 30;        % Tiempo de simulacion [s]
g = 9.81;       % Gravedad [m/s^2]
phi0 = 0;       % Angulo inicial del rail [rad]

disp('1. Rail fijo');
disp('2. Rail con giro libre');
modo = input('Elige el caso (1 o 2): ');

if m <= 0 || a <= 0
    error('La masa y a deben ser positivos.');
end

if modo == 1
    omega0 = 0;
elseif modo == 2
    omega0 = input('Velocidad angular inicial [rad/s] (recomendado 1rad/s): ');
else
    error('Debes elegir 1 o 2.');
end

% Momento angular inicial, que se conserva
ell = m*s0^2*omega0;

% Resolucion de las ecuaciones

% El vector u contiene:
% u(1) = s
% u(2) = ds/dt
% u(3) = phi

u0 = [s0; v0; phi0];

% Instantes en los que se obtiene la solucion
tiempos = linspace(0, tf, ceil(30*tf) + 1);

opciones = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

[t, u] = ode45(@(t,u) ecuaciones(t,u,m,a,g,ell), ...
              tiempos, u0, opciones);

s = u(:,1);
v = u(:,2);
phi = u(:,3);

% Velocidad angular obtenida del momento angular
if ell == 0
    omega = zeros(size(s));
else
    omega = ell./(m*s.^2);
end

% Coordenadas cartesianas

% s puede ser positiva o negativa: representa las dos ramas
x = s.*cos(phi);
y = s.*sin(phi);
z = a*s.^2;

% Energia y momento angular

Ec = 0.5*m*((1 + 4*a^2*s.^2).*v.^2 + s.^2.*omega.^2);

Ep = m*g*z;
E = Ec + Ep;

Lz = m*s.^2.*omega;

% Comprobacion de la conservacion de la energia
errorE = max(abs(E - E(1)));

fprintf('\nEnergia inicial: %.6f J\n', E(1));
fprintf('Desviacion maxima de la energia: %.3e J\n', errorE);

% Graficas

figure;

subplot(2,2,1);
plot(t, s, 'b', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('s [m]');
title('Posicion horizontal sobre el rail');
grid on;

subplot(2,2,2);
plot(t, omega, 'm', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('\omega [rad/s]');
title('Velocidad angular');
grid on;

subplot(2,2,3);
plot(t, Ec, 'b', t, Ep, 'r', t, E, 'k', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('Energia [J]');
legend('Cinetica', 'Potencial', 'Total');
title('Energia mecanica');
grid on;

subplot(2,2,4);
plot(t, Lz, 'g', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('L_z [kg m^2/s]');
title('Momento angular vertical');
grid on;

% Reaccion normal para el caso sin giro

if omega0 == 0
    N = m*(g + 2*a*v.^2)./sqrt(1 + 4*a^2*s.^2);

    figure;
    plot(t, N, 'b', 'LineWidth', 1.5);
    xlabel('Tiempo [s]');
    ylabel('N [N]');
    title('Reaccion normal del rail');
    grid on;
end

% Animacion del movimiento

figure;

% Puntos utilizados para dibujar el rail completo
R = max(0.5, 1.1*max(abs(s)));
sRail = linspace(-R, R, 200);
zRail = a*sRail.^2;

rail = plot3(sRail*cos(phi(1)), ...
             sRail*sin(phi(1)), zRail, 'b', 'LineWidth', 2);

hold on;

masa = plot3(x(1), y(1), z(1), 'ro', ...
             'MarkerFaceColor', 'r', 'MarkerSize', 10);

traza = plot3(x(1), y(1), z(1), 'r--');

% Eje vertical de giro
plot3([0 0], [0 0], [0 1.1*max(zRail)], 'k--');

xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
grid on;
axis equal;
view(40, 25);

xlim([-R R]);
ylim([-R R]);
zlim([0 1.1*max(zRail)]);

for k = 1:length(t)

    % Si se cierra la ventana, termina la animacion
    if ~isgraphics(rail)
        break;
    end

    % El rail gira un angulo phi
    set(rail, 'XData', sRail*cos(phi(k)), ...
              'YData', sRail*sin(phi(k)));

    % Posicion del cuerpo
    set(masa, 'XData', x(k), ...
              'YData', y(k), ...
              'ZData', z(k));

    % Trayectoria recorrida
    set(traza, 'XData', x(1:k), ...
               'YData', y(1:k), ...
               'ZData', z(1:k));

    title(sprintf('Movimiento del cuerpo: t = %.2f s', t(k)));

    drawnow;

    if k < length(t)
        pause(t(k+1) - t(k));
    end
end

% Funcion con las ecuaciones del movimiento

function du = ecuaciones(~, u, m, a, g, ell)

    s = u(1);
    v = u(2);

    % Si no hay momento angular, el rail no gira
    if ell == 0
        omega = 0;
    else
        omega = ell/(m*s^2);
    end

    % Aceleracion horizontal sobre el rail
    aceleracion = (s*omega^2 - 4*a^2*s*v^2 - 2*g*a*s) ...
                  /(1 + 4*a^2*s^2);

    du = [v; aceleracion; omega];

end