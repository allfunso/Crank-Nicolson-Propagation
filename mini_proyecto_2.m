% Propagación de ondas en espacio libre usando
% el método Crank-Nicolson
% Alfonso Castro Camino - A01733148

close all
clear

%% Método Crank-Nicolson

% Constantes del problema
lambda = 632.8e-9; % m
rangox = [-6e-3, 6e-3]; % m
nx = 511; % Impar para asegurar tener un punto en x=0
rangoz = [0, 1.5]; % m
nz = 200;

% Generar vectores de dominio en tiempo y propagación
x_ = linspace(rangox(1), rangox(2), nx);
z_ = linspace(rangoz(1), rangoz(2), nz);
dx = x_(2) - x_(1);
dz = z_(2) - z_(1);

% Matrices del método
k = 2 * pi / lambda; % Coeficiente difusión
kappa = - 1 / (2i * k);
alpha = kappa * dz / dx^2;

C = 2*(1+alpha)*eye(nx) + diag(-alpha*ones(nx-1, 1), 1) + diag(-alpha*ones(nx-1, 1), -1);
D = 2*(1-alpha)*eye(nx) + diag(alpha*ones(nx-1, 1), 1) + diag(alpha*ones(nx-1, 1), -1);

P = C \ D;
% f(n+1) = P * f(n)


%% 1.Experimento de Young

% Pantalla (onda incidente tiene amplitud unitaria)
h = 1e-3; % m
b = 0.2e-3; % m

% Condición inicial del campo
A = 1; % Amplitud
periodo = 2*b; % Periodo función de suavización
suav = sin(2*pi*x_/periodo).^2; % Función de suavización
U0 = A*double(abs(abs(x_)-h/2) <= b/2); % Función booleana pantalla incidente
U0 = U0 .* suav; % Campo inicial suavizado

figure(4)
plot(x_, U0)
hold on
%plot(x_, suav, ":")
hold off
title("Perfil de onda inicial")
xlabel("x (m)")

U = zeros(nx, nz);
U(:, 1) = U0;
Energ = zeros(1, nz);
Energ(1) = sum(abs(U0).^2 * dx);

for i = 1:nz-1
    U(:, i+1) = P*U(:, i);
    Energ(i+1) = sum(abs(U(:, i+1)).^2 * dx);
end

% Inciso a)
figure(1)
imagesc(abs(U))
title("Propagación")

% Inciso b)
figure(2)
E0 = Energ(1);
ERP = 100 * abs((E0 - Energ) / E0);
plot(z_, ERP)
title("Error relativo porcentual");
xlabel("z (m)")
ylabel("Error (%)")

% Inciso c)
figure(3)
distancias = 0:0.25:1.5;
hold on
for idx_z = floor(distancias/dz) + 1
    plot(x_, abs(U(:, idx_z)))
end
legend(string(distancias))


% Inciso d)
% Comprobar con método analítico
beta = 1/2 * k * b * sin(atan(x_/rangoz(2)));
gamma = 1/2 * k * h * sin(atan(x_/rangoz(2)));
I_f = (sin(beta)./beta).^2 .* cos(gamma);
%plot(x_, 0.3*abs(I_f))
delt_x = rangoz(2)*tan(lambda/h);
%xline(-3*delt_x:delt_x:3*delt_x, "LineStyle", "--", "Color", "r")
hold off

%% 2. Lente convergente a mitad de distancia
U = zeros(nx, nz);
U(:, 1) = U0;
Energ = zeros(1, nz);
Energ(1) = sum(abs(U0).^2 * dx);
nz_div = round(0.75/dz);

f = 0.75;
TL = exp(-1i * k / (2*f) * x_.^2);

for i = 1:nz-1%_div-1
    U(:, i+1) = P*U(:, i);
    if i == nz_div
        U(:, i+1) = U(:, i+1) .* TL.';
    end
    Energ(i+1) = sum(abs(U(:, i+1)).^2 * dx);
end

% Inciso a)
figure(4)
imagesc(abs(U))
title("Propagación")

% Inciso b)
figure(5)
E0 = Energ(1);
ERP = 100 * abs((E0 - Energ) / E0);
plot(z_, ERP)
title("Error relativo porcentual");
xlabel("z (m)")
ylabel("Error (%)")

% Inciso c)
figure(6)
distancias = 0:0.25:1.5;
hold on
for idx_z = floor(distancias/dz) + 1
    plot(x_, abs(U(:, idx_z)))
end
legend(string(distancias))
hold off
%% 3. Ejercicio libre

% Pantalla (onda incidente tiene amplitud 1)
h = 2e-3; % m
b = 0.4e-3; % m

% Condición inicial del campo
A = 1; % Amplitud
periodo = 2*b; % Periodo función de suavización
suav = cos(2*pi*x_/(0.5*h)).^2; % Función de suavización
asym = -0.33/h*x_ + 0.55; % Función de asimetrización
U0 = A*double(abs(abs(x_)-h/2) <= b/2); % Función booleana pantalla incidente
U0 = U0 .* asym .* suav; % Campo inicial con funciones aplicadas
U0(ceil(nx/2):ceil(nx/2)+4) = [0.25, 0.5, 0.55, 0.5, 0.25]; % Pequeña rendija en medio, ligeramente a la derecha

figure(7)
plot(x_, U0)
hold on
%plot(x_, suav, ":")
%plot(x_, asym, "--")
hold off
title("Perfil de onda inicial")
xlabel("x (m)")

U = zeros(nx, nz);
U(:, 1) = U0;
Energ = zeros(1, nz);
Energ(1) = sum(abs(U0).^2 * dx);
nz_div = round(0.9/dz);

f = 0.9;
TL = exp(-1i * k / (2*f) * x_.^2);

for i = 1:nz-1
    U(:, i+1) = P*U(:, i);
    if i == nz_div
        U(:, i+1) = U(:, i+1) .* TL.';
    end
    Energ(i+1) = sum(abs(U(:, i+1)).^2 * dx);
end

% Inciso a)
figure(8)
imagesc(abs(U))
title("Propagación")

% Inciso b)
figure(9)
E0 = Energ(1);
ERP = 100 * abs((E0 - Energ) / E0);
plot(z_, ERP)
title("Error relativo porcentual");
xlabel("z (m)")
ylabel("Error (%)")

% Inciso c)
figure(10)
distancias = [0, 0.5, 1, 1.25, 1.5];
hold on
for idx_z = floor(distancias/dz) + 1
    plot(x_, abs(U(:, idx_z)))
end
legend(string(distancias))
hold off
