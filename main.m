clear; close all;

%%% Физические константы
pc = struct;                  % pc - problem constants, константы задачи
pc.lambda1 = 0.6;             % Коэффициент теплопроводности воды, Вт / (м * K)
pc.c1 = 4180.6;               % Коэффициеент удельной теплоёмкости воды, Дж / (кг * К)
pc.rho1 = 1000;               % Плотность воды, кг/м^3
pc.a1_sq = pc.lambda1/...     % Коэффициент температуропроводности воды, м^2/с
    pc.c1/pc.rho1;            
pc.lambda2 = 2.33;            % Коэффициент теплопроводности льда, Вт / (м * K)
pc.c2 = 2110.0;               % Коэффициеент удельной теплоёмкости льда, Дж / (кг * К)
pc.rho2 = 916.7;              % Плотность льда, кг/м^3
pc.a2_sq = pc.lambda2/...     % Коэффициент температуропроводности льда, м^2/с
    pc.c2/pc.rho2;            
pc.qf = 330*1e3;              % Удельная теплота плавления льда, Дж / кг
%pc.rho = (rho1 + rho2)/2;     % Средняя плотность
%pc.L = 1;                     % Длина стержня, м
pc.Uf = 273.15;               % Температура фазового перехода, К

%%% Смешанные краевые условия
% Формат краевых условий:
% alpha00*u1(0, t) + alpha01*du1/dx(0, t) = g0(t) - для левого конца
% alpha10*u2(L, t) + alpha11*du2/dx(L, t) = g1(t) - для правого конца
bc = struct;                  % bc - boundary conditions
bc.alpha = zeros(2);
bc.alpha(1, 1) = 0;
bc.alpha(1, 2) = 1;
bc.g0 = @(t)(-0.5);
bc.alpha(2, 1) = 1;
bc.alpha(2, 2) = 0;
bc.g3 = @(t)(10 + 273.15);

%%% Параметры численного решения
Np = 100;            % Число узлов сетки для каждой фазы
%h = 1/N;             % Шаг по координате
tMax = 0.3*365.25*24*3600*0 + 3e6;        % Время, до которого необходимо моделировать, с
tau = 3600*24/24;     % Шаг по времени, с
%M = floor(tMax/tau); % Число шагов по времени

%%% Начальные условия
ic = struct;                 % ic - initial conditions
ic.s0 = 0;            % Начальные положения границы раздела сред, м
ic.s1 = 0.2;
ic.s2 = 0.7;
ic.s3 = 0.8;
ic.u1 = zeros(1, Np) + 273.15 + 10;
ic.u2 = zeros(1, Np) + 273.15 - 0;
ic.u3 = zeros(1, Np) + 273.15 + 10;

[U, X, T, s, t] = StefanProblemSolver(pc, bc, ic, Np, tau, tMax);
figure('DefaultAxesFontSize',15)%, 'windowState', 'maximized')
contourf(T, X, U - 273.15, 'LineColor', 'none', 'LevelStep', 0.5)
hold on
plot(t, s, '-w', 'LineWidth', 2)
hold off
xlabel("t, seconds")
ylabel("X, meters")
colormap(jet)
%caxis([-1 1])
hcb = colorbar;
hcb.Title.String = "Temperature, C";
%set(get(hcb,'Title'),'String','A Title')

