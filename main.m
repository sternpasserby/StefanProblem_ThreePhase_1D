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
bc.alpha = zeros(6, 2);
bc.alpha = [0 -pc.lambda1; 
            1 0
            1 0
            1 0
            1 0
            1 0];
bc.g0 = @(t)(0.05);
bc.g1 = @(t)(pc.Uf);
bc.g2 = @(t)(pc.Uf);
bc.g3 = @(t)(pc.Uf);
bc.g4 = @(t)(pc.Uf);
bc.g5 = @(t)(- 4.3 + 8*sin(2*pi*t/31556952 + pi/2) + 273.15); %% Внимание, здесь сдвиг по фазе на pi/2
%bc.g5 = @(t)(0 + 273.15);

%%% Параметры численного решения
Np = 500;            % Число узлов сетки для каждой фазы
%h = 1/N;             % Шаг по координате
tMax = 10*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24;     % Шаг по времени, с
%M = floor(tMax/tau); % Число шагов по времени

%%% Начальные условия
ic = struct;                 % ic - initial conditions
ic.s0 = 0;            % Начальные положения границы раздела сред, м
ic.s1 = 1;
ic.s2 = 9;
ic.s3 = 10;
ic.u1 = zeros(1, Np) + 273.15 + 0;
ic.u2 = zeros(1, Np) + 273.15 - 1;
ic.u3 = zeros(1, Np) + 273.15 + 1;

[U, X, T, s, t] = StefanProblemSolver(pc, bc, ic, Np, tau, tMax, 100, tau);
figure%('DefaultAxesFontSize',15)%, 'windowState', 'maximized')
subplot(5, 1, [2 5]);
contourf(T/3600/24, X, U - 273.15, 'LineColor', 'none', 'LevelStep', 0.5);
axis([-inf inf ic.s0 ic.s3])
hold on
plot(t/3600/24, s, '-w', 'LineWidth', 2)
hold off
xlabel("t, days")
ylabel("X, meters")
colormap(jet)
caxis([-12 8])
hcb = colorbar;
%hcb.Position = [0.9123 0.1101 0.03 0.6423];
hcb.Title.String = "T, C";
%hcb.Title.HorizontalAlignment = 'left';
%set(get(hcb,'Title'),'String','A Title')
%axis([2600 3450 79 80])

%subplot(5, 1, 1);
%hcb.Position(1) = ax.Position(1);
%plot(t/3600/24, bc.g5(t) - pc.Uf)
%axis([min(t)/3600/24 max(t)/3600/24 -inf inf])
%axis([2600 3450 -inf inf])
%ylabel("Temperature, C")

% plot(t/3600/24/365.25, s(2, :))
% xlabel("t, years")
% ylabel("X, meters")
% lg = legend("$s_2(t)$");
% lg.Interpreter = 'latex';

