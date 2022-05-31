% Скрипт для запуска численного решения трёхфазной задачи Стефана.
% Приведен пример задачи физических констант, граничных и начальных
% условий, визуализации результатов. Визуализация может занимать время
% порядка десятков секунд, даже если скрипт уже закончил работу.

clear; close all;

%%% Физические константы
pc = getPhysicalConstants("RealLife");

%%% Смешанные краевые условия
% Формат краевых условий:
% alpha00*u1(0, t) + alpha01*du1/dx(0, t) = g0(t) - для левого конца
% alpha10*u2(L, t) + alpha11*du2/dx(L, t) = g1(t) - для правого конца
bc = struct;                  % bc - boundary conditions
bc.alpha = [0 -pc.lambda1; 1 0];
bc.g0 = @(t)(52.6214/1000);
bc.g1 = @(t)(-4.3 + 8*sin(2*pi*t/31556952 + pi/2) + 273.15); %% Внимание, здесь сдвиг по фазе на pi/2

%%% Параметры численного решения
Np = [500 5000 500];            % Число узлов сетки для каждой фазы
tMax = 10*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24*365.25/365.25;        % Шаг по времени, с
tauSave = 3600*24*365.25/24;
NpBoundsSave = inf;

%%% Начальные условия
s = [0; 7; 78; 80];
x2 = linspace(s(2), s(3), Np(2));
u2 = -1*ones(1, length(x2)) + 273.15;
% Uf_adj = 273.15 - 7.43*1e-8*pc.rho2*9.81*( s(3) - s(2) );
% u2 = linspace(Uf_adj, pc.Uf, Np(2));  
accumRate = 0;
ic = struct('s', s, ...
            'dsdt', zeros(4, 1), ...
            'x1', linspace(s(1), s(2), Np(1)), ...
            'u1', 273.15 + zeros(Np(1), 1), ...
            'x2', x2, ...
            'u2', u2, ...
            'x3', linspace(s(3), s(4), Np(3)), ...
            'u3', 273.15 + 0*ones(Np(3), 1), ...
            'tInit', 0);
                                          
% %%% МОДЕЛИРОВАНИЕ ДЛЯ ТОЧКИ №64129 ИЗ РАЗРЕЖЕННОГО ГРИДА НАЧАЛЬНЫХ ДАННЫХ
% %%% Смешанные краевые условия
% % Формат краевых условий:
% % alpha00*u1(0, t) + alpha01*du1/dx(0, t) = g0(t) - для левого конца
% % alpha10*u2(L, t) + alpha11*du2/dx(L, t) = g1(t) - для правого конца
% bc = struct;                  % bc - boundary conditions
% bc.alpha = [0 -pc.lambda1; 1 0];
% bc.g0 = @(t)(52.6214/1000);
% bc.g1 = @(t)( -47.8797 + 0.082/(10*365.25*24*3600)*t + 14.5503*sin(2*pi*t/31556952) + 273.15);
% 
% %%% Параметры численного решения
% Np = [500 5000 500];            % Число узлов сетки для каждой фазы
% tMax = 100*365.25*24*3600;        % Время, до которого необходимо моделировать, с
% tau = 3600*24*365.25/12;        % Шаг по времени, с
% tauSave = 3600*24*365.25*1;
% NpBoundsSave = 100;
% 
% %%% Начальные условия
% s = [820; 3310 - 2490; 3310; 3310];
% x2 = linspace(s(2), s(3), Np(2));
% % WaisDivide = @(z)( -31.799 + 8.8595*1e-3*z - 9.4649*1e-6*z.^2 + 2.657*1e-9 * z.^3 );
% % u2 = flip( WaisDivide( ( x2-x2(1) )/( x2(end)-x2(1) )*3512 ) ) + 273.15;
% u2 = -2*ones(size(x2)) + pc.Uf;
% % Uf_adj = 273.15 - 7.43*1e-8*pc.rho2*9.81*( s(3) - s(2) );
% % u2 = linspace(Uf_adj, pc.Uf, Np(2));  
% accumRate = 85.577331466675000;
% ic = struct('s', s, ...
%             'dsdt', zeros(4, 1), ...
%             'x1', linspace(s(1), s(2), Np(1)), ...
%             'u1', 273.15 + zeros(Np(1), 1), ...
%             'x2', x2, ...
%             'u2', u2, ...
%             'x3', linspace(s(3), s(4), Np(3)), ...
%             'u3', 273.15 + 0*ones(Np(3), 1), ...
%             'tInit', 0);
        
% Если нет скомпилированного mex-файла, скомпилировать
if ~(isfile("mex_TDMA.mexw64") || isfile("mex_TDMA.mexa64"))
    mex -largeArrayDims mex_TDMA.cpp
end

[s, t, U, X, T] = StefanProblemSolver(pc, bc, ic, 'tau', tau, ...
                                                  'tauSave', tauSave, ...
                                                'tMax', tMax, ...
                                                'Np', Np,...
                                                'gridType', 'SigmoidBased', ...
                                                'NpSave', [100 1000 100], ...
                                                'accumRate', accumRate, ...
                                                'NpBoundsSave', NpBoundsSave);
 
%%% Построение графиков
% Общая диаграмма с температурами и границами фаз
figure
subplot(5, 1, 1)
plot(t/3600/24, bc.g1(t) - 273.15)
ylabel("T, C")
axis([-inf inf -inf inf])
subplot(5, 1, [2 5]);
contourf(repmat(T, length(X(:, 1)), 1)/3600/24, X(:, 1:end), U(:, 1:end) - 273.15, ...
    'LineColor', 'none', 'LevelStep', 0.5);
axis([-inf inf s(1, 1) max(s(4, :))])
hold on
plot(t/3600/24, s, '-w', 'LineWidth', 2)
hold off
xlabel("t, days")
ylabel("X, meters")
colormap(jet); 
caxis([-12 3])
%hcb = colorbar;
%hcb.Title.String = "T, C"; hcb.Title.Interpreter = 'latex'; hcb.TickLabelInterpreter = 'latex';

% Закон движения нижней кромки ледника
figure
plot(t/3600/24, s(2, :))
xlabel("time, days")
ylabel("s1, meters")

% Закон движения верхней кромки ледника и поверхности надледниковой воды
figure
plot(t/3600/24, s(3:4, :))
xlabel("time, days")
ylabel("s, meters")

% Изменение массы с поправкой на аккумуляцию
m = (s(2, :) - s(1, :))*pc.rho1 + ...
    (s(3, :) - s(2, :))*pc.rho2 + ...
    (s(4, :) - s(3, :))*pc.rho1 - accumRate/(365.25*24*3600)*t;
fprintf("Mass change: %e kg/m^2\n", m(end)-m(1));

figure
plot(t/3600/24, m)
xlabel("t, days")
ylabel("m, kg")

function savePlot(h, filename)
    savefig(h.Parent, filename);
    print(h.Parent, filename, '-dpng', '-r300');
    print(h.Parent, filename, '-depsc');
end

function pc = getPhysicalConstants(type)
    switch type
        case 'RealLife'
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
            pc.Uf = 273.15;               % Температура фазового перехода, К
        case 'AllOnes'
            pc = struct;
            pc.lambda1 = 1;
            pc.c1 = 1;
            pc.rho1 = 1;
            pc.a1_sq = 1;         
            pc.lambda2 = 1;
            pc.c2 = 1;
            pc.rho2 = 1;
            pc.a2_sq = 1;          
            pc.qf = 1;
            pc.Uf = 1;
        otherwise
            error(sprintf("Type %s is invalid.\nExpected type to match one of these values: %s, %s\n", ...
                type, "RealLife", "AllOnes"));
    end
    
end