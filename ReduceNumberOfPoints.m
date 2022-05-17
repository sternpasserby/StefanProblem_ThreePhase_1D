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
tMax = 4*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24*365.25/365.25;        % Шаг по времени, с
tauSave = 3600*24*365.25;

%%% Начальные условия
s = [0; 7; 80; 80];
x2 = linspace(s(2), s(3), Np(2));
u2 = -1*ones(1, length(x2)) + 273.15;
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

[s, t, U, X, T] = StefanProblemSolver(pc, bc, ic, 'tau', tau, ...
                                              'tauSave', tauSave, ...
                                              'tMax', tMax, ...
                                              'Np', Np,...
                                              'gridType', 'SigmoidBased', ...
                                              'NpSave', [100 1000 100], ...
                                              'accumRate', accumRate);
                                        
subplot(2, 2, 1)
setupPlot( plot(t, s(3:4, :), '--o') );
title( sprintf("Before, N = %d", length(t)) )
subplot(2, 2, 3)
setupPlot( plot(t, s(2, :), '--o') );

[tSp, sSp] = reduceNumOfPointsInS(t, s, 100);
subplot(2, 2, 2)
setupPlot( plot(tSp, sSp(3:4, :), '--o') );
title( sprintf("After, N = %d", length(tSp)) )
subplot(2, 2, 4)
setupPlot( plot(tSp, sSp(2, :), '--o') );

f = @() reduceNumOfPointsInS(t, s, 100);
time = timeit(f, 2);
fprintf("Original number of points: %d\n", length(t));
fprintf("     New number of points: %d\n", length(tSp));
fprintf("                    Ratio: %.2f\n", length(t)/length(tSp));
fprintf("       Time for reduction: %.2e sec\n", time)

function setupPlot(plotObj)
    ax = plotObj.Parent;
    ax.XLabel.String = 't, days';
    ax.YLabel.String = 's, meters';
    
    for i = 1:length(plotObj)
        plotObj(i).XData = plotObj(i).XData/( 3600*24 );
    end
end

function [tNew, sNew] = reduceNumOfPointsInS(t, s, newN)
    N = length(t);
    nPointsPerEvent_halved = 3; % Сколько точек выделить слева и справа на 
                                %   событие зарождения или исчезновения фазы

    % Поиск индексов, где происходит создание или вырождение фазы
    ds1 = s(2, :) - s(1, :);
    ds3 = s(4, :) - s(3, :);
    idPhase = zeros(1, N);
    k = 1;
    for i = 2:N
        if ( abs(ds1(i-1)) < 1e-6 && abs(ds1(i)) > 0 ) || ...
                ( abs(ds1(i-1)) > 0 && abs(ds1(i)) < 1e-6 ) || ...
                ( abs(ds3(i-1)) < 1e-6 && abs(ds3(i)) > 0 ) || ...
                ( abs(ds3(i-1)) > 0  && abs(ds3(i)) < 1e-6 )
            idI = max(i-nPointsPerEvent_halved, 1):min(i+nPointsPerEvent_halved, N);
            idPhase( k:k+length(idI)-1 ) = idI;
            k = k + length(idI);
        end
    end
    idPhase(k:end) = [];
    
    id = unique([ 1, idPhase, 1:ceil(N/newN):N, N ]);
    tNew = t(id);
    sNew = s(:, id);
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
