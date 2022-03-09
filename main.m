clear; close all;

%%% Физические константы
pc = getPhysicalConstants("RealLife");

%%% Смешанные краевые условия
% Формат краевых условий:
% alpha00*u1(0, t) + alpha01*du1/dx(0, t) = g0(t) - для левого конца
% alpha10*u2(L, t) + alpha11*du2/dx(L, t) = g1(t) - для правого конца
bc = struct;                  % bc - boundary conditions
bc.alpha = [0 -pc.lambda1; 1 0];
bc.g0 = @(t)(0.05);
bc.g1 = @(t)(- 4.3 + 8*sin(2*pi*t/31556952 + pi/2) + 273.15); %% Внимание, здесь сдвиг по фазе на pi/2

%%% Параметры численного решения
Np = 1000;            % Число узлов сетки для каждой фазы
tMax = 20*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24*14;     % Шаг по времени, с
tauSave = 3600*24*365.25/4;

%%% Начальные условия
% 1897000,-3000,820,0,3310,2490,52.621367149353,85.577331466675,1.6650680303574
ic = struct;
ic.s0 = 820;
ic.s1 = 3310 - 2490;
ic.s2 = 3310;
ic.s3 = 3310;
%ic.accumRate = 85.577331466675;
bc.g0 = @(t)(52.621367149353/1000);
ic.u1 = zeros(1, Np) + 273.15 + 0;
ic.u2 = zeros(1, Np) + 273.15 - 2;
ic.u3 = zeros(1, Np) + 273.15 + 1;

[s, t, U, X, T] = StefanProblemSolver(pc, bc, 'tau', 3600*24*30, ...
                                              'tauSave', 3600*24*30, ...
                                              'tMax', 20*365.25*24*3600);
% plot(s', '.')
% figure%('DefaultAxesFontSize',15)%, 'windowState', 'maximized')
% subplot(5, 1, [2 5]);
% contourf(T(:, 1:end)/3600/24, X(:, 1:end), U(:, 1:end) - 273.15, 'LineColor', 'none', 'LevelStep', 0.5);
% axis([-inf inf ic.s0 ic.s3])
% hold on
% plot(t/3600/24, s, '-w', 'LineWidth', 2)
% hold off
% xlabel("t, days")
% ylabel("X, meters")
% colormap(jet)
% caxis([-12 3])
% hcb = colorbar;
% hcb.Title.String = "T, C";

% Проверка закона сохранения массы
figure
% m = (s(2, :) - s(1, :))*pc.rho1 + (s(3, :) - s(2, :))*pc.rho2 + (s(4, :) - s(3, :))*pc.rho1;
% subplot(3, 1, 1)
% plot(t/3600/24, m, 'LineWidth', 2)
% xlabel("t, days")
% ylabel("m, kg")
% title("m(t)")
% set(gca, 'FontSize', 20)
subplot(2, 1, 1)
plot(t/3600/24, s, 'LineWidth', 2)
xlabel("t, days")
ylabel("X, meters")
axis([-inf inf min(s(4, :)) max(s(4, :))])
title("Координата поверхности ледника")
set(gca, 'FontSize', 20)
subplot(2, 1, 2)
plot(t/3600/24, s, 'LineWidth', 2)
xlabel("t, days")
ylabel("X, meters")
axis([-inf inf min(s(2, :)) max(s(2, :))])
title("Координата нижней кромки ледника")
set(gca, 'FontSize', 20)

% figure
% hold on
% tId = [1 400 800];
% legendArray = cell(1, length(tId));
% for i = 1:length(tId)
%     plot(X(101:200, tId(i)), U(101:200, tId(i)) - 273.15, 'LineWidth', 2);
%     legendArray(i) = {sprintf("t = %8.2f years", T(1, tId(i))/3600/24/365.25)};
% end
% hold off
% title("Распределение температуры льда в разные моменты времени.")
% xlabel("X, m");
% ylabel("T, C")
% legend(legendArray, 'location', 'best')
% axis([-inf inf -3 inf])
% set(gca, 'FontSize', 24)

% h = figure;
% axis tight manual
% filename = 'testnew51.gif';
% %vid = VideoWriter(filename);
% %open(vid);
% for i = 1:30:length(t)
%     bar([1;nan], [s(1, i), s(2, i) - s(1, i), s(3, i) - s(2, i), s(4, i) - s(3, i); nan(1,4)], 'stacked')
%     title(sprintf("time = %6.2f days", t(i)/3600/24));
%     drawnow
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0);
%     end
% end

%hcb.Position = [0.9123 0.1101 0.03 0.6423];
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
