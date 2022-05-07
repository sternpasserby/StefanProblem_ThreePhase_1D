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
%bc.g1 = @(t)(-4.3 + 8*sin(2*pi*t/31556952 + pi/2) + 273.15); %% Внимание, здесь сдвиг по фазе на pi/2
bc.g1 = @(t)( -47.8797 + 0.082/(10*365.25*24*3600)*t + 14.5503*sin(2*pi*t/31556952) + 273.15);

%%% Параметры численного решения
Np = [500 5000 500];            % Число узлов сетки для каждой фазы
tMax = 1000*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24*365.25/3;     % Шаг по времени, с
tauSave = 3600*24*365.25*10;

%%% Начальные условия
% 1897000,-3000,820,0,3310,2490,52.621367149353,85.577331466675,1.6650680303574
% ic = struct;
% ic.s0 = 820;
% ic.s1 = 3310 - 2490;
% ic.s2 = 3310;
% ic.s3 = 3310;
% %ic.accumRate = 85.577331466675;
% bc.g0 = @(t)(52.621367149353/1000);
% ic.u1 = zeros(1, Np) + 273.15 + 0;
% ic.u2 = zeros(1, Np) + 273.15 - 2;
% ic.u3 = zeros(1, Np) + 273.15 + 1;

WaisDivide = @(z)( -31.799 + 8.8595*1e-3*z - 9.4649*1e-6*z.^2 + 2.657*1e-9 * z.^3 );
% z = linspace(0, 4000, 10000); 
% plot(z, WaisDivide(z))
% Uf_adj = @(z)( - 7.43*1e-8*pc.rho2*9.81*z);
% hold on
% plot(z, Uf_adj(z));
% hold off

s = [820; 3310 - 2490; 3310; 3310];
x2 = linspace(s(2), s(3), Np(2));
%u2 = flip( WaisDivide( ( x2-x2(1) )/( x2(end)-x2(1) )*3512 ) ) + 273.15;
%u2 = -2*ones(size(x2)) + 273.15;
Uf_adj = (273.15 - 7.43*1e-8*pc.rho2*9.81*( s(3) - s(2) ))
u2 = linspace(Uf_adj, bc.g1(0), Np(2));
ic = struct('s', s, ...
            'dsdt', zeros(4, 1), ...
            'x1', linspace(s(1), s(2), Np(1)), ...
            'u1', 273.15 + zeros(Np(1), 1), ...
            'x2', x2, ...
            'u2', u2, ...
            'x3', linspace(s(3), s(4), Np(3)), ...
            'u3', 273.15 + ones(Np(3), 1), ...
            'tInit', 0);

[s, t, U, X, T] = StefanProblemSolver(pc, bc, ic, 'tau', tau, ...
                                              'tauSave', tauSave, ...
                                              'tMax', tMax, ...
                                              'Np', Np,...
                                              'gridType', 'SigmoidBased', ...
                                              'NpSave', [100 1000 100]);
plot(t, s)
 
% figure
% plot(X(:, 1), U(:, 1))
% hold on
% plot(X(:, end-1), U(:, end-1))
% hold off
% figure%('DefaultAxesFontSize',15)%, 'windowState', 'maximized')
% %subplot(5, 1, [2 5]);
% contourf(T(:, 1:end)/3600/24, X(:, 1:end), U(:, 1:end) - 273.15, 'LineColor', 'none', 'LevelStep', 0.5);
% axis([-inf inf s(1, 1) max(s(4, :))])
% hold on
% plot(t/3600/24, s, '-w', 'LineWidth', 2)
% hold off
% xlabel("t, days")
% ylabel("X, meters")
% colormap(jet)
% caxis([-12 3])
% hcb = colorbar;
% hcb.Title.String = "T, C";

%Проверка закона сохранения массы
% figure
% m = (s(2, :) - s(1, :))*pc.rho1 + (s(3, :) - s(2, :))*pc.rho2 + (s(4, :) - s(3, :))*pc.rho1;
% subplot(3, 1, 1)
% plot(t/3600/24, m, 'LineWidth', 2)
% xlabel("t, days")
% ylabel("m, kg")
% title("m(t)")
% set(gca, 'FontSize', 20)
% subplot(3, 1, 2)
% plot(t/3600/24, s, 'LineWidth', 2)
% xlabel("t, days")
% ylabel("X, meters")
% axis([-inf inf min(s(4, :)) max(s(4, :))])
% title("Координата поверхности ледника")
% %set(gca, 'FontSize', 20)
% subplot(3, 1, 3)
% plot(t/3600/24, s, 'LineWidth', 2)
% xlabel("t, days")
% ylabel("X, meters")
% axis([-inf inf min(s(2, :)) max(s(2, :))])
% title("Координата нижней кромки ледника")
%set(gca, 'FontSize', 20)

% figure
% hold on
% tId = [1 100 523];
% legendArray = cell(1, length(tId));
% for i = 1:length(tId)
%     plot(X(101:100+10000, tId(i)), U(101:100+10000, tId(i)) - 273.15, 'LineWidth', 2);
%     legendArray(i) = {sprintf("t = %8.2f years", T(1, tId(i))/3600/24/365.25)};
% end
% hold off
% title("Распределение температуры льда в разные моменты времени.")
% xlabel("X, m");
% ylabel("T, C")
% legend(legendArray, 'location', 'best')
% axis([-inf inf -inf inf])
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
