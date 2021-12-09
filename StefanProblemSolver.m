function [U, X, T, s, t] = StefanProblemSolver(pc, bc, ic, Np, tau, tMax, Np_save, tau_save)
%STEFANPROBLEMSOLVER Решатель трехфазной задачи Стефана
%   На вход подаются 3 структуры, число узлов сетки для каждой фазы Np,
%   размер временнОго шага tau и время моделирования tMax. Возвращает 3
%   матрицы U, X, T и вектор s.

%%% Распаковка структур
lambda1 = pc.lambda1;
c1 = pc.c1;
rho1 = pc.rho1;
a1_sq = pc.a1_sq;
lambda2 = pc.lambda2;
c2 = pc.c2;
rho2 = pc.rho2;
a2_sq = pc.a2_sq;
Uf = pc.Uf;
qf = pc.qf;
rho = (rho1 + rho2)/2;
alpha = bc.alpha;

% Задание характерных параметров для обезразмеривания
x0 = 1;
t0 = rho1*c1*x0*x0/lambda1;
U0 = Uf;
beta = qf*rho / (rho1*c1*U0);
kappa = c2*rho2/lambda2 * lambda1 / (c1*rho1); 

% Обезразмеривание исходных данных
g0 = @(t)(bc.g0(t*t0)/U0);
g1 = @(t)(bc.g1(t*t0)/U0);
g2 = @(t)(bc.g2(t*t0)/U0);
g3 = @(t)(bc.g3(t*t0)/U0);
g4 = @(t)(bc.g4(t*t0)/U0);
g5 = @(t)(bc.g5(t*t0)/U0);
u1 = ic.u1/U0;
u2 = ic.u2/U0;
u3 = ic.u3/U0;
tau = tau/t0;
tau_save = tau_save/t0;
tMax = tMax/t0;
Uf = Uf/U0;
alpha(:, 2) = bc.alpha(:, 2)/x0;

% Сохранение изначальных граничных условий для второй фазы
% alpha2_init = alpha(3:4, :);
% g2_init = g2;
% g3_init = g3;

N = Np - 1;
h = 1/N;             % Шаг по координате
M = round(tMax/tau); % Число шагов по времени
s0 = zeros(1, M + 1);
s1 = zeros(1, M + 1);
s2 = zeros(1, M + 1);
s3 = zeros(1, M + 1);
 t = zeros(1, M + 1);
s0(1) = ic.s0/x0;
s1(1) = ic.s1/x0;
s2(1) = ic.s2/x0;
s3(1) = ic.s3/x0;

nRows = min(3*Np, 3*Np_save);
nCols = min(M, round(tMax/tau_save) );
X = zeros(nRows, nCols + 1);
U = zeros(nRows, nCols + 1);
T = zeros(nRows, nCols + 1);

ksi = linspace(0, 1, Np)';
ksi_save = linspace(0, 1, nRows/3)';
A = sparse(Np);
b = zeros(Np, 1);

isUpperPhase = true;
isLowerPhase = true;
if s0(1) == s1(1)
    isLowerPhase = false;
end
if s3(1) == s2(1)
    isUpperPhase = true;
end

% Запись начальных условий в выходные массивы
if isLowerPhase
    x1 = s0(1) + ksi.*(s1(1) - s0(1));
    x1q = s0(1) + ksi_save.*(s1(1) - s0(1));
    u1q = interp1(x1, u1, x1q);
else
    x1q = ksi_save.*NaN;
    u1q = ksi_save.*NaN;
end
x2 = s1(1) + ksi.*(s2(1) - s1(1));
x2q = s1(1) + ksi_save.*(s2(1) - s1(1));
u2q = interp1(x2, u2, x2q);
if isUpperPhase
    x3 = s2(1) + ksi.*(s3(1) - s2(1));
    x3q = s2(1) + ksi_save.*(s3(1) - s2(1));
    u3q = interp1(x3, u3, x3q);
else
    x3q = ksi_save.*NaN;
    u3q = ksi_save.*NaN;
end
U(:, 1) = [u1q; u2q; u3q];
X(:, 1) = [x1q; x2q; x3q];

saveTime = tau_save;
saveId = 1;

time = tau;
tau0 = tau;     % Шаг по времени, заданный пользователем
n = 1;
dsMin = 0.01/x0;
dlMin = 1*1e-3/x0; % Нижняя граница толщины новой фазы

tic;
while time <= tMax
    u1_past = u1;
    u2_past = u2;
    u3_past = u3;
    
    % Вычисление скоростей движения границ
    C1 = 1/beta/( s1(n) - s0(n) );
    C2 = lambda2/(lambda1*beta*( s2(n) - s1(n) ));
    C3 = 1/beta/( s3(n) - s2(n) );
    ds0dt = 0;
    ds1dt = C2/(2*h)*(-3*u2_past(1) + 4*u2_past(2) - u2_past(3)) - ...
        C1/(2*h)*(3*u1_past(end) - 4*u1_past(end-1) + u1_past(end-2));
    ds2dt = C2/(2*h)*(-u2_past(end) + 4*u2_past(end-1) - 3*u2_past(end-2)) - ...
        C3/(2*h)*(-3*u3_past(1) + 4*u3_past(2) - u3_past(3));
    ds3dt = 0;
    
    % Интегрирование уравнений движения границ
    s0(n+1) = s0(n) + tau*ds0dt;
    s1(n+1) = s1(n) + tau*ds1dt;
    s2(n+1) = s2(n) + tau*ds2dt;
    s3(n+1) = s3(n) + tau*ds3dt;
    
    if abs(s2(n+1) - s2(n)) > dsMin && isUpperPhase
        tau = tau/2;
        continue;
    end
    if abs(s1(n+1)-s1(n)) > dsMin && isLowerPhase
        tau = tau/2;
        continue;
    end
    
    if ~isUpperPhase
        s2(n+1) = s3(n+1);
        ds2dt = ds3dt;
    end
    if ~isLowerPhase
        s1(n+1) = s0(n+1);
        ds1dt = ds0dt;
    end
    
    % Проседание льда из-за различной плотности льда и воды
    %dL = (s1(n+1) - s1(n))*( 1 - rho1/rho2 );
    %s1(n+1) = s1(n+1) + dL;
    %s2(n+1) = s2(n+1) + dL;
    %s3(n+1) = s3(n+1) + dL;
    
    % Вырождение верхней и нижней фаз
    if s2(n+1) >= s3(n+1) && isUpperPhase
        s2(n+1) = s3(n+1);
        ds2dt = ds3dt;
        %alpha(4, :) = alpha(6, :);
        %g3 = g5;
        isUpperPhase = false;
    end
    if s1(n+1) <= s0(n+1) && isLowerPhase
        s1(n+1) = s0(n+1);
        ds1dt = ds0dt;
        %alpha(3, :) = alpha(1, :);
        %g2 = g0;
        isLowerPhase = false;
    end
    
    % Получение распределения тепла для первой фазы (если она есть)
    if isLowerPhase
        [A, b] = getSysMat(u1_past, 1, tau, h, s1(n+1), s0(n+1), ds1dt, ds0dt, ...
           alpha(1:2, :), g0(time), g1(time));
        u1 = A \ b;
    end
    
    % Получение распределения тепла для второй фазы
    if isUpperPhase
        %sUpper = s2(n+1);
        %dsUpperdt = ds2dt;
        alphaUpper = alpha(4, :);
        gUpper = g3;
    else
        %sUpper = s3(n+1);
        %dsUpperdt = ds3dt;
        alphaUpper = alpha(6, :);
        gUpper = g5;
    end
    if isLowerPhase
        %sLower = s1(n+1);
        %dsLowerdt = ds1dt;
        alphaLower = alpha(3, :);
        gLower = g2;
    else
        %sLower = s0(n+1);
        %dsLowerdt = ds0dt;
        alphaLower = alpha(1, :);
        gLower = g0;
    end
    [A, b] = getSysMat(u2_past, kappa, tau, h, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
        [alphaLower; alphaUpper], gLower(time), gUpper(time));
    u2 = A \ b;
    %u2 = solveWithThomas(A, b);
    
    % Получение распределения тепла для третьей фазы
    if isUpperPhase
        [A, b] = getSysMat(u3_past, 1, tau, h, s3(n+1), s2(n+1), ds3dt, ds2dt, ...
           alpha(5:6, :), g4(time), g5(time));
        u3 = A \ b;
    end
        %u3 = solveWithThomas(A, b);
%     else
%         [A, b] = getSysMat(u2_past, kappa, tau, h, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
%            [alpha(3, :); alpha(6, :)], g2(time), g5(time));
%         u2 = A \ b;
%         %u2 = solveWithThomas(A, b);
%     end
    
    % Зарождение верхней фазы
    if u2(end) > Uf && u2(end-1) > Uf && ~isUpperPhase
        % Поиск номера последнего узла, который должен быть водой
        id = length(u2);
        for i = length(u2):-1:1
            if u2(i) < Uf
                id = i + 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = s1(n+1) + ksi.*(s2(n+1) - s1(n+1));
        dl = c2/qf*trapz(x(id:end), abs(u2(id:end) - Uf)*U0)/x0;
        
        if dl >= dlMin
            s2(n+1) = s3(n+1) - dl;
            u3 = ones(Np, 1)*Uf;
            u2(id:end) = Uf;
            u2 = interp1(x, u2, s1(n+1) + ksi.*(s2(n+1) - s1(n+1)));
            isUpperPhase = true;
        end
       
    end
    
    % Зарождение нижней фазы
    if u2(1) > Uf && u2(2) > Uf && ~isLowerPhase
        % Поиск номера последнего узла, который должен быть водой
        id = 1;
        for i = 1:length(u2)
            if u2(i) < Uf
                id = i - 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = s1(n+1) + ksi.*(s2(n+1) - s1(n+1));
        dl = c2/qf*trapz(x(1:id), abs(u2(1:id) - Uf)*U0)/x0;
        
        if dl >= dlMin
            s1(n+1) = s0(n+1) + dl;
            u1 = ones(Np, 1)*Uf;
            u2(1:id) = Uf;
            u2 = interp1(x, u2, s1(n+1) + ksi.*(s2(n+1) - s1(n+1)));
            isLowerPhase = true;
        end
    end
    
    % Запись результатов
    t(n + 1) = time;
    if (saveTime < time)
        fprintf("Progress: %4.2f%%\n", saveTime/tMax*100);
        saveTime = saveTime + tau_save;
        saveId = saveId + 1;
        
        x1 = s0(n+1) + ksi.*(s1(n+1) - s0(n+1));
        x2 = s1(n+1) + ksi.*(s2(n+1) - s1(n+1));
        x3 = s2(n+1) + ksi.*(s3(n+1) - s2(n+1));
        if isLowerPhase
            x1q = s0(n+1) + ksi_save.*(s1(n+1) - s0(n+1));
            u1q = interp1(x1, u1, x1q);
        else
            x1q = ksi_save.*NaN;
            u1q = ksi_save.*NaN;
        end
        x2q = s1(n+1) + ksi_save.*(s2(n+1) - s1(n+1));
        u2q = interp1(x2, u2, x2q);
        if isUpperPhase
            x3q = s2(n+1) + ksi_save.*(s3(n+1) - s2(n+1));
            u3q = interp1(x3, u3, x3q);
        else
            x3q = ksi_save.*NaN;
            u3q = ksi_save.*NaN;
        end
        
        %Nf = round(nRows/3);
%         if isUpperPhase
%             X(:, saveId) = [x1q; x2q; x3q];
%             U(:, saveId) = [interp1(x1, u1, x1q); interp1(x2, u2, x2q); interp1(x3, u3, x3q)];
%         else
%             X(:, saveId) = [ x1q; x2q; ksi_save.*NaN];
%             U(:, saveId) = [interp1(x1, u1, x1q); interp1(x2, u2, x2q); ksi_save.*NaN];
%         end
        U(:, saveId) = [u1q; u2q; u3q];
        X(:, saveId) = [x1q; x2q; x3q];
        T(:, saveId) = ones(nRows, 1)*time;
    end
    
    n = n+1;
    time = time + tau;
    tau = tau0;
end
s = [s0;s1;s2;s3];

% Масшабирование к исходной размерности
X = X*x0;
T = T*t0;
U = U*U0;
s = s*x0;
t = t*t0;

elapsedTime = toc;
fprintf("Elapsed time for Stefan Problem Solver: %4.2f sec.\n", elapsedTime);
end

