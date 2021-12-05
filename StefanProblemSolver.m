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
%alpha(1, 2) = bc.alpha(1, 2)/x0;
%alpha(2, 2) = bc.alpha(2, 2)/x0;

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
X = zeros(nRows, nCols);
U = zeros(nRows, nCols);
T = zeros(nRows, nCols);

ksi = linspace(0, 1, Np)';
ksi_save = linspace(0, 1, nRows/3)';
A = sparse(Np);
b = zeros(Np, 1);

saveTime = 0;
saveId = 0;

%htc = 119.7;
isThreePhase = true;

time = tau;
tau0 = tau;     % Шаг по времени, заданный пользователем
n = 1;
dsMin = 0.01/x0;

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
    
    if abs(s2(n+1) - s2(n)) > dsMin && isThreePhase
        tau = tau/2;
        continue;
    end
    
    % Проседание льда из-за различной плотности льда и воды
    %dL = (s1(n+1) - s1(n))*( 1 - rho1/rho2 );
    %s1(n+1) = s1(n+1) + dL;
    %s2(n+1) = s2(n+1) + dL;
    %s3(n+1) = s3(n+1) + dL;
   
    if s2(n+1) >= s3(n+1)
        s2(n+1) = s3(n+1);
        ds2dt = ds3dt;
        isThreePhase = false;
    end
    if s1(n+1) < s0(n+1)
        s1(n+1) = s1(n);
    end
    
    % Получение распределения тепла для первой фазы
    [A, b] = getSysMat(u1_past, 1, tau, h, s1(n+1), s0(n + 1), ds1dt, ds0dt, ...
       [alpha(1, :); alpha(2, :)], g0(time), g1(time));
%     [A, b] = getSysMat(u1_past, 1/a1_sq, tau, h, s1(n+1), s0(n + 1), ds1dt, ds0dt, ...
%        [alpha(1, :); 1 0], g0(n*tau), Uf);
    u1 = A \ b;
    
    % Получение распределения тепла для второй фазы
    if isThreePhase
        [A, b] = getSysMat(u2_past, kappa, tau, h, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
           [alpha(3, :); alpha(4, :)], g2(time), g3(time));
%         [A, b] = getSysMat(u2_past, 1/a2_sq, tau, h, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
%            [1 0; 1 0], Uf, Uf);
        u2 = A \ b;
        %u2 = solveWithThomas(A, b);

        % Получение распределения тепла для третьей фазы
        [A, b] = getSysMat(u3_past, 1, tau, h, s3(n+1), s2(n+1), ds3dt, ds2dt, ...
           [alpha(5, :); alpha(6, :)], g4(time), g5(time));
%         [A, b] = getSysMat(u3_past, 1/a1_sq, tau, h, s3(n+1), s2(n+1), ds3dt, ds2dt, ...
%            [1 0; alpha(2, :)], Uf, g3(n*tau));
        u3 = A \ b;
        %u3 = solveWithThomas(A, b);
    else
        %ds2dt = C2/(2*h)*(-u2_past(end) + 4*u2_past(end-1) - 3*u2_past(end-2));
        [A, b] = getSysMat(u2_past, kappa, tau, h, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
           [alpha(3, :); alpha(6, :)], g2(time), g5(time));
%         [A, b] = getSysMat(u2_past, 1/a2_sq, tau, h, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
%             [1 0; alpha(2, :)], Uf, g3(n*tau));
        u2 = A \ b;
        %u2 = solveWithThomas(A, b);
        
        u3 = u2;
    end
    
    if u2(end) > Uf && u2(end-1) > Uf && ~isThreePhase
        % Поиск номера последнего узла, который должен быть водой
        id = length(u2);
        for s = length(u2):-1:1
            if u2(s) < Uf
                id = s + 1;
                break;
            end
        end
        
        x = s1(n+1) + ksi.*(s2(n+1) - s1(n+1));
        dl = c2/qf*trapz(x(id:end), abs(u2(id:end) - Uf)*U0)/x0;
        s2(n+1) = s3(n+1) - dl;
        u3 = ones(Np, 1)*Uf;
        u2 = interp1(x, u2, s1(n+1) + ksi.*(s2(n+1) - s1(n+1)));
        
        %plot(s1(n+1) + ksi.*(s2(n+1) - s1(n+1)), u2, s2(n+1) + ksi.*(s3(n+1) - s2(n+1)), u3)
        isThreePhase = true;
    end
    
    % Запись результатов
    t(n + 1) = time;
    if (saveTime < time)
        fprintf("Progress: %4.2f%%\n", saveTime/tMax*100);
        saveTime = saveTime + tau_save;
        saveId = saveId + 1;
        
        x1q = s0(n+1) + ksi_save.*(s1(n+1) - s0(n));
        x2q = s1(n+1) + ksi_save.*(s2(n+1) - s1(n+1));
        x1 = s0(n+1) + ksi.*(s1(n+1) - s0(n));
        x2 = s1(n+1) + ksi.*(s2(n+1) - s1(n+1));
        if isThreePhase
            x3q = s2(n+1) + ksi_save.*(s3(n+1) - s2(n+1));
            x3 = s2(n+1) + ksi.*(s3(n+1) - s2(n+1));
            X(:, saveId) = [x1q; x2q; x3q];
            U(:, saveId) = [interp1(x1, u1, x1q); interp1(x2, u2, x2q); interp1(x3, u3, x3q)];
        else
            X(:, saveId) = [ x1q; x2q; ksi_save.*NaN];
            U(:, saveId) = [interp1(x1, u1, x1q); interp1(x2, u2, x2q); ksi_save.*NaN];
        end
        T(:, saveId) = ones(nRows, 1)*time;
        %U(:, saveId) = [u1; u2; u3];
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

