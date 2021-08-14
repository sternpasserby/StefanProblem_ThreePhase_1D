function [U, X, T, s, t] = StefanProblemSolver(pc, bc, ic, Np, tau, tMax)
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
g0 = bc.g0;
g3 = bc.g3;
u1 = ic.u1;
u2 = ic.u2;
u3 = ic.u3;

N = Np - 1;
h = 1/N;             % Шаг по координате
M = floor(tMax/tau); % Число шагов по времени
s0 = zeros(1, M + 1);
s1 = zeros(1, M + 1);
s2 = zeros(1, M + 1);
s3 = zeros(1, M + 1);
 t = zeros(1, M + 1);
s1(1) = ic.s1;
s2(1) = ic.s2;
s3(1) = ic.s3;

ksi = linspace(0, 1, Np)';
A = zeros(Np);
b = zeros(Np, 1);

%maxNRows = 500;                 % Максимальное число строк для выходных матрицы
maxNCols = 500;                 % Максимальное число столбцой для выходных матриц
%nRows = min(2*Np, maxNRows);
nRows = 3*Np;
nCols = min(M, maxNCols);
X = zeros(nRows, nCols);
U = zeros(nRows, nCols);
T = zeros(nRows, nCols);

saveTime = 0;
saveStep = tMax/nCols;
saveId = 0;

tic;
for n = 1:M
    u1_past = u1;
    u2_past = u2;
    u3_past = u3;
    
    % Вычисление скоростей движения границ
    C1 = lambda2/qf/rho/( s2(n) - s1(n) );
    C2 = lambda1/qf/rho/( s1(n) - s0(n) );
    ds0dt = 0;
    ds1dt = C1/(2*h)*(-3*u2_past(1) + 4*u2_past(2) - u2_past(3)) - ...
        C2/(2*h)*(3*u1_past(end) - 4*u1_past(end-1) + u1_past(end-2));
    C2 = lambda1/qf/rho/( s3(n) - s2(n) );
    ds2dt = C1/(2*h)*(3*u2_past(end) - 4*u2_past(end-1) + u2_past(end-2)) - ...
        C2/(2*h)*(-3*u3_past(1) + 4*u3_past(2) - u3_past(3));
    ds3dt = 0;
    
    % Интегрирование уравнений движения границ
    s0(n+1) = s0(n) + tau*ds0dt;
    s1(n+1) = s1(n) + tau*ds1dt;
    s2(n+1) = s2(n) + tau*ds2dt;
    s3(n+1) = s3(n) + tau*ds3dt;
    
    % Проседание льда из-за различной плотности льда и воды
    dL = (s1(n+1) - s1(n))*( 1 - rho1/rho2 );
    s1(n+1) = s1(n+1) + dL;
    s2(n+1) = s2(n+1) + dL;
    s3(n+1) = s3(n+1) + dL;
    
    % Получение распределения тепла для первой фазы
    [A, b] = getSysMat(u1_past, 1/a1_sq, tau, h, s1(n+1), s0(n + 1), ds1dt, ds0dt, ...
        [alpha(1, :); 1 0], Uf, g0(n*tau));
    u1 = A \ b;
    
    % Получение распределения тепла для второй фазы
    [A, b] = getSysMat(u2_past, 1/a2_sq, tau, h, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
        [1 0; 1 0], Uf, Uf);
    u2 = A \ b;
    
    % Получение распределения тепла для третьей фазы
    [A, b] = getSysMat(u3_past, 1/a1_sq, tau, h, s3(n+1), s2(n+1), ds3dt, ds2dt, ...
        [1 0; alpha(2, :)], g3(n*tau), Uf);
    u3 = A \ b;
    
    % Запись результатов
    t(n + 1) = tau*n;
    if (saveTime < t(n+1))
        saveTime = saveTime + saveStep;
        saveId = saveId + 1;
        X(:, saveId) = [ s0(n) + ksi.*(s1(n+1) - s0(n)); s1(n+1) + ksi.*(s2(n+1) - s1(n+1));...
            s2(n+1) + ksi.*(s3(n+1) - s2(n+1))];
        T(:, saveId) = ones(nRows, 1)*tau*n;
        U(:, saveId) = [u1; u2; u3];
    end
end
s = [s0;s1;s2;s3];
elapsedTime = toc;
fprintf("Elapsed time for Stefan Problem Solver: %4.2f sec.\n", elapsedTime);
end

