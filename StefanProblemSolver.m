function [s, t, U, X, T] = StefanProblemSolver(pc, bc, varargin)
%STEFANPROBLEMSOLVER Решатель трехфазной задачи Стефана
%   На вход подаются 3 структуры, число узлов сетки для каждой фазы Np,
%   размер временнОго шага tau и время моделирования tMax. Возвращает 3
%   матрицы U, X, T и вектор s.

%%% ПАРСИНГ ВХОДНЫХ ПАРАМЕТРОВ
% Задание значений входных параметров по умолчанию
defaultIc = struct('s', [-130; -129; 1100; 1100], ...
                   'dsdt', zeros(4, 1), ...
                   'x1', linspace(-130, -129, 100), ...
                   'u1', 273.15 + zeros(100, 1), ...
                   'x2', linspace(-129, 1100, 1000), ...
                   'u2', 273.15 - 2*ones(1000, 1), ...
                   'x3', ones(100, 1)*1100, ...
                   'u3', 273.15 + ones(100, 1), ...
                   'tInit', 0);
defaultAccumRate = 85.5;
defaultTau = 3600*24*14;
defaultTMax = 20*365.25*24*3600; 
defaultChLength = defaultIc.s(3) - defaultIc.s(2);
defaultChTemperature = pc.Uf;
defaultGridType = 'NonUniform';
defaultNpSave = [100 200 100];
defaultTauSave = defaultTau;
defaultDs = 0.01;
defaultMinNewPhaseThickness = 1*1e-3;

% Настройка объекта parserObj типа InputParser
parserObj = inputParser;
parserObj.StructExpand = false;
addRequired(parserObj, 'pc');
addRequired(parserObj, 'bc');
addOptional(parserObj, 'ic', defaultIc);
addParameter(parserObj, 'accumRate', defaultAccumRate);
addParameter(parserObj, 'tau', defaultTau);
addParameter(parserObj, 'tMax', defaultTMax);
addParameter(parserObj, 'chLength', defaultChLength);
addParameter(parserObj, 'chTemperature', defaultChTemperature);
addParameter(parserObj, 'gridType', defaultGridType);
addParameter(parserObj, 'NpSave', defaultNpSave);
addParameter(parserObj, 'tauSave', defaultTauSave);
addParameter(parserObj, 'minDs', defaultDs);
addParameter(parserObj, 'minNewPhaseThickness', defaultMinNewPhaseThickness);

% Распаковка структур
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
  alpha = bc.alpha;
rho = (rho1 + rho2)/2;

% Получение значений параметров от объекта parserObj типа InputParser
parse(parserObj, pc, bc, varargin{:});
accumRate = parserObj.Results.accumRate;
ic = parserObj.Results.ic;
tau = parserObj.Results.tau;
tMax = parserObj.Results.tMax;
x0 = parserObj.Results.chLength;
U0 = parserObj.Results.chTemperature;
gridType = parserObj.Results.gridType;
NpSave = parserObj.Results.NpSave;
tauSave = parserObj.Results.tauSave;
minDs = parserObj.Results.minDs;
minNewPhaseThickness = parserObj.Results.minNewPhaseThickness;
         
% Вычисление безразмерных параметров
t0 = rho1*c1*x0*x0/lambda1;
beta = qf*rho / (rho1*c1*U0);
kappa = c2*rho2/lambda2 * lambda1 / (c1*rho1); 

% Обезразмеривание исходных данных
g0 = @(t)(bc.g0(t*t0)/U0);
g1 = @(t)(bc.g1(t*t0)/U0);
u1 = ic.u1/U0;
u2 = ic.u2/U0;
u3 = ic.u3/U0;
% x1 = ic.x1/x0;
% x2 = ic.x2/x0;
% x3 = ic.x3/x0;
tInit = ic.tInit/t0;
tau = tau/t0;
tauSave = tauSave/t0;
tMax = tMax/t0;
Uf = Uf/U0;
alpha(:, 2) = bc.alpha(:, 2)/x0;
minDs = minDs/x0;
minNewPhaseThickness = minNewPhaseThickness/x0;
s0 = ic.s(1)/x0;

numOfTimeSteps = ceil(tMax/tau);
Np = [length(u1) length(u2) length(u3)];
ksiNew = getGrid(Np(1));
if ~(ic.s(1)==ic.s(2))
    ksi = ( ic.x1 - ic.s(1) ) ./ ( ic.s(2) - ic.s(1) );
    %u1 = csInterp(ksi, u1, ksiNew, [alpha(1, :); 1 0], g0(tInit), Uf);
    u1 = interp1(ksi, u1, ksiNew, 'linear', 'extrap')';
else
    ksi = ic.x1.*NaN;
    u1 = u1.*NaN;
end
ph1 = struct('ksi', ksiNew, ...
             'u', u1, ...
             's', zeros(1, numOfTimeSteps + 1), ...
             'dsdt', 0, ...
             'exists', ~(ic.s(1)==ic.s(2)) );
ph1.s(1) = ic.s(2)/x0;
ksi = ( ic.x2 - ic.s(2) ) ./ ( ic.s(3) - ic.s(2) );
ksiNew = getGrid(Np(2));
ph2 = struct('ksi', ksiNew, ...
             'u', interp1(ksi, u2, ksiNew, 'linear', 'extrap')', ... %'u', csInterp(ksi, u2, ksiNew, [1 0; 1 0], Uf, Uf), ...
             's', zeros(1, numOfTimeSteps + 1), ...
             'dsdt', 0, ...
             'exists', true );
ph2.s(1) = ic.s(3)/x0;
ksi = ( ic.x3 - ic.s(3) ) ./ ( ic.s(4) - ic.s(3) );
ksiNew = getGrid(Np(3));
if ~(ic.s(3)==ic.s(4))
    ksi = ( ic.x1 - ic.s(1) ) ./ ( ic.s(2) - ic.s(1) );
    %u3 = csInterp(ksi, u3, ksiNew, [1 0; alpha(2, :)], Uf, g1(tInit));
    u3 = interp1(ksi, u3, ksiNew, 'linear', 'extrap')';
else
    ksi = ic.x3.*NaN;
    u3 = u3.*NaN;
end
ph3 = struct('ksi', ksiNew, ...
             'u', u3, ...
             's', zeros(1, numOfTimeSteps + 1), ...
             'dsdt', 0, ...
             'exists', ~(ic.s(3)==ic.s(4)) );
ph3.s(1) = ic.s(4)/x0;

t = zeros(1, numOfTimeSteps + 1);

nRows = sum(NpSave);
nCols = min(numOfTimeSteps, ceil(tMax/tauSave) ) + 1;
X = zeros(nRows, nCols);
U = zeros(nRows, nCols);
T = zeros(nRows, nCols);

A1 = sparse(Np(1));
A2 = sparse(Np(2));
A3 = sparse(Np(3));
b1 = zeros(Np(1), 1);
b2 = zeros(Np(2), 1);
b3 = zeros(Np(3), 1);

% Запись начальных условий в выходные массивы
ksiSave1 = getGrid(NpSave(1));
if ph1.exists
    %u1q = csInterp(ph1.ksi, ph1.u, ksiSave1, [alpha(1, :); 1 0], g0(tInit), Uf);
    u1q = interp1(ph1.ksi, ph1.u, ksiSave1, 'linear', 'extrap')';
    x1q = s0 + ksiSave1*( ph1.s(1) - s0 );
else
    x1q = ksiSave1.*NaN;
    u1q = ksiSave1'.*NaN;
end
ksiSave2 = getGrid(NpSave(2));
%u2q = csInterp(ph2.ksi, ph2.u, ksiSave2, [1 0; 1 0], Uf, Uf);
u2q = interp1(ph2.ksi, ph2.u, ksiSave2, 'linear', 'extrap')';
x2q = ph1.s(1) + ksiSave2*( ph2.s(1) - ph1.s(1) );
ksiSave3 = getGrid(NpSave(3));
if ph3.exists
    %u3q = csInterp(ph3.ksi, ph3.u, ksiSave3, [1 0; alpha(2, :)], Uf, g1(tInit));
    u3q = interp1(ph3.ksi, ph3.u, ksiSave3, 'linear', 'extrap')';
    x3q = ph2.s(1) + ksiSave3*( ph3.s(1) - ph2.s(1) );
else
    x3q = ksiSave3.*NaN;
    u3q = ksiSave3'.*NaN;
end
U(:, 1) = [u1q; u2q; u3q];
X(:, 1) = [x1q'; x2q'; x3q'];

saveTime = tauSave;
saveId = 1;

time = tInit;
tau0 = tau;     % Шаг по времени, заданный пользователем
n = 1;

tic;
while time <= tMax
    u1_past = ph1.u;
    u2_past = ph2.u;
    u3_past = ph3.u;
    
    % Вычисление скоростей движения границ
    C1 = 1/beta/( ph1.s(n) - s0 );
    C2 = lambda2/(lambda1*beta*( ph2.s(n) - ph1.s(n) ));
    C3 = 1/beta/( ph3.s(n) - ph2.s(n) );
    h0 = ph2.ksi(2) - ph2.ksi(1);
    h1 = ph2.ksi(3) - ph2.ksi(2);
    h_Npm1 = ph1.ksi( Np(1) ) - ph1.ksi( Np(1)-1 );
    h_Npm2 = ph1.ksi( Np(1)-1 ) - ph1.ksi( Np(1)-2 );
    ph1.dsdt = C2/( h0*h1*(h0+h1) )*( ( u2_past(2) - u2_past(1) )*(h0+h1)^2 - (u2_past(3)-u2_past(1))*h0^2 ) - ...
        C1/( h_Npm1*h_Npm2*(h_Npm1+h_Npm2) )*...
        ( ( u1_past(Np(1)) - u1_past(Np(1)-1) )*(h_Npm2+h_Npm1)^2 + (u1_past(Np(1)-2) - u1_past(Np(1)))*h_Npm1^2 );
    h_Npm1 = ph2.ksi( Np(2) ) - ph2.ksi( Np(2)-1 );
    h_Npm2 = ph2.ksi( Np(2)-1 ) - ph2.ksi( Np(2)-2 );
    h0 = ph3.ksi(2) - ph3.ksi(1);
    h1 = ph3.ksi(3) - ph3.ksi(2);
    ph2.dsdt = C2/( h_Npm1*h_Npm2*(h_Npm1+h_Npm2) )*...
        ( ( u2_past(Np(2)) - u2_past(Np(2)-1) )*(h_Npm2+h_Npm1)^2 + (u2_past(Np(2)-2)-u2_past(Np(2)))*h_Npm1^2 ) - ...
        C3/( h0*h1*(h0+h1) )*( ( u3_past(2) - u3_past(1) )*(h0+h1)^2 - (u3_past(3)-u3_past(1))*h0^2 );
    ph3.dsdt = 0;
    
    % Интегрирование уравнений движения границ
    ph1.s(n+1) = ph1.s(n) + tau*ph1.dsdt;
    ph2.s(n+1) = ph2.s(n) + tau*ph2.dsdt;
    ph3.s(n+1) = ph3.s(n) + tau*ph3.dsdt;
    
    if (abs( ph2.s(n+1)-ph2.s(n) ) > minDs && ph2.exists) || ...
            (abs( ph1.s(n+1)-ph1.s(n) ) > minDs && ph1.exists)
        tau = tau/2;
        continue;
    end
    
    %fprintf("%10.4e %10.4e %10.4e %10.4e\n", tau, ph1.dsdt, ph2.dsdt, ph3.dsdt);
    
    if ~ph3.exists
        ph2.s(n+1) = ph3.s(n+1);
        ph2.dsdt = ph3.dsdt;
    end
    if ~ph1.exists
        ph1.s(n+1) = s0;
        ph1.dsdt = 0;
    end
    
    % Вырождение верхней и нижней фаз
    if ph2.s(n+1) >= ph3.s(n+1) && ph3.exists
        ph2.s(n+1) = ph3.s(n+1);
        ph2.dsdt = ph3.dsdt;
        ph3.exists = false;
    end
    if ph1.s(n+1) <= s0 && ph1.exists
        ph1.s(n+1) = s0;
        ph1.dsdt = 0;
        ph1.exists = false;
    end
    time = time + tau;
    
    Uf_adj = (273.15 - 7.43*1e-8*rho2*9.81*( ph2.s(n+1)-ph1.s(n+1) )*x0)/U0;
    
    % Получение распределения тепла для первой фазы (если она есть)
    if ph1.exists
        [A1, b1] = getSysMat(u1_past, 1, tau, ph1.ksi, ph1.s(n+1), s0, ph1.dsdt, 0, ...
           [alpha(1, :); 1 0], g0(time), Uf_adj);
        ph1.u = solveWithBackslash(A1, b1);
    end
    
    % Получение распределения тепла для второй фазы
    if ph3.exists
        alphaUpper = [1 0];
        gUpper = Uf;
    else
        alphaUpper = alpha(2, :);
        gUpper = g1(time);
    end
    if ph1.exists
        alphaLower = [1 0];
        gLower = Uf_adj;
    else
        alphaLower = alpha(1, :);
        gLower = g0(time);
    end
    [A2, b2] = getSysMat(u2_past, kappa, tau, ph2.ksi, ph2.s(n+1), ph1.s(n+1), ph2.dsdt, ph1.dsdt, ...
        [alphaLower; alphaUpper], gLower, gUpper);
    ph2.u = solveWithBackslash(A2, b2);
    
    % Получение распределения тепла для третьей фазы
    if ph3.exists
        [A3, b3] = getSysMat(u3_past, 1, tau, ph3.ksi, ph3.s(n+1), ph2.s(n+1), ph3.dsdt, ph2.dsdt, ...
           [1 0; alpha(2, :)], Uf, g1(time));
        ph3.u = solveWithBackslash(A3, b3);
    end
    
    % Зарождение верхней фазы
    if ph2.u(end) > Uf && ph2.u(end-1) > Uf && ~ph3.exists
        % Поиск номера последнего узла, который должен быть водой
        id = Np(2);
        for i = Np(2):-1:1
            if ph2.u(i) < Uf
                id = i + 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = ph1.s(n+1) + ph2.ksi.*(ph2.s(n+1) - ph1.s(n+1));
        dl = c2*rho2/qf/rho1*trapz(x(id:end), abs(ph2.u(id:end) - Uf)*U0);
        
        if dl >= minNewPhaseThickness
            ph2.s(n+1) = ph3.s(n+1) - dl;
            ph3.u = ones(Np(3), 1)*Uf;
            ph2.u(id:end) = Uf;
            %ph2.u = csInterp(x, ph2.u, ph1.s(n+1) + ph2.ksi.*(ph2.s(n+1) - ph1.s(n+1)), [1 0; 1 0], Uf_adj, Uf);
            ph2.u = interp1(x, ph2.u, ph1.s(n+1) + ph2.ksi.*(ph2.s(n+1) - ph1.s(n+1)), 'linear', 'extrap')';
            ph3.exists = true;
        end
       
    end
    
    % Зарождение нижней фазы
    if ph2.u(1) > Uf_adj && ph2.u(2) > Uf_adj && ~ph1.exists
        % Поиск номера последнего узла, который должен быть водой
        id = 1;
        for i = 1:Np(2)
            if pf2.u(i) < Uf_adj
                id = i - 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = ps1.s(n+1) + ph2.ksi.*(ph2.s(n+1) - ph1.s(n+1));
        dl = c2*rho2/qf/rho1*trapz(x(1:id), abs(ph2.u(1:id) - Uf_adj)*U0);
        
        if dl >= minNewPhaseThickness
            ph1.s(n+1) = s0 + dl;
            ph1.u = ones(Np(1), 1)*Uf_adj;
            ph2.u(1:id) = Uf_adj;
            ph2.u = interp1(x, ph2.u, ph1.s1(n+1) + ph2.ksi.*(ph2.s2(n+1) - ph1.s1(n+1)), 'linear', 'extrap')';
            ph1.exists = true;
        end
    end
    
    % Проседание льда из-за различной плотности льда и воды
    dL = (ph1.s(n+1) - ph1.s(n) - (ph2.s(n+1) - ph2.s(n)) )*( 1 - rho1/rho2 );
    ph2.s(n+1) = ph2.s(n+1) + dL;
    ph3.s(n+1) = ph3.s(n+1) + dL;
    
    % Аккумуляция
    dL = accumRate/(365.25*24*3600)*tau*t0/rho2/x0;
    ph2.s(n+1) = ph2.s(n+1) + dL;
    ph3.s(n+1) = ph2.s(n+1) + dL;
    
    % Запись результатов
    t(n + 1) = time;
    if (saveTime < time)
        fprintf("Progress: %4.2f%%\n", saveTime/tMax*100);
        saveTime = saveTime + tauSave;
        saveId = saveId + 1;
        
        %ksiSave1 = getGrid(NpSave(1));
        if ph1.exists
            %u1q = csInterp(ph1.ksi, ph1.u, ksiSave1, [alpha(1, :); 1 0], g0(time), Uf);
            u1q = interp1(ph1.ksi, ph1.u, ksiSave1, 'linear', 'extrap')';
            x1q = s0 + ksiSave1*( ph1.s(1) - s0 );
        else
            x1q = ksiSave1.*NaN;
            u1q = ksiSave1'.*NaN;
        end
        ksiSave2 = getGrid(NpSave(2));
        %u2q = csInterp(ph2.ksi, ph2.u, ksiSave2, [1 0; 1 0], Uf, Uf);
        u2q = interp1(ph2.ksi, ph2.u, ksiSave2, 'linear', 'extrap')';
        x2q = ph1.s(1) + ksiSave2*( ph2.s(1) - ph1.s(1) );
        ksiSave3 = getGrid(NpSave(3));
        if ph3.exists
            %u3q = csInterp(ph3.ksi, ph3.u, ksiSave3, [1 0; alpha(2, :)], Uf, g1(time));
            u3q = interp1(ph3.ksi, ph3.u, ksiSave3, 'linear', 'extrap')';
            x3q = ph2.s(1) + ksiSave3*( ph3.s(1) - ph2.s(1) );
        else
            x3q = ksiSave3.*NaN;
            u3q = ksiSave3'.*NaN;
        end
        U(:, saveId) = [u1q; u2q; u3q];
        X(:, saveId) = [x1q'; x2q'; x3q'];
        T(:, saveId) = ones(nRows, 1)*time;
    end
    
    n = n+1;
    tau = tau0;
end
s = [ ones(1, length(t))*s0; ph1.s; ph2.s; ph3.s];

% Масшабирование к исходной размерности
X = X*x0;
T = T*t0;
U = U*U0;
s = s*x0;
t = t*t0;

elapsedTime = toc;
fprintf("Elapsed time for Stefan Problem Solver: %4.2f sec.\n", elapsedTime);
end

function xNew = getGrid(Np)
    X = [0 1/4 3/4 1];
    H = 0.1;
    ppF = struct();
    ppF.form = 'pp';
    ppF.breaks = X;
    c1 = ( (1-H)/3 + H )*(X(2)-X(1));
    ppF.coefs = [ (1-H)/3/(X(2)-X(1))^2 (H-1)/(X(2)-X(1)) 1 0; ...
                   0 0 H c1; ...
                   (1-H)/3/(X(4)-X(3))^2 0 H c1 + H*(X(3)-X(2))];
    ppF.pieces = 3;
    ppF.order = 4;
    ppF.dim = 1;
    
    der_ppF = fnder(ppF, 1);
    
    y = linspace( ppval(ppF, 0), ppval(ppF, 1), Np );
    xNew = linspace(0, 1, Np);
    for i = 1:10
        xNew = xNew - (ppval(ppF, xNew)-y)./(ppval(der_ppF, xNew));
        %fprintf("Error: %6.2e\n", max(abs(xNew-temp)));
        %temp = xNew;
    end
    
    if ~issorted(xNew)
        error("xNew is not monotonically increasing!")
    end
    if ~all(xNew >= 0 & xNew <= 1)
        error("Some of xNew elements lie outside of [0, 1] domain!");
    end
end

function yNew = csInterp(x, y, xq, alpha, g0, g1)
    conds = [1 1];
    if alpha(1, 2) == 0
        conds(1) = 0;
        C0 = g0/alpha(1, 1);
    else
        conds(1) = 1;
        C0 = (g0 - alpha(1, 1)*y(1)) / alpha(1, 2);
    end
    if alpha(2, 2) == 0
        conds(2) = 0;
        C1 = g1/alpha(2, 1);
    else
        conds(2) = 1;
        C1 = (g1 - alpha(2, 1)*y(end)) / alpha(2, 2);
    end
    pp = csape(x, [C0; y; C1]', conds);
    yNew = ppval(pp, xq)';
end
