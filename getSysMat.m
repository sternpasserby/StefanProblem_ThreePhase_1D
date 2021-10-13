function [A, b] = getSysMat(u_past, kappa, tau, h, s_j, s_jm1, ...
    ds_jdt, ds_jm1dt, alpha, g_jm1, g_j)
%GETSYSTEMMATRIX Вычисление матрицы системы для решения задачи
%Стефана методом конечных разностей с применением метода выпрямления
%фронта. Граничные условия учитываются.
%   u_past - распределение температур на предыдущем временном шаге

C1 = kappa*(s_j - s_jm1)^2/tau;
Np = length(u_past);
A = zeros(Np);
b = zeros(Np, 1);

C3 = alpha(1, 2)/(2*h);
A(1, 1) = alpha(1, 1)*(s_j - s_jm1) - 3*C3;
A(1, 2) = 4*C3;
A(1, 3) = -C3;
b(1) = g_jm1*(s_j - s_jm1);
for i = 2:Np-1
    C2 = kappa*(s_j - s_jm1)*( (1 - (i-1)*h)*ds_jm1dt + (i-1)*h*ds_jdt)/(2*h);
    A(i, i - 1) = 1/h/h - C2;
    A(i, i) = -(2/h/h + C1);
    A(i, i + 1) = 1/h/h + C2;
    b(i) = -C1*u_past(i);
end
C4 = alpha(2, 2)/(2*h);
A(end, end - 2) = C4;
A(end, end - 1) = -4*C4;
A(end, end) = 3*C4 + alpha(2, 1)*(s_j - s_jm1);
b(end) = g_j*(s_j - s_jm1);
end

