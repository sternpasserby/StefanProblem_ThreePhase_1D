function [A, b] = getSysMat(u_past, kappa, tau, h, s_j, s_jm1, ...
    ds_jdt, ds_jm1dt)
%GETSYSTEMMATRIX Вычисление матрицы системы для решения трёхфазной задачи
%Стефана методом конечных разностей с применением метода выпрямления
%фронта. Граничные условия не учитываются.
%   u_past - распределение температур на предыдущем временном шаге

C1 = kappa*(s_j - s_jm1)^2/tau;
Np = length(u_past);
A = zeros(Np);
b = zeros(Np, 1);
for i = 2:Np-1
    C2 = kappa*(s_j - s_jm1)*( (1 - (i-1)*h)*ds_jm1dt + (i-1)*h*ds_jdt)/(2*h);
    A(i, i) = -(2/h/h + C1);
    A(i, i - 1) = 1/h/h - C2;
    A(i, i + 1) = 1/h/h + C2;
    b(i) = -C1*u_past(i);
end

end

