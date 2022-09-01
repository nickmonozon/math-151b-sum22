%% Runge-Kutta Method for Systems of 2 Differential Equations

%% Input information
a = 0;          % left endpoint
b = 1;          % right endpoint
h = 0.2;        % stepsize
N = (b-a)/h;    % number of subintervals
alpha1 = 1;     % initial condition 1
alpha2 = 1;     % initial condition 2

f1 = @(t,u1,u2) 3*u1 + 2*u2 - (2*t^2 + 1)*exp(2*t);
f2 = @(t,u1,u2) 4*u1 + u2 + (t^2 + 2*t - 4)*exp(2*t);

% Exact solutions
u1 = @(t) 1/3*exp(5*t) - 1/3*exp(-t) + exp(2*t);
u2 = @(t) 1/3*exp(5*t) + 2/3*exp(-t) + t^2*exp(2*t);

%% Performing the first method

% Starting values
t = a;
w1 = alpha1;
w2 = alpha2;

% Output header and starting iteration
fprintf('t_i \t w_1i \t\t u_1i \t\t |w_1i - u_1i| \t w_2i \t\t u_2i \t\t |w_2i - u_2i| \n')            
fprintf('%.1f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \n',t,w1,u1(t), ...
    abs(w1-u1(t)),w2,u2(t), abs(w2-u2(t)))

for i=1:N

    k(1,1) = h * f1(t, w1, w2);
    k(1,2) = h * f2(t, w1, w2);

    k(2,1) = h * f1(t + h/2, w1 + k(1,1)/2, w2 + k(1,2)/2);
    k(2,2) = h * f2(t + h/2, w1 + k(1,1)/2, w2 + k(1,2)/2);

    k(3,1) = h * f1(t + h/2, w1 + k(2,1)/2, w2 + k(2,2)/2);
    k(3,2) = h * f2(t + h/2, w1 + k(2,1)/2, w2 + k(2,2)/2);

    k(4,1) = h * f1(t + h, w1 + k(3,1), w2 + k(3,2));
    k(4,2) = h * f2(t + h, w1 + k(3,1), w2 + k(3,2));

    w1 = w1 + (k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))/6;
    w2 = w2 + (k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))/6;

    t = a + i*h;

    fprintf('%.1f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \n',t,w1,u1(t), ...
    abs(w1-u1(t)),w2,u2(t), abs(w2-u2(t)))
end