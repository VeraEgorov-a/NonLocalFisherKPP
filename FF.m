% Front-fixing method for non local Fisher-KPP problem
% 2-step transformation + explicit FDM
% Vera Egorova. september, 2024
function [X,W,H,G,t] = FF(scheme, M,T, mu, h0, u0, alpha2, f, J, K,dt)
% Inputs:
% scheme: Numerical scheme indicator (0 = Lax, 1 = Forward-Time Backward-Space (FTB-FS), 2 = Forward-Time Centered-Space (FTCS)).
% M: Number of spatial grid points.
% T: Final simulation time.
% mu: Parameter controlling the boundary dynamics.
% h0: Initial length of the domain.
% u0: Function handle for the initial population distribution \( u_0(x) \).
% alpha2: Diffusion coefficient.
% f: Function handle for the growth term \( f(u) \).
% J: Kernel function for nonlocal diffusion.
% K: Cumulative distribution function of the kernel.
% dt(optional): Time step size (default is computed based on stability criteria).
% 
% Outputs:
% X: Evolving spatial grid at each time step.
% W: Solution matrix, where each row corresponds to the population density at a specific time.
% H: Upper boundary values \( h(t) \).
% G: Lower boundary values \( g(t) \).
% t: Time vector.
h = 1/M; % spatial step-size
if nargin < 11
    k = h/mu;
else
    k = dt;
end
t = 0:k:T;
N = length(t);
r = linspace(0,1,M+1); % spatial grid

ht  = h0; gt = -h0;
lt = ht-gt;

%inverse transformation
x = r*lt+gt;
w = u0(x);

% initial conditions
L = zeros(N,1);
G = zeros(N,1);
W = zeros(N,M+1);
X = zeros(N,M+1);
L(1) = lt; G(1) = gt;
W(1,:) = w;
X(1,:) = x;


quadr = @(w) h/3*(w(1)+  4*sum(w(2:2:end-1)) + 2*sum(w(3:2:end-2)) +w(end));
fprintf('FF starts: dr = %.2e, dt = %.2e, T = %d \n\n',h,k, T);
Start = cputime;
for n = 1:N-1
    % interior nodes
    I1 = quadr(w.*(1 - K(lt*r)));
    I2 = quadr(w.*K(lt*(r-1)));
    adv_coef =-mu*((1-r)*I1 - r*I2);
    
    qn = w.*(1- K(lt*r) + K(lt*(r-1)));
    pn = w.*(1-K(lt*r));
    adv_coef1 = mu*(r*quadr(qn) - quadr(pn));
    s = (adv_coef>=0);
    C = adv_coef1.*k/h;

    for j = 2:M
        I3 = quadr(w.*J(lt*(r(j)-r)));
        switch scheme
            case 0 % Lax
                W(n+1,j) = w(j) +0.5*C(j)*(w(j+1)-w(j-1))+0.5*C(j)^2*(w(j-1)-2*w(j)+w(j+1))...
                    +k*alpha2*(lt*I3 - w(j))+k*f(w(j));
                %+0.5*k*(f(w(j+1)) + f(w(j-1)));                 
                     
            case 1 % FTB-FS
                W(n+1,j) = w(j) + k*adv_coef(j)*...
                    (s(j)*(w(j+1)-w(j))/h+... % upwind, if coef>0
                    (1-s(j))*(w(j)-w(j-1))/h)+...%downwind, if coef<0
                    k*alpha2*(lt*I3 - w(j))+k*f(w(j));
            case 2 % FTCS
                W(n+1,j) = w(j) + k*adv_coef(j)*(w(j+1)-w(j-1))/(2*h)+... % central
                    k*alpha2*(lt*I3 - w(j))+k*f(w(j));
        end
    end
    % compute l^{n+1}:
    lt = lt*(1 + k*mu*quadr(w.*(1-K(lt*r) + K(lt*(r-1)))));
    gt = gt-k*mu*lt*( quadr(w.*(1-K(lt*r))));
    L(n+1) = lt;
    G(n+1) = gt;
    w = W(n+1,:);
    X(n+1,:) = r*lt+gt;
end
H = G+L;
Finish = cputime - Start;
fprintf("done!, T = %d, CPU time: %.2f s\n", T, Finish)
end



