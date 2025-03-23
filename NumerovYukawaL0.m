%% FIND STATES AND ENERGIES W/NUMEROV FOR L=0
%% 1. FORMATO
set(groot,'defaultAxesFontSize',18)                    
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% 2. PARAMETERS, VARIABLES, AND SHOOTING METHOD TO FIND NUMBER OF STATES AND ENERGIES
clear;clc;close all;
% Space discretization
dr = 0.01;

% Limits
r0 = 0; rf = 500;
R = r0:dr:rf; nr = length(R);

% Parameters relevant to the problem
    % {In this next lines, I declared the vector containing the lambda and 
    % alpha vectors to be used when obtaining the number of states. The
    % algorithm will find all states and energies associated with all of 
    % the values contained in this vectos, so for example if you write
    % lambda_vec = 3:5
    % a_vec = 0.5:0.1:1
    % the algorithm will save the number of states and values of energy 
    % of each combination of these values.
    % lambda = 3, a = 0.5
    % lambda = 3, a = 0.6
    % ...
    % lambda = 4, a = 0.5
    % ...}% 

lambda_vec = 2:3;    % Coupling constant
a_vec = 0.5:0.1:0.6;       % Screening parameter
l = 0;             % Angular momentum   [DO NOT CHANGE]
    
% Vectors lenght
na = length(a_vec);     nlambda = length(lambda_vec);

% Constants hbar and mass
h = 1;  m = 1;  

% Matrices used to store the values and number of energies
E = zeros(na, nlambda, 10);
EN = zeros(na, nlambda);

% Shooting method to find En(lambda, a, l)
tic     % Start timer

% Array to store values of the critical values of alpha
acrit = zeros(size(lambda_vec));
for lambda = lambda_vec

    for a = a_vec

        % Yukawa potential
        V = @(r) -lambda*exp(-a*r)./r;
        En = 0;     % Number of states
        
        % Quantum number b
        for n = 0:length(E)
            
            % Starting parameters for the numerical integration
            if n == 0
                E0 = -30;
                de = 1;
                enb = E0-1;
            else
                E0 = -order(E(a_vec==a, lambda_vec==lambda, n));
                de = 1;
                enb =E(a_vec==a, lambda_vec==lambda, n);
            end

            % Stop condition, U(0) = 0;
            dig = 0;

            U0 = 100;
            it = 0;
            itLim = 250; % Maximum number of iterations
            
            while abs(U0)>1e-10
                
                % Effective Potential
                Veff = @(r) V(r) + l*(l+1)*h^2./(2*m*r.^2);

                % q Function for Numerov's Method
                q = @(r) (2*m/h^2)*(E0 - Veff(r));
                
                % Wave function
                U = zeros(1, nr);

                % Boundary conditions
                U(end) = 0;
                U(end-1) = 1e-4;
                
                % Reverse Numerov Integration
                for i = nr-2:-1:1
                    ri = R(i);      % r_{i}
                    ri1 = R(i+1);   % r_{i+1}
                    ri2 = R(i+2);   % r_{i+2}
                    
                    Ci2 = (1+(dr^2/12)*q(ri2))*U(i+2);    % C_{i+2}   
                    Ci1 = (2-(5*dr^2/6)*q(ri1))*U(i+1);   % C_{i+1}
                    Ci =  (1+(dr^2/12)*q(ri));           % C_{i}
                    
                    U(i) = (Ci1-Ci2)/Ci;    % Wave function value at r_i
                    if i <= 1 && l == 0     % Approximation for r=0
                        ri3 = R(i+3);
                        U(i) = -2*U(i+1) + U(i+2) +  (1/12)*dr^2*(13*q(ri1)*U(i+1) - 2*q(ri2)*U(i+2) + q(ri3)*U(i+3));
                    % Not used
                    elseif i<=1 && l > 0
                        ri3 = R(i+3);

                        A2 = (2^-(l+1))*(4*l^4 + 28*l^3 - l^2 - 217*l + 240);
                        A3 = (3^-(l+2))*(2*l^4 + 41*l^3 + 232*l^2 + 301*l - 540);
                        B1 = 2*l^2 + 15*l - 5;
                        B2 = (2^(1-l))*(4*l^2 + 24*l- 25);
                        B3 = (3^-(l+1))*(6*l^2 + 9*l - 15);
                        Q = 2*l^4 + 5*l^3 - -68*l^2 + 85*l - 60;
                        a2 = A2/Q;
                        a3 = A3/Q;
                        b1 = B1/Q;
                        b2 = B2/Q;
                        b3 = B3/Q;
                        
                        U(i) = U(i+1) + a2*U(i+2) + a3*U(i+3) + dr^2*(b1*q(ri1)*U(i+1) + b2*q(ri2)*U(i+2) + b3*q(ri3)*U(i+3));
                        % U(i) = -U(i+1) - a2*U(i+2) - a3*U(i+3) + dr^2*(b1*q(ri1)*U(i+1) + b2*q(ri2)*U(i+2) + b3*q(ri3)*U(i+3));
                    end
                    
                    % If wave function gets to big, normalize it
                    if U(i) > 1e100
                        U = U/sqrt(sum(abs(U).^2*dr));
                    end
                end
                Un = U/sqrt(sum(abs(U).^2*dr));
                
                % Uncomment to show a graph of the wave function being
                % integrated

                % plot(R, Un)
                % title("n = " + string(n) + "        E = "+  string(E0))
                % xlabel("$\alpha$ = " + string(a) + "  $\lambda = $" + string(lambda))
                % xlim([-0.1, rf])
                % drawnow
                % pause(1/60)
                
               

                U0 = Un(1);         % Get value of the wave function at r=0
                rts = raiz(U);      % Get number of roots
                it = it+1;  
                eb = E0;            % Energy of the this iteration

                %{Here the value of energy is adjusted. The conditions
                % were set after studying the behavior of the wave function
                % depending if the energy values is bigger or smaller than
                % actual value of energy for the particular quantum number
                % n}%
                if (rts == n && E0 > enb) || rts > n+1 || E0>0 || (E0==0 && rts==n+1)
                    E0 = E0 - de + de*0.1;
                    de = de*0.1;
                    dig = dig + 1;
                else
                    E0 = E0 + de;
                end

                % Values smaller than 1e15 are considered 0
                if E0 > -1e-15
                    E0 = 0;
                end
                if it > itLim
                    break
                end
            end
            
            % If the value of energy is too small, then there are no more
            % bound states
            if abs(E0) <= 1e-15
                break
            end

            % Store value of energy
            E(a_vec==a, lambda_vec==lambda, n+1) = E0; 
            En = En + 1;    % If another state is found, increase number
                            % of states
        end
        EN(a_vec==a, lambda_vec==lambda) = En; % Save number of states
    end
end
toc
%% 3. ENERGIES TABLE
clc
fig = uifigure;
uit = uitable(fig, "Data",EN);
set(uit, 'columnname', "𝜆 = " + string(lambda_vec))
set(uit, 'rowname', "𝛼 = " + string(a_vec))
s = uistyle("Interpreter", "latex");
addStyle(uit, s);
%% ENERGIES VALUES

alpha = 0.5;
lambda = 3;
ENERGIES = getEnergies(alpha, lambda, E, a_vec, lambda_vec)
















% Functions
function UM = fdo(a, lambda, E, l)

dr = 0.01;
r0 = 0; rf = 400;
R = r0:dr:rf; nr = length(R);

nE = length(E);

UM = zeros(nE, nr);

V = @(r) -lambda*exp(-a*r)./r;
m = 1; h = 1;

for n = 1:nE
    
    % Potencial efectivo
    Veff = @(r) V(r) + l*(l+1)*h^2./(2*m*r.^2);
    % Función 'q' para Numerov
    q = @(r) (2*m/h^2)*(E(n) - Veff(r));

    % Funcion de onda
    U = zeros(1, nr);

    % Condiciones iniciales
    U(end) = 0;
    U(end-1) = 1e-4;

    % Algoritmo de Numerov pa tras
    for i = nr-2:-1:1
        ri = R(i);      % r_{i}
        ri1 = R(i+1);   % r_{i+1}
        ri2 = R(i+2);   % r_{i+2}

        Ci2 = (1+(dr^2/12)*q(ri2))*U(i+2);    % C_{i+2}
        Ci1 = (2-(5*dr^2/6)*q(ri1))*U(i+1);   % C_{i+1}
        Ci =  (1+(dr^2/12)*q(ri));           % C_{i}

        U(i) = (Ci1-Ci2)/Ci;
        if i <= 1 && l == 0
            ri3 = R(i+3);
            U(i) = -2*U(i+1) + U(i+2) +  (1/12)*dr^2*(13*q(ri1)*U(i+1) - 2*q(ri2)*U(i+2) + q(ri3)*U(i+3));
        elseif i<=1 && l > 0
            ri3 = R(i+3);

            A2 = (2^-(l+1))*(4*l^4 + 28*l^3 - l^2 - 217*l + 240);
            A3 = (3^-(l+2))*(2*l^4 + 41*l^3 + 232*l^2 + 301*l - 540);
            B1 = 2*l^2 + 15*l - 5;
            B2 = (2^(1-l))*(4*l^2 + 24*l- 25);
            B3 = (3^-(l+1))*(6*l^2 + 9*l - 15);
            Q = 2*l^4 + 5*l^3 - -68*l^2 + 85*l - 60;
            a2 = A2/Q;
            a3 = A3/Q;
            b1 = B1/Q;
            b2 = B2/Q;
            b3 = B3/Q;

            U(i) = U(i+1) + a2*U(i+2) + a3*U(i+3) + dr^2*(b1*q(ri1)*U(i+1) + b2*q(ri2)*U(i+2) + b3*q(ri3)*U(i+3));
            % U(i) = -U(i+1) - a2*U(i+2) - a3*U(i+3) + dr^2*(b1*q(ri1)*U(i+1) + b2*q(ri2)*U(i+2) + b3*q(ri3)*U(i+3));
            
        end
        

        if abs(U(i)) > 1e10
            U = U/sqrt(sum(abs(U).^2*dr));
        end
        % if  mod(i, 100)==0
        % plot(R, U)
        % drawnow
        % pause(1/144)
        % end
    end
    
    

    nrm = 1/sqrt(sum(abs(U.^2)*dr));

    Un = nrm*U;
    
    
    UM(n, :) = Un;
end
end

function ENERGIAS = getEnergies(alpha, lambda, E, a_vec, lambda_vec)
e1  = E(a_vec==alpha, lambda_vec==lambda, :);
e2 = e1(e1~=0);

if isempty(e2)
    fprintf("This configuration has no bound states.\n")
    ENERGIAS = 0;
    return
end

ENERGIAS = reshape(e2, 1, []);
fprintf("This configuration has "+ string(length(ENERGIAS)) + " bound states.\n")


end

function nE0 = order(E0)

E0 = abs(E0);
ord = 0;

while E0 < 1

    E0 = E0*10;
    ord = ord + 1;

end

E0 = E0 - mod(E0, 1);

nE0 = (E0)/(10^ord);

end


function zeros = raiz(y)

zeros = 0;

for i = 1:length(y)-1

    if sign(y(i+1)) == 0
        continue
    end

    if sign(y(i)) == 0

        if sign(y(i-1)) ~= sign(y(i+1))
            zeros = zeros + 1;

        end
        continue

    end

    if sign(y(i)) ~= sign(y(i+1))
        zeros = zeros + 1;
       
    end

end
end

function dV = derOrd4(V, t)

N = length(V);
dx = abs(t(1)-t(2));

D = zeros(N, N);

val = [1/12, -2/3, 0, 2/3, -1/12];
% val = [-1/2, 0, 1/2];

nv = length(val);
nm = floor(nv/2);

for i = -nm:nm

    D = D + diag(ones(N-abs(i), 1)*val(i+nm+1), i);

end

r1 = zeros(1, length(D));
rf = r1;

r1(1) = -25/12;
r1(2) = 4;
r1(3) = -3;
r1(4) = 4/3;
r1(5) = -1/4;

rf(end) = 25/12;
rf(end-1) = -4;
rf(end-2) = 3;
rf(end-3) = -4/3;
rf(end-4) = 1/4;

r2 = [0, r1];
r2(end) = [];

rf1 = [rf, 0];
rf1(1) = [];

D(1, :) = r1;
D(2, :) = r2;
D(end, :) = rf;
D(end-1, :) = rf1;

dV = (D*V')/dx;

end
