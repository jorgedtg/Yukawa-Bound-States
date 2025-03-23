%% PARAMETERS, VARIABLES, AND SHOOTING METHOD TO FIND NUMBER OF STATES AND ENERGIES
clear;clc;close all;


dx = 0.01; x0 = 1e-5; xf = 500;
X = x0:dx:xf;   nx = length(X);

hb = 1;  m = 1; 

a = 0.3; lambda = 1;    % INPUT

l = 0;

V = @(x) -lambda*exp(-a*x)./x;


E = zeros(1, 100);
en = 1;
de = 1;
it = 1;

figure;
for n = 0:length(E)
    if n == 0
        E0 = -30;
        de = 1;
    else
        E0 = -order(E(n));
        de = 1;
    end
    dig = 0;
    it = 0;

    Veff = @(x) V(x) + l*(l+1)*hb^2./(2*m*x.^2);

    while dig < 10
        h = dx;
        
        U = zeros(1, nx);
        z = zeros(1, nx);

        U(end) = 0;
        z(end) = -1e-100;

        
        q = @(x) (2*m/hb^2)*(E0 - Veff(x));
        fy = @(x, u, z) z;
        fz = @(x, u, z) -q(x)*u;
        
        
        broke = false;
        % figure;
        for i = nx:-1:2

            k1y=fy(X(i), U(i), z(i));
            k1z=fz(X(i), U(i), z(i));
            k2y=fy(X(i)-h/2, U(i)-h/2*k1y, z(i)-h/2*k1z);
            k2z=fz(X(i)-h/2, U(i)-h/2*k1y, z(i)-h/2*k1z);
            k3y=fy(X(i)-h/2, U(i)-h/2*k2y,z(i)-h/2*k2z);
            k3z=fz(X(i)-h/2, U(i)-h/2*k2y,z(i)-h/2*k2z);
            k4y=fy(X(i)-h, U(i)-h*k3y, z(i)-h*k3z);
            k4z=fz(X(i)-h, U(i)-h*k3y, z(i)-h*k3z);
            U(i-1)=U(i)-h/6*(k1y+2*k2y+2*k3y+k4y);
            z(i-1)=z(i)-h/6*(k1z+2*k2z+2*k3z+k4z);
            if U(i-1) > 1e10
                
                CU = 1/sqrt(sum(abs(U).^2*dx));
                U = CU*U;
                Cz = 1/sqrt(sum(abs(z).^2*dx));
                z = Cz*z;
            
            end
        end
        
        % plot(X, U)
        % title("n = " + string(n) + "        E = "+  string(E0))
        % xlim([0, 40])
        % xlabel(string(it))
        % drawnow
        % pause(1/60)
       
        rts = raiz(U);

        if rts ~= n && ~broke
            E0 = E0 - de + de*0.1;
            de = de*0.1;
            dig = dig + 1;
        else
            E0 = E0 + de;
        end

        it  = it + 1;
    end
    if E0 > -1e-6
        break
    end
    E(n+1) = E0;
end














function nE0 = order(E0)

E0 = abs(E0);
ord = 0;

while E0 < 1

    E0 = E0*10;
    ord = ord + 1;

end

E0 = E0 - mod(E0, 1);

nE0 = E0/(10^ord);

end





function zeros = raiz(y)

zeros = 0;

for i = 1:length(y)-1

    if sign(y(i+1)) == 0
        continue
    end

    if sign(y(i)) == 0
        if i == 1
            continue
        end
        if sign(y(i-1)) ~= sign(y(i+1))
            zeros = zeros + 1;

        end
        continue

    end

    if sign(y(i)) ~= sign(y(i+1))
        % y(i)
        % y(i+1)
        zeros = zeros + 1;
        % pause(3)
        
    end

end
end

