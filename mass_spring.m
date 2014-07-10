Ti = 0; % start time
Tf = 10; % final time 10 seconds
m1 = 4;
m2 = 1;
m3 = 3;
Er = (4/3)*1000;

k1 = 12;
k2 = 12;
k3 = 12;
gamma = 0.2;

delt = .01; % step size
period = pi;

F = @(time) 10*sin(period*time);

n = (Tf - Ti)/delt + 1;

x = zeros(n,4);  
X_temp = zeros(1,3);

x(1, :) = [0,0,0, gamma];
x(2, :) = [0,0,0, gamma];  % boundary conditions, x2 is found by x2 = x1 +delt*x1'


% start solving using 2nd order centered differences, explicit solve for F,
% implicit solve for FC

for i = 2:(n-1)
       
       DE = @(x_next) [m1*(x_next(1) - 2*x(i,1) + x(i-1,1)) + (delt^2)*(k1*x(i,1)-k2*(x(i,2)-x(i,1)));
                            m2*(x_next(2) - 2*x(i,2) + x(i-1,2)) + (delt^2)*(k2*(x(i,2)-x(i,1))-F((i-1)*delt) + fcol(x_next(2), x_next(3) , gamma, Er));
                            m3*(x_next(3) - 2*x(i,3) + x(i-1,3)) + (delt^2)*(k3*x(i,3) + fcol(x_next(2), x_next(3) , gamma, Er))];
                        
        options = optimset('TolFun',10^(-10), 'display', 'off');
       x(i+1,1:3) = fsolve(DE, [x(i,1), x(i,2), x(i,3)], options);
       
       x(i+1, 4) = (gamma - x(i+1,2) - x(i+1,3));
end


time = linspace(0,10,n);

figure(1);
plot(time, x(:,1), time, x(:,2), time, x(:,3));
legend('X_1', 'X_2', 'X_3');
xlabel('time');
ylabel('position');

figure(2);
plot(time, x(:,4));
legend('distance between mass 2 and 3');
xlabel('time');
ylabel('distance');
       

                        
        
    
