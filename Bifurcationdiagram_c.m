clear; clc;
a = 0.2; b = 0.2; global c; 
crange = 1:0.05:5;       % Range for parameter c
k = 0; tspan = 0:0.1:500; % Time interval for solving Rossler system
xmax = [];                  % A matrix for storing the sorted value of x1
for c = crange 
    f = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)]; 
    x0 = [2 -3 2];                  % initial condition for Rossler system
    k = k + 1; 
    [t,x] = ode45(f,tspan,x0);    % call ode() to solve Rossler system
    count = find(t>100);          % find all the t_values which is >10
    x = x(count,:); 
    j = 1; 
    n = length(x(:,1));      % find the length of vector x1(x in our problem)
    for i=2 : n-1 
        % check for the min value in 1st column of sol matrix
        if (x(i-1,1)+eps) < x(i,1) && x(i,1) > (x(i+1,1)+eps)
            xmax(k,j)=x(i,1); % Sorting the values of x1 in increasing order
            j=j+1; 
        end 
    end
    % generating bifurcation map by plotting j-1 element of kth row each time
    if j>1 
        plot(c,xmax(k,1:j-1),'k.','MarkerSize',1); 
    end 
    hold on; 
    index(k)=j-1; 
end 
xlabel('c'); 
ylabel('x max'); 
title('Bifurcation for c'); 