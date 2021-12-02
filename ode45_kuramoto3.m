clc, clf, clear, close all


 
N=5;   %Number of points

K=1;    %Coupling constant

steps = 0.01; %Step of interation

tmin=0;
tmax=30;
allt = tmin:steps:tmax;

 
theta0=2*pi*rand(2*N+1,1); %Random starting point for 2N+1 numbers of point
 
omegaTemp = 0.5*ones(N,1);

Omega = [omegaTemp; -omegaTemp; 1];
 %Random starting omega for 2N+1 numbers of point

Adj = zeros(2*N+1);
for iAdj = 1:2*N+1 %Adjcency matrix for this system
    for jAdj = 1:2*N+1
        if iAdj <= N && jAdj <= N
            Adj(iAdj,jAdj)= 1;
        end
        if iAdj > N && iAdj <= 2*N && jAdj > N && jAdj <= 2*N
            Adj(iAdj,jAdj) = 1;
        end
        if iAdj == N && jAdj == 2*N+1
            Adj(iAdj,jAdj) = 1;
        end
        if jAdj == N && iAdj == 2*N+1
            Adj(iAdj,jAdj) = 1;
        end
        if iAdj == 2*N && jAdj == 2*N+1
            Adj(iAdj,jAdj) = 1;
        end
        if jAdj == 2*N && iAdj == 2*N+1
            Adj(iAdj,jAdj) = 1;
        end
        if iAdj == jAdj
            Adj(iAdj,jAdj) = 0;
        end
    end
end


%LHS function
tspan = [tmin tmax];

% MATLAB can create a new function funname by going
%funname = @(x1,x2,x3) x1 + x2*x3;

ode_kuramoto_for_solver= @(t,theta) ode_kuramoto(t,theta,K,N,Omega,Adj);
sol = ode45(ode_kuramoto_for_solver,tspan,theta0); %ode45 function
alltheta = deval(sol,allt);



 %plotgraph function

        for i=1:numel(allt)
         s=linspace(0,2*pi,100);%Plot original circles
         originalx=cos(s);
         originaly=sin(s);
 
         xaxis=cos(alltheta(1:N-1,i));%Plot the moving points 
         yaxis=sin(alltheta(1:N-1,i));%in axis, from 1 to N
         
         xaxis2=cos(alltheta(N,i));%Plot the moving points 
         yaxis2=sin(alltheta(N,i));%in axis, node N
         
         xaxis3=cos(alltheta(N+1:2*N-1,i));%Plot the moving points 
         yaxis3=sin(alltheta(N+1:2*N-1,i));%in axis, from N+1 to 2N
                  
         xaxis4=cos(alltheta(2*N,i));%Plot the moving points 
         yaxis4=sin(alltheta(2*N,i));%in axis, node 2N
         
         xaxis5=cos(alltheta(2*N+1,i));%Plot the moving points 
         yaxis5=sin(alltheta(2*N+1,i));%in axis, from 1 to 2N+1
         
         clf
         hold on
         plot(originalx,originaly,'k');
         plot(xaxis,yaxis,'p','color','r'); %1 to N
         
         plot(xaxis2,yaxis2,'+','color','#D95319'); %N

         plot(xaxis3,yaxis3,'o','color','b'); %N+1 to 2N

         plot(xaxis4,yaxis4,'x','color','#0072BD'); %2N

         plot(xaxis5,yaxis5,'d','color','k'); %1 to 2N+1
         
         axis([-3 3 -3 3]);
         drawnow limitrate

        end




function dthetadt = ode_kuramoto(~,theta,K,N,omega,Adj)  %Original kuramoto function (RHS function)
                       

             
     dthetadt = zeros(size(theta));
     %    dthetadt = ode(t,theta);
 
     for j=1:2*N+1  %For node 1 to (2N+1) using adjcency matrix
         dthetadt(j) = omega(j);
         for innerj = 1:2*N+1
            dthetadt(j) = dthetadt(j) +  K/N*Adj(j,innerj)*sin(theta(innerj)-theta(j));
         end
     end

end
