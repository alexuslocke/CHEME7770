function [T,Y] = call_osc()
 tspan = [0 30];
 y1_0 = 1;
 y2_0 = 10;
 [T,Y] = ode15s(@osc,tspan,[y1_0 y2_0]);
 plot(T,Y(:,1),'o',T,Y(:,2),'-o')
 xlabel('time'); ylabel('concetration');
end
