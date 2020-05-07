function dydt = osc(t,y)
 dydt = zeros(2,1); % this creates an empty column
 ror = 1;
rr = 100;
roa = 100;
ra = 5000;
da = 30;
 %vector that you can fill with your two derivatives:
 dydt(1) = -da*y(1) + (roa+ra*y(1)^2)/(1+y(1)^2+y(2)^2);
 dydt(2) = -y(2) + (ror+rr*y(1)^2)/(1+y(1)^2);
 %In this case, y(1) is y1 and y(2) is y2, and dydt(1)
 %is dy1/dt and dydt(2) is dy2/dt.
end
