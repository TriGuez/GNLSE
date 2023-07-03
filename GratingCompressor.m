function [betas] = GratingCompressor(N, theta_i, lambda)
format long g
c = 299792458;
d = 1/(N*1e3);
theta_i = theta_i*pi/180;
L = 1;
w0 = 2*pi*c/lambda;
m = 1;
theta_d = asin(m*(lambda/d)-sin(theta_i));

beta2 = (-1/c)*(m^2/w0^3)*((2*pi*c)/d)^2*(1/(cos(theta_d)^3));
f = ((1/w0)+((m/w0^2)*((2*pi*c)/d)*(sin(theta_d)/(cos(theta_d)^2))));
beta3 = -3*beta2*f;
beta4 = beta2*(12*f^2+3*(m^2/w0^4)*((2*pi*c)/d)^2*(1/(cos(theta_d)^4)));
betas = [beta2 beta3 beta4];



% beta2 = ((-(-1)^2*lambda*lambda*lambda*L)./(2*pi*c*c*d*d))*(1-((lambda/d)-sin(theta_i))^2)^(-3/2);
