%Chapter 2
%Problem #1:
% dx=10x10^-6; L=5x10^-3; 
% Sample Number:M=L/dx=500
% Nyquist Frequency: f=1/2dx= 5x10^4 cycles/m
% Frequency Sample Interval: df=1/L= 200cycles/m
% Range of coordinates in spatial domain:
%       x=[-L/2:1/L:dx:L/2-dx]=[-2.5*10^-3,2.5*10^-3]
% Range of coordinates in frequency domain:
%       fx=[-1/2dx:1/L:1/L:1/2dx-1/L]=[-5*10^4,5*10^4]

%Problem #2: 
%a) g(x,y)=circ(sqrt(x^2+y^2)/w),w=1 mm;
%   1) Effective Bandwidth: B= 5/w = 5 cycles/mm
%   2) Maximum Sample Interval (dx): dx = w/10 = 0.1 mm
%   3) Side Length with 256 Samples: L = M*dx = 25.6 mm
%b) g(x,y)=exp(-(x^2+y^2)/w^2),w = 1 mm;
%   1) Effective Bandwidth: B= 0.79/w = 5 cycles/mm
%        1 = integral_0^(2π) integral_0^B π w^2 exp(-r^2 π / w^2) r dr dθ
%        B = sqrt((log(π) - log(π - 1))/π) = 0.626 cycles/mm
%   2) Maximum Sample Interval (dx): dx = 1/2B = 0.72 mm
%   3) Side Length with 256 Samples: L = M*dx = 185 mm

%Problem #3:
% g(x,y)=exp(-pi(x^2+y^2)/w^2)
% is equal to g(x,y)=0.01 at 1.21w; Solved graphically.

%Problem 4
%a) Bandwidth of sinc(x/w)*sinc(y/w);
%      B=w; dx=1/2B=1/2w
%b) Bandwidth of sinc^2(x/w)*sinc^2(y/w);
%      B=w/2; dx=1/2B=1/w

%Problem 5:
%a) g1(x,y)=del(x/d)*del(y/d), g2(x,y)=circ(sqrt(x^2+y^2)/d)
%      g1(x,y)*g2(x,y); Side Length: L = M*dx=M/2B=6d
%b     g2(x,y)*g2(x,y); Side Length: L = M*dx=M/2B=4d;
