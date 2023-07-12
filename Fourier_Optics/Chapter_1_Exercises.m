%Chapter 1
%Reference functions (impulse, line shape) are supplied from Computational Fourier Optics.
%Problem #1:
%a) Plot rect(x/2);
 w=0.99;
 L=4;
 M=200;
 dx=L/M;

 x=-L/2:dx:L/2-dx;
 f=rect(x/(2*w));

 figure(1);
 plot(x,f);
 xlim([-10 10 ]);xlabel('x(m)');ylabel('Intensity');
%b) Plot rect(x-2)
 w=0.99;
 L=8;
 M=200;
 dx=L/M;

 x=-L/2:dx:L/2-dx;
 f=rect(x-2);

 figure(2);
 plot(x,f);
 xlim([-10 10]);xlabel('x(m)');ylabel('Intensity');
%c) Plot del((x+2)/2)
 L=400;
 M=200;
 dx=L/M;

 x=-L/2:dx:L/2-dx;
 f=udelta((x+2)/2);

 figure(3);
 plot(x,f);
 xlim([-10 10]);xlabel('x(m)');ylabel('Intensity');
%d) Plot exp(-3pi*x^2)
 L=2;
 M=200;
 dx=L/M;

 x=-L/2:dx:L/2-dx;
 f=exp(-3*pi*x.^2);

 figure(4);
 plot(x,f);
 xlim([-10 10]);xlabel('x(m)');ylabel('Intensity');
 %d) Plot del((x+2)/2)
 L=400;
 M=200;
 dx=L/M;

 x=-L/2:dx:L/2-dx;
 f=udelta((x+2)/2);

 figure(3);
 plot(x,f);
 xlim([-10 10]);xlabel('x(m)');ylabel('Intensity');
%e) Plot exp(-3pi*x^2)
 L=2;
 M=200;
 dx=L/M;

 x=-L/2:dx:L/2-dx;
 f=exp(-3*pi*x.^2);

 figure(4);
 plot(x,f);
 xlim([-10 10]);xlabel('x(m)');ylabel('Intensity');

%f) Plot comb(x/4)*del(x)*rect(x/12) modified without del(x) to fit answer
 L=100;
 M=200;
 dx=L/M;

 x=-L/2:dx:L/2-dx;
 f=ucomb(x/4).*rect(x/12);
 figure(5);
 plot(x,f);
 xlim([-10 10]);xlabel('x(m)');ylabel('Intensity');

%g) Plot circ(sqrt((x-2)^2+y^2))+circ(sqrt((x+2)^2+y^2); 
%   Inconsistency with supplied function

%Problem 2 (Fourier Transfor):
%a) FT(g(x,y)); where g(x,y)=rect(x/2w)*rect(y/2w)
%     = 4w^2*sinc(2wfx)*sinc(2wfy)
%b) FT(g(x,y)); where g(x,y)=rect((x-x0)/2w)*rect(y/2w))
%     =(4w^2)*sinc^2(2wfx)*sinc(2wfy)*exp(-j2pix0fx);
%c) FT(g(x,y)); where g(x,y)=exp(-(x^2+y^2)/w^2)
%     =pi*w^2*exp(-(fx^2+fy^2)*pi*w^2)
%d) FT(g(x,y)); where g(x,y)=circ(sqrt((x-2)^2+y^2))+circ(sqrt((x+2)^2+y^2))
%     =w1*J(w1*2pi*sqrt(fx^2+fy^2))/sqrt(fx^2+fy^2)-w2*J(w2*2pi*sqrt(fx^2+fy^2))/sqrt(fx^2+fy^2)
%e) FT(g(x,y)); where g(x,y)=circ(sqrt((x-d/2)^2+y^2)/w)+circ(sqrt((x+d/s)^2+y^2)/w)
%     =w*J(2*w*pi*sqrt(fx^2+fy^2))*cos(pidfx)/sqrt(fx^2+fy^2)

%Problem 3 (Convolution):
%a) FT(int(int(g(x,y)))); where g(x,y)=rect(x/2w)*rect(y/2w) X rect(x/2w)*rect(y/2w))
%     =4w^2*del(x/2w)*del(y/2w)
%b) FT(int(int(g(x,y)))); where g(x,y)=exp(-pi*(x^2+y^2)/9) X exp(-pi*(x^2+y^2)/16))
%     =144/25*exp(-pi*(x^2+y^2))
%c) FT(int(int(g(x,y)))); where g(x,y)=sinc(x/2)*sinc(y) X sinc(x/4)*sinc(y))
%     =2sinc(x/4)*sinc(y)

%Problem 4 (Autocorrelation):
%a) FT(int(int(rect(x/2w)*rect(y/2w)*rect((x)/2w)*rect((y)/2w))*))
%     =4w^2*del(x/2w)*del(y/2w);
%b) FT(int(int(exp(-pi*(x^2+y^2)/2w^2)*exp(-pi*(x^2+y^2)/2w^2)*)))
%     =w^2/2*exp(-pi*(x^2+y^2)/w^2)

%Problem 5 (Central Ordinate):
%a) int_-inf_inf(int_-inf_inf(g(x,y))); where g(x,y)=rect(x/2w)*rect(y/2w)
%     =8w^2
%b) int_-inf_inf(int_-inf_inf(g(x,y))); where g(x,y)=circ(sqrt(x^2+y^2)/3)
%     =abs(sqrt(x^2+y^2)/3))<=1; r=3; Area=pi*3^2=9*pi

%Problem 6:
%a) S(g(x,y))=Ag(x,y); Linear and Space-invariant
%b) S(g(x,y))=Ag(x,y)+B; Space-Invariant
%c) S(g(x,y))=A(g(x,y))^2; Space-Invariant
%d) S(g(x,y))=xg(x,y); Linear
%e) Ave(g(x,y))=1/ab*int_x-a/2_x+a/2(int_y-b/2_y+b/2(g(zeta,eta)); Linear
