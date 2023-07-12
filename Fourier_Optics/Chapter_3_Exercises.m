n%Chapter 3
%Problem 1
%a) 
   % function[out]=tri(x)
   %   t=1-abs(x);
   %   mask=abs(x)<=1;
   %   out=t.*mask;
   % end
%b)
  b=0.1;
  l=2;
  M=200;
%c)
  figure(1);
  x=-l/2:b:l/2-b;
  y=tri(x);
  plot(x,y);
  xlabel('x(m)');ylabel('Intensity');
%d) Computational fourier transform of a triangle function is sinc.
  fx=-1/(2*b):1/l:1/(2*b)-1/l;
  f0=fftshift(y);
  F0=fft(f0)*b;
  F=fftshift(F0);
%e)
  F_an=sinc(fx).^2;
%f)
  figure(2);
  hold on;
  plot(fx,abs(F),fx,abs(F_an),':');
  xlabel('x(m)');ylabel('Intensity');
%Problem 2:
%a)
  wa=0.3;
  wb=0.2;
  L=2;
  M=200;
  dx=L/M;
  x=[-L/2:dx:L/2-dx];
  fa=exp(-pi*(x.^2)/wa^2);
  fb=exp(-pi*(x.^2)/wb^2);

  fx=[-1/(2*dx):1/L:1/(2*dx)-1/L];
  Fa_an=wa*wb/sqrt(wa^2+wb^2)*exp(-pi*x.^2/(wa^2+wb^2));
  figure(3);
  plot(x,fa,x,fb,'--');
  xlabel('x(m)');ylabel('Intensity');legend('Gaussian A','Gaussian B')

  Fa=fft(fa);
  Fb=fft(fb);
  F0=Fa.*Fb;
  f0=ifft(F0)*dx;
  f=fftshift(f0);
  
  figure(4);
  plot(x,f,x,abs(Fa_an));title('Convolution');
  xlabel('x(m)');ylabel('Intensity'); legend('Computational', 'Analytic');
%Problem 3:
%a) 
%   function[out]=circ(r);
%     out=abs(r)<=1;
%   end
%b) 
 r=0.015;
 l=0.2;
 M=200;
 dx=l/M;
 x=-l/2:dx:l/2-dx;
 [X,Y]=meshgrid(x,x);
 f=circ(sqrt(X.^2+Y.^2)/r);
 figure(5);
 imagesc(x,x,nthroot(f,3));
 colormap('gray'); axis square; axis xy; xlabel('x(m)');ylabel('y(m)');
