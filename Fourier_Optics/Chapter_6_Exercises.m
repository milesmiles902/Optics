%Chapter 6
%Problem #1:
%Equation 6.6 - alpha<=lambda(1/2dx-B1)
%ztan(alpha)=(L-D2)/2=(L-2B1*lambda*z)/2
%Assume small angle approximation: z*alpha=(L-2B1*lambda*z)/2

%Problem #2:
    z=2000;
    alpha=5.0e-10;
    zf=1002;
    deg=pi/180;
    alpha=5.0e-5;
    theta=45*deg;
    
    L1=0.5;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    y1=x1;
    
    lambda=0.5*10^-6;
    k=2*pi/lambda;
    w=0.051;
    z=2000;
    
    [X1,Y1]=meshgrid(x1,y1);
    u1=rect(X1/(2*w)).*rect(Y1/(2*w));
    %[u1]=tilt(u1,L1,lambda,alpha,theta);
    [u1]=focus(u1,L1,lambda,zf)
    I1=abs(u1.^2);
    
    figure(1)
    imagesc(x1,y1,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('z=0m')
    
    u2=propTF(u1,L1,lambda,z);
    
    x2=x1;
    y2=y1;
    I2=abs(u2).^2;
    
    figure(2);
    imagesc(x2,y2,I2);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title(['z= ',num2str(z),' m'])
    
    figure(3);
    plot(x2,I2(M/2+1,:));
    xlabel('x (m)'); ylabel("irradiance");
    title(['z= ',num2str(z),' m'])
    
    figure(4)
    plot(x2,abs(u2(M/2+1,:)));
    xlabel('x(m)');ylabel('Magnitude');
    title(['z= ', num2str(z), ' m']);
    
    figure(5);
    plot(x2,unwrap(angle(u2(M/2+1,:))));
    xlabel('x (m)'); ylabel('Phase (rad)');
    title(['z= ',num2str(z),' m'])

%Problem #3:
%a/b) 
    z=0.25;
    zf=2000;
    
    L1=0.025;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    y1=x1;
    
    lambda=0.5*10^-6;
    k=2*pi/lambda;
    w=0.0125;
    lz=L1*z;
    
    [X1,Y1]=meshgrid(x1,y1);
    u1=(w^2/lz).*(jinc(w/lz*sqrt(X1.^2+Y1.^2)));
    [u1]=focus(u1,L1,lambda,zf)
    I1=abs(u1.^2);
    
    figure(6)
    imagesc(x1,y1,nthroot(I1,3));
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('z=0m')
    
    u2=propFF(u1,L1,lambda,z);
    
    x2=x1;
    y2=y1;
    I2=abs(u2).^2;
    
    figure(7);
    imagesc(x2,y2,I2);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title(['z= ',num2str(z),' m'])
    
    figure(8);
    plot(x2,I2(M/2+1,:));
    xlabel('x (m)'); ylabel("irradiance");
    title(['z= ',num2str(z),' m'])
    
    figure(9)
    plot(x2,abs(u2(M/2+1,:)));
    xlabel('x(m)');ylabel('Magnitude');
    title(['z= ', num2str(z), ' m']);
    
    figure(10);
    plot(x2,unwrap(angle(u2(M/2+1,:))));
    xlabel('x (m)'); ylabel('Phase (rad)');
    title(['z= ',num2str(z),' m'])
%c) The irradiance in the focal plane of a cylindrical lens is possibly described by a Fraunhofer pattern, but the coefficient of equation (6.18) requires the dimensions to be 2-D.

%Problem 4: Inconsistent with propFF, but does produce a similar irradiance function.
    %a) Generate the Airy pattern of Figure 6.7
    L=0.25;               %Side length(m)
    M=250;                 %Samples
    dx=L/M;                %Sample Interval
    x=-L/2:dx:L/2-dx; y=x; %Coord
    [X,Y]=meshgrid(x,y);
    
    w=1e-1;                %x half-width
    lambda=0.5e-6;        %Wavelength
    z=10;                    %Propagation Distance
    k=2*pi/lambda;          %Wavenumber
    lz=lambda*z;
    
    %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)
    r=sqrt(X.^2+Y.^2)
    I2=(1/lz)^2.*((w^2.*jinc(2*pi*w.*r/lz)./(w.*r/lz))).^2;
    figure(11);
    imagesc(x,y,nthroot(I2,10));
    xlabel('x(m)');ylabel('y(m)');
    colormap('gray');
    axis square;
    axis xy;
%Problem 5
%ta(r)=1/2[1+cos(k/2f*r^2)]circ(r/w);
%a) ta(x^2+y^2)=1/2[1+cos(k/2f*(x^2+y^2)]circ(sqrt(x^2+y^2)/w); The transmitted
% field function will periodicially oscillate about 1/2 in the x-and-y
% directions.
%b) Fraunhofer Criterion:
%   Sample Interval: dF<=1/(lz2|f_max|)
%   Spatial Domain Criterion: dx>lz/L
%   Oversampling: dx>lz/L;
%   Fresnel Criterion:
%   Sample Interval: dF<=1/(lz2|f_max|)
%   Undersampled Impulse Response: dx>lz/L
%   Oversampled Impulse response: dx<lz/L
%   The sampling criterion for this plate would be dF<=1/(lz2|1.5|)
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=10;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+cos(k/(2*f).*(X.^2+Y.^2))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(12)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=10m')
%c) The sampling criterion for this plate is described as critically sampled being
%   dx=L/M=10^-4; lz/L=10^-4
%d) The propagation sampling at z=10 m would be equivalent to f=10
%descriptions in Chapter 6 of Computational Fourier Optics. When the
%relationship of dx=L/M meets an equality, then there are no differences
%between the transfer and impulse response functions.
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=5;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+cos(k/(2*f).*(X.^2+Y.^2))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(13)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=5m')
    
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=8;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+cos(k/(2*f).*(X.^2+Y.^2))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(14)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=8m')
    
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=10;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+cos(k/(2*f).*(X.^2+Y.^2))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(15)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=10m')
    
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=12;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+cos(k/(2*f).*(X.^2+Y.^2))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(16)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=12m')

%f)
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=12;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+sign(cos(k/(2*f).*(X.^2+Y.^2)))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(17)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=12m')
%g)
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=10;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+cos(k/(2*f).*(X.^2+Y.^2))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(18)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=10m')
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=10;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+sign(cos(k/(2*f).*(X.^2+Y.^2)))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(19)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=10m')
% Efficiency is described as the ratio of output/input. The zone plate with
% sign() and without sign() shows equal amounts of efficiency for their periodically
% oscillating functions.
%h) 
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=10;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    u1t=u1.*1/2.*(1+sign(cos(k/(2*f).*(X.^2+Y.^2)))).*circ(sqrt(X.^2+Y.^2)/w)
    I1=abs(u1t.^2);
    
    figure(20)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=10m')
    
    L=50/1000;
    M=500'
    lambda=0.5/10^6;
    k=2*pi/lambda;
    w=6.25/1000;
    f=10;
    dx=L/M;
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    r=sqrt(X.^2+Y.^2);
    
    u1=circ(r/w);
    [u1f]=focus(u1,L,lambda,f)
    I1=abs(u1f.^2);
    
    figure(21)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray'); xlabel('x(m)');ylabel('y(m)');
    title('f=10m')

%Problem 6:
  lambda=0.5e-6;
  f=0.5;
  P=1e-4;
  D1=1e-3;
  L1=1e-2;
  M=1000;
  dx1=L1/M;
  x1=-L1/2:dx1:L1/2-dx1;
  [X1,Y1]=meshgrid(x1,x1);

  %construct
  fc=fft(fftshift(ucomb(x1/P)));
  fr=fft(fftshift(rect(x1/(P/2))));
  ux=ifftshift(ifft(fc.*fr));
  u1=repmat(ux,M,1);
  u1=u1.*rect(X1/D1).*rect(Y1/D1);

  [u2,L2]=propFF(u1,L1,lambda,f);
  dx2=L2/M;
  x2=-L/2:dx2:L2/2; y2=x2;
  I2=abs(u2).^2;

  figure(22);
  imagesc(x2,y2,I2);
  axis square; axis xy;
  colormap('gray'); xlabel('x(m)');ylabel('y(m)');
  title('f=0.5 m');

  M=1001;
  dx1=L1/M;
  x1=-L1/2:dx1:L1/2-dx1;
  [X1,Y1]=meshgrid(x1,x1);

  %construct
  fc=fft(fftshift(ucomb(x1/P)));
  fr=fft(fftshift(rect(x1/(P/2))));
  ux=ifftshift(ifft(fc.*fr));
  u1=repmat(ux,M,1);
  u1=u1.*rect(X1/D1).*rect(Y1/D1);

  [u2,L2]=propFF(u1,L1,lambda,f);
  dx2=L2/M;
  x2=-L/2:dx2:L2/2; y2=x2;
  I2=abs(u2).^2;

  figure(23);
  imagesc(x2,y2,I2);
  axis square; axis xy;
  colormap('gray'); xlabel('x(m)');ylabel('y(m)');
  title('f=0.5 m');

  M=2000;
  dx1=L1/M;
  x1=-L1/2:dx1:L1/2-dx1;
  [X1,Y1]=meshgrid(x1,x1);

  %construct
  fc=fft(fftshift(ucomb(x1/P)));
  fr=fft(fftshift(rect(x1/(P/2))));
  ux=ifftshift(ifft(fc.*fr));
  u1=repmat(ux,M,1);
  u1=u1.*rect(X1/D1).*rect(Y1/D1);

  [u2,L2]=propFF(u1,L1,lambda,f);
  dx2=L2/M;
  x2=-L/2:dx2:L2/2; y2=x2;
  I2=abs(u2).^2;

  figure(24);
  imagesc(x2,y2,I2);
  axis square; axis xy;
  colormap('gray'); xlabel('x(m)');ylabel('y(m)');
  title('f=0.5 m');

 % The observation plane is effected by changing the number of samples (M).
 % When samples are incremented (+1,+2,...) the resolution is worse, until
 % M=2000, which is double the original base, M=4000, and M=8000m, etc.

%Problem 7:
%a)
    lambda=0.5e-6;
    f=0.5;
    P=1e-4;
    D1=1.02e-3;
    L1=1e-2;
    M=500;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    [X1,Y1]=meshgrid(x1,x1);
    m=pi/4;

    u1=exp(i*m*cos(2*pi.*X1./P)).*rect(X1/D1).*rect(Y1/D1);
    
    %Fraunhofer Pattern
    [u2,L2]=propFF(u1,L1,lambda,f);
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2; y2=x2;
    I2=abs(u2).^2;
    
    figure(25)
    imagesc(x2,y2,nthroot(I2,3))
    xlabel('x(m)');ylabel('y(m)');
    figure(26)
    plot(x2,I2(M/2+1,:));
    xlabel('x(m)');ylabel('Irradiance');
    
    [X2,Y2]=meshgrid(x2,y2);
    lf=lambda*f;
    u2a=(1/lf)*D1^2/2*sinc(D1/lf*Y2).*(sinc(D1/lf*X2)-1/2*sinc(D1/lf*(X2+lf/P))-1/2*sinc(D1/lf*(X2-lf/P)));
    I2a=abs(u2a).^2;
    
    figure(27);
    imagesc(x2,y2,nthroot(I2a,3))
    xlabel('x(m)');ylabel('y(m)');
    figure(28);
    plot(x2,I2a(M/2+1,:));
    xlabel('x(m)');ylabel('y(m)');
%b)
    lambda=0.5e-6;
    f=0.5;
    P=1e-4;
    D1=1.02e-3;
    L1=1e-2;
    M=500;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    [X1,Y1]=meshgrid(x1,x1);
    m=pi/4;

    u1=exp(i*m.*rect(X1./(P/2)).*ucomb(X1./P)./P).*rect(X1/D1).*rect(Y1/D1);
    
    %Fraunhofer Pattern
    [u2,L2]=propFF(u1,L1,lambda,f);
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2; y2=x2;
    I2=abs(u2).^2;
    
    figure(29)
    imagesc(x2,y2,nthroot(I2,3))
    xlabel('x(m)');ylabel('y(m)');
    figure(30)
    plot(x2,I2(M/2+1,:))
    xlabel('x(m)');ylabel('Irradiance');
    [X2,Y2]=meshgrid(x2,y2);
    lf=lambda*f;
    u2a=(1/lf)*D1^2/2*sinc(D1/lf*Y2).*(sinc(D1/lf*X2)-1/2*sinc(D1/lf*(X2+lf/P))-1/2*sinc(D1/lf*(X2-lf/P)));
    I2a=abs(u2a).^2;
    
    figure(31);
    imagesc(x2,y2,nthroot(I2a,3))
    xlabel('x(m)');ylabel('y(m)');
    figure(32);
    plot(x2,I2a(M/2+1,:));
    xlabel('x(m)');ylabel('Irradiance (W/m)');
%c)
    lambda=0.5e-6;
    f=0.5;
    P=1e-4;
    D1=1.02e-3;
    L1=1e-2;
    M=500;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    [X1,Y1]=meshgrid(x1,x1);
    m=pi/1;

    u1=exp(i*m.*udelta(X1./(P/2)).*ucomb(X1./P)./P).*rect(X1/D1).*rect(Y1/D1);
    
    %Fraunhofer Pattern
    [u2,L2]=propFF(u1,L1,lambda,f);
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2; y2=x2;
    I2=abs(u2).^2;
    
    figure(33)
    imagesc(x2,y2,nthroot(I2,3))
    xlabel('x(m)');ylabel('y(m)');
    figure(34)
    plot(x2,I2(M/2+1,:))
    xlabel('x(m)');ylabel('Irradiance (W/m)');
    [X2,Y2]=meshgrid(x2,y2);
    lf=lambda*f;
    u2a=(1/lf)*D1^2/2*sinc(D1/lf*Y2).*(sinc(D1/lf*X2)-1/2*sinc(D1/lf*(X2+lf/P))-1/2*sinc(D1/lf*(X2-lf/P)));
    I2a=abs(u2a).^2;
    
    figure(35);
    imagesc(x2,y2,nthroot(I2a,3))
    xlabel('x(m)');ylabel('y(m)');
    figure(36);
    plot(x2,I2a(M/2+1,:));
    xlabel('x(m)');ylabel('Irradiance (W/m)');
% Problem 8: %%%Unsure why the output display does not correspond to answer
% key. 
%a) 
    lambda=0.5e-6;
    f=0.5;
    P=2e-6;
    D1=1e-4;
    L1=0.8e-3;
    M=3200;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;

    u1=udelta(x1./(P/2)).*ucomb(x1./P)./P.*rect(x1/D1);
    
    %Fraunhofer Pattern
    [u2,L2]=propFF(u1,L1,lambda,f);
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2; 
    I2=abs(u2).^2;
    
    figure(37)
    plot(x2,abs(u2));
    xlabel('x(m)');ylabel('Intensity (W/m)');