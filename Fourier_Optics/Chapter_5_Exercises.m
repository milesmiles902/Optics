%Chapter 5: Exercises
%Problem 1a)
  %Side Length dX=lz/L; L=lz/dx; L=0.5;
  %Sample Interval dX=lz/L=1/1000(1/m); 
  %Nyquist Frequency 1/2dX=500(1/s);
%Problem 1b)
  %Source Effective Bandwidth B1=5/w=5/0.05m=100(1/m);
%Problem 1c)
  %Fresnel Number N=w^2/lz=.05m^2/(0.5/10^6*1000)=5;
  %The Fresnel Region is best described by a Fresnel Number less than 1.
  %This is within the Fresnel Region (>>1)
%Problem 1d)
  %Simulation of Fresnel propagation

  z=500;
  L1=0.5;
  M=500;
  dx1=L1/M;
  x1=-L1/2:dx1:L1/2-dx1;
  y1=x1;

  lambda=0.5*10^-6;
  k=2*pi/lambda;
  w=0.05;

  [X1,Y1]=meshgrid(x1,y1);
  u1=rect(X1/(2*w)).*rect(Y1/(2*w));
  I1=abs(u1.^2);

  figure(1);
  imagesc(x1,y1,I1);
  axis square; axis xy;
  colormap('gray'); xlabel('x(m)');ylabel('y(m)');
  title('z=0 m');

  u2=propIR(u1,L1,lambda,z);
  x2=x1;
  y2=y1;
  I2=abs(u2).^2;

  figure(2);
  imagesc(x2,y2,I2);
  axis square; axis xy;
  colormap('gray');xlabel('x(m)');ylabel('y(m)');
  title(['z= ',num2str(z),' m']);
  
  figure(3);
  plot(x2,I2(M/2+1,:));
  xlabel('x (m)');ylabel('Irradiance');
  title(['z= ',num2str(z),' m']);

  figure(4)
  plot(x2,abs(u2(M/2+1,:)));
  xlabel('x (m)');ylabel('Magnitude');
  title(['z= ',num2str(z), ' m']);

  figure(5);
  plot(x2,unwrap(angle(u2(M/2+1,:))));
  xlabel('x (m)');ylabel('Phase(rad)');
  title(['z= ',num2str(z),' m'])

  z=1000;
  L1=0.5;
  M=500;
  dx1=L1/M;
  x1=-L1/2:dx1:L1/2-dx1;
  y1=x1;

  lambda=0.5*10^-6;
  k=2*pi/lambda;
  w=0.05;

  [X1,Y1]=meshgrid(x1,y1);
  u1=rect(X1/(2*w)).*rect(Y1/(2*w));
  I1=abs(u1.^2);

  figure(6);
  imagesc(x1,y1,I1);
  axis square; axis xy;
  colormap('gray'); xlabel('x(m)');ylabel('y(m)');
  title('z=0 m');

  u2=propIR(u1,L1,lambda,z);
  x2=x1;
  y2=y1;
  I2=abs(u2).^2;

  figure(7);
  imagesc(x2,y2,I2);
  axis square; axis xy;
  colormap('gray');xlabel('x(m)');ylabel('y(m)');
  title(['z= ',num2str(z),' m']);
  
  figure(8);
  plot(x2,I2(M/2+1,:));
  xlabel('x (m)');ylabel('Irradiance');
  title(['z= ',num2str(z),' m']);

  figure(9)
  plot(x2,abs(u2(M/2+1,:)));
  xlabel('x (m)');ylabel('Magnitude');
  title(['z= ',num2str(z), ' m']);

  figure(10);
  plot(x2,unwrap(angle(u2(M/2+1,:))));
  xlabel('x (m)');ylabel('Phase(rad)');
  title(['z= ',num2str(z),' m'])

  z=2000;
  L1=0.5;
  M=500;
  dx1=L1/M;
  x1=-L1/2:dx1:L1/2-dx1;
  y1=x1;

  lambda=0.5*10^-6;
  k=2*pi/lambda;
  w=0.05;

  [X1,Y1]=meshgrid(x1,y1);
  u1=rect(X1/(2*w)).*rect(Y1/(2*w));
  I1=abs(u1.^2);

  figure(11);
  imagesc(x1,y1,I1);
  axis square; axis xy;
  colormap('gray'); xlabel('x(m)');ylabel('y(m)');
  title('z=0 m');

  u2=propTF(u1,L1,lambda,z);
  x2=x1;
  y2=y1;
  I2=abs(u2).^2;

  figure(12);
  imagesc(x2,y2,I2);
  axis square; axis xy;
  colormap('gray');xlabel('x(m)');ylabel('y(m)');
  title(['z= ',num2str(z),' m']);
  
  figure(13);
  plot(x2,I2(M/2+1,:));
  xlabel('x (m)');ylabel('Irradiance');
  title(['z= ',num2str(z),' m']);

  figure(14)
  plot(x2,abs(u2(M/2+1,:)));
  xlabel('x (m)');ylabel('Magnitude');
  title(['z= ',num2str(z), ' m']);

  figure(15);
  plot(x2,unwrap(angle(u2(M/2+1,:))));
  xlabel('x (m)');ylabel('Phase(rad)');
  title(['z= ',num2str(z),' m'])


%Problem 2a)
  L=0.002;                 %Side length(m)
  M=500;                 %Samples
  dx=L/M;                %Sample Interval
  x=-L/2:dx:L/2-dx; y=x; %Coord
  [X,Y]=meshgrid(x,y);

  wx=1e-4;                %x half-width
  wy=5e-5;                %y half-width
  lambda=0.633e-6;        %Wavelength
  z=0.005;                  %Propagation Distance
  k=2*pi/lambda;          %Wavenumber
  lz=lambda*z;
  %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)

  u1=rect(X/(2*wx)).*rect(Y/(2*wy));
  u2=propTF(u1,L,lambda,z);
 
  I=abs(u2).^2;
  figure(1);
  imagesc(x,y,nthroot(I,3));
  xlabel('x(m)');ylabel('y(m)');
  colormap('gray');
  axis square;
  axis xy;
  figure(2);
  plot(x,I(M/2+1,:))
  xlabel('x(m)');ylabel("Irradiance")

%Problem 2b)

  L=0.025;               %Side length(m)
  M=500;                 %Samples
  dx=L/M;                %Sample Interval
  x=-L/2:dx:L/2-dx; y=x; %Coord
  [X,Y]=meshgrid(x,y);

  w1=1e-3;                %x half-width
  w2=2e-4;                %y half-width
  lambda=0.633e-6;        %Wavelength
  z=5;                  %Propagation Distance
  k=2*pi/lambda;          %Wavenumber
  lz=lambda*z;

  %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)
  r=sqrt(X.^2+Y.^2)
  
  %u1=exp(i*k*z).*exp(i*k*(r.^2)/2/z)/(i*lz).*((w1^2.*jinc(2*pi*w1.*r/lz)./(w1.*r/lz))-(w2^2.*jinc(2*pi*w2.*r/lz)./(w2.*r/lz)));
  u1=circ(r./w1)-circ(r./w2);
  u2=propIR(u1,L,lambda,z);
  I=abs(u2).^2;
  figure(3);
  imagesc(x,y,nthroot(I,20));
  xlabel('x(m)');ylabel('y(m)');
  colormap('gray');
  axis square;
  axis xy;
  figure(16);
  plot(x,I(M/2+1,:))
  xlabel('x(m)');ylabel("Irradiance")
%Problem 2c)
  L=0.025;               %Side length(m)
  M=500;                 %Samples
  dx=L/M;                %Sample Interval
  x=-L/2:dx:L/2-dx; y=x; %Coord
  [X,Y]=meshgrid(x,y);

  w=1e-3;                %x half-width
  s=4e-3;                %y half-width
  lambda=0.633e-6;        %Wavelength
  z=5;                    %Propagation Distance
  k=2*pi/lambda;          %Wavenumber
  lz=lambda*z;

  %Irradiance
  %I(x,y)=|U(x,y)|^2=2exp(jkz)*exp(jkz(x^2+y^2)/lz)*exp(w^2*jinc(2pi*r/lz))*cos(pi*s*x/lz)/jlz
  u1=circ(r./w).*cos(pi.*s.*x/lz)
  u2=propIR(u1,L,lambda,z);
  r=sqrt(X.^2+Y.^2)
  I=abs(u2).^2
  figure(5);
  imagesc(x,y,nthroot(I,5));
  xlabel('x(m)');ylabel('y(m)');
  colormap('gray');
  axis square;
  axis xy;
  figure(17);
  plot(x,I(M/2+1,:))
  xlabel('x(m)');ylabel("Irradiance")

% Critical Sample - dX=lz/L; 
% Problem 2a) dx=0.633x10^-12; dx=1.266x10^-11; dx=0.633x10^-11
%             Calculated dx=wx/500=2x10-7
%             | Observation plane limited. Small "z" or lambda
% Problem 2b) dx=1.26x10^-4; dx=5.06x10^-4; dx=1.26x10^-3
%             Calculated dx=wx/500=2x10-6
%             | Source Bandwidth limited. Large "z" or lambda
% Problem 2c) dx=1.26x10^-4; dx=5.06x10^-4; dx=1.26x10^-3
%             Calculated dx=wx/500=2x10-6
%             | Observation plane limited, but closer to critical sampling.
% Problem 3a)
% sqr_beam propagation example
%
% sqr_beam propagation example
%
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
  u2=propTF(u1,L,lambda,z);
  Pu=sum(sum(u2*dx,2)*dx,1)
  I1=abs(u1.^2);
  PI=sum(sum(I1*dx,2)*dx,1)
  figure(18);
  imagesc(x1,y1,I1);
  axis square; axis xy;
  colormap('gray');xlabel('x (m)');ylabel('y (m)');
  title('z= 0 m');

% My power calculation was two-orders of magnitude lower than the book.
% The propIR function produced differing source and observation plane
% powers.
% Problem 4a)
  %{
  function[u2]=TFRS(u1,L,lambda,z);
  %propagation - Rayleigh-Sommerfeld diffraction transfer function
  %assumes same x and y side lengths and
  %uniform sampling
  %u1 - source plane field
  %L - source and observation plane side length
  %lambda - wavelength
  %z - propagation distance
  %u2 - observation plane field

  [M,N]=size(u1);              %get input field array size
  dx=L/M;                      %sample interval
  k=2*pi/lambda;               %wavenumber

  fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %Frequency Coordinates
  [FX,FY]=meshgrid(fx,fx);

  H=exp(j*k*z*sqrt(1-(lambda*FX).^2+(lambda*FY).^2)); %Transfer function
  H=fftshift(H);                       %shift trans func
  U1=fft2(fftshift(u1));               %shift, fft, src field
  U2=H.*U1;                            %multiply
  u2=ifftshift(ifft2(U2));             %inv fft, center obs field
  end
  
  function[u2]=IRRS(u1,L,lambda,z);
  %propagation - Rayleigh-Sommerfeld diffraction impulse response approach
  %assumes same x and y side lengths and
  %uniform sampling
  %u1 - source plane field
  %L - source and observation plane side length
  %lambda - wavelength
  %z - propagation distance
  %u2 - observation plane field

  [M,N]=size(u1);              %get input field array size
  dx=L/M;                      %sample interval
  k=2*pi/lambda;               %wavenumber

  x=-L/2:dx:L/2-dx; %Frequency Coordinates
  [X,Y]=meshgrid(x,x);

  h=z/(j*lambda)*exp(j*k*sqrt(X.^2+Y.^2))/(X.^2+Y.^2); %impulse function
  H=fft2(fftshift(h))*dx^2;                       %shift trans func
  U1=fft2(fftshift(u1));               %shift, fft, src field
  U2=H.*U1;                            %multiply
  u2=ifftshift(ifft2(U2));             %inv fft, center obs field
  end
  %}
% Problem 4b)
% The routines were saved, and processed below.
% sqr_beam propagation example

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
  u1=rect(X1./(2*w)).*rect(Y1./(2*w));
  u2=IRRS(u1,L1,lambda,z);
  I1=abs(u2.^2);
  figure(19)
  imagesc(x1,y1,I1);
  axis square; axis xy;
  colormap('gray');xlabel('x (m)');ylabel('y (m)');
  title('z= 2000 m');
  figure(20);
  plot(x1,I1(M/2+1,:))
  xlabel('x(m)');ylabel("Irradiance")
% There are differences between Rayleigh-Sommerfeld and Fresnel results,
% specifically, intensity.

%Problem 4c: Fresnel Numbers dX=lz/L;
% z=1000; dx=0.002;lambda*z/L=0.5*10^-6*2000/0.5=0.001 -- Transfer Func
% z=2000; dx=0.002;lambda*z/L=0.5*10^-6*2000/0.5=0.002 -- Transfer/Impulse
% z=4000; dx=0.002;lambda*z/L=0.5*10^-6*2000/0.5=0.004 -- Impulse Func
% z=20000; dx=0.002;lambda*z/L=0.5*10^-6*2000/0.5=0.02 -- Impulse Func

%Problem 5a:Gaussian Irradiance I(x,y,z)=I(w0/w(z))^2*exp(-2r^2/w(z)^2)
%            Beam Irradiance w(z)=w0*sqrt(1+(z/zR)^2)
%            Rayleigh Range  zR=pi*w0^2/lambda
%            Optical field   U0(x,y)=A*exp(-r^2/w0^2)
%            w0=1mm,lambda=0.633um,A=1V/m,I0=|A0|^2/(2*eta) W/m^2, eta=377 Ohms.
    L=15/1000; % Side length
    M=250;     % Sample Size
    dx=L/M;    % Sample Divisions
    x=-L/2:dx:L/2-dx;
    y=x;
    [X,Y]=meshgrid(x,y);
    
    lambda=0.633*10^-6 %Wavelength
    k=2*pi/lambda;     %Wavevector
    z=10 ;              %Propagation Distance
    
    A=1;
    w0=1/1000;
    U=A*exp(-(X.^2+Y.^2)/w0^2);
    I1=abs(U).^2;
    figure(21)
    imagesc(x,y,I1);
    axis square; axis xy;
    colormap('gray');xlabel('x (m)');ylabel('y (m)');
    
    eta=377;
    I0=abs(A).^2/(2*eta)
    zR=pi*w0^2/lambda;
    wZ=w0*sqrt(1+(z/zR)^2)
    I2=I0*(w0/wZ)^2*exp(-2*(X.^2+Y.^2)/wZ^2)
    figure(22)
    imagesc(x,y,I2);
    axis square; axis xy;
    colormap('gray');xlabel('x (m)');ylabel('y (m)');
%Problem 5b: Propagation distance for critical sampling is closest to 1m.
%Problem 5c: Fraunhofer irradiance expression I=abs(FFT(FFT(U0(x,y)))).^2
%Problem 6: H(fx,fy,z)=H(fx,fy,z1)H(fx,fy,z2-z1)...H(fx,fy,zn-zn-1)
% sqr_beam propagation example
%}
    L1=0.5;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    y1=x1;
    
    lambda=0.5*10^-6;
    k=2*pi/lambda;
    w=0.051;
    z=20000;
    
    [X1,Y1]=meshgrid(x1,y1);
    u1=rect(X1/(2*w)).*rect(Y1/(2*w));
    
    w=0.011;
    u2=propTF(u1,L1,lambda,2000);
    u3=propTF(u1,L1,lambda,4000);
    u4=propTF(u1,L1,lambda,6000);
    u5=propTF(u1,L1,lambda,8000);
    u6=propTF(u1,L1,lambda,10000);
    u7=propTF(u1,L1,lambda,12000);
    u8=propTF(u1,L1,lambda,14000);
    u9=propTF(u1,L1,lambda,16000);
    u10=propTF(u1,L1,lambda,18000);
    u11=propTF(u1,L1,lambda,20000);
    I2=abs(u2).^2
    I3=abs(u3).^2
    I4=abs(u4).^2
    I5=abs(u5).^2
    I6=abs(u6).^2
    I7=abs(u7).^2
    I8=abs(u8).^2
    I9=abs(u9).^2
    I10=abs(u10).^2
    I11=abs(u11).^2
    IT=I2*I3*I4*I5*I6*I7*I8*I9*I10*I11;
    figure(23)
    imagesc(x1,y1,nthroot(IT,3));
    xlabel('x (m)');ylabel('y (m)');
    T2=propTF(u1,L1,lambda,z);
    Irr=abs(T2).^2;
    figure(24)
    imagesc(x1,y1,nthroot(Irr,3))
    axis square; axis xy;
    colormap('gray');xlabel('x(m)');ylabel('y(m)');
    title(['z= ',num2str(z),' m']);
%The results should be the same from a split-step result and single Transfer Function Propagation.
% Problem 7a: Aperature : L1= 2mm; z= 5m; M=500;
  L=2/10;              %Side length(m)
  M=500;                 %Samples
  dx=L/M;                %Sample Interval
  x=-L/2:dx:L/2-dx; y=x; %Coord
  [X,Y]=meshgrid(x,y);

  wx=1*10^-4;                %x half-width
  wy=5*10^-5;                %y half-width
  lambda=0.633e-6;        %Wavelength
  z=5;                    %Propagation Distance
  k=2*pi/lambda;          %Wavenumber
  lz=lambda*z;

  %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)
  I=(4*wx*wy/lz)^2.*sinc(2*wx.*X/lz).^2.*sinc(2*wy.*Y/lz).^2;
  figure(25);
  plot(x,I(M/2+1,:))
  xlabel('x(m)');ylabel("Irradiance")
  xlim([-0.1 0.1])
% Problem 7b: Aperature : L= 2.5cm; z= 50m; M=500;
  L=0.25;               %Side length(m)
  M=500;                 %Samples
  dx=L/M;                %Sample Interval
  x=-L/2:dx:L/2-dx; y=x; %Coord
  [X,Y]=meshgrid(x,y);

  w1=1*10^-3;                %x half-width
  w2=2*10^-4;                %y half-width
  lambda=0.633*10^-6;        %Wavelength
  z=50;                    %Propagation Distance
  k=2*pi/lambda;          %Wavenumber
  lz=lambda*z;

  %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)
  r=sqrt(X.^2+Y.^2)
  I2=(1/lz)^2.*((w1^2.*jinc(2*pi*w1.*r/lz)./(w1.*r/lz))-(w2^2.*jinc(2*pi*w2.*r/lz)./(w2.*r/lz))).^2;
  figure(26);
  plot(x,nthroot(I2(M/2+1,:),2))
  %imagesc(x,y,nthroot(I2,2))
  xlabel('x(m)');ylabel("Irradiance")
  xlim([-0.05 0.05])
% Problem 7c: Aperature : L= 2.5cm; z= 50m; M=500; Slightly inconsistent.
% Also, reference problem 4.4 was solved incorrectly too.
  L=2.5/100;               %Side length(m)
  M=500;                 %Samples
  dx=L/M;                %Sample Interval
  x=-L/2:dx:L/2-dx; y=x; %Coord
  [X,Y]=meshgrid(x,y);

  w=1*10^3;                %x half-width
  s=4*10^-3;                %y half-width
  lambda=0.633*10^-6;        %Wavelength
  z=50;                    %Propagation Distance
  k=2*pi/lambda;          %Wavenumber
  lz=lambda*z;

  %Irradiance
  %I(x,y)=|U(x,y)|^2=2exp(jkz)*exp(jkz(x^2+y^2)/lz)*exp(w^2*jinc(2pi*r/lz))*cos(pi*s*x/lz)/jlz
  r=sqrt(X.^2+Y.^2)
  I3=(2/lz)^2.*((w^2.*jinc(2*pi*w.*r/lz)./(w.*r/lz))).^2.*cos(pi*s.*x/lz).^2;
  
  figure(27);
  plot(x,I3(M/2+1,:))
  xlabel('x(m)');ylabel("Irradiance")
  xlim([-0.01 0.01])
%Problem 8
%Eqn 2.11: B=5/w; dx=w/10; Square aperature, square or rectangle function.
%X-and-Y coordinates; 10x10 Samples=100 Samples
%Problem 9: Fresnel Two Step Propagator
    L1=0.5;
    L2=0.2;
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
    u2=prop2step(u1,L1,L2,lambda,z);
    I2=abs(u2).^2;
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2;
    y2=x2;
    
    figure(28);
    imagesc(x2,y2,nthroot(I2,3));
    xlabel('x (m)');ylabel('y (m)');
    figure(29);
    plot(x2,I2(M/2+1,:));
    xlabel('x (m)');ylabel('Irradiance');
    
    L1=0.5;
    L2=0.4;
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
    u2=prop2step(u1,L1,L2,lambda,z);
    I2=abs(u2).^2;
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2;
    y2=x2;
    
    figure(30);
    imagesc(x2,y2,nthroot(I2,3));
    xlabel('x (m)');ylabel('y (m)');
    figure(31);
    plot(x2,I2(M/2+1,:));
    xlabel('x (m)');ylabel('Irradiance');
    L1=0.5;
    L2=0.6;
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
    u2=prop2step(u1,L1,L2,lambda,z);
    I2=abs(u2).^2;
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2;
    y2=x2;
    
    figure(32);
    imagesc(x2,y2,nthroot(I2,3));
    xlabel('x (m)');ylabel('y (m)');
    figure(33);
    plot(x2,I2(M/2+1,:));
    xlabel('x (m)');ylabel('Irradiance');
%}
%The apparent features are shown in Figure 2 with over-sampled data.
%Problem 9b
    L1=0.5;
    L2=10;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    y1=x1;
    
    lambda=0.5*10^-6;
    k=2*pi/lambda;
    w=0.051;
    z=20000;
    
    [X1,Y1]=meshgrid(x1,y1);
    u1=rect(X1/(2*w)).*rect(Y1/(2*w));
    u2=prop2step(u1,L1,L2,lambda,z);
    I2=abs(u2).^2;
    dx2=L2/M;
    x2=-L2/2:dx2:L2/2-dx2;
    y2=x2;
    
    figure(34);
    imagesc(x2,y2,nthroot(I2,3));
    xlabel('x (m)');ylabel('Irradiance');
    figure(35);
    plot(x2,I2(M/2+1,:));
    xlabel('x (m)');ylabel('Irradiance');
 %  When L2 is large, there is oversampling, and as L2
 %  is reduced, low intensity generates choppy, pixelated data that is undersampled.