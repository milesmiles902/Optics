%Chapter 9:
%Problem 1
%a) S(v)=1/dv*rect((v-v0)/dv); tau=1/dv

    lambda0=650e-9;
    c=3e8;
    k0=2*pi/lambda0;
    nu0=c/lambda0;
    
    N=51;
    delnu=2e9;
    b=delnu/(2*sqrt(log(2)));
    dnu=4*delnu/N;
    
    L1=50e-3;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    x1=fftshift(x1);
    [X1,Y1]=meshgrid(x1,x1);
    
    w=1e-3;
    dels=5e-3;
    deld=5e-2;
    f=0.25;
    lf=lambda0*f;
    
    I2=zeros(M);
    for n=1:N
        nu=(n-(N+1)/2)*dnu+nu0;
        S=1/delnu*rect((nu-nu0)/delnu);
        k=2*pi*nu/c;
        u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
        u2=1/lf*(fft2(u))*dx1^2;
        I2=I2+S*(abs(u2).^2)*dnu;
    end
    
    I2=ifftshift(I2);
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(1);
    imagesc(x2,y2,I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy; colormap('gray');
    
    figure(2);
    plot(x2,I2(M/2+1,:));
    xlabel('x(m)');ylabel('Irradiance')

%b)|gam(t)|=intS(v)*exp(-j2pi*nu*tau)dv
%          =sinc(dnu*tau)
%c) The simulation result shows similar results to analytic.
%Problem 2)
%a)
    lambda0=650e-9;
    c=3e8;
    k0=2*pi/lambda0;
    nu0=c/lambda0;
    
    N=51;
    delnu=2e9;
    b=delnu/(2*sqrt(log(2)));
    dnu=4*delnu/N;
    
    L1=50e-3;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    x1=fftshift(x1);
    [X1,Y1]=meshgrid(x1,x1);
    
    w=1e-3;
    dels=5e-3;
    deld=5e-2;
    f=0.25;
    lf=lambda0*f;
    
    I2=zeros(M);
    n=1
    nu=(n-(N+1)/2)*dnu+nu0;
    S=1/(sqrt(pi)*b)*exp(-(nu-nu0)^2/b^2);
    k=2*pi*nu/c;
    u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
    u2=1/lf*(fft2(u))*dx1^2;
    I2=I2+S*(abs(u2).^2)*dnu;
  
    I2=ifftshift(I2);
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(3);
    imagesc(x2,y2,I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy; colormap('gray');
    
    figure(4);
    plot(x2,I2(M/2,:));
    xlabel('x(m)');ylabel('Irradiance')
%b) The optical path difference limits intensity and symmetry of
    %irradiance profiles.
    deld=0;
    f=0.25;
    lf=lambda0*f;
    
    I2=zeros(M);
    n=1
    nu=(n-(N+1)/2)*dnu+nu0;
    S=1/(sqrt(pi)*b)*exp(-(nu-nu0)^2/b^2);
    k=2*pi*nu/c;
    u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
    u2=1/lf*(fft2(u))*dx1^2;
    I2=I2+S*(abs(u2).^2)*dnu;
  
    I2=ifftshift(I2);
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(3);
    imagesc(x2,y2,I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy; colormap('gray');
    
    figure(4);
    plot(x2,I2(M/2,:));
    xlabel('x(m)');ylabel('Irradiance')

    deld=lambda0/4;
    f=0.25;
    lf=lambda0*f;
    
    I2=zeros(M);
    n=1
    nu=(n-(N+1)/2)*dnu+nu0;
    S=1/(sqrt(pi)*b)*exp(-(nu-nu0)^2/b^2);
    k=2*pi*nu/c;
    u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
    u2=1/lf*(fft2(u))*dx1^2;
    I2=I2+S*(abs(u2).^2)*dnu;
  
    I2=ifftshift(I2);
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(5);
    imagesc(x2,y2,I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy; colormap('gray');
    
    figure(6);
    plot(x2,I2(M/2,:));
    xlabel('x(m)');ylabel('Irradiance')


    deld=lambda0/2;
    f=0.25;
    lf=lambda0*f;
    
    I2=zeros(M);
    n=1
    nu=(n-(N+1)/2)*dnu+nu0;
    S=1/(sqrt(pi)*b)*exp(-(nu-nu0)^2/b^2);
    k=2*pi*nu/c;
    u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
    u2=1/lf*(fft2(u))*dx1^2;
    I2=I2+S*(abs(u2).^2)*dnu;
  
    I2=ifftshift(I2);
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(7);
    imagesc(x2,y2,I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy; colormap('gray');
    
    figure(8);
    plot(x2,I2(M/2,:));
    xlabel('x(m)');ylabel('Irradiance')


    deld=lambda0*3/4;
    f=0.25;
    lf=lambda0*f;
    
    I2=zeros(M);
    n=1
    nu=(n-(N+1)/2)*dnu+nu0;
    S=1/(sqrt(pi)*b)*exp(-(nu-nu0)^2/b^2);
    k=2*pi*nu/c;
    u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
    u2=1/lf*(fft2(u))*dx1^2;
    I2=I2+S*(abs(u2).^2)*dnu;
  
    I2=ifftshift(I2);
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(9);
    imagesc(x2,y2,I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy; colormap('gray');
    
    figure(10);
    plot(x2,I2(M/2,:));
    xlabel('x(m)');ylabel('Irradiance')
%Problem 3:
%Visibility V=2A1A2gam(dD/c)/(A1^2+A2^2)
%a) V=0.75*0.20*gam(5*10^-2 m/(2.998*10^8 m/s))/(0.75^2+0.20^2);
%    =0.25*exp(-(pidnu*1.7*10^-10/(2sqrtln2))^2)*exp(-i2piv0(1.7*10^-10))
%    =0.25*exp(-(pi*2e9*1.7*10^-10/(2sqrtln2))^2)*exp(-i2pi*4.61*10^14*(1.7*10^-10))
%    =0.13*exp(156740pi*j)=0.13
%b) 
    lambda0=650e-9;
    c=3e8;
    k0=2*pi/lambda0;
    nu0=c/lambda0;
    
    N=51;
    delnu=2e9;
    b=delnu/(2*sqrt(log(2)));
    dnu=4*delnu/N;
    
    L1=50e-3;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1;
    x1=fftshift(x1);
    [X1,Y1]=meshgrid(x1,x1);
    
    w=1e-3;
    dels=5e-3;
    deld=5e-2;
    f=0.25;
    lf=lambda0*f;

    A1=0.75;
    A2=0.20;
    
    I2=zeros(M);
    for n=1:N
        nu=(n-(N+1)/2)*dnu+nu0;
        S=1/(sqrt(pi)*b)*exp(-(nu-nu0)^2/b^2);
        k=2*pi*nu/c;
        u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
        u2=1/lf*(fft2(u))*dx1^2;
        I2=I2+S*(abs(u2).^2)*dnu;
    end
    
    I2=ifftshift(I2);
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(11);
    imagesc(x2,y2,0.13.*I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy; colormap('gray');
    
    figure(12);
    plot(x2,0.13.*I2(M/2,:));
    xlabel('x(m)');ylabel('Irradiance')

%c: Analytical Expression:
%     u=V*circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
%     where V=2A1A2gam(dD/c)/(A1^2+A2^2
%   See figures of part b
%Problem 4: pc_spatial listed in the book is disfunctional.
%Problem 5: pc_spatial listed in the book is disfunctional.  
%Problem 6: 
%a) I=I0/del^2*exp(-2(x^2+y^2)/w0^2*del^2);
%   del=[1+(2x/kw0^2)^2(1+2*w0^2/lcr^2))^1/2
%   U1=exp(-x^2+y^2/w0^2)
    lambda=1.06e-9;
    L1=50e-2;
    M=250;
    dx1=L1/M;
    x1=-L1/2:dx1:L1/2-dx1; 
    x1=fftshift(x1);
    [X1,Y1]=meshgrid(x1,x1);
    w0=2e-3;
    N=100;
    Lcr=1e-2;
    z=1000;
    k=2*pi/lambda;
    
    I0=10;
    del=sqrt((1+(2*z/(k*w0^2)))^2*(1+2*w0^2/Lcr^2))
    I2=I0/del^2*exp(-2*(X1.^2+Y1.^2)/(w0^2*del^2));

    I2=ifftshift(I2)/N;
    x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
    y2=x2;
    
    figure(13);
    imagesc(x2,y2,I2);
    xlabel('x(m)');ylabel('y(m)');
    axis square; axis xy;
    colormap('gray');
    
    figure(14);
    plot(x2,I2(M/2+1,:));
    xlabel('x(m)');ylabel('Irradiance');