%Chapter: 4
%Problem Set - #4a
    L=0.2;               %Side length(m)
    M=250;                 %Samples
    dx=L/M;                %Sample Interval
    x=-L/2:dx:L/2-dx; y=x; %Coord
    [X,Y]=meshgrid(x,y);
    
    wx=1e-4;                %x half-width
    wy=5e-5;                %y half-width
    lambda=0.633e-6;        %Wavelength
    z=5;                    %Propagation Distance
    k=2*pi/lambda;          %Wavenumber
    lz=lambda*z;
    
    %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)
    I=(4*wx*wy/lz)^2.*sinc(2*wx.*X/lz).^2.*sinc(2*wy.*Y/lz).^2;
    figure(1);
    imagesc(x,y,nthroot(I,3));
    xlabel('x(m)');ylabel('y(m)');
    colormap('gray');
    axis square;
    axis xy;
    figure(2);
    plot(x,I(M/2+1,:))
    xlabel('x(m)');ylabel("Irradiance")

%Problem Set - #4b
    L=0.2;               %Side length(m)
    M=250;                 %Samples
    dx=L/M;                %Sample Interval
    x=-L/2:dx:L/2-dx; y=x; %Coord
    [X,Y]=meshgrid(x,y);
    
    w1=1e-3;                %x half-width
    w2=2e-4;                %y half-width
    lambda=0.633e-6;        %Wavelength
    z=50;                    %Propagation Distance
    k=2*pi/lambda;          %Wavenumber
    lz=lambda*z;
    
    %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)
    r=sqrt(X.^2+Y.^2)
    I2=(1/lz)^2.*((w1^2.*jinc(2*pi*w1.*r/lz)./(w1.*r/lz))-(w2^2.*jinc(2*pi*w2.*r/lz)./(w2.*r/lz))).^2;
    figure(3);
    imagesc(x,y,nthroot(I2,20));
    xlabel('x(m)');ylabel('y(m)');
    colormap('gray');
    axis square;
    axis xy;
    figure(4);
    plot(x,I2(M/2+1,:))
    xlabel('x(m)');ylabel("Irradiance")
    
    %Problem Set - #4c
    L=0.2;               %Side length(m)
    M=250;                 %Samples
    dx=L/M;                %Sample Interval
    x=-L/2:dx:L/2-dx; y=x; %Coord
    [X,Y]=meshgrid(x,y);
    
    w1=1e-3;                %x half-width
    w2=2e-4;                %y half-width
    lambda=0.633e-6;        %Wavelength
    z=50;                    %Propagation Distance
    k=2*pi/lambda;          %Wavenumber
    lz=lambda*z;
    
    %Irradiance I(x,y)=|U(x,y)|^2=exp(jkz)*exp(jk(x^2+y^2)/2/z)*(4wx*wy)*sinc(2*wx.*X/lambda/z)sinc(2*wy.*Y/lambda/z)
    r=sqrt(X.^2+Y.^2)
    I2=(1/lz)^2.*((w1^2.*jinc(2*pi*w1.*r/lz)./(w1.*r/lz))-(w2^2.*jinc(2*pi*w2.*r/lz)./(w2.*r/lz))).^2;
    figure(3);
    imagesc(x,y,nthroot(I2,20));
    xlabel('x(m)');ylabel('y(m)');
    colormap('gray');
    axis square;
    axis xy;
    figure(4);
    plot(x,I2(M/2+1,:))
    xlabel('x(m)');ylabel("Irradiance")

%Problem Set - #4c
    L=0.1;               %Side length(m)
    M=250;                 %Samples
    dx=L/M;                %Sample Interval
    x=-L/2:dx:L/2-dx; y=x; %Coord
    [X,Y]=meshgrid(x,y);
    
    w=1e-3;                %x half-width
    s=4e-3;                %y half-width
    lambda=0.633e-6;        %Wavelength
    z=50;                    %Propagation Distance
    k=2*pi/lambda;          %Wavenumber
    lz=lambda*z;
    
    %Irradiance
    %I(x,y)=|U(x,y)|^2=2exp(jkz)*exp(jkz(x^2+y^2)/lz)*exp(w^2*jinc(2pi*r/lz))*cos(pi*s*x/lz)/jlz
    r=sqrt(X.^2+Y.^2)
    I3=(2/lz)^2.*((w^2.*jinc(2*pi*w.*r/lz)./(w1.*r/lz))).^2.*cos(pi*s.*x/lz).^2;
    figure(5);
    imagesc(x,y,nthroot(I3,20));
    xlabel('x(m)');ylabel('y(m)');
    colormap('gray');
    axis square;
    axis xy;
    figure(6);
    plot(x,I2(M/2+1,:))
    xlabel('x(m)');ylabel("Irradiance")


%Problem Set - #5
    L=0.1;               %Side length(m)
    M=250;                 %Samples
    dx=L/M;                %Sample Interval
    x=-L/2:dx:L/2-dx; y=x; %Coord
    [X,Y]=meshgrid(x,y);
    
    w=1e-3;                %x half-width
    s=4e-3;                %y half-width
    lambda=0.633e-6;        %Wavelength
    z=50;                    %Propagation Distance
    k=2*pi/lambda;          %Wavenumber
    lz=lambda*z;
    A1=1;                  %Amplitude Irradiance
    A2=0.4                 %Amplitude Irradiance
    %Irradiance
    %I(x,y)=|U(x,y)|^2=2exp(jkz)*exp(jkz(x^2+y^2)/lz)*exp(w^2*jinc(2pi*r/lz))*cos(pi*s*x/lz)/jlz
    r=sqrt(X.^2+Y.^2)
    I4=((A1^2+A2^2)/lz^2).*((w^2.*jinc(2*pi*w.*r/lz)./(w.*r/lz))).^2.*(1+((2*A1*A2)/(A1^2+A2^2)).*cos(pi*s.*x/lz));
    figure(7);
    imagesc(x,y,nthroot(I4,20));
    xlabel('x(m)');ylabel('y(m)');
    colormap('gray');
    axis square;
    axis xy;
    figure(8);
    plot(x,I2(M/2+1,:))
    xlabel('x(m)');ylabel("Irradiance")
