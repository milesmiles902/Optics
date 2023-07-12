%Chapter 7

%Problem 1:
    %a) A thin lens, r=10 mm, f=50 mm, z=200 mm;
    r=10*10^-3;
    f=50*10^-3;
    z=200*10^-3;
    
    %The image distance (1/z2=1/f-1/z1) is: 66 mm.
    %The transverse magnification (M=y2/y1=-z2/z1) is: -1/3
    
    %b) 
    % The exit pupil diameter (Dxp=lens diameter) is 2*radius=20 mm
    % The axial image point (zxp) is z2= 66 mm.
    % While the paraxial working number (f/#) or f-number is zxp/Dxp=3.33.
    
    %c)The approximate image distance such that z1>>f is 1/z2=1/f-1/z1.
    
    % Problem 2:
    % Suppose Dxp=0.5in, zxp=4 inch, and lambda=0.587 um.
    %a) The f-number(f/#) is zxp/Dxp is: 8.
    %   The coherent cutoff frequency (2f=2w/lz=1/lf) is: 3x10^13 1/s=212cyc/mm
    %   The incoherent cutoff frequency (2f=2w/lz=1/lf) is: 3x10^13 1/s
    %   cycle/mm
    
    %b) The sample interval requirement for an ideal image in a
    % diffraction-limited simulation is dU<=lf/2: 2.348 um
    
    %c) Suppose side length (L) is 1mm. The sampling requirements are L<=Mlf/2
    %   So then, the sampling requirement (M) would be: 426 samples

% Problem 3:
    % Assume zxp=50 mm, lamdba=0.5 um, half width wx = 1 mm, wy = 0.5 mm, L =
    % 1mm, M = 250;
    %a) The coherent cutoff frequency (f=w/lz) in x-direction: 40000 cyc/m
    %   The coherent cutoff frequency (f=w/lz) in y-direction: 20000 cyc/m
    %   The incoherent cutoff frequency (f=2w/lz) in x-direction: 80000 cyc/m
    %   The incoherent cutoff frequency (f=2w/lz) in y-direction: 40000 cyc/m
    
    %b) 
    A=imread('USAF1951B20','png');
    [M,N]=size(A);
    A=flipud(A);
    Ig=single(A);
    Ig=Ig./max(max(Ig));
    ug=sqrt(Ig(:,1:208,1)).*exp(j*2*pi*rand(M));
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du;v=u;
    
    figure(1);
    imagesc(u,v,Ig);
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    lambda=0.5*10^-6;
    wxp=1e-3;
    wyp=0.5e-3;
    zxp=50e-3;
    fx0=wxp/(lambda*zxp);
    fy0=wxp/(lambda*zxp);
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fv=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fv);
    H=meshgrid(rect(-lambda*zxp*fu/(2*wxp)),rect(-lambda*zxp*fv/(2*wyp)));
    
    figure(2);
    surf(fu,fv,H.*.99);
    camlight left; lighting phong;
    colormap('gray');
    shading interp;
    ylabel('fu (cyc/m)'); xlabel('fv (cyc/m)');
    
    H=fftshift(H);
    Gg=fft2(fftshift(ug));
    Gi=Gg(:,1:208,1).*H;
    ui=ifftshift(ifft2(Gi));
    Ii=(abs(ui)).^2;
    
    figure(3);
    imagesc(u,v,nthroot(Ii,2));
    colormap('gray');xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    figure(4);
    vvalue=-0.8e-4;
    vindex=round(vvalue/du+(M/2+1));
    plot(u,Ii(vindex,:),u,Ig(vindex,1:208,1),':');
    xlabel('u(m)');ylabel('Irradiance');
    
    %c)
    A=imread('USAF1951B20','png');
    [M,N]=size(A);
    A=flipud(A);
    Ig=single(A);
    Ig=Ig./max(max(Ig));
    
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du;v=u;
    
    lambda=0.5*10^-6;
    wxp=1e-3;
    wyp=0.5e-3;
    zxp=50e-3;
    fx0=wxp/(lambda*zxp);
    fy0=wxp/(lambda*zxp);
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fv=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fv);
    
    [Fu,Fv]=meshgrid(fu,fv);
    H=circ(sqrt(Fu.^2+Fv.^2)/sqrt(fx0.^2+fy0.^2));
    OTF=ifft2(abs(fft2(fftshift(H))).^2);
    OTF=OTF/OTF(1,1);
    figure(5);
    surf(fu,fv,fftshift(abs(OTF)));
    camlight left; lighting phong;
    colormap('gray');
    shading interp;
    ylabel('fu(cyc/m)');xlabel('fv (cyc/m)');zlabel('Intensity');
    
    Gg=fft2(fftshift(Ig));
    Gi=Gg(:,1:208,3).*OTF;
    Ii=ifftshift(ifft2(Gi));
    Ii=real(Ii);mask=Ii>=0;Ii=mask.*Ii;
    
    figure(6);
    imagesc(u,v,nthroot(Ii,2));
    colormap('gray');xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    figure(7);
    vvalue=0.2e-4;
    vindex=round(vvalue/du+(M/2+1));
    plot(u,Ii(vindex,:),u,Ig(vindex,1:208,3),':');
    xlabel('u(m)');ylabel('Irradiance');

    %d) The spatial resolution in the x vs y directions is different for
    % coherent imaging due to slit width changes.

%Problem 4
    % The axial image point (zxp) is: 50 mm.
    % The wavelength (lambda) is 0.5 uM.
    % The annular outer radius (w0) is: 1 mm.
    % The side length (L) is: 1 mm
    % The Number of Samples (M): 250
    
    %a) The coherent cutoff frequency (f=w/lz): 40000 cyc/m
    %   The incoherent cutoff frequency (f=2w/lz): 80000 cyc/m
    %b)  
    A=imread('USAF1951B20','png');
    [M,N]=size(A);
    A=flipud(A);
    Ig=single(A);
    Ig=Ig./max(max(Ig));
    ug=sqrt(Ig(:,1:208,1));
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du;v=u;
    
    figure(8);
    imagesc(u,v,Ig);
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    lambda=0.5*10^-6;
    wxpo=1*10^-3;
    wxpi=5*10^-4;
    zxp=50*10^-3;
    f0i=wxpi/(lambda*zxp);
    f0o=wxpo/(lambda*zxp);
    fuo=-1/(2*du):1/L:1/(2*du)-(1/L);
    fui=-1/(2*du):1/L:1/(2*du)-(1/L);
    fvi=fui;
    fvo=fuo;
    [Fui,Fvi]=meshgrid(fui,fvi);
    [Fuo,Fvo]=meshgrid(fuo,fvo);
    H=circ(sqrt(Fuo.^2+Fvo.^2)/f0o)-circ(sqrt(Fui.^2+Fvi.^2)/f0i);
    
    figure(9);
    surf(fuo,fvo,H.*.99);
    camlight left; lighting phong;
    colormap('gray');
    shading interp;
    ylabel('fu (cyc/m)'); xlabel('fv (cyc/m)');zlabel('Intensity');
    
    H=fftshift(H);
    Gg=fft2(fftshift(ug));
    Gi=Gg.*H(1:208,1:208);
    ui=ifftshift(ifft2(Gi));
    Ii=(abs(ui)).^2;
    
    figure(10);
    imagesc(u,v,nthroot(Ii,2));
    colormap('gray');xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    figure(11);
    vvalue=-0.8e-4;
    vindex=round(vvalue/du+(M/2+1));
    plot(u,Ii(vindex,:),u,Ig(vindex,1:208,1),':');
    xlabel('u(m)');ylabel('Irradiance');

  %b)
  A=imread('USAF1951B20','png');
    [M,N]=size(A);
    A=flipud(A);
    Ig=single(A);
    Ig=Ig./max(max(Ig));
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du;v=u;
    
    figure(12);
    imagesc(u,v,Ig);
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    lambda=0.5*10^-6;
    wxpo=1*10^-3;
    wxpi=5*10^-4;
    zxp=50*10^-3;
    f0i=wxpi/(lambda*zxp);
    f0o=wxpo/(lambda*zxp);
    fuo=-1/(2*du):1/L:1/(2*du)-(1/L);
    fui=-1/(2*du):1/L:1/(2*du)-(1/L);
    fvi=fui;
    fvo=fuo;
    [Fui,Fvi]=meshgrid(fui,fvi);
    [Fuo,Fvo]=meshgrid(fuo,fvo);
    H=circ(sqrt(Fuo.^2+Fvo.^2)/f0o)-circ(sqrt(Fui.^2+Fvi.^2)/f0i);
  
    OTF=ifft2(abs(fft2(fftshift(H))).^2);
    OTF=OTF/OTF(1,1);
    figure(13);
    surf(fu,fv,fftshift(abs(OTF)));
    camlight left; lighting phong;
    colormap('gray');
    shading interp;
    ylabel('fu(cyc/m)');xlabel('fv (cyc/m)');zlabel('Intensity');
    
    Gg=fft2(fftshift(Ig));
    Gi=Gg(:,1:208,3).*OTF;
    Ii=ifftshift(ifft2(Gi));
    Ii=real(Ii);mask=Ii>=0;Ii=mask.*Ii;
    
    figure(14);
    imagesc(u,v,nthroot(Ii,2));
    colormap('gray');xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    figure(15);
    vvalue=0.2e-4;
    vindex=round(vvalue/du+(M/2+1));
    plot(u,Ii(vindex,:),u,Ig(vindex,1:208,3),':');
    xlabel('u(m)');ylabel('Irradiance');

%d) Losing lower frequency DC components effects the resolution of
%   fine-grained imaging

%Problem 5
%a) 
    A=imread('USAF1951B20','png');
    [M,N]=size(A);
    A=flipud(A);
    Ig=single(A);
    Ig=Ig./max(max(Ig));
    ug=sqrt(Ig(:,1:208,1)).*exp(j*rand(M));
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du;v=u;
    
    figure(16);
    imagesc(u,v,Ig);
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    lambda=0.5*10^-6;
    wxp=3*10^-3;
    zxp=50*10^-3;
    f0=wxp/(lambda*zxp);
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fv=fu;
    [Fu,Fv]=meshgrid(fu,fv);
    H=circ(sqrt(Fu.^2+Fv.^2)/f0);
    
    figure(17);
    surf(fu,fv,H.*.99);
    camlight left; lighting phong;
    colormap('gray');
    shading interp;
    ylabel('fu (cyc/m)'); xlabel('fv (cyc/m)');zlabel('Intensity');
    
    H=fftshift(H);
    Gg=fft2(fftshift(ug));
    Gi=Gg.*H(1:208,1:208);
    ui=ifftshift(ifft2(Gi));
    Ii=(abs(ui)).^2;
    
    figure(18);
    imagesc(u,v,nthroot(Ii,2));
    colormap('gray');xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    figure(19);
    vvalue=-0.8e-4;
    vindex=round(vvalue/du+(M/2+1));
    plot(u,Ii(vindex,:),u,Ig(vindex,1:208,1),':');
    xlabel('u(m)');ylabel('Irradiance');

    % As the 2pi multiplier changes, the irradiated image is subject to
    % increased or decreased resolution. The 2pi multiplier represents phase angle.

    %b) Aperature size (wxp) alters the intensity of the irradiated image i.e.
    %   more or less light.
%Problem 6:
    % Consider the coherent image transfer function
    % H(fu,fv)=circ(sqrt(fu^2+fv^2)/f0)
    %a) The impulse response : h(u,v)=FT(H(fu,fv))=(J1(2pif0*sqrt(u^2+v^2))/sqrt(u^2+v^2))^2
    %b) The incoherent point spread function:
    % h(u,v)=FT(H(fu,fv))=FT(J1(2piw/lz*sqrt(fu^2+fv^2))/sqrt(fu^2+fv^2))^2);
    % Therefore, h(u,v)=circ(sqrt(u^2+v^2)/(w/lz))
    %c) The full width of the central lobe of an Airy pattern PSF
    % represents:2D=2*1.22*lz/w=2.44*lf#
%Problem 7:
    % The Raleigh criterion distance: d=1.22*l*f#
    %a) The calculated value of Raleigh criterion: d=1.22*0.5*10-6*20=12.2um;
    %b) The number of samples required for the intervalue (du) to span
    %distance (d) would be: du=lf#/2= 5 samples.
    %c) Ideal image frame consisting of two point sources seperated by
    %distance (d): 
    A=imread('USAF1951B20','png');
    [M,N]=size(A);
    A=flipud(A);
    S=10;
    Ig=zeros(M)
    Ig(M/2+1,M/2+1-round(S/2))=1;
    Ig(M/2+1,M/2+1+round(S/2))=1;
    ug=sqrt(Ig(:,1:208,1));
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du;v=u;
    
    figure(20);
    imagesc(u,v,Ig);
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;

    lambda=0.5*10^-6;
    wxp=1*10^-3;
    zxp=50*10^-3;
    f0=wxp/(lambda*zxp);
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fv=fu;
    [Fu,Fv]=meshgrid(fu,fv);
    H=circ(sqrt(Fu.^2+Fv.^2)/f0);
    
    figure(21);
    surf(fu,fv,H.*.99);
    camlight left; lighting phong;
    colormap('gray');
    shading interp;
    ylabel('fu (cyc/m)'); xlabel('fv (cyc/m)');zlabel('Intensity')
    
    H=fftshift(H);
    Gg=fft2(fftshift(ug));
    Gi=Gg.*H(1:208,1:208);
    ui=ifftshift(ifft2(Gi));
    Ii=(abs(ui)).^2;
    
    figure(22);
    imagesc(u,v,nthroot(Ii,2));
    colormap('gray');xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    figure(23);
    vvalue=-0.8e-4;
    vindex=round(vvalue/du+(M/2+1));
    plot(u,Ii(vindex,:),u,Ig(vindex,1:208,1),':');
    xlabel('u(m)');ylabel('Irradiance');
%Problem 8:
%a) Altitude: 700 km, Diameter = 1 m, lambda = 0.5 um
%   F#obj=z/D=700,000 m/1m=700000
%b) Incoherent Cutoff Frequency: f=w/lz= 1.43 1/m
%   Rayleigh Resolution Distnace: S=1.22l*f# = 0.427 m
%c) Maximum Sample Interval: dU<=lf#/2=0.175m
%   Object plane side length (1024x1024): L<=Mlf#/2=3.58 m
%Problem 9:
%a) Convert the USAF image to a phase object
    A=imread('USAF1951B20','png');
    [M,N]=size(A);
    A=flipud(A);
    Ig=single(A);
    Ig=Ig./max(Ig);
    theta=sqrt(Ig)*pi/100;
    ug=exp(j*theta);
    L=0.3e-3;
    du=L/M;
    u=-L/2:du:L/2-du;v=u;
    
    figure(24);
    imagesc(u,v,Ig(:,1:208,1));
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    lambda=0.5*10^-6;
    wxp=6.25e-3;
    zxp=125e-3;
    f0=wxp/(lambda*zxp);
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fv=fu;
    [Fu,Fv]=meshgrid(fu,fv);
    H=circ(sqrt(Fu.^2+Fv.^2)/f0);
    
    figure(25);
    surf(fu,fv,H.*.99);
    camlight left; lighting phong;
    colormap('gray');
    shading interp;
    ylabel('fu (cyc/m)'); xlabel('fv (cyc/m)');zlabel('Intensity')
    
    H=fftshift(H);
    Gg=fft2(fftshift(ug));
    Gi=Gg(:,1:208,1).*H;
    ui=ifftshift(ifft2(Gi));
    Ii=(abs(ui)).^2;
    
    figure(26);
    imagesc(u,v,nthroot(Ii,2));
    colormap('gray');xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
    
    figure(27);
    vvalue=-0.8e-4;
    vindex=round(vvalue/du+(M/2+1));
    plot(u,Ii(vindex,:),u,Ig(vindex,1:208,1),':');
    xlabel('u(m)');ylabel('Irradiance');
%b) Examine the square magnitude:
    figure(28)
    imagesc(u,v,abs(ug(:,1:208,1)).^2);
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
%c) See problem b
%d)  
    figure(29)
    imagesc(u,v,abs(ug(:,1:208,1)*exp(j*3*pi/2)).^2)
    colormap('gray'); xlabel('u(m)');ylabel('v(m)');
    axis square;
    axis xy;
%e) The analytic results show slight change of the overall intensity.
%Problem 10:
%a) 
    w1x=0.2;
    w1y=1;
    w2x=0.5;
    w2y=0.2;
    w3x=0.5;
    w3y=0.2;
    
    Lx=10;
    Ly=10;
    M=500;
    dx=Lx/M;
    dy=Ly/M;
    x=-Lx/2:dx:Lx/2-dx;
    y=-Ly/2:dy:Ly/2-dy;
    [X,Y]=meshgrid(x,y);
    g1=rect((X+2)/(2*w1x)).*rect(Y/(2*w1y))+rect((X+1.3)/(2*w2x)).*rect((Y-0.8)/(2*w2y))+rect((X+1.3)/(2*w3x)).*rect(Y/(2*w3y));
    figure(30);
    imagesc(x,y,g1);
    colormap('gray');
    axis square;
    axis xy;
    xlabel('x(m)');ylabel('y(m)');
%b)
    z=2000;

    lambda=0.5*10^-6;
    k=2*pi/lambda;
    w=0.5;
    z=2000;
    zf=1000;

    u1=rect(X/(2*w)).*rect(Y/(2*w));
    u2=propFF(u1,Lx,lambda,z);
    [u3]=focus(u2,Lx,lambda,z);
    I1=abs(u3.^2);
    figure(31)
    g1=g1+I1./1000;
    imshow(g1);
    colormap('gray');
    axis square;
    axis xy;
    xlabel('x(m)');ylabel('y(m)');
%c,d,e,f are incomplete.
