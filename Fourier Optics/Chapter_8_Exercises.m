%Chapter 8
%Problem 2
%a) (u,v)=(0,0), w_040=1,w_{d,131,222,220,311}=0;
%{
    L=1e-2;
    M=250;
    dx=L/M;
    x=-L/2:dx:L/2-dx; y=x;
    [X,Y]=meshgrid(x,y);
    u0=0;v0=0;
    wd=0;w040=1;w131=0;w222=0;w220=0;w311=0;
    w=seidel_5(u0,v0,X,Y,wd,w040,w131,w222,w220,w311);
    P=circ(sqrt(X.^2+Y.^2));
    mask=(P==0);
    w(mask)=NaN;
    
    figure(1);
    surfc(x,y,w);
    camlight left; lighting phong;
    colormap('gray');shading interp;
    xlabel('x');ylabel('y');
%b) 
    L=1e-2;
    M=250;
    dx=L/M;
    x=-L/2:dx:L/2-dx; y=x;
    [X,Y]=meshgrid(x,y);
    u0=0;v0=0;
    wd=-(X.^2+Y.^2);w040=1;w131=0;w222=0;w220=0;w311=0;
    w=seidel_5(u0,v0,X,Y,wd,w040,w131,w222,w220,w311);
    P=circ(sqrt(X.^2+Y.^2));
    mask=(P==0);
    w(mask)=NaN;
    
    figure(1);
    surfc(x,y,w);
    camlight left; lighting phong;
    colormap('gray');shading interp;
    xlabel('x');ylabel('y');
% Problem 3: dD=8f^2*w_{d}; dD=+0.5 creates the narrowest peak, with a
% w_{d} value of 2.5x10^-6 m.
    dD1=-1e-3;
    dD2=-0.5e-3;
    dD3=0e-3;
    dD4=0.5e-3;
    dD5=1e-3;

    M=1024;
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du; v=u;
    
    lambda=0.55e-6;
    k=2*pi/lambda;
    Dxp=20e-3;wxp=Dxp/2;
    zxp=100e-3;
    fnum=zxp/(2*wxp);
    lz=lambda*zxp;
    twof0=1/(lambda*fnum);
    
    u0=0;v0=0;
    
    wd=dD1/(8*fnum^2);
    w040=4.963*lambda;
    w131=2.637*lambda;
    w222=9.025*lambda;
    w220=7.536*lambda;
    w311=0.157*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fu);
    
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    
    figure(1);
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF');xlim([-1e-4 1e-4]);
    hold on;
    
    wd=dD2/(8*fnum^2);
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF'); xlim([-1e-4 1e-4]);
    
    wd=dD3/(8*fnum^2);
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF');xlim([-1e-4 1e-4]);
   
    wd=dD4/(8*fnum^2);
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF');xlim([-1e-4 1e-4]);
    
    wd=dD5/(8*fnum^2);
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF');xlim([-1e-4 1e-4]);
    legend('dD1','dD2','dD3','dD4','dD5');
% Problem 4: fnum=2.5; f=200 mm; D=40 mm; h=3.5mm; lambda=0.633 uM;
    M=1024;
    L=0.1e-3;
    du=L/M;
    u=-L/2:du:L/2-du; v=u;
    
    lambda=0.633e-6;
    k=2*pi/lambda;
    Dxp=40e-3;wxp=Dxp/2;
    zxp=100e-3;
    fnum=2.5;
    lz=lambda*zxp;
    twof0=1/(lambda*fnum);
    
    u0a=0;v0a=0;
    u0b=0;v0b=1;
    u0c=0.707;v0c=0.707;
    
    wd=0;
    w040=0;
    w131=-1.3792*lambda;
    w222=0.4815*lambda;
    w220=0;
    w311=0;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fu);
    
    Wa=seidel_5(u0a,v0a,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    Ha=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*Wa);
    h2a=abs(ifftshift(ifft2(fftshift(Ha)))).^2;

    Wb=seidel_5(u0b,v0b,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    Hb=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*Wb);
    h2b=abs(ifftshift(ifft2(fftshift(Hb)))).^2;
 
    Wc=seidel_5(u0c,v0c,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    Hc=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*Wc);
    h2c=abs(ifftshift(ifft2(fftshift(Hc)))).^2;

    figure(1);
    plot(u,h2a(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    hold all;
    plot(u,h2b(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    plot(u,h2c(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    legend('u0a','u0b','u0c')
    figure(2);
    plot(u,h2a(:,M/2+1));xlabel('v(m)');ylabel('PSF');
    hold all;
    plot(u,h2b(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    plot(u,h2c(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    legend('v0a','v0b','v0c');
%Problem 5: Critically satisfied means dW/dx < L/4f 
% = 2W_{d}x+4W_{040}((x^2+y^2)x)+W_{131}u(3x^2+y^2)+2W_{222}u^2*x+2W_220*u^2*x+W_{311}*u^3;
% W_d=19.4
    M=1024;
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du; v=u;
    
    lambda=0.55e-6;
    k=2*pi/lambda;
    Dxp=20e-3;wxp=Dxp/2;
    zxp=100e-3;
    fnum=zxp/(2*wxp);
    lz=lambda*zxp;
    twof0=1/(lambda*fnum);
    
    u0=1;v0=0;
    
    wd=19.4*lambda;
    w040=4.963*lambda;
    w131=2.637*lambda;
    w222=9.025*lambda;
    w220=7.536*lambda;
    w311=0.157*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fu);
    
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    figure(1);
    imagesc(u,v,angle(H));
    xlabel('u(m)');ylabel('v(m)');colormap('gray');
    
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    
    figure(2);
    imagesc(u,v,nthroot(h2,2));
    xlabel('u(m)');ylabel('v(m)');colormap('gray');
    figure(3);
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    figure(4);
    plot(u,h2(:,M/2+1));xlabel('v(m)');ylabel('PSF');
    
    MTF=fft2(fftshift(h2));
    MTF=abs(MTF/MTF(1,1));
    MTF=ifftshift(MTF);
    
    MTF_an=(2/pi)*(acos(fu/twof0)-(fu/twof0).*sqrt(1-(fu/twof0).^2));
    MTF_an=MTF_an.*rect(fu/(2*twof0));
    
    figure(5);
    plot(fu,MTF(M/2+1,:),fu,MTF(:,M/2+1),':',fu,MTF_an,'--');
    axis([0 150000 0 1]);
    legend('u MTF','v MTF','diff limit');
    xlabel('f(cyc/m)');ylabel('Modulation');

%Problem 6:
%a) 
    function[w]=zernike_37(X,Y,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,n32,n33,n34,n35,n36,n37)
      p=sqrt(X.^2+Y.^2);
      A=atan2(Y,X);
      c1=7.70185740
      c2=0.00000000*4^(1/2)*(p)*COS(A);
      c3=1.04836665*4^(1/2)*(p)*SIN(A);
      c4=4.95277338*3^(1/2)*(2*p^2-1);
      c5=0.00000000*6^(1/2)*(p^2)*SIN(2*A);
      c6=-1.83668298*6^(1/2)*(p^2)*COS(2*A);
      c7=0.38477648*8^(1/2)*(3*p^3-2*p)*SIN(A);
      c8=0.00000000*8^(1/2)*(3*p^3-2*p)*COS(A);
      c9=-0.00829929*8^(1/2)*(p^3)*SIN(3*A);
      c10=0.00000000*8^(1/2)*(p^3)*COS(3*A);
      c11=0.39476169*5^(1/2)*(6*p^4-3*p^2+1);
      c12=-0.00582307*10^(1/2)*(4*p^4-3*p^2)*COS(2*A);
      c13=0.00000000*10^(1/2)*(4*p^4-3*p^2)*SIN(2*A);
      c14=0.00007007*10^(1/2)*(p^4)*COS(4*A);
      c15=0.00000000*10^(1/2)*(p^4)*SIN(4*A);
      c16=0.00000000*12^(1/2)*(10*p^5-12*p^3+3*p)*COS(A);
      c17=0.00782021*12^(1/2)*(10*p^5-12*p^3+3*p)*SIN(A);
      c18=0.00000000*12^(1/2)*(5*p^5-4*p^3)*COS(3*A);
      c19=-0.00010086*12^(1/2)*(5*p^5-4*p^3)*SIN(3*A);
      c20=0.00000000*12^(1/2)*(p/5)*COS(5*A);
      c21=0.00000048*12^(1/2)*(p^5)*SIN(5*A);
      c22=0.00233053*7^(1/2)*(12*p^6-30*p^4+12*p^2-1)
      c23=0.00000000*14^(1/2)*(15*p^6-20*p^4+6*p^2)*SIN(2*A);
      c24=-0.00010805*14^(1/2)*(15*p^6-20*p^4+6*p^2)*COS(2*A)
      c25=0.00000000*14^(1/2)*(6*p^6-5*p^4)*SIN(4*A);
      c26=0.00000126*14^(1/2)*SIN(6*A);
      c27=0.00000000*14^(1/2)*(p^6)*SIN(6*A);
      c28=0.00000006*14^(1/2)*(p^6)*COS(6*A);
      c29=0.00008738*16^(1/2)*(35*p^7-60*p^5+30*p^3-4*p)*SIN(A);
      c30=0.00000000*16^(1/2)*(35*p^7-60*p^5+30*p^3-4*p)*COS(A);
      c31=-0.00000163*16^(1/2)*(21*p^7-30*p^5+10*p^3)*SIN(3*A);
      c32=0.00000000*16^(1/2)*(21*p^7-30*p^5+10*p^3)*COS(3*A);
      c33=0.00000000*16^(1/2)*(7*p^7-6*p^5)*SIN(5*A);
      c34=0.00000000*16^(1/2)*(7*p^7-6*p^5)*COS(5*A);
      c35=0.00000002*16^(1/2)*(p^7)*SIN(7*A);
      c36=0.00000000*16^(1/2)*(p^7)*COS(7*A);
      c37=0.00001672*9^(1/2)*(70*p^8-140*p^6+90*p^5-20*p^2+1)
      w=c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12+c13+c14+c15+c16+c17+c18+c19+c20+c21+c22+c23+c24+c25+c26+c27+c28+c29+c30+c31+c32+c33+C34+c35+c36+c37;
    end
%b) Incomplete
    M=1024;
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du; v=u;
    
    lambda=0.55e-6;
    k=2*pi/lambda;
    Dxp=20e-3;wxp=Dxp/2;
    zxp=100e-3;
    fnum=zxp/(2*wxp);
    lz=lambda*zxp;
    twof0=1/(lambda*fnum);
    
    u0=0;v0=1;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fu);
    
    W=zernike_37(u0,v0*lambda);
    
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    figure(1);
    imagesc(u,v,angle(H));
    xlabel('u(m)');ylabel('v(m)');colormap('gray');
    
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    
    figure(2);
    imagesc(u,v,nthroot(h2,2));
    xlabel('u(m)');ylabel('v(m)');colormap('gray');
    figure(3);
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    figure(4);
    plot(u,h2(:,M/2+1));xlabel('v(m)');ylabel('PSF');
    
    MTF=fft2(fftshift(h2));
    MTF=abs(MTF/MTF(1,1));
    MTF=ifftshift(MTF);
    
    MTF_an=(2/pi)*(acos(fu/twof0)-(fu/twof0).*sqrt(1-(fu/twof0).^2));
    MTF_an=MTF_an.*rect(fu/(2*twof0));
    
    figure(5);
    plot(fu,MTF(M/2+1,:),fu,MTF(:,M/2+1),':',fu,MTF_an,'--');
    axis([0 150000 0 1]);
    legend('u MTF','v MTF','diff limit');
    xlabel('f(cyc/m)');ylabel('Modulation');

%Problem 7
%a)
    M=250;
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du; v=u;
    
    lambda=0.5*10^-6;
    k=2*pi/lambda;
    wxp=2.5e-3;
    zxp=100e-3;
    fnum=zxp/(2*wxp);
    
    twof0=1/(lambda*fnum);
    fN=1/(2*du);
    
    wd=1*lambda;
    w040=0*lambda;
    w131=0*lambda;
    w222=0*lambda;
    w220=0*lambda;
    w311=0*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fu=fftshift(fu);
    [Fu,Fv]=meshgrid(fu,fu);
    
    I=zeros(M);
    
    %for u0=[-0.7:-0.7/3:0.7]
    %    for v0=[-0.7:0.7/3:0.7]
             W=seidel_5(u,v,-2*lambda*fnum*Fu,-2*lambda*fnum*Fv,wd,w040,w131,w222,w220,w311);
             H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum).*exp(-j*k*W);
             h2=abs(ifftshift(ifft2(H))).^2;
             h2=circshift(h2,[round(v*M/2),round(u*M/2)]);
             I=I+h2;
    %    end
    %end

    figure(1)
    imagesc(u,v,nthroot(I,1));
    xlabel('u(m)');ylabel('v(m)');
    colormap('gray'); axis square; axis xy;

   
    wd=0*lambda;
    w040=1*lambda;
    w131=0*lambda;
    w222=0*lambda;
    w220=0*lambda;
    w311=0*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fu=fftshift(fu);
    [Fu,Fv]=meshgrid(fu,fu);
    
    I=zeros(M);
    
    %for u0=[-0.7:-0.7/3:0.7]
    %    for v0=[-0.7:0.7/3:0.7]
             W=seidel_5(u,v,-2*lambda*fnum*Fu,-2*lambda*fnum*Fv,wd,w040,w131,w222,w220,w311);
             H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum).*exp(-j*k*W);
             h2=abs(ifftshift(ifft2(H))).^2;
             h2=circshift(h2,[round(v*M/2),round(u*M/2)]);
             I=I+h2;
    %    end
    %end
%b)
    figure(2)
    imagesc(u,v,nthroot(I,1));
    xlabel('u(m)');ylabel('v(m)');
    colormap('gray'); axis square; axis xy;

    wd=0*lambda;
    w040=0*lambda;
    w131=2*lambda;
    w222=0*lambda;
    w220=0*lambda;
    w311=0*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fu=fftshift(fu);
    [Fu,Fv]=meshgrid(fu,fu);
    
    I=zeros(M);
    
    %for u0=[-0.7:-0.7/3:0.7]
    %    for v0=[-0.7:0.7/3:0.7]
             W=seidel_5(u,v,-2*lambda*fnum*Fu,-2*lambda*fnum*Fv,wd,w040,w131,w222,w220,w311);
             H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum).*exp(-j*k*W);
             h2=abs(ifftshift(ifft2(H))).^2;
             h2=circshift(h2,[round(v*M/2),round(u*M/2)]);
             I=I+h2;
    %    end
    %end
%c)
    figure(3)
    imagesc(u,v,nthroot(I,1));
    xlabel('u(m)');ylabel('v(m)');
    colormap('gray'); axis square; axis xy;

     wd=0*lambda;
    w040=0*lambda;
    w131=0*lambda;
    w222=0*lambda;
    w220=1*lambda;
    w311=0*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fu=fftshift(fu);
    [Fu,Fv]=meshgrid(fu,fu);
    
    I=zeros(M);
    
    %for u0=[-0.7:-0.7/3:0.7]
    %    for v0=[-0.7:0.7/3:0.7]
             W=seidel_5(u,v,-2*lambda*fnum*Fu,-2*lambda*fnum*Fv,wd,w040,w131,w222,w220,w311);
             H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum).*exp(-j*k*W);
             h2=abs(ifftshift(ifft2(H))).^2;
             h2=circshift(h2,[round(v*M/2),round(u*M/2)]);
             I=I+h2;
    %    end
    %end
%d)
    figure(4)
    imagesc(u,v,nthroot(I,1));
    xlabel('u(m)');ylabel('v(m)');
    colormap('gray'); axis square; axis xy;

    wd=0*lambda;
    w040=0*lambda;
    w131=0*lambda;
    w222=0*lambda;
    w220=1*lambda;
    w311=0*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fu=fftshift(fu);
    [Fu,Fv]=meshgrid(fu,fu);
    
    I=zeros(M);
    
    %for u0=[-0.7:-0.7/3:0.7]
    %    for v0=[-0.7:0.7/3:0.7]
             W=seidel_5(u,v,-2*lambda*fnum*Fu,-2*lambda*fnum*Fv,wd,w040,w131,w222,w220,w311);
             H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum).*exp(-j*k*W);
             h2=abs(ifftshift(ifft2(H))).^2;
             h2=circshift(h2,[round(v*M/2),round(u*M/2)]);
             I=I+h2;
    %    end
    %end
%e)
    figure(5)
    imagesc(u,v,nthroot(I,1));
    xlabel('u(m)');ylabel('v(m)');
    colormap('gray'); axis square; axis xy;

     wd=0*lambda;
    w040=0*lambda;
    w131=0*lambda;
    w222=0*lambda;
    w220=0*lambda;
    w311=3*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    fu=fftshift(fu);
    [Fu,Fv]=meshgrid(fu,fu);
    
    I=zeros(M);
    
    for u0=[-0.7:-0.7/3:0.7]
        for v0=[-0.7:0.7/3:0.7]
             W=seidel_5(u0,v0,-2*lambda*fnum*Fu,-2*lambda*fnum*Fv,wd,w040,w131,w222,w220,w311);
             H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum).*exp(-j*k*W);
             h2=abs(ifftshift(ifft2(H))).^2;
             h2=circshift(h2,[round(v*M/2),round(u*M/2)]);
             I=I+h2;
        end
    end
%f)
    figure(6)
    imagesc(u,v,nthroot(I,1));
    xlabel('u(m)');ylabel('v(m)');
    colormap('gray'); axis square; axis xy;

%Problem 8
%a) 
    M=1024;
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du; v=u;
    
    lambda=0.55e-6;
    k=2*pi/lambda;
    Dxp=20e-3;wxp=Dxp/2;
    zxp=100e-3;
    fnum=zxp/(2*wxp);
    lz=lambda*zxp;
    twof0=1/(lambda*fnum);
    
    u0=0;v0=0;
  
    wd=0*lambda;
    w040=0*lambda;
    w131=0*lambda;
    w222=0*lambda;
    w220=0*lambda;
    w311=0*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fu);
    
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    figure(1);
    imagesc(u,v,angle(H));
    xlabel('u(m)');ylabel('v(m)');colormap('gray');
    
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    Iu=max(max(h2));

    figure(2);
    imagesc(u,v,nthroot(h2,2));
    xlabel('u(m)');ylabel('v(m)');colormap('gray');
    figure(3);
    plot(u,h2(M/2+1,:));xlabel('u(m)');ylabel('PSF');
    figure(4);
    plot(u,h2(:,M/2+1));xlabel('v(m)');ylabel('PSF');
    
    MTF=fft2(fftshift(h2));
    MTF=abs(MTF/MTF(1,1));
    MTF=ifftshift(MTF);
    
    MTF_an=(2/pi)*(acos(fu/twof0)-(fu/twof0).*sqrt(1-(fu/twof0).^2));
    MTF_an=MTF_an.*rect(fu/(2*twof0));
    
    figure(5);
    plot(fu,MTF(M/2+1,:),fu,MTF(:,M/2+1),':',fu,MTF_an,'--');
    axis([0 150000 0 1]);
    legend('u MTF','v MTF','diff limit');
    xlabel('f(cyc/m)');ylabel('Modulation');
    %}
%b/c) 
    M=1024;
    L=1e-3;
    du=L/M;
    u=-L/2:du:L/2-du; v=u;
    
    lambda=0.55e-6;
    k=2*pi/lambda;
    Dxp=20e-3;wxp=Dxp/2;
    zxp=100e-3;
    fnum=zxp/(2*wxp);
    lz=lambda*zxp;
    twof0=1/(lambda*fnum);
    coef=[0,0.1,0.25,0.5,0.75,1];
    u0=0;v0=0;

    wd=0*lambda;
    w040=coef.*lambda;
    w131=0*lambda;
    w222=0*lambda;
    w220=0*lambda;
    w311=0*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fu);
    Iu=[1:6];
    for i =1:numel(coef)
        W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040(i),w131,w222,w220,w311);
        
        H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
        h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    
        Iu(i)=max(max(h2));
    end
    
    figure(1);
    plot(coef,Iu(:)/Iu(1));
    %axis([0 150000 0 1]);
    xlabel('W040 Coefficient');ylabel('Strehl Ratio(S)');
%d)  
    hold on;
    k=1/lambda*10000;
    sigma=(1/12+1/180).*w040;
    plot(coef,exp(-k*k*sigma));
%e) 
    wd=0*lambda;
    w040=4.963*lambda;
    w131=2.637*lambda;
    w222=9.025*lambda;
    w220=7.536*lambda;
    w311=0.157*lambda;
    
    fu=-1/(2*du):1/L:1/(2*du)-(1/L);
    [Fu,Fv]=meshgrid(fu,fu);
    
    W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,wd,w040,w131,w222,w220,w311);
    
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);

    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    Ia=max(max(h2));
    S=Ia/Iu(1);