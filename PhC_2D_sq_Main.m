clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Plotting parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

band=1;
Field=1;
Epsilon=0;

AAbs=1;               %% Plot abs(E)
RReal=0;              %% Plot real(E)
IImag=0;              %% Plot imag(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TM=0;
TE=1;

Nx=32;         % number of points on the x grid % has to be a power of 2 (32,64,128,256,512,...)
Ny=32;         % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...)
NGx=10;        % number of harmonics % has to be 2 times -1 smaller than x
NGy=11;        % number of harmonics % has to be 2 times -1 smaller than y

Nkx=15;        % number of points on the k space for the dispersion
Nky=Nkx;       % number of points on the k space for the dispersion

nmodes=5;      % number of solutions asked
Np=2;          % number of period to plot for the Field

n1 =1;         %% optical index material 1
n2 = sqrt(12); %% optical index material 2

NormUnits=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Building of the index Geometry %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NormUnits==1
  L=1;
elseif NormUnits==0
  L=1e-6;  
end

Lx=L;
Ly=L;
x=linspace(-Lx/2, Lx/2, Nx);
y=linspace(-Ly/2, Ly/2, Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);

[XX,YY] = meshgrid(x,y);
a=0.38;
idx = (XX.^2 + YY.^2) < (a*L)^2;

eps = idx*n2^2 + (1-idx)*n1^2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Building Epsilon in Fourier space %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGx = 2*floor(NGx/2);           %% round to lower even number
NGy = 2*floor(NGy/2);           %% round to lower even number

Gamma=1./eps;
Gammak = fftshift(fft2(Gamma))*dx*dy/Lx/Ly;
Gammak = Gammak(Ny/2-NGy+1:Ny/2+NGy+1 , Nx/2-NGx+1:Nx/2+NGx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Reciprocal lattice vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gx = (-NGx/2:NGx/2)'*2*pi/Lx;
Gy = (-NGy/2:NGy/2)'*2*pi/Ly;

NGx=length(Gx);
NGy=length(Gy);
NG=NGx*NGy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Building of k-space vector %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kx=linspace( 0 , pi/Lx , Nkx);
ky=linspace( 0 , pi/Ly , Nky);

k=[
kx'                    kx'*0
ky'*0+kx(end)          ky'   
sort(kx,'descend')'    sort(ky,'descend')'
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (TE==1 && TM==1) || (TE==0 && TM==0)
  display('Error: Select "TM" or "TE"')
  break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Building first part of Hamitonian that is not depending on k %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HHH=zeros(NGy,NGx,NGy,NGx);

for ix=1:NGx
for jx=1:NGx
    for iy=1:NGy
    for jy=1:NGy
        HHH(iy,ix,jy,jx) = Gammak(iy-jy+NGy,ix-jx+NGx );
    end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k(:,1))
  
  [psi,f0]=PhC2D_sq_PWE_f(x,y,Gx,Gy,k(i,:),HHH,nmodes,TE,TM);
  
  E(:,:,:,i)=psi;
  
  if NormUnits==1
    FF(:,i) = f0 * Lx / (2*pi);
  elseif NormUnits==0
    FF(:,i) = f0 * c / (2*pi) *1e-12;     % Convert in THz
    lambda(:,i)=2*pi./f0*1e6;             % Convert in wavelength (um)
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if AAbs==1
  EE=abs(E);
end
if RReal==1
  EE=real(E);
end
if IImag==1
  EE=imag(E);
end

if NormUnits==0
  x=x*1e6;
  y=y*1e6;
  Lx=Lx*1e6;
  Ly=Ly*1e6;
  k=k*1e-6;
end

if Field==1
    
    if TM==1 && TE==0
      figure('position',[100 50 1000 1000],'name','Ez');
    elseif TE==1 && TM==0
      figure('position',[100 50 1000 1000],'name','Exy');
    end
    colormap(jet)
    
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,1+3*ii)
                hold on
                pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,1) )
                contour( (x+(i-1/2)*Lx), y+(j-1/2)*Ly ,abs(eps),1,'linewidth',2,'linecolor','w')
            end       
        end
        shading flat
        
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('\Gamma: w=' , num2str(FF(1+ii,1), '%.2f') ))
          xlabel('x (norm. units)')
          ylabel('y (norm. units)')
        elseif NormUnits==0 
          title(strcat('\Gamma: \lambda=' , num2str(lambda(1+ii,1), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 Np*Lx])
        ylim([0 Np*Ly])
    end
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,2+3*ii)
                hold on
                pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,1*Nkx) )
                contour( (x+(i-1/2)*Lx), y+(j-1/2)*Ly ,abs(eps),1,'linewidth',2,'linecolor','w')
            end
        end
        shading flat
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('X: w=' , num2str(FF(1+ii,length(k)/3), '%.2f') ))
          xlabel('x (norm. units)')
          ylabel('y (norm. units)')
        elseif NormUnits==0 
          title(strcat('X: \lambda=' , num2str(lambda(1+ii,length(k)/3), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 Np*Lx])
        ylim([0 Np*Ly])
    end
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,3+3*ii)
                hold on
                pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,2*Nkx) )
                contour( (x+(i-1/2)*Lx), y+(j-1/2)*Ly ,abs(eps),1,'linewidth',2,'linecolor','w')
            end
        end
        shading flat
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('M: w=' , num2str(FF(1+ii,length(k)*2/3), '%.2f') ))
          xlabel('x (norm. units)')
          ylabel('y (norm. units)')
        elseif NormUnits==0 
          title(strcat('M: \lambda=' , num2str(lambda(1+ii,length(k)*2/3), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 Np*Lx])
        ylim([0 Np*Ly])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Epsilon==1
  
  figure('position',[1100 50 500 400])
  subplot(111)
  hold on
     
  for i=1:Np
    for j=1:Np
        pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , real(eps) )
    end
  end
  shading flat
  colormap(jet)
  c=colorbar;
  title(c,'Epsilon')
  if NormUnits==1
    xlabel('x (norm. units)')
    ylabel('y (norm. units)')
  elseif NormUnits==0 
    xlabel('x (um)')
    ylabel('y (um)')
  end
  %axis equal

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if band==1
    
    figure('position',[50 570 900 450])
    
    subplot(1,2,1)
    hold on;%grid on;
    
    plot(0:length(k)-1,real(FF(1:nmodes,:))','o-')
    
    yscale=get(gca,'ylim');
    xlim([0 length(k)-1])
    
    plot( [1/3*length(k)    1/3*length(k)] , yscale , 'k')
    plot( [2/3*length(k)    2/3*length(k)] , yscale , 'k')
    plot( [3/3*length(k)    3/3*length(k)] , yscale , 'k')
    
    text(0/3*length(k) , -0.05*yscale(2) , ' \Gamma')
    text(1/3*length(k) , -0.05*yscale(2) , ' X'     )
    text(2/3*length(k) , -0.05*yscale(2) , ' M'     )
    text(3/3*length(k) , -0.05*yscale(2) , ' \Gamma')
    %xlabel('k')
    set(gca,'xticklabel',[])
        
    if NormUnits==1
      ylabel('w (2\pi/Ltot)')
    elseif NormUnits==0 
      ylabel('f (THz)')
    end
    title(strcat('R/a=',num2str(a),'; n1=',num2str(n1,'%.2f'),'; n2=',num2str(n2,'%.2f')  ))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(1,2,2)
    hold on;grid on;
    LW=2;
    
    for i=1:nmodes
      plot3( k(:,1)*Lx/pi , k(:,2)*Ly/pi , real(FF(i,:))','o-')
    end
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    xlabel('kx (pi/Lx)')
    ylabel('ky (pi/Ly)')
    if NormUnits==1
      zlabel('w (2\pi/Ltot)')
    elseif NormUnits==0 
      zlabel('f (THz)')
    end
    xlim([-1 1]*max(k(:,1))*Lx/pi)
    ylim([-1 1]*max(k(:,2))*Ly/pi)
    view(-30,15)
    
    plot3( [-1 1] , +[1 1] , [0 0] ,'r', 'linewidth',LW )
    plot3( [-1 1] , -[1 1] , [0 0] ,'r', 'linewidth',LW )
    plot3( +[1 1] , [-1 1] , [0 0] ,'r', 'linewidth',LW )
    plot3( -[1 1] , [-1 1] , [0 0] ,'r', 'linewidth',LW )
    
    plot3( [-1 1] , +[1 1] , [1 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    plot3( [-1 1] , -[1 1] , [1 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    plot3( +[1 1] , [-1 1] , [1 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    plot3( -[1 1] , [-1 1] , [1 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    
    plot3( +[1 1] , +[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    plot3( -[1 1] , -[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    plot3( +[1 1] , -[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    plot3( -[1 1] , +[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'r', 'linewidth',LW )
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%