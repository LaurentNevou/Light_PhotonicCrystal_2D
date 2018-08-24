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
Epsilon=1;

AAbs=1;               %% Plot abs(E)
RReal=0;              %% Plot real(E)
IImag=0;              %% Plot imag(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TM=1;
TE=0;

Nx=64;         % number of points on the x grid % has to be a power of 2 (32,64,128,256,512,...)
Ny=64;         % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...)
NGx=20;        % number of harmonics % has to be 2 times -1 smaller than x
NGy=16;        % number of harmonics % has to be 2 times -1 smaller than y

Nkx=15;        % number of points on the k space for the dispersion

nmodes=5;      % number of solutions asked
Np=2;          % number of period to plot for the Field

n1 =1;         %% optical index material 1
n2 = sqrt(12); %% optical index material 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Building of the index Geometry %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx=52.4e-6;              % map in x direction
Ly=20e-6;                % map in y direction

x=linspace(-Lx/2, Lx/2, Nx);
y=linspace(-Ly/2, Ly/2, Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);

[XX,YY] = meshgrid(x,y);


idx1=(abs(XX)< 39.4e-6/2) .* (abs(YY)< 7.5e-6/2) ;
idx2=(abs(XX)< 1) .* (abs(YY)< 2.5e-6/2) ;
idx3=(abs(YY)> 18e-6/2) ;

idx = idx1|idx2;

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

if Nkx==1
    kx=0;
else
    kx=linspace( 0 , pi/Lx , Nkx)';
end

k=[kx kx*0];

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
  FF(:,i) = f0 * c / (2*pi) *1e-12;     % Convert in THz
  lambda(:,i)=2*pi./f0*1e6;             % Convert in wavelength (um)
  
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

x=x*1e6;
y=y*1e6;
Lx=Lx*1e6;
Ly=Ly*1e6;
k=k*1e-6;

if Field==1
    
    if TM==1 && TE==0
      figure('position',[100 50 1000 1000],'name','Ez');
    elseif TE==1 && TM==0
      figure('position',[100 50 1000 1000],'name','Exy');
    end
    colormap(jet)
    
    for ii=0:nmodes-1
        for i=1:Np
          subplot(nmodes,2,1+2*ii)
          hold on
          pcolor( x+(i-1/2)*Lx , y+Ly/2 , EE(:,:,ii+1,1) )
          contour( (x+(i-1/2)*Lx), y+Ly/2 ,abs(eps),1,'linewidth',2,'linecolor','w')
        end
        shading flat
        
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        
        title(strcat('\Gamma: \lambda=' , num2str(lambda(1+ii,1), '%.2f') , 'um' ))
        xlabel('x (um)')
        ylabel('y (um)')
        xlim([0 Np*Lx])
        ylim([0 Ly])
    end
    for ii=0:nmodes-1
        for i=1:Np
          subplot(nmodes,2,2+2*ii)
          hold on
          pcolor( x+(i-1/2)*Lx , y+Ly/2 , EE(:,:,ii+1,end) )
          contour( (x+(i-1/2)*Lx), y+Ly/2 ,abs(eps),1,'linewidth',2,'linecolor','w')
        end
        shading flat
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        
        title(strcat('X: \lambda=' , num2str(lambda(1+ii,end), '%.2f') , 'um' ))
        xlabel('x (um)')
        ylabel('y (um)')
        xlim([0 Np*Lx])
        ylim([0 Ly])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Epsilon==1
  
  figure('position',[1100 50 600 300])
  subplot(111)
  hold on
  
  for i=1:Np
    pcolor( x+(i-1/2)*Lx , y+Ly/2 , real(eps) )
  end
  
  shading flat
  colormap(jet)
  c=colorbar;
  title(c,'Epsilon')
  xlabel('x (um)')
  ylabel('y (um)')
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if band==1
    
    figure('position',[50 570 500 450])
    
    subplot(111)
    hold on;%grid on;
    
    plot(k(:,1)*Lx/pi,real(FF(1:nmodes,:))','o-')
    
    yscale=get(gca,'ylim');
    xlim([0 1])
    xlabel('k')
    ylabel('f (THz)')
    title(strcat(' n1=',num2str(n1,'%.2f'),'; n2=',num2str(n2,'%.2f')  ))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%