
%% FDTD 3D Maxwell's equations
%%%% Problem Description

%%%% Media:
%%%%       Vaccum  
%%%% Boundaries: 
%%%%            ABC:Absorbing Boundary Condition
%%%%            Liao's 4th Order
%%%%            Non-regular Grid

%%%% by: Esteban Jiménez Mejía
%%%% resteban.jimenez@gmail.com

clear all
close all
clc


Mu_Zero=4*pi*1e-7;
E_Zero=8.854187817e-12;
c=1/sqrt(Mu_Zero*E_Zero);


%% minimum delta values (to be in accordance with the paper)
deltax=1;
deltay=1;
deltaz=1;


%%%% Total distances in meters
Lx=500;
Ly=1000;
Lz=2000;


%%%% the space discretization, however is can be non-uniformly distributed
x=0:deltax:Lx;      %% space discretization
y=0:deltay:Ly;      %% space discretization
z=0:deltaz:Lz;      %% space discretization


x=[0:deltax:5 Geometric_Mesh(5+deltax,245,10,deltax,2) 245+deltax:deltax:255 Geometric_Mesh(255+deltax,495,10,deltax,2) 495+deltax:deltax:Lx];
y=[0:deltay:5 Geometric_Mesh(5+deltay,95,10,deltay,2) 95+deltay:deltay:105 Geometric_Mesh(105+deltay,995,10,deltay,2) 995+deltay:deltay:Ly];
z=[0:deltaz:30 Geometric_Mesh(30+deltaz,1995,10,deltaz,2) 1995+deltaz:deltaz:Lz]; 



%%%%%% space discretization visualization
[xx1,yy1]=meshgrid(x,y);
[xx2,zz2]=meshgrid(x,z);


figure(30)
mesh(xx1,yy1,zeros(size(xx1)))
view(0,90)
xlabel('X [m]','Fontsize',14)
ylabel('Y [m]','Fontsize',14)
set(gca,'Fontsize',14)


figure(31)
mesh(xx2,zz2,zeros(size(xx2)))
view(0,90)
xlabel('X [m]','Fontsize',14)
ylabel('Z [m]','Fontsize',14)
set(gca,'Fontsize',14)


NDK=length(x)-1;
NDL=length(y)-1;
NDZ=length(z)-1;

dx=x(2:NDK+1)-x(1:NDK);     %%%the delta X values can be irregular sized
dy=y(2:NDL+1)-y(1:NDL);     %%%the delta Y values can be irregular sized
dz=z(2:NDZ+1)-z(1:NDZ);     %%%the delta Z values can be irregular sized

dt=0.9*(1/(c*sqrt((1/(min(dx)^2))+(1/(min(dy)^2))+(1/(min(dz)^2)))));                %% minimum time discretization (Courant's factor)

tsim=0:dt:2000*dt;

Ex=zeros(NDK,NDL+1,NDZ+1);  %% 3D Matrix
Ey=zeros(NDK+1,NDL,NDZ+1);  %% 3D Matrix
Ez=zeros(NDK+1,NDL+1,NDZ);  %% 3D Matrix

Hx=zeros(NDK+1,NDL,NDZ);  %% 3D Matrix
Hy=zeros(NDK,NDL+1,NDZ);  %% 3D Matrix
Hz=zeros(NDK,NDL,NDZ+1);  %% 3D Matrix


%%%%% Displacement Currents
Jx=zeros(NDK,NDL+1,NDZ+1);  %% 3D Matrix
Jy=zeros(NDK+1,NDL,NDZ+1);  %% 3D Matrix
Jz=zeros(NDK+1,NDL+1,NDZ);  %% 3D Matrix

%%%%% Delta Spaces for Electric Fields
DX=((x(2:NDK+1)-x(1:NDK)));
DY=((y(2:NDL+1)-y(1:NDL)));
DZ=((z(2:NDZ+1)-z(1:NDZ)));


%%%%% Deltas for Hx calculation
DZ_Ey=zeros(NDK+1,NDL,NDZ);
for i=1:length(DZ)
    DZ_Ey(:,:,i)=ones(NDK+1,NDL).*DZ(i);
end

DY_Ez=zeros(NDK+1,NDL,NDZ);
for i=1:length(DY)
    DY_Ez(:,i,:)=ones(NDK+1,NDZ).*DY(i);
end

%%%%% Deltas for Hy calculation

DX_Ez=zeros(NDK,NDL+1,NDZ);
for i=1:length(DX)
    DX_Ez(i,:,:)=ones(NDL+1,NDZ).*DX(i);
end

DZ_Ex=zeros(NDK,NDL+1,NDZ);
for i=1:length(DZ)
    DZ_Ex(:,:,i)=ones(NDK,NDL+1).*DZ(i);
end

%%%%% Deltas for Hz calculation
DY_Ex=zeros(NDK,NDL,NDZ+1);
for i=1:length(DY)
    DY_Ex(:,i,:)=ones(NDK,NDZ+1).*DY(i);
end

DX_Ey=zeros(NDK,NDL,NDZ+1);
for i=1:length(DX)
    DX_Ey(i,:,:)=ones(NDL,NDZ+1).*DX(i);
end

%%%%% Delta Spaces for Magnetic Fields
x_H=x(1:NDK)+0.5*DX;
y_H=y(1:NDL)+0.5*DY;
z_H=z(1:NDZ)+0.5*DZ;

DX_H=x_H(2:NDK)-x_H(1:NDK-1);
DY_H=y_H(2:NDL)-y_H(1:NDL-1);
DZ_H=z_H(2:NDZ)-z_H(1:NDZ-1);


%%%%% Deltas for Ex calculation
DZ_Hy=zeros(NDK,NDL-1,NDZ-1);
for i=1:length(DZ_H)
    DZ_Hy(:,:,i)=ones(NDK,NDL-1).*DZ_H(i);
end

DY_HZ=zeros(NDK,NDL-1,NDZ-1);
for i=1:length(DY_H)

   DY_HZ(:,i,:)=ones(NDK,NDZ-1).*DY_H(i);
end

%%%%% Deltas for Ey calculation
DX_Hz=zeros(NDK-1,NDL,NDZ-1);
for i=1:length(DX_H)
    DX_Hz(i,:,:)=ones(NDL,NDZ-1).*DX_H(i);
end

DZ_Hx=zeros(NDK-1,NDL,NDZ-1);
for i=1:length(DZ_H)
    DZ_Hx(:,:,i)=ones(NDK-1,NDL).*DZ_H(i);
end


%%%%% Deltas for Ez calculation
DX_Hy=zeros(NDK-1,NDL-1,NDZ);
for i=1:length(DX_H)
    DX_Hy(i,:,:)=ones(NDL-1,NDZ).*DX_H(i);
end


DY_Hx=zeros(NDK-1,NDL-1,NDZ);
for i=1:length(DY_H)
    DY_Hx(:,i,:)=ones(NDK-1,NDZ).*DY_H(i);
end

%%% Conductive Material Grid Components
sigma_mx=zeros(size(Ex));  %% 3D Matrix
sigma_my=zeros(size(Ey));  %% 3D Matrix
sigma_mz=zeros(size(Ez));  %% 3D Matrix

%%% Permitivity Material Grid Components
E_Zero_mx=E_Zero*ones(size(Ex));  %% 3D Matrix
E_Zero_my=E_Zero*ones(size(Ey));  %% 3D Matrix
E_Zero_mz=E_Zero*ones(size(Ez));  %% 3D Matrix

%%% Permeability Material Grid Components
Mu_Zero_mx=Mu_Zero*ones(size(Hx));  %% 3D Matrix
Mu_Zero_my=Mu_Zero*ones(size(Hy));  %% 3D Matrix
Mu_Zero_mz=Mu_Zero*ones(size(Hz));  %% 3D Matrix




%% Thin wire Location & Inclusion (coordinates of the initial wire in meters)
%%%%% vertical wire
xw1_ini=250;
yw1_ini=100;
zw1_ini=11;

xw1_end=250;
yw1_end=100;
zw1_end=1600;

iw1_ini=find(x==xw1_ini);       %x location of z-directed voltage source
jw1_ini=find(y==yw1_ini);       %y location of z-directed voltage source
kw1_ini=find(z==zw1_ini);       %z location of z-directed voltage source (initial z cell for the wire)

iw1_end=find(x==xw1_end);       %x location of z-directed voltage source
jw1_end=find(y==yw1_end);       %y location of z-directed voltage source
kw1_end=find(z==zw1_end);       %z location of z-directed voltage source (initial z cell for the wire)

kw1_end=181;

%% Source Conditions
V=5e6;
rise_time=1e-6;
a=0.3;
m=V/rise_time;        
       for ii = 1 : length(tsim),
             if tsim(ii) <= rise_time   
                    ramp(ii) = m*tsim(ii);
             else
                    ramp(ii) = V;
             end
        end 

        I1=(2/(120*pi))*atan(pi./(2*log(sqrt(3e8^2*tsim.^2)/a)));
        I0=2*conv(diff(ramp),I1);
        I =real(I0(1:length(tsim)));
% % %
%%%%%%%%%% Heidler 
% % % % % Io1=9.9e3;
% % % % % Nu_1=0.845;
% % % % % n=2;
% % % % % Tao_1=0.072e-6;
% % % % % Tao_2=5e-6;
% % % % % Io2=7.5e3;
% % % % % Tao_3=100e-6;
% % % % % Tao_4=6e-6;
% % % % % 
% % % % % I=(Io1/Nu_1)*(((tsim/Tao_1).^n)./(1+(tsim/Tao_1).^n)).*exp(-tsim/Tao_2) + Io2*(exp(-tsim/Tao_3)-exp(-tsim/Tao_4));
% % 

   %%% Active Source Plotting
    figure(1)
    plot(tsim*1e6,I/1e3,'linewidth',2)
    title({'Source Excitation','Current along z-axis'},'Fontsize',14)
    xlabel('Time [\mus]','Fontsize',14)
    ylabel('Current [kA]','Fontsize',14)
    set(gca,'Fontsize',14)
    
    is=iw1_ini;
    js=jw1_ini;
    ksource=kw1_ini-1;       %% one cell less than the first cell of the wire
    
%%%%% Internal resistance of the source
Rs=50e12;       %%% [resistance in ohms]

    %%% Z directed Electric Field
    Cez1=(1+((dt*DZ(ksource))/(2*E_Zero*DX_H(is-1)*DY_H(js-1)*Rs)));
    Cez2=(1-((dt*DZ(ksource))/(2*E_Zero*DX_H(is-1)*DY_H(js-1)*Rs)));
    
    Vsz=I*Rs;

    %%% Active Source Plotting
    figure(2)
    plot(tsim*1e6,Vsz/1e3,'linewidth',2)
    title({'Source Excitation','Gaussian Voltage along z-axis'},'Fontsize',14)
    xlabel('Time [ns]','Fontsize',14)
    ylabel('Volts [kV]','Fontsize',14)
    set(gca,'Fontsize',14)
    
     

%% Thin Wires Electric Fields 
    %%%Vertical thin wire inlcusion

    WlNDZ=kw1_end-kw1_ini;       %%% length in cells units of wire
    a2=0.3;       %%% conductor radii lead-wire
    m_factor2=1.471/log(DX_H(iw1_ini-1)/a2);       %%%% cálculo del factor para modificación de los parámetros circundantes del medio
    
    Ez(iw1_ini,jw1_ini,kw1_ini:kw1_end-1)=0; %%% internal electric field equal to zero

    E_Zero_mx(iw1_ini-1:iw1_ini,jw1_ini,kw1_ini:kw1_end)=(m_factor2*E_Zero);
    E_Zero_my(iw1_ini,jw1_ini-1:jw1_ini,kw1_ini:kw1_end)=(m_factor2*E_Zero);
% 
    Mu_Zero_mx(iw1_ini,jw1_ini-1:jw1_ini,kw1_ini:kw1_end)=(Mu_Zero/m_factor2);
    Mu_Zero_my(iw1_ini-1:iw1_ini,jw1_ini,kw1_ini:kw1_end)=(Mu_Zero/m_factor2);


%% Conductive Soil Inclusion
sigma_mx(:,:,1:ksource)=1e16;
sigma_my(:,:,1:ksource)=1e16;
sigma_mz(:,:,1:ksource-1)=1e16;


%%%%% Variables Initialization
    Vini=zeros(1,length(tsim));

    
    Iini=zeros(1,length(tsim));
    Imiddle=zeros(1,length(tsim));
    Iend=zeros(1,length(tsim));

    I1=zeros(1,length(tsim));
    I2=zeros(1,length(tsim));
    I3=zeros(1,length(tsim));
    I4=zeros(1,length(tsim));
    I5=zeros(1,length(tsim));
   
    %%%%% Electric Fields Calculation
    %%%% Position 1
    Exp1=zeros(1,length(tsim));
    Eyp1=zeros(1,length(tsim));
    Ezp1=zeros(1,length(tsim));

    %%%% Position 2
    Exp2=zeros(1,length(tsim));
    Eyp2=zeros(1,length(tsim));
    Ezp2=zeros(1,length(tsim));

    %%%% Position 3
    Exp3=zeros(1,length(tsim));
    Eyp3=zeros(1,length(tsim));
    Ezp3=zeros(1,length(tsim));
   
    %%%% Position 4
    Exp4=zeros(1,length(tsim));
    Eyp4=zeros(1,length(tsim));
    Ezp4=zeros(1,length(tsim));
    
    
%%%%% Current Calculation Determination
IalZ=zeros(1,WlNDZ);

 Ez_past=Ez(is,js,ksource);
 Vsource_t=zeros(1,length(tsim));

%%%% for LIAO´s boundaries
Ex_y2=zeros(NDK,NDZ+1,5); 
Ex_y3=zeros(NDK,NDZ+1,5); 
Ex_y4=zeros(NDK,NDZ+1,5);
Ex_y5=zeros(NDK,NDZ+1,5);

Ex_yNDL_3=zeros(NDK,NDZ+1,5);
Ex_yNDL_2=zeros(NDK,NDZ+1,5);
Ex_yNDL_1=zeros(NDK,NDZ+1,5);
Ex_yNDL=zeros(NDK,NDZ+1,5);

Ex_z2=zeros(NDK,NDL+1,5); 
Ex_z3=zeros(NDK,NDL+1,5); 
Ex_z4=zeros(NDK,NDL+1,5); 
Ex_z5=zeros(NDK,NDL+1,5); 

Ex_zNDZ_3=zeros(NDK,NDL+1,5); 
Ex_zNDZ_2=zeros(NDK,NDL+1,5); 
Ex_zNDZ_1=zeros(NDK,NDL+1,5); 
Ex_zNDZ=zeros(NDK,NDL+1,5); 

Ey_x2=zeros(NDL,NDZ+1,5);
Ey_x3=zeros(NDL,NDZ+1,5);
Ey_x4=zeros(NDL,NDZ+1,5);
Ey_x5=zeros(NDL,NDZ+1,5);

Ey_xNDK_3=zeros(NDL,NDZ+1,5);
Ey_xNDK_2=zeros(NDL,NDZ+1,5);
Ey_xNDK_1=zeros(NDL,NDZ+1,5);
Ey_xNDK=zeros(NDL,NDZ+1,5);

Ey_z2=zeros(NDK+1,NDL,5);
Ey_z3=zeros(NDK+1,NDL,5);
Ey_z4=zeros(NDK+1,NDL,5);
Ey_z5=zeros(NDK+1,NDL,5);

Ey_zNDZ_3=zeros(NDK+1,NDL,5);
Ey_zNDZ_2=zeros(NDK+1,NDL,5);
Ey_zNDZ_1=zeros(NDK+1,NDL,5);
Ey_zNDZ=zeros(NDK+1,NDL,5);

Ez_x2=zeros(NDL+1,NDZ,5);
Ez_x3=zeros(NDL+1,NDZ,5);
Ez_x4=zeros(NDL+1,NDZ,5);
Ez_x5=zeros(NDL+1,NDZ,5);

Ez_xNDK_3=zeros(NDL+1,NDZ,5);
Ez_xNDK_2=zeros(NDL+1,NDZ,5);
Ez_xNDK_1=zeros(NDL+1,NDZ,5);
Ez_xNDK=zeros(NDL+1,NDZ,5);

Ez_y2=zeros(NDK+1,NDZ,5);
Ez_y3=zeros(NDK+1,NDZ,5);
Ez_y4=zeros(NDK+1,NDZ,5);
Ez_y5=zeros(NDK+1,NDZ,5);

Ez_yNDL_3=zeros(NDK+1,NDZ,5);
Ez_yNDL_2=zeros(NDK+1,NDZ,5);
Ez_yNDL_1=zeros(NDK+1,NDZ,5);
Ez_yNDL=zeros(NDK+1,NDZ,5);

Jz(iw1_ini,jw1_ini,kw1_ini:kw1_end-1)=0; %%% initialization of the Displacement Current

Rind=1*DZ(kw1_ini:kw1_end-1);
Lind=3e-6*DZ(kw1_ini:kw1_end-1);

C1jz(1,1,:)=-(((Rind/2)-(Lind/dt))./((Rind/2)+(Lind/dt)));
C2jz(1,1,:)=(DZ(kw1_ini:kw1_end-1)./((DX_H(iw1_ini-1)).*(DY_H(jw1_ini-1)))).*(1./((Rind/2)+(Lind/dt)));

for t=1:length(tsim)
  
    
    %%%% RL load Updating
    Jz(iw1_ini,jw1_ini,kw1_ini:kw1_end-1)=C1jz(1,1,:).*Jz(iw1_ini,jw1_ini,kw1_ini:kw1_end-1)+C2jz.*Ez(iw1_ini,jw1_ini,kw1_ini:kw1_end-1);
%     
%% INIT Electromagnetic Fields Calculation ---------------
     %%% Magnetic Field Calculation
    Hx(:,:,:)=Hx(:,:,:)+(dt./Mu_Zero_mx).*(1./DZ_Ey).*(Ey(:,:,2:NDZ+1)-Ey(:,:,1:NDZ))-(dt./Mu_Zero_mx).*(1./DY_Ez).*(Ez(:,2:NDL+1,:)-Ez(:,1:NDL,:));
    Hy(:,:,:)=Hy(:,:,:)+(dt./Mu_Zero_my).*(1./DX_Ez).*(Ez(2:NDK+1,:,:)-Ez(1:NDK,:,:))-(dt./Mu_Zero_my).*(1./DZ_Ex).*(Ex(:,:,2:NDZ+1)-Ex(:,:,1:NDZ));
    Hz(:,:,:)=Hz(:,:,:)+(dt./Mu_Zero_mz).*(1./DY_Ex).*(Ex(:,2:NDL+1,:)-Ex(:,1:NDL,:))-(dt./Mu_Zero_mz).*(1./DX_Ey).*(Ey(2:NDK+1,:,:)-Ey(1:NDK,:,:));
    
    Bx=Mu_Zero.*Hx;
    By=Mu_Zero.*Hy;
    Bz=Mu_Zero.*Hz;
    
    %%% Electric Field Calculation
    Ex(:,2:NDL,2:NDZ)=((2*E_Zero_mx(:,2:NDL,2:NDZ)-dt*sigma_mx(:,2:NDL,2:NDZ))./(2*E_Zero_mx(:,2:NDL,2:NDZ)+dt*sigma_mx(:,2:NDL,2:NDZ))).*Ex(:,2:NDL,2:NDZ)+(2*dt./(2*E_Zero_mx(:,2:NDL,2:NDZ)+dt*sigma_mx(:,2:NDL,2:NDZ))).*(1./DY_HZ).*(Hz(:,2:NDL,2:NDZ)-Hz(:,1:NDL-1,2:NDZ))-(2*dt./(2*E_Zero_mx(:,2:NDL,2:NDZ)+dt*sigma_mx(:,2:NDL,2:NDZ))).*(1./DZ_Hy).*(Hy(:,2:NDL,2:NDZ)-Hy(:,2:NDL,1:NDZ-1))-(2*dt./(2*E_Zero_mx(:,2:NDL,2:NDZ)+dt*sigma_mx(:,2:NDL,2:NDZ))).*Jx(:,2:NDL,2:NDZ);
    Ey(2:NDK,:,2:NDZ)=((2*E_Zero_my(2:NDK,:,2:NDZ)-dt*sigma_my(2:NDK,:,2:NDZ))./(2*E_Zero_my(2:NDK,:,2:NDZ)+dt*sigma_my(2:NDK,:,2:NDZ))).*Ey(2:NDK,:,2:NDZ)+(2*dt./(2*E_Zero_my(2:NDK,:,2:NDZ)+dt*sigma_my(2:NDK,:,2:NDZ))).*(1./DZ_Hx).*(Hx(2:NDK,:,2:NDZ)-Hx(2:NDK,:,1:NDZ-1))-(2*dt./(2*E_Zero_my(2:NDK,:,2:NDZ)+dt*sigma_my(2:NDK,:,2:NDZ))).*(1./DX_Hz).*(Hz(2:NDK,:,2:NDZ)-Hz(1:NDK-1,:,2:NDZ))-(2*dt./(2*E_Zero_my(2:NDK,:,2:NDZ)+dt*sigma_my(2:NDK,:,2:NDZ))).*Jy(2:NDK,:,2:NDZ);
    Ez(2:NDK,2:NDL,:)=((2*E_Zero_mz(2:NDK,2:NDL,:)-dt*sigma_mz(2:NDK,2:NDL,:))./(2*E_Zero_mz(2:NDK,2:NDL,:)+dt*sigma_mz(2:NDK,2:NDL,:))).*Ez(2:NDK,2:NDL,:)+(2*dt./(2*E_Zero_mz(2:NDK,2:NDL,:)+dt*sigma_mz(2:NDK,2:NDL,:))).*(1./DX_Hy).*(Hy(2:NDK,2:NDL,:)-Hy(1:NDK-1,2:NDL,:))-(2*dt./(2*E_Zero_mz(2:NDK,2:NDL,:)+dt*sigma_mz(2:NDK,2:NDL,:))).*(1./DY_Hx).*(Hx(2:NDK,2:NDL,:)-Hx(2:NDK,1:NDL-1,:))-(2*dt./(2*E_Zero_mz(2:NDK,2:NDL,:)+dt*sigma_mz(2:NDK,2:NDL,:))).*Jz(2:NDK,2:NDL,:);
        
 %% END Electromagnetic Fields Calculation ---------------
 
 
 
 %% INIT Source Updating Calculation ---------------
 
    %%% Source Updating   
    Ez(is,js,ksource)=(Cez2/Cez1)*Ez_past+(dt/(Cez1*Mu_Zero*E_Zero*DX_H(is-1)))*(By(is,js,ksource)-By(is-1,js,ksource))-(dt/(Cez1*Mu_Zero*E_Zero*DY_H(js-1)))*(Bx(is,js,ksource)-Bx(is,js-1,ksource))-(dt/(Cez1*E_Zero*DX_H(is-1)*DY_H(js-1)*Rs))*Vsz(t);
    Ez_past=Ez(is,js,ksource);       %%%source Ez at the last time step to take it into account for the next step
    
    %%% Vsource Calculation just with the source cell.
    Vsource_t(t)=-Ez(is,js,ksource)*DZ(ksource);
 
% % %  %% Forcing internal electric fields to be zero
% % %     Ez(iw1_ini,jw1_ini,kw1_ini:kw1_end-1)=0; %%% internal electric field equal to zero
    
    
    %%%% Electric Fields Measurement
    %%%50m from the channel base
    Exp1(t)=Ex(iw1_ini,jw1_ini+3,20);
    Eyp1(t)=Ey(iw1_ini,jw1_ini+3,20);
    Ezp1(t)=Ez(iw1_ini,jw1_ini+3,20);

% %     %%%% 100m from the channel Base
    Exp2(t)=Ex(iw1_ini,jw1_ini+6,20);
    Eyp2(t)=Ey(iw1_ini,jw1_ini+6,20);
    Ezp2(t)=Ez(iw1_ini,jw1_ini+6,20);
% %     
% % %%%% 151m from the channel Base
    Exp3(t)=Ex(iw1_ini,jw1_ini+9,20);
    Eyp3(t)=Ey(iw1_ini,jw1_ini+9,20);
    Ezp3(t)=Ez(iw1_ini,jw1_ini+9,20);
% %     
% %     %%%% 201m from the channel Base
    Exp4(t)=Ex(iw1_ini,jw1_ini+12,20);
    Eyp4(t)=Ey(iw1_ini,jw1_ini+12,20);
    Ezp4(t)=Ez(iw1_ini,jw1_ini+12,20);
% %     
 
 %% currents and voltage calculation
   Iini(t)=(DX_H(is-1)./Mu_Zero_my(is,js,ksource)).*(By(is,js,ksource)-By(is-1,js,ksource))+(DY_H(js-1)./Mu_Zero_mx(is,js,ksource)).*(Bx(is,js-1,ksource)-Bx(is,js,ksource));
      
   %%%% Assessment of the current along the z-wire
   IalZ(1,:)=(DX_H(is-1)./Mu_Zero_my(is,js,kw1_ini:kw1_end-1)).*(By(is,js,kw1_ini:kw1_end-1)-By(is-1,js,kw1_ini:kw1_end-1))+(DY_H(js-1)./Mu_Zero_mx(is,js,kw1_ini:kw1_end-1)).*(Bx(is,js-1,kw1_ini:kw1_end-1)-Bx(is,js,kw1_ini:kw1_end-1));
   
%    Imiddle(t)=IalZ(round(length(IalZ)/2));
   Imiddle(t)=(DX_H(is-1)./Mu_Zero_my(is,js,ksource+4)).*(By(is,js,ksource+4)-By(is-1,js,ksource+4))+(DY_H(js-1)./Mu_Zero_mx(is,js,ksource+4)).*(Bx(is,js-1,ksource+4)-Bx(is,js,ksource+4));
   Iend(t)=IalZ(length(IalZ));
   
   
   [~,I]=min(abs(z-(2+6)));
   I1(t)=IalZ(find(z==z(I))-6);
   
   [~,I]=min(abs(z-(20+6)));
   I2(t)=IalZ(find(z==z(I))-6);
   
   [~,I]=min(abs(z-(60+6)));
   I3(t)=IalZ(find(z==z(I))-6);
   
   [~,I]=min(abs(z-(100+6)));
   I4(t)=IalZ(find(z==z(I))-6);
   
   [~,I]=min(abs(z-(500+6)));
   I5(t)=IalZ(find(z==z(I))-6);
   
   %%%% Voltage at the beginnig  and end of the Line
    Ez_ini(1,:)=Ez(iw1_ini,jw1_ini,ksource:kw1_end-1);
    Vini(t)=-sum(Ez_ini.*dz(ksource:kw1_end-1));
    
    
 
 %%
%%%% For Ex Boundaries Refreshing
%%%%%%% this operation is saving the past time values at each space
%%%%%%% positions
          
           Ex_y2(:,:,5)=Ex_y2(:,:,4);
           Ex_y3(:,:,5)=Ex_y3(:,:,4);
           Ex_y4(:,:,5)=Ex_y4(:,:,4);
           Ex_y5(:,:,5)=Ex_y5(:,:,4);

           Ex_y2(:,:,4)=Ex_y2(:,:,3);
           Ex_y3(:,:,4)=Ex_y3(:,:,3);
           Ex_y4(:,:,4)=Ex_y4(:,:,3);
           Ex_y5(:,:,4)=Ex_y5(:,:,3);

           Ex_y2(:,:,3)=Ex_y2(:,:,2);
           Ex_y3(:,:,3)=Ex_y3(:,:,2);
           Ex_y4(:,:,3)=Ex_y4(:,:,2);
           Ex_y5(:,:,3)=Ex_y5(:,:,2);

           Ex_y2(:,:,2)=Ex_y2(:,:,1);
           Ex_y3(:,:,2)=Ex_y3(:,:,1);
           Ex_y4(:,:,2)=Ex_y4(:,:,1);
           Ex_y5(:,:,2)=Ex_y5(:,:,1);

           Ex_y2(:,:,1)=Ex(:,2,:);
           Ex_y3(:,:,1)=Ex(:,3,:);
           Ex_y4(:,:,1)=Ex(:,4,:);
           Ex_y5(:,:,1)=Ex(:,5,:);


           Ex_yNDL(:,:,5)=Ex_yNDL(:,:,4);
           Ex_yNDL_1(:,:,5)=Ex_yNDL_1(:,:,4);
           Ex_yNDL_2(:,:,5)=Ex_yNDL_2(:,:,4);
           Ex_yNDL_3(:,:,5)=Ex_yNDL_3(:,:,4);
           
           Ex_yNDL(:,:,4)=Ex_yNDL(:,:,3);
           Ex_yNDL_1(:,:,4)=Ex_yNDL_1(:,:,3);
           Ex_yNDL_2(:,:,4)=Ex_yNDL_2(:,:,3);
           Ex_yNDL_3(:,:,4)=Ex_yNDL_3(:,:,3);
           
           Ex_yNDL(:,:,3)=Ex_yNDL(:,:,2);
           Ex_yNDL_1(:,:,3)=Ex_yNDL_1(:,:,2);
           Ex_yNDL_2(:,:,3)=Ex_yNDL_2(:,:,2);
           Ex_yNDL_3(:,:,3)=Ex_yNDL_3(:,:,2);

           Ex_yNDL(:,:,2)=Ex_yNDL(:,:,1);
           Ex_yNDL_1(:,:,2)=Ex_yNDL_1(:,:,1);
           Ex_yNDL_2(:,:,2)=Ex_yNDL_2(:,:,1);
           Ex_yNDL_3(:,:,2)=Ex_yNDL_3(:,:,1);

           Ex_yNDL(:,:,1)=Ex(:,NDL,:);
           Ex_yNDL_1(:,:,1)=Ex(:,NDL-1,:);
           Ex_yNDL_2(:,:,1)=Ex(:,NDL-2,:);
           Ex_yNDL_3(:,:,1)=Ex(:,NDL-3,:);
        
        
            Ex_z2(:,:,5)=Ex_z2(:,:,4);
            Ex_z3(:,:,5)=Ex_z3(:,:,4);
            Ex_z4(:,:,5)=Ex_z4(:,:,4);
            Ex_z5(:,:,5)=Ex_z5(:,:,4);
           
           
            Ex_z2(:,:,4)=Ex_z2(:,:,3);
            Ex_z3(:,:,4)=Ex_z3(:,:,3);
            Ex_z4(:,:,4)=Ex_z4(:,:,3);
            Ex_z5(:,:,4)=Ex_z5(:,:,3);
           
            Ex_z2(:,:,3)=Ex_z2(:,:,2);
            Ex_z3(:,:,3)=Ex_z3(:,:,2);
            Ex_z4(:,:,3)=Ex_z4(:,:,2);
            Ex_z5(:,:,3)=Ex_z5(:,:,2);
            
            Ex_z2(:,:,2)=Ex_z2(:,:,1);
            Ex_z3(:,:,2)=Ex_z3(:,:,1);
            Ex_z4(:,:,2)=Ex_z4(:,:,1);
            Ex_z5(:,:,2)=Ex_z5(:,:,1);

            Ex_z2(:,:,1)=Ex(:,:,2);
            Ex_z3(:,:,1)=Ex(:,:,3);
            Ex_z4(:,:,1)=Ex(:,:,4);
            Ex_z5(:,:,1)=Ex(:,:,5);
        

            
            Ex_zNDZ_3(:,:,5)=Ex_zNDZ_3(:,:,4);
            Ex_zNDZ_2(:,:,5)=Ex_zNDZ_2(:,:,4);
            Ex_zNDZ_1(:,:,5)=Ex_zNDZ_1(:,:,4);
            Ex_zNDZ(:,:,5)=Ex_zNDZ(:,:,4);
            
            Ex_zNDZ_3(:,:,4)=Ex_zNDZ_3(:,:,3);
            Ex_zNDZ_2(:,:,4)=Ex_zNDZ_2(:,:,3);
            Ex_zNDZ_1(:,:,4)=Ex_zNDZ_1(:,:,3);
            Ex_zNDZ(:,:,4)=Ex_zNDZ(:,:,3);
            
            Ex_zNDZ_3(:,:,3)=Ex_zNDZ_3(:,:,2);
            Ex_zNDZ_2(:,:,3)=Ex_zNDZ_2(:,:,2);
            Ex_zNDZ_1(:,:,3)=Ex_zNDZ_1(:,:,2);
            Ex_zNDZ(:,:,3)=Ex_zNDZ(:,:,2);
            
            Ex_zNDZ_3(:,:,2)=Ex_zNDZ_3(:,:,1);
            Ex_zNDZ_2(:,:,2)=Ex_zNDZ_2(:,:,1);
            Ex_zNDZ_1(:,:,2)=Ex_zNDZ_1(:,:,1);
            Ex_zNDZ(:,:,2)=Ex_zNDZ(:,:,1);

            Ex_zNDZ_3(:,:,1)=Ex(:,:,NDZ-3);
            Ex_zNDZ_2(:,:,1)=Ex(:,:,NDZ-2);
            Ex_zNDZ_1(:,:,1)=Ex(:,:,NDZ-1);
            Ex_zNDZ(:,:,1)=Ex(:,:,NDZ);
        
        %%%% For Ey Boundaries Refreshing
            Ey_x2(:,:,5)=Ey_x2(:,:,4); 
            Ey_x3(:,:,5)=Ey_x3(:,:,4); 
            Ey_x4(:,:,5)=Ey_x4(:,:,4); 
            Ey_x5(:,:,5)=Ey_x5(:,:,4); 
        
        
            Ey_x2(:,:,4)=Ey_x2(:,:,3); 
            Ey_x3(:,:,4)=Ey_x3(:,:,3); 
            Ey_x4(:,:,4)=Ey_x4(:,:,3); 
            Ey_x5(:,:,4)=Ey_x5(:,:,3); 
        
            Ey_x2(:,:,3)=Ey_x2(:,:,2); 
            Ey_x3(:,:,3)=Ey_x3(:,:,2); 
            Ey_x4(:,:,3)=Ey_x4(:,:,2); 
            Ey_x5(:,:,3)=Ey_x5(:,:,2); 
  
            Ey_x2(:,:,2)=Ey_x2(:,:,1); 
            Ey_x3(:,:,2)=Ey_x3(:,:,1); 
            Ey_x4(:,:,2)=Ey_x4(:,:,1); 
            Ey_x5(:,:,2)=Ey_x5(:,:,1); 

        
            Ey_x2(:,:,1)=Ey(2,:,:); 
            Ey_x3(:,:,1)=Ey(3,:,:); 
            Ey_x4(:,:,1)=Ey(4,:,:); 
            Ey_x5(:,:,1)=Ey(5,:,:); 

            Ey_xNDK_3(:,:,5)=Ey_xNDK_3(:,:,4);
            Ey_xNDK_2(:,:,5)=Ey_xNDK_2(:,:,4);
            Ey_xNDK_1(:,:,5)=Ey_xNDK_1(:,:,4);
            Ey_xNDK(:,:,5)=Ey_xNDK(:,:,4);
            
            Ey_xNDK_3(:,:,4)=Ey_xNDK_3(:,:,3);
            Ey_xNDK_2(:,:,4)=Ey_xNDK_2(:,:,3);
            Ey_xNDK_1(:,:,4)=Ey_xNDK_1(:,:,3);
            Ey_xNDK(:,:,4)=Ey_xNDK(:,:,3);
            
            Ey_xNDK_3(:,:,3)=Ey_xNDK_3(:,:,2);
            Ey_xNDK_2(:,:,3)=Ey_xNDK_2(:,:,2);
            Ey_xNDK_1(:,:,3)=Ey_xNDK_1(:,:,2);
            Ey_xNDK(:,:,3)=Ey_xNDK(:,:,2);
            
            Ey_xNDK_3(:,:,2)=Ey_xNDK_3(:,:,1);
            Ey_xNDK_2(:,:,2)=Ey_xNDK_2(:,:,1);
            Ey_xNDK_1(:,:,2)=Ey_xNDK_1(:,:,1);
            Ey_xNDK(:,:,2)=Ey_xNDK(:,:,1);

            
            Ey_xNDK_3(:,:,1)=Ey(NDK-3,:,:);
            Ey_xNDK_2(:,:,1)=Ey(NDK-2,:,:);
            Ey_xNDK_1(:,:,1)=Ey(NDK-1,:,:);
            Ey_xNDK(:,:,1)=Ey(NDK,:,:);


            Ey_z2(:,:,5)=Ey_z2(:,:,4);
            Ey_z3(:,:,5)=Ey_z3(:,:,4);
            Ey_z4(:,:,5)=Ey_z4(:,:,4);
            Ey_z5(:,:,5)=Ey_z5(:,:,4);
            
            Ey_z2(:,:,4)=Ey_z2(:,:,3);
            Ey_z3(:,:,4)=Ey_z3(:,:,3);
            Ey_z4(:,:,4)=Ey_z4(:,:,3);
            Ey_z5(:,:,4)=Ey_z5(:,:,3);
            
            Ey_z2(:,:,3)=Ey_z2(:,:,2);
            Ey_z3(:,:,3)=Ey_z3(:,:,2);
            Ey_z4(:,:,3)=Ey_z4(:,:,2);
            Ey_z5(:,:,3)=Ey_z5(:,:,2);
            
            Ey_z2(:,:,2)=Ey_z2(:,:,1);
            Ey_z3(:,:,2)=Ey_z3(:,:,1);
            Ey_z4(:,:,2)=Ey_z4(:,:,1);
            Ey_z5(:,:,2)=Ey_z5(:,:,1);

            Ey_z2(:,:,1)=Ey(:,:,2);
            Ey_z3(:,:,1)=Ey(:,:,3);
            Ey_z4(:,:,1)=Ey(:,:,4);
            Ey_z5(:,:,1)=Ey(:,:,5);

            
            Ey_zNDZ_3(:,:,5)=Ey_zNDZ_3(:,:,4);
            Ey_zNDZ_2(:,:,5)=Ey_zNDZ_2(:,:,4);
            Ey_zNDZ_1(:,:,5)=Ey_zNDZ_1(:,:,4);
            Ey_zNDZ(:,:,5)=Ey_zNDZ(:,:,4);
            
            Ey_zNDZ_3(:,:,4)=Ey_zNDZ_3(:,:,3);
            Ey_zNDZ_2(:,:,4)=Ey_zNDZ_2(:,:,3);
            Ey_zNDZ_1(:,:,4)=Ey_zNDZ_1(:,:,3);
            Ey_zNDZ(:,:,4)=Ey_zNDZ(:,:,3);
                        
            Ey_zNDZ_3(:,:,3)=Ey_zNDZ_3(:,:,2);
            Ey_zNDZ_2(:,:,3)=Ey_zNDZ_2(:,:,2);
            Ey_zNDZ_1(:,:,3)=Ey_zNDZ_1(:,:,2);
            Ey_zNDZ(:,:,3)=Ey_zNDZ(:,:,2);
            
            Ey_zNDZ_3(:,:,2)=Ey_zNDZ_3(:,:,1);
            Ey_zNDZ_2(:,:,2)=Ey_zNDZ_2(:,:,1);
            Ey_zNDZ_1(:,:,2)=Ey_zNDZ_1(:,:,1);
            Ey_zNDZ(:,:,2)=Ey_zNDZ(:,:,1);

            Ey_zNDZ_3(:,:,1)=Ey(:,:,NDZ-3);
            Ey_zNDZ_2(:,:,1)=Ey(:,:,NDZ-2);
            Ey_zNDZ_1(:,:,1)=Ey(:,:,NDZ-1);
            Ey_zNDZ(:,:,1)=Ey(:,:,NDZ);
        
        %%%% For Ez Boundaries Refreshing
        Ez_x2(:,:,5)=Ez_x2(:,:,4);
        Ez_x3(:,:,5)=Ez_x3(:,:,4);
        Ez_x4(:,:,5)=Ez_x4(:,:,4);
        Ez_x5(:,:,5)=Ez_x5(:,:,4);
        
        Ez_x2(:,:,4)=Ez_x2(:,:,3);
        Ez_x3(:,:,4)=Ez_x3(:,:,3);
        Ez_x4(:,:,4)=Ez_x4(:,:,3);
        Ez_x5(:,:,4)=Ez_x5(:,:,3);
        
        Ez_x2(:,:,3)=Ez_x2(:,:,2);
        Ez_x3(:,:,3)=Ez_x3(:,:,2);
        Ez_x4(:,:,3)=Ez_x4(:,:,2);
        Ez_x5(:,:,3)=Ez_x5(:,:,2);

        Ez_x2(:,:,2)=Ez_x2(:,:,1);
        Ez_x3(:,:,2)=Ez_x3(:,:,1);
        Ez_x4(:,:,2)=Ez_x4(:,:,1);
        Ez_x5(:,:,2)=Ez_x5(:,:,1);
        
        Ez_x2(:,:,1)=Ez(2,:,:);
        Ez_x3(:,:,1)=Ez(3,:,:);
        Ez_x4(:,:,1)=Ez(4,:,:);
        Ez_x5(:,:,1)=Ez(5,:,:);
        

        Ez_xNDK_3(:,:,5)=Ez_xNDK_3(:,:,4);
        Ez_xNDK_2(:,:,5)=Ez_xNDK_2(:,:,4);
        Ez_xNDK_1(:,:,5)=Ez_xNDK_1(:,:,4);
        Ez_xNDK(:,:,5)=Ez_xNDK(:,:,4);
                
        Ez_xNDK_3(:,:,4)=Ez_xNDK_3(:,:,3);
        Ez_xNDK_2(:,:,4)=Ez_xNDK_2(:,:,3);
        Ez_xNDK_1(:,:,4)=Ez_xNDK_1(:,:,3);
        Ez_xNDK(:,:,4)=Ez_xNDK(:,:,3);
        
        Ez_xNDK_3(:,:,3)=Ez_xNDK_3(:,:,2);
        Ez_xNDK_2(:,:,3)=Ez_xNDK_2(:,:,2);
        Ez_xNDK_1(:,:,3)=Ez_xNDK_1(:,:,2);
        Ez_xNDK(:,:,3)=Ez_xNDK(:,:,2);
        
        Ez_xNDK_3(:,:,2)=Ez_xNDK_3(:,:,1);
        Ez_xNDK_2(:,:,2)=Ez_xNDK_2(:,:,1);
        Ez_xNDK_1(:,:,2)=Ez_xNDK_1(:,:,1);
        Ez_xNDK(:,:,2)=Ez_xNDK(:,:,1);

        Ez_xNDK_3(:,:,1)=Ez(NDK-3,:,:);
        Ez_xNDK_2(:,:,1)=Ez(NDK-2,:,:);
        Ez_xNDK_1(:,:,1)=Ez(NDK-1,:,:);
        Ez_xNDK(:,:,1)=Ez(NDK,:,:);
        

        
        Ez_y2(:,:,5)=Ez_y2(:,:,4); 
        Ez_y3(:,:,5)=Ez_y3(:,:,4); 
        Ez_y4(:,:,5)=Ez_y4(:,:,4);
        Ez_y5(:,:,5)=Ez_y5(:,:,4);
        
        Ez_y2(:,:,4)=Ez_y2(:,:,3); 
        Ez_y3(:,:,4)=Ez_y3(:,:,3); 
        Ez_y4(:,:,4)=Ez_y4(:,:,3);
        Ez_y5(:,:,4)=Ez_y5(:,:,3);
        
        Ez_y2(:,:,3)=Ez_y2(:,:,2); 
        Ez_y3(:,:,3)=Ez_y3(:,:,2); 
        Ez_y4(:,:,3)=Ez_y4(:,:,2); 
        Ez_y5(:,:,3)=Ez_y5(:,:,2); 
        
        Ez_y2(:,:,2)=Ez_y2(:,:,1); 
        Ez_y3(:,:,2)=Ez_y3(:,:,1); 
        Ez_y4(:,:,2)=Ez_y4(:,:,1); 
        Ez_y5(:,:,2)=Ez_y5(:,:,1); 

        
        Ez_y2(:,:,1)=Ez(:,2,:); 
        Ez_y3(:,:,1)=Ez(:,3,:); 
        Ez_y4(:,:,1)=Ez(:,4,:); 
        Ez_y5(:,:,1)=Ez(:,5,:); 

        
        Ez_yNDL_3(:,:,5)=Ez_yNDL_3(:,:,4);      
        Ez_yNDL_2(:,:,5)=Ez_yNDL_2(:,:,4); 
        Ez_yNDL_1(:,:,5)=Ez_yNDL_1(:,:,4); 
        Ez_yNDL(:,:,5)=Ez_yNDL(:,:,4);
        
        Ez_yNDL_3(:,:,4)=Ez_yNDL_3(:,:,3);      
        Ez_yNDL_2(:,:,4)=Ez_yNDL_2(:,:,3); 
        Ez_yNDL_1(:,:,4)=Ez_yNDL_1(:,:,3); 
        Ez_yNDL(:,:,4)=Ez_yNDL(:,:,3); 
        
        Ez_yNDL_3(:,:,3)=Ez_yNDL_3(:,:,2); 
        Ez_yNDL_2(:,:,3)=Ez_yNDL_2(:,:,2); 
        Ez_yNDL_1(:,:,3)=Ez_yNDL_1(:,:,2); 
        Ez_yNDL(:,:,3)=Ez_yNDL(:,:,2); 
        
        Ez_yNDL_3(:,:,2)=Ez_yNDL_3(:,:,1); 
        Ez_yNDL_2(:,:,2)=Ez_yNDL_2(:,:,1); 
        Ez_yNDL_1(:,:,2)=Ez_yNDL_1(:,:,1); 
        Ez_yNDL(:,:,2)=Ez_yNDL(:,:,1); 
        
        Ez_yNDL_3(:,:,1)=Ez(:,NDL-3,:); 
        Ez_yNDL_2(:,:,1)=Ez(:,NDL-2,:); 
        Ez_yNDL_1(:,:,1)=Ez(:,NDL-1,:); 
        Ez_yNDL(:,:,1)=Ez(:,NDL,:); 
            

%%% ===== Refresh the LIAO 3rd order boundary conditions
if t>4    
%%%for x direction

Ex(:,1,:)=4*Ex_y2(:,:,2)-6*Ex_y3(:,:,3)+4*Ex_y4(:,:,4)-Ex_y5(:,:,5);
Ex(:,NDL+1,:)=4*Ex_yNDL(:,:,2)-6*Ex_yNDL_1(:,:,3)+4*Ex_yNDL_2(:,:,4)-Ex_yNDL_3(:,:,5);

Ex(:,:,1)=4*Ex_z2(:,:,2)-6*Ex_z3(:,:,3)+4*Ex_z4(:,:,4)-Ex_z5(:,:,5);
Ex(:,:,NDZ+1)=4*Ex_zNDZ(:,:,2)-6*Ex_zNDZ_1(:,:,3)+4*Ex_zNDZ_2(:,:,4)-Ex_zNDZ_3(:,:,5);

%%%for y direction   
   
Ey(1,:,:)=4*Ey_x2(:,:,2)-6*Ey_x3(:,:,3)+4*Ey_x4(:,:,4)-Ey_x5(:,:,5);
Ey(NDK+1,:,:)=4*Ey_xNDK(:,:,2)-6*Ey_xNDK_1(:,:,3)+4*Ey_xNDK_2(:,:,4)-Ey_xNDK_3(:,:,5);

Ey(:,:,1)=4*Ey_z2(:,:,2)-6*Ey_z3(:,:,3)+4*Ey_z4(:,:,4)-Ey_z5(:,:,5);
Ey(:,:,NDZ+1)=4*Ey_zNDZ(:,:,2)-6*Ey_zNDZ_1(:,:,3)+4*Ey_zNDZ_2(:,:,4)-Ey_zNDZ_3(:,:,5);

%     %%% for z direction
Ez(1,:,:)=4*Ez_x2(:,:,2)-6*Ez_x3(:,:,3)+4*Ez_x4(:,:,4)-Ez_x5(:,:,5);
Ez(NDK+1,:,:)=4*Ez_xNDK(:,:,2)-6*Ez_xNDK_1(:,:,3)+4*Ez_xNDK_2(:,:,4)-Ez_xNDK_3(:,:,5);

Ez(:,1,:)=4*Ez_y2(:,:,2)-6*Ez_y3(:,:,3)+4*Ez_y4(:,:,4)-Ez_y5(:,:,5);
Ez(:,NDL+1,:)=4*Ez_yNDL(:,:,2)-6*Ez_yNDL_1(:,:,3)+4*Ez_yNDL_2(:,:,4)-Ez_yNDL_3(:,:,5);

end


    
%%%% Wave Propagation visualization
Ep1(:,:)=Ex(1:NDK,round(jw1_ini/2),:);
Ep2(:,:)=Ey(:,round(jw1_ini/2),:);
Ep3(:,:)=Ez(:,round(jw1_ini/2),1:NDZ);  %%% at middle of conductor


    figure(2)
    subplot(2,3,1)
    plot(Ez(:,jw1_ini,17),'linewidth',2)
    hold on
    plot(Ez(:,35,6),'linewidth',2)
    plot(Ez(:,35,5),'r','linewidth',2)
    hold off
    title('Ez')
    
    subplot(2,3,2)
    surf(Ep1)
    view(-90,90)
    shading flat
    title('Ex','Fontsize',14)
    xlabel('X','Fontsize',14)
    ylabel('Z','Fontsize',14)
    set(gca,'Fontsize',14)
    
    
    subplot(2,3,3)
    plot(Ex(:,jw1_ini,17),'linewidth',2)
    hold on
    plot(Ex(:,35,15),'r--','linewidth',2)
    plot(Ex(:,35,6),'linewidth',2)
    plot(Ex(:,35,5),'r','linewidth',2)
    hold off
    title('Ex')
    
    subplot(2,3,4)
    surf(Ep2)
    view(-90,90)
    shading flat
    title('Ey','Fontsize',14)
    xlabel('X','Fontsize',14)
    ylabel('Z','Fontsize',14)
    set(gca,'Fontsize',14)
    
    subplot(2,3,5)
    surf(Ep3)
    shading flat
    view(0,90)
     view(-90,90)
    title('Ez','Fontsize',14)
    xlabel('X','Fontsize',14)
    zlabel('Z','Fontsize',14)
    
    subplot(2,3,6)
    plot(IalZ,'linewidth',2)
    title('Current Along the wire','Fontsize',14)
    xlabel('Length of Wire','Fontsize',14)
    ylabel('Current [A]','Fontsize',14)
    set(gca,'Fontsize',14)
    
end



figure(3)
plot(tsim*1e6,Vini/1e3,'k--','linewidth',2)
hold on
plot(tsim*1e6,Vsource_t/1e3,'b-.','linewidth',2)
plot(tsim*1e6,Vsz/1e3,'b-.','linewidth',2)
title({'Voltages at Line ends'},'Fontsize',14)
ylabel('Voltage [V]','Fontsize',14)
xlabel('time [\mus]','Fontsize',14)
legend({'Vini Channel','Vsource'},'Fontsize',14)
set(gca,'Fontsize',14)



figure(4)
figure
plot(tsim*1e6,Iini/1e3,'k--','linewidth',2)
hold on
plot(tsim*1e6,Imiddle/1e3,'k--','linewidth',2)
title('Current at the Base Channel','Fontsize',14)
ylabel('Current [kA]','Fontsize',14)
xlabel('time [\mus]','Fontsize',14)
set(gca,'Fontsize',14)


figure(5)
plot(tsim*1e6,Iini/1e3,'k--','linewidth',2)
hold on
plot(tsim*1e6,I1/1e3,'r--','linewidth',2)
plot(tsim*1e6,I2/1e3,'r--','linewidth',2)
plot(tsim*1e6,I3/1e3,'r--','linewidth',2)
plot(tsim*1e6,I4/1e3,'r--','linewidth',2)
plot(tsim*1e6,I5/1e3,'r--','linewidth',2)
title('Currents along the Channel','Fontsize',14)
ylabel('Current [kA]','Fontsize',14)
xlabel('time [\mus]','Fontsize',14)
set(gca,'Fontsize',14)


figure(6)
subplot(2,2,1)
plot(tsim*1e6,Exp1/1e3,'k--','linewidth',2)
hold on
plot(tsim*1e6,Exp2/1e3,'k--','linewidth',2)
plot(tsim*1e6,Exp3/1e3,'k--','linewidth',2)
plot(tsim*1e6,Exp4/1e3,'k--','linewidth',2)
title('Ex','Fontsize',14)
ylabel('Electric Field [kV/m]','Fontsize',14)
xlabel('time [\mus]','Fontsize',14)
set(gca,'Fontsize',14)

subplot(2,2,2)
plot(tsim*1e6,Eyp1/1e3,'k--','linewidth',2)
hold on
plot(tsim*1e6,Eyp2/1e3,'k--','linewidth',2)
plot(tsim*1e6,Eyp3/1e3,'k--','linewidth',2)
plot(tsim*1e6,Eyp4/1e3,'k--','linewidth',2)
title('Ey','Fontsize',14)
ylabel('Electric Field [kV/m]','Fontsize',14)
xlabel('time [\mus]','Fontsize',14)
set(gca,'Fontsize',14)

subplot(2,2,3:4)
plot(tsim*1e6,Ezp1/1e3,'k--','linewidth',2)
hold on
plot(tsim*1e6,Ezp2/1e3,'k--','linewidth',2)
plot(tsim*1e6,Ezp3/1e3,'k--','linewidth',2)
plot(tsim*1e6,Ezp4/1e3,'k--','linewidth',2)
title('Ez','Fontsize',14)
ylabel('Electric Field [kV/m]','Fontsize',14)
xlabel('time [\mus]','Fontsize',14)
set(gca,'Fontsize',14)





