clc;
clear all;
z=3;
rho=1.21;
%raggio palare massimo
R=1.43;
%raggio del mozzo: attenzione, non deve essere troppo piccolo!
Rm=0.20;
%numero di elementi discreti
N=14;
%coefficienti del profilo scelto trovati con test xfoil
Cl=0.95;
Cd=0.02;
%velocità del vento indisturbata
v_inf=12;
alfa=6; %angolo di attacco, deg
Omega=342; %velocità di rotazione, rpm
MAXITER=100;    %numero massimo di iterazioni
%unità di misura adeguate
alfa=alfa*(pi/180); %rad
Omega=2*pi*Omega/60; %rad/s

r(1)=Rm;
for i=2:1:(N+1)
    r(i)=r(i-1)+R/N;
end
for i=1:1:N+1
    lambdar(i)=Omega*r(i)/v_inf;
    a(i)=0.330;
    aprimo(i)=(1-3*a(i))/(4*a(i)-1);
    iter(i)=0;
    
    while iter(i)<MAXITER
        theta(i)=atan((1-a(i))/((1+aprimo(i))*lambdar(i)));
        sigma(i)=4*(1-cos(theta(i)))/Cl;
        c(i)=2*pi*r(i)*sigma(i)/z;
        
        Cx(i)=Cl*cos(theta(i))+Cd*sin(theta(i));
        Cu(i)=Cl*sin(theta(i))-Cd*cos(theta(i));
    
        anew(i)=((4*sin(theta(i))*sin(theta(i)))/(sigma(i).*Cx(i))+1)^(-1);
        aprimonew(i)=((4*sin(theta(i))*cos(theta(i))/(sigma(i).*Cu(i)))-1)^(-1);
    
        err=max([abs(anew(i)-a(i)), abs(aprimonew(i)-aprimo(i))]);
    
            if err>0.0001        
                a(i)=anew(i);
                aprimo(i)=aprimonew(i);
                else break
            end
    iter(i)=iter(i)+1;
    end
end
% disp([r', lambdar', a', aprimo', rad2deg((theta))', sigma', anew', aprimonew', iter', c']);

figure(1)
hold on
grid minor
xlabel('r [m]');    
ylabel('[deg]');
scatter(r, rad2deg(theta),'x','r');
thetadeg=rad2deg((2/3)*atan(1./lambdar));
phi=rad2deg(theta-alfa);
plot(r, thetadeg,'b');
plot(r, phi,'g');
legend('theta calcolo iterativo', 'theta con formula pratica','phi');
hold off

figure (2)
hold on
grid minor
xlabel('r [m]');
ylabel('[-]');
plot(r, a);
plot(r, aprimo,'.-');
plot(r, sigma,'.');
plot(r, c, 'g');
legend('a','aprimo','sigma locale','corda');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcolo dei triangoli di velocità per ogni elemento di pala
u=(1+aprimo).*Omega.*r;     % velocità tangenziale del vento
v=(1-a).*v_inf.*ones(1,length(u));  %velocità assoluta del vento
w=sqrt(u.^2+v.^2); % velocità relativa del vento

figure (3)
hold on
grid minor
xlabel('r [m]');
ylabel('[m/s]');
plot(r, u);
plot(r, v);
plot(r, w);
legend('u=vel. tang,','v=vel. assol.','w=vel. rel.');
hold off

%calcolo del fattore di Prandtl - perdite # finito di pale
f=(z*(R+Rm-r))./(2*r.*sin(theta));
F=(2/pi)*acos((exp(1)).^(-f));
F=round(F, 3);

% calcolo delle forze per ogni elemento discreto della pala
dFu=0.5*rho*F.*Cu.*c.*w.^2; 
dFx=0.5*rho*F.*Cx.*c.*w.^2;
dM=r.*dFu;

figure (4)
hold on
grid minor
xlabel('r [m]');
ylabel('[-]');
title('forze adimensionalizzate');
plot(r, dFu./max(dFu));
plot(r, dFx./max(dFx));
plot(r, dM./max(dM));
plot(r, F,'.');
legend('delta Fu,','delta Fx','delta M','F=fattore di Prandtl');
hold off

%quantità per 1 pala
Fu=trapz(r, dFu);   %forza tangenziale complessiva
Fx=trapz(r, dFx);   %forza assiale complessiva
M=trapz(r, dM);     %coppia complessiva
P=z*Omega*M;        %potenza totale (z pale) della turbina
Cp=P/(0.5*rho*v_inf^3*pi*R^2);  %cifra di potenza teorico globale

format shortG
disp('       Fx [N]     Ftang [N]      M [N*m]      P [W]       Cp[-]');
disp([Fx, Fu, M, P, round(Cp, 3)]);

disp(['          r        lambdar         a         aprimo        theta        sigma          c']);
disp([r', round(lambdar, 3)', round(a, 3)', round(aprimo, 3)',...
    round(rad2deg((theta)), 1)', round(sigma, 2)', round(c, 3)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apre il file con le coordinate adimensionalizzate 
fid1=fopen('S823_coordinates.txt','r');
if fid1==-1
    disp('File open failed');
else
    while feof(fid1)==0
       G = cell2mat(textscan(fid1, '%f %f', 'headerlines', 1));
       x_adim=G(:,1);       
       y_adim=G(:,2); 
    end
    closeresult=fclose(fid1);
    if closeresult==0
        disp('File closed successful');
    else
        disp('File closed failed');
    end
end


figure(5)
axis auto
for j=1:1:length(c)
    %calcola le coordinate in [m] a partire dalle coordinate adimensionali
   x=c(j).*x_adim;
   y=c(j).*y_adim;
   %salva le coordinate dei punti che formano il profilo in una struttura
   %sistema di riferimento: la corda di ogni profilo inizia (x(1); y(1))
   %nell'origine (0:0)
   C_0(j).id=j;   C_0(j).xData=x;   C_0(j).yData=y; C_0(j).zData=ones(length(x), 1)*r(j);
   
   %calcola centro di massa del profilo alare
   poligono(j)=polyshape(x, y);
   [x_c(j), y_c(j)]=centroid(poligono(j));
   %trasla ogni profilo in modo che il suo centro di massa cada
   %nell'origine degli assi cartesiani 
   for k=1:length(x)
      x(k)=x(k)-x_c(j);
      y(k)=y(k)-y_c(j);   
   end
   %disegna i profili in 3D
%    profilo=plot3(x, y, ones(length(x), 1)*r(j), 'k');
   profilo=fill3(x, y, ones(length(x), 1)*r(j), ones(length(x), 1)*r(j));
   %ruota i profili di theta [deg] attorno al loro centro di massa (portato
   %in (0;0);
   rotate(profilo, [0; 0; 1], (-1).*(phi(j)), [0; 0; 0]);

   C_cdm(j).id=j;   
   C_cdm(j).xData=x;  
   C_cdm(j).yData=y; 
   C_cdm(j).zData=C_0(j).zData;
   hold on
    fid(j) = fopen(sprintf('profilo_pala_%d.txt',j), 'w');
   fprintf(fid(j), '%6.2f %6.2f %6.2f\r\n', [1000*C_cdm(j).xData, 1000*C_cdm(j).yData, zeros(length(x), 1) ]');
   fclose(fid(j))
end
hold off


 

