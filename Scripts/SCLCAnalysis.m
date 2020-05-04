% SCLC analysis script which sets up a device in the usual manner, does a
% JV and analyses the JV given certain J-V power laws in regions of the
% curve. It then plots mobility for the relavent regions and conductivitiy
% for Ohmic regions. Charge density plots are made for the ohmic and M-G
% regions
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code

% Create parameters objects for Spiro/MAPI/TiO2 and PEDOT/MAPI/PCBM devices
prompt = ('Please select the input file you wish to analyse');
FileName = string(uigetfile('/Users/Will/Documents/MATLAB/GitHub/Driftfusion/Input_files/.csv'));

%define the parameter(s) you wish to vary 
% Phi_left = [-3.8 -3.9 -4.0];      
% Phi_right = Phi_left;
% 
% paramvar = Phi_left;

Ncat=  [1e0 1e15 1e18 1e21];
Nani = Ncat;

paramvar = Ncat;

% Ncatstr = cell(1,size(Ncat,2)) ;
% Ncatstr(:) = [NaN] ;
% NcatstrMG = cell(1,size(Ncat,2)) ;
% NcatstrMG(:) = [NaN] ;
% NcatstrOhm = cell(1,size(Ncat,2)) ;
% NcatstrOhm(:) = [NaN] ;

for k=1:size(paramvar,2)

par = pc(FileName);

%Define your parameter to vary. 

par.Ncat = paramvar(k);
par.Nani = paramvar(k);
% par.Phi_left = paramvar(k);             
% par.Phi_right = paramvar(k);              

Paramvarstr{k} = num2str(paramvar(k),'%10.3e\n');
% par.mue = mue(k);
% par.muh = mue(k);
% Ncatstr{k} = num2str(mue(k),'%10.3e\n');

% Find equilibrium solutions
soleq{1,k} = Paramvarstr{k};
soleq{2,k} = equilibrate(par);


% Perform a current voltage scan with frozen ions to 100V
JV{1,k} = Paramvarstr{k};
JV{2,k} = doJV(soleq{2,k}.ion, 1, 1000, 0, 0, 0, 50, 1);


%% ANALYSIS %%

   % Call dfana to obtain band energies and QFLs for this worksapce           
[Jf,jf,xf] = dfana.calcJ(JV{2,k}.dk.f);
[Jr,jr,xr] = dfana.calcJ(JV{2,k}.dk.r);
            
 Jtotf = Jf.tot(:,end);
 Jtotr = Jr.tot(:,end);
 Jtot = [Jtotf; Jtotr];  
 Vappf = (dfana.calcVapp(JV{2,k}.dk.f))';
 Vappr = (dfana.calcVapp(JV{2,k}.dk.r))';
 Vapp = [Vappf; Vappr];
 gradJV = gradient(log(Jtot))./gradient(log(Vapp));
 
 MaxGrad = max(gradJV);
%  MaxGrad = sprintf('Max Grad is %.4f',MaxGrad);
%     disp(MaxGrad)

MG= NaN(size(Vapp));               %Prealocate MG and Ohm yo have the same size as Vapp but to have NaN where the analysis condition is not satisfied
Ohm = NaN(size(Vapp));

count1 = 1;
count2 = 1;
 
 for i=1:size(gradJV)
     
     if gradJV(i)>1.6 && gradJV(i)<2.3
         
         MG(i,1) = ((Jtot(i,:))./((Vapp(i,:).^2))).*(8./9).*((par.d.^3)./(par.epp.*par.epp0.*par.e));                   
         %NcatstrMG{count1,1} = Ncatstr{k};
         count1= count1+1;
         
     elseif gradJV(i)<1.3 && gradJV(i)>0.7
         
         Ohm(i,1) = (Jtot(i,:))./(Vapp(i,:));
         %NcatstrOhm{count2,1} = Ncatstr{k};
         count2= count2+1;  
         
     end 
 end 
 
% plot the current voltage curve
dfplot.JV(JV{2,k},1)
dfplot.JVSCLC(JV{2,k})

% figure(4)
% hold on
% 
% figure(10)
% hold on
 
MGtot(:,k) = MG;
Ohmtot(:,k) = Ohm;  
 
figure(200)
plot(Vapp, MGtot(:,k),'.-','LineWidth',1,'MarkerSize',10)%,'DisplayName',NcatstrMG) 
set(gca,'FontSize',16)
xlabel('Vapp')
ylabel('MG Mobility [cm^2 V^-1 s^-1]')
title('Mott-Gurney Mobilitiy (set at 20)') 
legend('show')
grid on
hold on 
 
figure(201)
plot(Vapp,Ohmtot(:,k),'.-','LineWidth',1,'MarkerSize',10)%,'DisplayName', NcatstrOhm)
set(gca,'FontSize',16)
xlabel('Vapp')
ylabel('Ohmic Conductivity [J V^-1]')
title('Ohmic Conductivity')
legend('show')
grid on
hold on 

mu_MG(k) = nanmean(MG);
G(k)= nanmean(Ohm);



end 
