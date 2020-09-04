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
% Phi_left = [-5.4 -5.375 -5.35 -5.325 -5.3 -5.25 -5.2 -5.15];      
% Phi_right = Phi_left;
% 
% paramvar = Phi_left;

% d = [10e-4 5e-4];
% paramvar = d;
% density = 200/1e-5;

scanrate = [0.01 0.1 1 10 100];
paramvar = scanrate;


% V_prebias = [0.1];% 1];
% V_prebias = [-6 -5 -4 -3 -2 -1 -0.8 -0.6 -0.4 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 -0.025 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.6 0.8 1 2 3 4 5 6];
%V_prebias = [-0.5 -0.1 0 0.1 0.5];
% paramvar= V_prebias;

% Ncat = [1e13,1e15,1e16,1e17,1e18,1e19]; 
% Ncat=  [1e5,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20];
% Nani = Ncat;
% paramvar = Ncat;

% Ncatstr = cell(1,size(Ncat,2)) ;
% Ncatstr(:) = [NaN] ;
% NcatstrMG = cell(1,size(Ncat,2)) ;
% NcatstrMG(:) = [NaN] ;
% NcatstrOhm = cell(1,size(Ncat,2)) ;
% NcatstrOhm(:) = [NaN] ;

for k=1:size(paramvar,2)

    identifier = 0;
    identifier2 = 0;
    
par = pc(FileName);

%Define your parameter to vary. 

% par.Ncat = paramvar(k);
% par.Nani = paramvar(k);
% par.Phi_left = paramvar(k);             
% par.Phi_right = paramvar(k);
% par.d = paramvar(k);
% layerpoints = [20000 10000];
% par.layer_points = layerpoints(k);
%par.layer_points = round(min([1000 density*paramvar(k)]));
par.SRHset = 0;
par = refresh_device(par);


Paramvarstr{k} = num2str(paramvar(k),'%10.3e\n');
% par.mue = mue(k);
% par.muh = mue(k);
% Ncatstr{k} = num2str(mue(k),'%10.3e\n');

% Find equilibrium solutions
soleq{1,k} = Paramvarstr{k};
soleq{2,k} = equilibrate(par);

% sol_prebias{1,k} = Paramvarstr{k};

Conductance_eq(k) = 1./(trapz(par.xx,1./(par.e.*(soleq{2,k}.ion.u(end,:,2).*par.mue + soleq{2,k}.ion.u(end,:,3).*par.muh))));

JV{1,k} = Paramvarstr{k};
JV{2,k} = doJV(soleq{2,k}.ion, paramvar(k), 500, 0, 1, 0, 30, 1);

% if paramvar(k) == 0 
%     figure(2)
%     hold on 
%     JV{2,k} = doJV(soleq{2,k}.ion, 1, 2500, 0, 0, 0, 30, 1);
%     sol_prebias{2,k} = soleq{2,k}.ion;
%     Conductance_prebias(k) = Conductance_eq(k);
%     Conductance_prebias0V(k) = Conductance_eq(k);
% else  
%     sol_prebias{2,k} = jumptoV(soleq{2,k}.ion, paramvar(k), 500, 1, 0, 1, 0);
%     Conductance_prebias(k) =  (1./(trapz(par.xx,1./(par.e.*(sol_prebias{2,k}.u(end,:,2).*par.mue + sol_prebias{2,k}.u(end,:,3).*par.muh)))));
%     figure(2)
%     hold on 
%     scantoJV{1,k} = doJV(sol_prebias{2,k},0.1,500,0,0,paramvar(k),0,1);
%     % prebias_0V{1,k} = stabilize(scantoJV{1,k}.dk.f);
%     Conductance_prebias0V(k) =  (1./(trapz(par.xx,1./(par.e.*(scantoJV{1,k}.dk.f.u(end,:,2).*par.mue + scantoJV{1,k}.dk.f.u(end,:,3).*par.muh)))));
%     % Perform a current voltage scan with frozen ions to 100V
%     JV{2,k} = doJV(scantoJV{1,k}.dk.f, 1, 2500, 0, 0, 0, 30, 1);
% end

%% ANALYSIS %%
% Call dfana to obtain band energies and QFLs for this worksapce           
[Jf,jf,xf] = dfana.calcJ(JV{2,k}.dk.f);
[Jr,jr,xr] = dfana.calcJ(JV{2,k}.dk.r);
            
 Jtotf = Jf.tot(:,end);
 Jtotr = Jr.tot(:,end);
 Jtot = Jtotf;
%  Jtot = [Jtotf; Jtotr]; 
 lgJtot = log(abs(Jtot));
%  lgJtot([1 end],:) = [];
 Vappf = (dfana.calcVapp(JV{2,k}.dk.f))';
 Vappr = (dfana.calcVapp(JV{2,k}.dk.r))';
 Vapp = Vappf; 
 %  Vapp = [Vappf; Vappr];
 lgVapp = log(abs(Vapp));
%  lgVapp([1 end],:) = [];
 gradJV = [];
 gradJV = gradient(log(Jtot))./gradient(log(Vapp));
 
 MaxGrad = max(gradJV);
 MaxGrad_loc = find(gradJV==MaxGrad);
%  MaxGrad = sprintf('Max Grad is %.4f',MaxGrad);
%     disp(MaxGrad)

MG= NaN(size(Vapp));               %Prealocate MG and Ohm yo have the same size as Vapp but to have NaN where the analysis condition is not satisfied
Ohm = NaN(size(Vapp));
Conductance_V = NaN(size(Vapp));
% 
count1 = 1;
count2 = 1;
index = 1;
%  
 for i=1:size(gradJV)
     
     if gradJV(i)>1.65 && gradJV(i)<2.35
         
         MG(i,1) = ((Jtot(i,:))./((Vapp(i,:).^2))).*(8./9).*(((par.d).^3)./(par.epp.*par.epp0.*par.e));                   
         %NcatstrMG{count1,1} = Ncatstr{k};
%          count1= count1+1;
     end    
     if gradJV(i)< 2.2 && gradJV(i)>0.7
         
         Ohm(i,1) = ((Jtot(i,:))./(Vapp(i,:)));
         Conductance_V(i,1) = 1./(trapz(par.xx,1./(par.e.*(JV{2,k}.dk.f.u(i,:,2).*par.mue + JV{2,k}.dk.f.u(i,:,3).*par.muh))));       
         %NcatstrOhm{count2,1} = Ncatstr{k};
%          count2= count2+1;  
         
     end 
 end

MGtot(:,k) = MG;
Ohmtot(:,k) = Ohm;  
ConductanceVtot(:,k) = Conductance_V;
 
MG_muexp1(k)=((Jtot(MaxGrad_loc,:))./((Vapp(MaxGrad_loc,:).^2))).*(8./9).*(((par.d).^3)./(par.epp.*par.epp0.*par.e));
 
for t=1:(size(gradJV))
     
     if gradJV(t)>1.6 && gradJV(t)<2.4 && gradJV(t)>gradJV(t-1)
        
        identifier = 1;
         
        MGregionlgJ(count1) = lgJtot(t);
        MGregionlgV(count1) = lgVapp(t);
        count1 = count1 + 1; 
        %calcualte the mobilities from the gradients of the grad and create
        %strings of the equations to add as labels to the graphs .
     end 
     
     if gradJV(t)>2.4 && gradJV(t)>gradJV(t-1)
         
         TrapregionV(count2) = Vapp(t);
         count2 = count2 +1;
         identifier2 = 1;
         %Calculate the trap density in regions where the gradient exceeds
         %2.5, identify the first point in this region and use this. 
     end    
     
end

        if identifier == 1

            logfit = fit(MGregionlgV',MGregionlgJ','poly1');
            fitcoefficients{1,k} = logfit.p1;
            fitcoefficients{2,k} = logfit.p2; 
            Mulogfit(k) = (8.*(par.d).^3).*exp(fitcoefficients{2,k})./(9.*par.epp.*par.epp0.*par.e);
            mstring = num2str(fitcoefficients{1,k});
            intstring = num2str(fitcoefficients{2,k});
            eqlabel{k} = (['log(J) = ' mstring 'log(V) + ' intstring]);
            
        elseif identifier == 0 
            
            fitcoefficients{1,k} = NaN;
            fitcoefficients{2,k} = NaN; 
            Mulogfit(k) = NaN;
            mstring = num2str(NaN);
            intstring = num2str(NaN);
            eqlabel{k} = (['log(J) = ' mstring 'log(V) + ' intstring]);

        end
       
        if identifier2 == 1
            Vons(k) = TrapregionV(1);
            TrapDensity(k) = (2.*par.epp.*par.epp0.*Vons(k))./((par.d).^2) ;
            index = index + 1;
        elseif identifier2 == 0 
            Vons(k) = NaN;
            TrapDensity(k) = NaN;
        end 
        
       

% plot the current voltage curve
dfplot.JV(JV{2,k},1)
hold on 
dfplot.JVSCLC(JV{2,k})

figure(200)
plot(Vapp, MGtot(:,k),'.-','LineWidth',1,'MarkerSize',10)%,'DisplayName',NcatstrMG) 
set(gca,'FontSize',16)
xlabel('Vapp')
ylabel('MG Mobility [cm^2 V^-1 s^-1]')
title('Mott-Gurney Mobility (set at 20)') 
legend('show')
grid on
hold on 
 
figure(201)
plot(Vapp,Ohmtot(:,k),'.-','LineWidth',1,'MarkerSize',10)%,'DisplayName', NcatstrOhm)
hold on
set(gca,'FontSize',16)
plot(Vapp,ConductanceVtot(:,k),'--','LineWidth',1,'MarkerSize',10)
xlabel('Vapp')
ylabel('Ohmic Conductance [S cm^-2]')
title('Ohmic Conductance')
legend('show')
grid on
hold on 

y=fitcoefficients{1,k}.*lgVapp + fitcoefficients{2,k};

figure(210)
plot(lgVapp, lgJtot,'.','LineWidth',1,'MarkerSize',10);
hold on 
plot(lgVapp,y,'-','LineWidth',1,'MarkerSize',10);
legend('show')
xlabel('log(Vapp)')
ylabel('log(J)')
title('Typical analysis of experimental SCLC data')
grid on
hold on 

VonsLabel = Vons(imag(Vons)==0);

figure(998)
loglog(Vapp,Jtot,'-','LineWidth',1)
%xline(VonsLabel(index),'--',{'Onset Voltage'})
xlabel('log(Vapp)')
ylabel('log(J)')
title('logJ vs logV')
hold on

mu_MG(k) = nanmean(MG);
G(k)= nanmean(Ohm);


end

