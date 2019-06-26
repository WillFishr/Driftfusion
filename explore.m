classdef explore
    
    methods (Static)
        
        function parexsol = explore2par(par_base, parnames, parvalues, JVpnts, JVstatsswitch, Vocstableswitch, SCLCswitch)
            % EXPLOREPAR is used to explore 2 different parameters using a parallel pool.
            % The code is likely to require modification for individual parameters
            % owing to possible dependencies.
            % PAR_BASE is the base parameter set
            % PARNAMES is a cell array with the parameter names in - check these
            % carefully to avoid heartache later
            % PARVALUES is matrix with the parameter value ranges e.g.
            % JVSTATSSWITCH = calculate JV statistic e.g. Jsc, Voc etc.
            % VOCSTABLESWITCH = Run find Voc for stablised solution
            % SCLCswitch
            
            tic
            disp('Starting parameter exploration');
            disp(['Parameter 1: ', parnames(1)]);
            disp(['Parameter 2: ', parnames(2)]);
            parval1 = cell2mat(parvalues(1));
            parval2 = cell2mat(parvalues(2));
            str1 = char(parnames(1));
            str2 = char(parnames(2));
            
            j = 1;
            
            parfor i = 1:length(parval1)
                
                par = par_base;
                par.Ana = 0;                    
                par = exploreparhelper(par, str1, parval1(i));
                                            %This region can be used to set
                                            %conditions for the parameter
                                            %exploration
                par.PhiC = par.PhiA;        %Set paramters to co-vary here
                
                if JVstatsswitch == 1
                    Voc_f = zeros(1, length(parval2));
                    Voc_r = zeros(1, length(parval2));
                    Jsc_f = zeros(1, length(parval2));
                    Jsc_r = zeros(1, length(parval2));
                    mpp_f = zeros(1, length(parval2));
                    mpp_r = zeros(1, length(parval2));
                    FF_f = zeros(1, length(parval2));
                    FF_r = zeros(1, length(parval2));
                end
                
                if Vocstableswitch
                    Voc_stable = zeros(1, length(parval2));
                    PLint = zeros(1, length(parval2));
                end
                
                if SCLCswitch
                    SCLCVapp = zeros(1, JVpnts);
                    SCLCJ = zeros(length(parval2), JVpnts);
                    SCLCgradJV = zeros(length(parval2), JVpnts);
                    SCLCmu_MG = zeros(length(parval2), JVpnts);
                    SCLCmax_gradJV = zeros(1, length(parval2));
                    SCLCmu_MG_maxgrad = zeros(1, length(parval2));
                   % SCLCmu_MG_GRAD2 = zeros(length(parval2), JVpnts);
                    SCLC_delta_mue = zeros(length(parval2), JVpnts); 
                    SCLC_delta_muh = zeros(length(parval2), JVpnts);
                    
                end
                
                for j = 1:length(parval2)
                    
                    runN = (i-1)*length(parval2) + j;
                    disp(['Run no. ', num2str(runN), ', taun = ', num2str(parval1(i)), ', E0 = ', num2str(parval2(j))]);
                    
                    par = exploreparhelper(par, str2, parval2(j));
                    
                    par.dev = pc.builddev(par);
                    
                    soleq = equilibrate(par);
                    % JV = doJV(soleq.i_sr, 50e-3, 100, 1, 1e-10, 0, 1.5, 2);
                    JV = doJV(soleq.eq, 50e-3, 100, 1, 0, 0, 10, 1);
                    SCLCstats = SCLCana(JV.dk.f, 2);
                    
                    if SCLCswitch
                        SCLCVapp(j,:) = SCLCstats.Vapp;
                        SCLCJ(j,:) = SCLCstats.J;
                        SCLCgradJV(j,:) = SCLCstats.gradJV';
                        SCLCmu_MG(j,:) = SCLCstats.mu_MG';
                        SCLCmax_gradJV(j) = SCLCstats.max_gradJV;
                        SCLCmu_MG_maxgrad(j) = SCLCstats.mu_MG_maxgrad;
                        %SCLCmu_MG_GRAD2(j) = SCLCstats.mu_MG_GRAD2;
                        SCLC_delta_mue(j,:) = SCLCstats.delta_mue';
                        SCLC_delta_muh(j,:) = SCLCstats.delta_muh';
                        
                    end
                    
                    if JVstatsswitch == 1
                        Voc_f(j) = JV.stats.Voc_f;
                        Voc_r(j) = JV.stats.Voc_r;
                        Jsc_f(j) = JV.stats.Jsc_f;
                        Jsc_r(j) = JV.stats.Jsc_r;
                        mpp_f(j) = JV.stats.mpp_f;
                        mpp_r(j) = JV.stats.mpp_r;
                        FF_f(j) = JV.stats.FF_f;
                        FF_r(j) = JV.stats.FF_r;
                    end
                    
                    if Vocstableswitch == 1
                        % For PL
                        [sol_Voc, Voc] = findVoc(soleq.i_sr, 1e-6, Voc_f(j), (Voc_f(j)+0.1))
                        Voc_stable(j) = Voc;
                        PLint(j) = sol_Voc.PLint(end);
                    end
                    
                end
                
                if JVstatsswitch == 1
                    A(i,:) = Voc_f;
                    B(i,:) = Voc_r;
                    C(i,:) = Jsc_f;
                    D(i,:) = Jsc_r;
                    E(i,:) = mpp_f;
                    F(i,:) = mpp_r;
                    G(i,:) = FF_f;
                    H(i,:) = FF_r;
                end
                
                if Vocstableswitch == 1
                    J(i,:) = Voc_stable;
                    K(i,:) = PLint;
                end
                
                if SCLCswitch
                    AA(:,:,i) = SCLCVapp;
                    BB(:,:,i) = SCLCJ;
                    CC(:,:,i) = SCLCgradJV;
                    DD(:,:,i) = SCLCmu_MG;
                    EE(i,:) = SCLCmax_gradJV;
                    FF(i,:) = SCLCmu_MG_maxgrad;
                    %GG(:,:,i) = SCLCmu_MG_GRAD2;
                    HH(:,:,i) = SCLC_delta_mue;
                    II(:,:,i) = SCLC_delta_muh;
                    
                end
                
            end
            
            if JVstatsswitch == 1
                parexsol.stats.Voc_f = A;
                parexsol.stats.Voc_r = B;
                parexsol.stats.Jsc_f = C;
                parexsol.stats.Jsc_r = D;
                parexsol.stats.mpp_f = E;
                parexsol.stats.mpp_r = F;
                parexsol.stats.FF_f = G;
                parexsol.stats.FF_r = H;
            end
            
            if Vocstableswitch == 1
                parexsol.stats.Voc_stable = J;
                parexsol.stats.PLint = K;
            end
            
            if SCLCswitch
                parexsol.SCLCstats.Vapp = AA;
                parexsol.SCLCstats.J = BB;
                parexsol.SCLCstats.gradJV = CC; 
                parexsol.SCLCstats.mu_MG = DD;
                parexsol.SCLCstats.max_gradJV = EE;
                parexsol.SCLCstats.mu_MG_maxgrad = FF;
               % parexsol.SCLCstats.mu_MG_GRAD2 = GG; 
                parexsol.SCLCstats.delta_mue = HH;
                parexsol.SCLCstats.delta_muh = II; 
            end
            
            parexsol.parnames = parnames;
            parexsol.parvalues = parvalues;
            parexsol.parval1 = parval1;
            parexsol.parval2 = parval2;
            parexsol.par_base = par_base;
            
            toc
        end
        
        function par = exploreparhelper(par, parname, parvalue)
            % takes parameter set and sets parname to parvalue- workaround for parallel
            % computing loops
            eval(['par.',parname,'=parvalue']);
        end
        
        function plotPL(parexsol)
            
            figure(3000)
            s1 = surf(parexsol.parval1, parexsol.parval2, parexsol.stats.PLint);
            ylabel(parexsol.parnames(1))
            xlabel(parexsol.parnames(2))
            set(s1,'YScale','log');
            zlabel('PL intensity [cm-2s-1]')
            shading interp
            colorbar
            
        end
        
        function plotVoc(parexsol)
            
            offset = parexsol.parval2-parexsol.par_base.IP(1);
            
            figure(3001)
            surf(offset, parexsol.parval1, parexsol.stats.Voc_f)
            s1 = gca;
            %ylabel('Ion density [cm-3]')
            ylabel('p-type SRH time constant [s]')
            xlabel('\Phi_A')%/p-type VB-Fermi level offset [eV]')
            zlabel('Voc F scan [V]')
            xlim([offset(1), offset(end)]);
            ylim([parexsol.parval1(1), parexsol.parval1(end)])
            set(s1,'YScale','log');
            shading interp
            colorbar
            %caxis([0.75, 0.95])
            
            figure(3002)
            surf(offset, parexsol.parval1, parexsol.stats.Voc_r)
            s1 = gca;
            %ylabel('Ion density [cm-3]')
            ylabel('p-type SRH time constant [s]')
            xlabel('p-type VB-Fermi level offset [eV]')
            zlabel('Voc R scan [V]')
            xlim([offset(1), offset(end)]);
            ylim([parexsol.parval1(1), parexsol.parval1(end)])
            set(s1,'YScale','log');
            shading interp
            colorbar
            %caxis([0.75, 0.95])
            
        end
        
        function plotVocstable(parexsol)
            
            offset = parexsol.parval2-parexsol.par_base.IP(1);
            
            figure(3001)
            surf(offset, parexsol.parval1, parexsol.stats.Voc_stable)
            s1 = gca;
            ylabel('Mobile ion density [cm-3]')
            %ylabel('p-type SRH time constant [s]')
            xlabel('p-type VB-Fermi level offset [eV]')
            zlabel('Voc F scan [V]')
            xlim([offset(1), offset(end)]);
            ylim([parexsol.parval1(1), parexsol.parval1(end)])
            set(s1,'YScale','log');
            shading interp
            colorbar
            caxis([0.75, 0.95])
            %caxis([1.05, 1.15])
        end
        
        
        function plotJscF(parexsol)
            
            offset = parexsol.parval2-parexsol.par_base.IP(1);
            
            figure(3001)
            surf(offset, parexsol.parval1, parexsol.stats.Jsc_f)
            s1 = gca;
            %ylabel('Mobile ion density [cm-3]')
            ylabel('p-type SRH time constant [s]')
            xlabel('p-type VB-Fermi level offset [eV]')
            zlabel('Jsc F scan [Acm-2]')
            xlim([offset(1), offset(end)]);
            ylim([parexsol.parval1(1), parexsol.parval1(end)])
            set(s1,'YScale','log');
            shading interp
            colorbar
            %caxis([0.75, 0.95])
            %caxis([1.05, 1.15])
        end
        
        function plotdeltamu(parexsol)
            
            figure(005)
            surf(parexsol.parval1, parexsol.parval2, squeeze(parexsol.SCLCstats.delta_mue(:,1,:)))
            s1 = gca;
            title('mu MG - mu input')
            xlim([parexsol.parval1(1,1) parexsol.parval1(1,end)])
            ylim([parexsol.parval2(1,1) parexsol.parval2(1,end)])
            xlabel(parexsol.parnames{1,1})
            ylabel(parexsol.parnames{1,2})
            zlabel('Calculated mobility [cm2V-1s-1]')
%             set(s1,'XScale','log');
            set(s1,'YScale','log');
%             set(s1,'ZScale','log');
            shading interp
            cb = colorbar();
%             cb.Ruler.Scale = 'log';
            cb.Ruler.MinorTick = 'on';
%            caxis([1e-10, 1])            
%            xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
%             
%             figure(3004)
%             semilogy(parexsol.parval2-parexsol.par_base.IPi, parexsol.SCLCstats.mu_MG_maxgrad(1,:), '-o');
%             xlabel('Fermi-IP offset [eV]')
%             ylabel('Calculated mobility [cm2V-1s-1]')
%             ylim([0,2.2])
%             xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
        end
        
        
        function plotmu_MG(parexsol)
            
            figure(007)
            surf(parexsol.parval1, parexsol.parval2, squeeze(parexsol.SCLCstats.mu_MG(:,1,:)))
            s1 = gca;
            title('mu MG')
            xlim([parexsol.parval1(1,1) parexsol.parval1(1,end)])
            ylim([parexsol.parval2(1,1) parexsol.parval2(1,end)])
            xlabel(parexsol.parnames{1,1})
            ylabel(parexsol.parnames{1,2})
            zlabel('Calculated mobility [cm2V-1s-1]')
%             set(s1,'XScale','log');
            set(s1,'YScale','log');
%             set(s1,'ZScale','log');
            shading interp
            cb = colorbar();
%             cb.Ruler.Scale = 'log';
            cb.Ruler.MinorTick = 'on';
%            caxis([1e-10, 1])            
%            xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
%             
%             figure(3004)
%             semilogy(parexsol.parval2-parexsol.par_base.IPi, parexsol.SCLCstats.mu_MG_maxgrad(1,:), '-o');
%             xlabel('Fermi-IP offset [eV]')
%             ylabel('Calculated mobility [cm2V-1s-1]')
%             ylim([0,2.2])
%             xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
        end
        
         function plotmu_MG_MaxGrad(parexsol)
            
            figure(006)
            surf(parexsol.parval1, parexsol.parval2, parexsol.SCLCstats.mu_MG_maxgrad)
            s1 = gca;
            title('mu MG at MaxGrad')
            xlim([parexsol.parval1(1,1) parexsol.parval1(1,end)])
            ylim([parexsol.parval2(1,1) parexsol.parval2(1,end)])
            xlabel(parexsol.parnames{1,1})
            ylabel(parexsol.parnames{1,2})
            zlabel('Calculated mobility [cm2V-1s-1]')
%             set(s1,'XScale','log');
            set(s1,'YScale','log');
%             set(s1,'ZScale','log');
            shading interp
            cb = colorbar();
%             cb.Ruler.Scale = 'log';
            cb.Ruler.MinorTick = 'on';
%            caxis([1e-10, 1])            
%            xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
%             
%             figure(3004)
%             semilogy(parexsol.parval2-parexsol.par_base.IPi, parexsol.SCLCstats.mu_MG_maxgrad(1,:), '-o');
%             xlabel('Fermi-IP offset [eV]')
%             ylabel('Calculated mobility [cm2V-1s-1]')
%             ylim([0,2.2])
%             xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
        end
        
        function plot_MaxGrad(parexsol)
            
            figure(008)
            surf(parexsol.parval1, parexsol.parval2,  parexsol.SCLCstats.max_gradJV);
            s1 = gca;
            title('MaxGrad')
            xlim([parexsol.parval1(1,1) parexsol.parval1(1,end)])
            ylim([parexsol.parval2(1,1) parexsol.parval2(1,end)])
            xlabel(parexsol.parnames{1,1})
            ylabel(parexsol.parnames{1,2})
            zlabel('Max Gradient of log(J)/log(V)')
%             set(s1,'XScale','log');
            set(s1,'YScale','log');
%             set(s1,'ZScale','log');
            shading interp
            cb = colorbar();
%             cb.Ruler.Scale = 'log';
            cb.Ruler.MinorTick = 'on';
%            caxis([1e-10, 1])            
%            xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
%             
%             figure(3004)
%             semilogy(parexsol.parval2-parexsol.par_base.IPi, parexsol.SCLCstats.mu_MG_maxgrad(1,:), '-o');
%             xlabel('Fermi-IP offset [eV]')
%             ylabel('Calculated mobility [cm2V-1s-1]')
%             ylim([0,2.2])
%             xlim([parexsol.parval2(1)-parexsol.par_base.IPi, parexsol.parval2(end)-parexsol.par_base.IPi])
        end
        
    end
    
end
