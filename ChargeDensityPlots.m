 function ChargeDensityPlots(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]

            % tarr is a time time array for the time you wish to plot
            if length(varargin) == 1
                sol = varargin{1};
                tarr = sol.t(end);
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 2
                sol = varargin{1};
                tarr = varargin{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 3
                sol = varargin{1};
                tarr = varargin{2};
                xrange = varargin{3};
                pointtype = 't';
            end

            % Call dfana to obtain band energies and QFLs
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);

            xnm = x*1e7;    % x in nm for plotting

            for i = 1:length(tarr)
                % find the tarr
                p1 = find(sol.t <= tarr(i));
                p1 = p1(end);

               
                FH1 = figure(5);
                % Final Charge Densities
                
                plot(xnm, n(p1, :), xnm, p(p1, :));
                hold on

                % Ionic space charge density
                
                plot(xnm, a(p1,:))%-par.dev.Nani);
                hold on
                
                plot(xnm, c(p1,:))%-par.dev.Ncat);
                
            end

            figure(5)
            grid on 
            ylabel('Carrier density [cm-3]')
            legend('\itn', '\itp','Anions','Cations')
            xlabel('Position [nm]')
            xlim([xrange(1), xrange(2)]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            

        end