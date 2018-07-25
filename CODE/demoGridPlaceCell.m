%    This file is part of the Grid Cell Hexagonal Fourier Transform (GCHFT)
%    software for simulating grid fields and other spatial-domain signals
%    based on an Hexagonal Fourier Transform.
%
%    The GCHFT is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    The GCHFT software is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this software.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2018 Ulises Rodriguez Dominguez and Jeremy B. Caplan.
%    -------------------------------------------------------------------------
%    -------------------------------------------------------------------------


%------------------ Simulation of place fields with a  --------------------
%------------------ linear combination of sets of      --------------------
%------------------ phase-aligned grid fields at       --------------------
%------------------ different scales and weighted by   --------------------
%------------------ a Gaussian spatial window.         --------------------



clear -all;

% define parameters
R          = 32;%128                      % <-- length of hexagonal axes.
N          = 3*R*R;                       % <-- number of hexagonal elements on the hexagonal
                                          %     region of support.
theta      = 0;                           % <-- orientation of the grid in radians.
sigma      = 0;                           % <-- standard deviation for noise to add to each frequency
                                          %     point in hexagonal coordinates (drawn from Normal(0,sigma)).

Ox        = 0;%32;                        % <-- offset in horizontal axis in cartesian coordinates.
Oy        = 0;%32*sqrt(3)/2;              % <-- offset in vertical axis in cartesian coordinates.
fixedPhase  = 0;                          % <-- 1 for fixed phase (in which case the phase is
                                          %     derived from [Ox,Oy]) and 0 for random 
                                          %     local + non-local phase.

fPowSet = [1, 2, 3, 4, 5, 6, 7];          % <-- set of frequency powers to set the frequency
                                          %     of grid fields according to the rule
                                          %     sqrt(2)^ fPowSet(index) rounded to the nearest
                                          %     integer, which is in cycles across the 
                                          %     hexagonal region of support.

s       = 40;                             % <-- scale parameter for the size of the place field.
range  = 200;%256;                        % <-- range in cartesian coordinates in which to vary
                                          %     the local and non-local part of the phases for
                                          %     the offset (phase) clusters that will drive the
                                          %     place fields.
nplaceFields = 1;                         % <-- number of place fields to obtain (every place field
                                          %     takes some time to be computed on processors with around
                                          %     2 Ghz).
nphaseClusters = 4;                       % <-- number of offset (phase) clusters that will drive
                                          %     each place field (default to 4 based on the findings
                                          %     from Hayman and Jeffery (2008).


% definition of a red-green-black map to plot negative-positive-zero
% values respectively in the frequency domain.
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
redGreenBlackMap = [redColorMap; greenColorMap; zeros(1, 256)]';


cartesianCoordsX = zeros(N,1);            % memory for visualization in cartesian coordinates.
cartesianCoordsY = zeros(N,1);            % memory for visualization in cartesian coordinates.
wx               = zeros(N,1);            % memory to visualize the contribution of each grid
                                          % field for the place field formation after applying
                                          % a Gaussian weighting.
            
outFolder = 'out_place/';                 % <-- output folder.



for nplacecell=1:nplaceFields
    placeField       = zeros(N,1); % memory to visualize the current nplacecell placefield.
    
    % set the phases for the grid fields with a local part (r0_loc) and a
    % non-local part (r0)
    if( fixedPhase == 0 )
        r0     = range*rand(1,2) - 0.5*range;
        r0_loc = 0.1*range*rand(nphaseClusters,2) - 0.5*0.1*range;
        % print to file the matrix with the phases for reproducibility
        filePhases = fopen(strcat(outFolder,'filePhases.txt'),'a+');
        for j = 1:size(r0,1)
            fprintf(filePhases, '%2.2f\t', r0(j,:));
            fprintf(filePhases,'\n');
        end
        fprintf(filePhases,'-------\n');
        for j = 1:size(r0_loc,1)
            fprintf(filePhases, '%2.2f\t', r0_loc(j,:));
            fprintf(filePhases,'\n');
        end
        fprintf(filePhases,'------------------------------------\n');
        
        fclose(filePhases);
    end
    
    
    for nphasecluster=1:nphaseClusters
        
        % process the current frequency-varying grid-cell set for the current
        % offset (phase) cluster with ID nphasecluster.
        for fPow=1:size(fPowSet,2)
            omega = round( sqrt(2)^(fPowSet(fPow)) ); %<-- follow sqrt(2) increase rule for the frequency.
            
            S          = N*ones(6,1);         % <-- scaling factors for amplitude of each frequency point.
            
            if( fixedPhase == 1)
                ox = Ox; oy = Oy;
            else
                ox = r0(1,1)+r0_loc(nphasecluster,1);
                oy = r0(1,2)+r0_loc(nphasecluster,2);
            end
            
            % obtain grid field over the hexagonal region of support [spatial domain, freq. domain]
            [x,X] = gridCell(ox,oy,theta,omega,false,S,R,sigma);
            
            % map to cartesian coordinates
            cartesianCoordsX(:,1) = 0;
            cartesianCoordsY(:,1) = 0;
            wx(:,1)               = 0;
            ind = 1;
            for r1=(-R):(0+R-1)
               for r2=(-R):(1*R-1)
                  if( r1-r2>=-R && r1-r2<R) %  if indexes are inside hexagon in spatial domain
                     cartesianCoordsX(ind)   = r1 - 0.5*r2;
                     cartesianCoordsY(ind)   = 0.5*sqrt(3.0)*r2;
                     % weight by Gaussian here (at same phase (offset) as current grid module)
                     dist_2 = (r1 - ox)^2 + (r2 - oy)^2 - (r1 - ox)*(r2 - oy);
                     w = exp(-dist_2 / (2*s*s));
                     wx(ind,1) = w;
                     placeField(ind) = placeField(ind) + w*x(ind,1);
                     ind = ind + 1;
                  end
               end
            end

            minFreq = min(min(real(X)),min(imag(X)));
            maxFreq = max(max(real(X)),max(imag(X)));
            
            % plot the real part of the frequency domain of the grid cell in cartesian coordinates
            fig1 = figure('Visible','off');
            scatter(cartesianCoordsX,cartesianCoordsY,16.5,real(X),'o','filled');

            % Apply the colormap.
            colormap(redGreenBlackMap);caxis([-N, N]);%colormap(jet);caxis([minFreq, maxFreq]);
            axis image off;
            set(gca,'LooseInset',get(gca,'TightInset'));
            % separate each place field and each phase grid cluster ID in appropriate directories.
            figTitle = outFolder;
            figTitle = strcat( figTitle, strcat( strcat(num2str(nplacecell),'/'), strcat(num2str(nphasecluster),'/') ) );
            if( ~ exist(figTitle,'dir') )
               mkdir(figTitle);
            end
            figTitle = strcat(figTitle,'GC_FreqHIDTFT_R_');figTitle = strcat(figTitle, num2str(R));
            figTitle = strcat(figTitle,'_REAL_');
            figTitle = strcat(figTitle,'_fpow_');figTitle = strcat(figTitle,num2str(fPow));
            figTitle = strcat(figTitle,'.png');
            saveas(fig1, figTitle,'png');

            % plot the imaginary part of the frequency domain of the grid cell in cartesian coordinates
            fig2 = figure('Visible','off');
            scatter(cartesianCoordsX,cartesianCoordsY,16.5,imag(X),'o','filled');
            colormap(redGreenBlackMap);caxis([-N, N]);
            axis image off;
            set(gca,'LooseInset',get(gca,'TightInset'));
            % separate each place field and each phase grid cluster ID in appropriate directories.
            figTitle = outFolder;
            figTitle = strcat( figTitle, strcat( strcat(num2str(nplacecell),'/'), strcat(num2str(nphasecluster),'/') ) );
            if( ~ exist(figTitle,'dir') )
               mkdir(figTitle);
            end
            figTitle = strcat(figTitle,'GC_Freq_R_');figTitle = strcat(figTitle, num2str(R));
            figTitle = strcat(figTitle,'_IMAG_');
            figTitle = strcat(figTitle,'_fpow_');figTitle = strcat(figTitle,num2str(fPow));
            figTitle = strcat(figTitle,'.png');
            saveas(fig2, figTitle,'png');

            minX = min(real(x));maxX = max(real(x));

            % plot the spatial domain of the grid cell (the grid field) in cartesian coordinates
            fig3_1 = figure('Visible','off');
            scatter(cartesianCoordsX,cartesianCoordsY,16.5,real(x),'o','filled');
            colormap(jet);caxis([minX, maxX]);
            axis image off;
            set(gca,'LooseInset',get(gca,'TightInset'));
            % separate each place field and each phase grid cluster ID in appropriate directories.
            figTitle = outFolder;
            figTitle = strcat( figTitle, strcat( strcat(num2str(nplacecell),'/'), strcat(num2str(nphasecluster),'/') ) );
            if( ~ exist(figTitle,'dir') )
               mkdir(figTitle);
            end
            figTitle = strcat(figTitle,'GC_HIDTFT_R_');figTitle = strcat(figTitle, num2str(R));
            figTitle = strcat(figTitle,'_REAL_');figTitle = strcat(figTitle,'_fpow_');
            figTitle = strcat(figTitle,num2str(fPow));figTitle = strcat(figTitle,'.png');
            saveas(fig3_1, figTitle,'png');


            % plot the grid field Gaussian weight to visualize the
            % contribution of the current grid cell
            x_nonneg = x - minX;
            wx(:,1) = wx(:,1) .* x_nonneg(:,1);
            minX = min(real(wx));
            maxX = max(real(wx));
            fig3_2 = figure('Visible','off');
            scatter(cartesianCoordsX,cartesianCoordsY,16.5,real(wx),'o','filled');
            colormap(jet);caxis([minX, maxX]);
            axis image off;
            set(gca,'LooseInset',get(gca,'TightInset'));
            % separate each place field and each phase grid cluster ID in appropriate directories.
            figTitle = outFolder;
            figTitle = strcat( figTitle, strcat( strcat(num2str(nplacecell),'/'), strcat(num2str(nphasecluster),'/') ) );
            if( ~ exist(figTitle,'dir') )
               mkdir(figTitle);
            end
            figTitle = strcat(figTitle,'GC_W_HIDTFT_R_');figTitle = strcat(figTitle, num2str(R));
            figTitle = strcat(figTitle,'_REAL_');figTitle = strcat(figTitle,'_fpow_');
            figTitle = strcat(figTitle,num2str(fPow));figTitle = strcat(figTitle,'.png');
            saveas(fig3_2, figTitle,'png');
            
            % force some memory to be freed
            clear x_nonneg;
            clear x;
            clear X;
            clear wx;
            
        end
    end



   % plot the resulting place fields
   %----------------------------------------------------------------
   
   maxRePC_Cont    = max(max(real(placeField(:,:) )));
   minRePC_Cont    = min(min(real(placeField(:,:) )));
   PC_Control_Re = (real(placeField(:,:) )-minRePC_Cont)/(maxRePC_Cont-minRePC_Cont);
   
   % first plot the untresholded place field
   figPC1 = figure('Visible','off');
   scatter(cartesianCoordsX,cartesianCoordsY,16.5,real(PC_Control_Re),'o','filled');
   colormap(jet);
   axis image off;
   set(gca,'LooseInset',get(gca,'TightInset'));
    % separate each place field ID in appropriate directories.
    figTitle = outFolder;
    figTitle = strcat( figTitle, strcat(num2str(nplacecell),'/') );
    if( ~ exist(figTitle,'dir') )
       mkdir(figTitle);
    end
    figTitle = strcat(figTitle,'PC_HEXGC_R_');
    figTitle = strcat(figTitle, num2str(R));
    figTitle = strcat(figTitle,'_s_');figTitle = strcat(figTitle,num2str(s));
    figTitle = strcat(figTitle,'_ID_');figTitle = strcat(figTitle,num2str(nplacecell));
    figTitle = strcat(figTitle,'.png');
    saveas(figPC1, figTitle,'png');
   
   
    % then plot the tresholded place field
    PC_Control_Re( PC_Control_Re< (2/3) ) = 0;
    figPC1_2 = figure('Visible','off');
    scatter(cartesianCoordsX,cartesianCoordsY,16.5,real(PC_Control_Re),'o','filled');
    colormap(jet);
    axis image off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    % separate each place field ID in appropriate directories.
    figTitle = outFolder;
    figTitle = strcat( figTitle, strcat(num2str(nplacecell),'/') );
    if( ~ exist(figTitle,'dir') )
       mkdir(figTitle);
    end
    figTitle = strcat(figTitle,'PCThresh_HEXGC_R_');
    figTitle = strcat(figTitle, num2str(R));
    figTitle = strcat(figTitle,'_s_');figTitle = strcat(figTitle,num2str(s));
    figTitle = strcat(figTitle,'_ID_');figTitle = strcat(figTitle,num2str(nplacecell));
    figTitle = strcat(figTitle,'.png');
    saveas(figPC1_2, figTitle,'png');
    
end

