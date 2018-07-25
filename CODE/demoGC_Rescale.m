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


%------------------ Simulation of grid field rescaling --------------------
%------------------ along a single direction by moving --------------------
%------------------ the grid-cell frequency points     --------------------
%------------------ closer or farther away from the    --------------------
%------------------ origin.                            --------------------



clear -all;

% define parameters
R          = 32;%128;                % <-- length of hexagonal axis.
N          = R*(3*R);                % <-- number of hexagonal elements on the hexagonal
                                     %     region of support.
theta      = 0;                      % <-- orientation of the grid in radians.
sigma      = 0;                      % <-- standard deviation for noise to add to each frequency
                                     %     point in hexagonal coordinates (drawn from Normal(0,sigma)).

Ox        = 0;                       % <-- cartesian horizontal offset to set the grid field phase.
Oy        = 0;                       % <-- cartesian vertical offset to set the grid field phase.

outFolder = 'out_rescale/';          % <-- output folder for the figures.

omega = 4;                           % <-- frequency in cycles across the 
                                     %     hexagonal region of support (should be integer).

fS = ones(6,2,2);                    % <-- set of scaling factors to scale the position of
                                     %     each frequency point.
% vertical doubling of frequency
%   (scaling each frequency point in hexagonal coordinates)
fS(1,1,1) =   1 ;  fS(1,2,1) =  0;
fS(2,1,1) =  0.5;  fS(2,2,1) =  2;
fS(3,1,1) =  1.5;  fS(3,2,1) =  2;
fS(4,1,1) =  -1 ;  fS(4,2,1) =  0;
fS(5,1,1) = -0.5;  fS(5,2,1) = -2;
fS(6,1,1) = -1.5;  fS(6,2,1) = -2;
% horizontal doubling of frequency
%   (scaling each frequency point in hexagonal coordinates)
fS(1,1,2) =   2 ;  fS(1,2,2) =  0;
fS(2,1,2) = -0.5;  fS(2,2,2) =  1;
fS(3,1,2) =  1.5;  fS(3,2,2) =  1;
fS(4,1,2) =  -2 ;  fS(4,2,2) =  0;
fS(5,1,2) =  0.5;  fS(5,2,2) = -1;
fS(6,1,2) = -1.5;  fS(6,2,2) = -1;



for fsIndex=1:size(fS,3)
    S          = N*ones(6,1);          % <-- scaling factors for each frequency point.

    % obtain grid field over the hexagonal region of support    
    [x,X] = gridCell(Ox,Oy,theta,omega,fS(:,:,fsIndex),S,R,sigma);

    % map to cartesian coordinates
    cartesianCoordsX = zeros(N,1);
    cartesianCoordsY = zeros(N,1);
    ind = 1;
    for r1=(-R):(0+R-1)
       for r2=(-R):(1*R-1)
          if( r1-r2>=-R && r1-r2<R) %  if indexes are inside hexagonal region
                                    %  in spatial domain
             cartesianCoordsX(ind)   = r1 - 0.5*r2;
             cartesianCoordsY(ind)   = 0.5*sqrt(3.0)*r2;
             ind = ind + 1;
          end
       end
    end

    minFreq = min(min(real(X)),min(imag(X)));
    maxFreq = max(max(real(X)),max(imag(X)));

    % plot the frequency domain in cartesian coordinates
    % -------------------------------------------------------------------------

    % plot the real frequency domain
    fig1 = figure('Visible','off');
    scatter(cartesianCoordsX,cartesianCoordsY,16.5,real(X),'o','filled');
    greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
    redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
    redGreenBlackMap = [redColorMap; greenColorMap; zeros(1, 256)]';
    colormap(redGreenBlackMap);
    caxis([-N, N]);
    axis image off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    figTitle = strcat(outFolder,'GC_Freq_R_');figTitle = strcat(figTitle, num2str(R));
    figTitle = strcat(figTitle,strcat('_REAL_fsID_',num2str(fsIndex)));
    figTitle = strcat(figTitle,'.png');
    saveas(fig1, figTitle,'png');

    % plot the imaginary frequency domain
    fig2 = figure('Visible','off');
    scatter(cartesianCoordsX,cartesianCoordsY,16.5,imag(X),'o','filled');
    colormap(redGreenBlackMap);
    caxis([-N, N]);
    axis image off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    figTitle = strcat(outFolder,'GC_Freq_R_');figTitle = strcat(figTitle, num2str(R));
    figTitle = strcat(figTitle,strcat('_IMAG_fsID_',num2str(fsIndex)));
    figTitle = strcat(figTitle,'.png');
    saveas(fig2, figTitle,'png');

    minX = min(real(x));
    maxX = max(real(x));

    % plot the spatial domain (grid field) in cartesian coordinates
    % -------------------------------------------------------------------------

    fig3_1 = figure('Visible','off');
    scatter(cartesianCoordsX,cartesianCoordsY,16.5,real(x),'o','filled');
    colormap(jet);
    caxis([minX, maxX]);%caxis([0, maxX]);
    axis image off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    figTitle = strcat(outFolder,'GC_HIDFT_R_');figTitle = strcat(figTitle, num2str(R));
    figTitle = strcat(figTitle,strcat('_REAL_fsID_',num2str(fsIndex)));
    figTitle = strcat(figTitle,'.png');
    saveas(fig3_1, figTitle,'png');

end

%--------------------------------------------------------------------------
