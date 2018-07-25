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


%------------------ Simulation of a non-regular grid  --------------------
%------------------ field by adding Gaussian noise to --------------------
%------------------ the position of the frequency     --------------------
%------------------ coordinates. A regular grid field --------------------
%------------------ can also be obtained by setting   --------------------
%------------------ sigma = 0.                        --------------------



clear -all;

% define parameters
R          = 32;%128;                % <-- length of hexagonal axis.
N          = R*(3*R);                % <-- number of hexagonal elements on the hexagonal
                                     %     region of support.
theta      = 0;                      % <-- orientation of the grid in radians.
sigma      = 0.5;                    % <-- standard deviation for noise to add to each frequency
                                     %     point in hexagonal coordinates (drawn from Normal(0,sigma)).

Ox        = 0;                       % <-- cartesian horizontal offset to set the grid field phase.
Oy        = 0;                       % <-- cartesian vertical offset to set the grid field phase.

outFolder = 'out_noise/';          % <-- output folder for the figures.
%outFolder = 'out_regular/';        % <-- uncomment this for a regular grid field
                                   %     with no noise and set sigma = 0.

freqPow = [1,4,7];                 % <-- set of frequency powers to set the frequency
                                   %     of grid fields according to the rule
                                   %     sqrt(2)^ freqPow(index) rounded to the nearest
                                   %     integer, which is in cycles across the 
                                   %     hexagonal region of support.



for index=1:size(freqPow,2)
    S          = N*ones(6,1);          % <-- scaling factors for each frequency point.
    omega      = round(sqrt(2)^freqPow(1,index)); % <-- frequency in cycles across the
                                                  %     hexagonal region of support
                                                  %     (should be integer).

    % obtain grid field over the hexagonal region of support    
    [x,X] = gridCell(Ox,Oy,theta,omega,false,S,R,sigma);

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
    figTitle = strcat(figTitle,strcat('_REAL_omega_',num2str(omega)));
    figTitle = strcat(figTitle,strcat('_sigma_',num2str(sigma)));
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
    figTitle = strcat(figTitle,strcat('_IMAG_omega_',num2str(omega)));
    figTitle = strcat(figTitle,strcat('_sigma_',num2str(sigma)));
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
    figTitle = strcat(figTitle,strcat('_REAL_omega_',num2str(omega)));
    figTitle = strcat(figTitle,strcat('_sigma_',num2str(sigma)));
    figTitle = strcat(figTitle,'.png');
    saveas(fig3_1, figTitle,'png');

end

%--------------------------------------------------------------------------
