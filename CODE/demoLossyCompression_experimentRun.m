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


%------------------ Simulation of a population of grid --------------------
%------------------ cells acting as filters at certain --------------------
%------------------ frequency bands and orientations.  --------------------
%------------------ The input 2D signal is either a    --------------------
%------------------ natural image or a uniform white   --------------------
%------------------ noise image.                       --------------------
%------------------ Natural images can be obtained from--------------------
%------------------ http://cvcl.mit.edu/database.htm   --------------------
%------------------ for which the following should be  --------------------
%------------------ cited:

%Oliva, A., & Torralba, A. (2001). Modeling the shape of the scene: A holistic representation of the spatial envelope. International Journal of Computer Vision, 42 (3), 145-175.



clear -all;

expBaseIndex    = 1;                      % <-- index for exponential base selection in the array
                                          %     expbase. Set this index between 1 to 5.

expbase    = [sqrt(2),2,2*sqrt(2),exp(1),4]; % exponential base for the case of exponential increase
                                             % in frequency, in which the remaining linear case
                                             % is based on, in order to delimit a maximum high-frequency
                                             % band and to use the same number of grid cells (sets of
                                             % 6 frequency elements).


% number of frequency bands to consider for the linear
%   sampling with fixed bandwidth (according to exp. base)
nModules   = 6;
% formula  nModules = floor(  log(R) / log(expbase)  )       [ and subtract 1 if the result is exact to avoid problems at the edge ]
% nModules = 13 for expbase=sqrt(2)
% nModules = 6 for expbase=2
% nModules = 4 for expbase=2*sqrt(2)
% nModules = 4 for expbase=e
% nModules = 3 for expbase=4


isotropicFreqVariance = 5;

% an angle for each of the nModules modules
angles     = [0, 5, 30,35, 0,5, 30,35, 30,35, 30,35, 30];
delt       = 10;


% definition of a red-green-black map to plot negative-positive-zero
% values respectively in the frequency domain.
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
redGreenBlackMap = [redColorMap; greenColorMap; zeros(1, 256)]';

% FOR COMMAND LINE LINUX: matlab -nodisplay -nodesktop -r "run /path/file.m"

for set_ind=1:2
for file_ind=10:59

%fileName = strcat(strcat('DATASET/set0',strcat(num2str(set_ind),'/')),strcat(num2str(file_ind),'.jpg'));
fileName = 'image.jpg';%'whitenoise';     % <-- file name for the image to process,
                                          %     which should be a square image and
                                          %     in a format that MATLAB can read as
                                          %     a colour image with 3 channels (3D
                                          %     array) such as a colour JPG image.
                                          %     One channel from this colour image
                                          %     is then used to work in grayscale.
                                          %     If set to 'whitenoise', an grayscale
                                          %     image with uniform white noise is
                                          %     used instead.

outFolder = 'out_lossy_compression/';     % <-- output folder.
outFileResults = 'RESULTS_expbase01.txt'; 

% ------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------

if( strcmp(fileName, 'whitenoise') )
    fileId    = fileName;
    outSubFold= fileId;
    f_path = strcat(outFolder,strcat(strcat('set0',num2str(set_ind)),strcat('_',outSubFold)));
    mkdir(f_path);
    R      = 128;
    N      = 3*R*R;
    S      = N*ones(6,1);                 % <-- scaling factors for amplitude of each frequency point.
    % generate uniform white noise image vector
    fhex   = rand(N,1);
    %      map between -1 and 1
    fhex   = 2*fhex - 1.0;
    
    % obtain indices Ind for the white noise image fhex, its HDFFT version fhexP,
    %   and the cartesian coordinates (cartX and cartY) in which to visualize the results.
    [fhex,fhexP,cartX,cartY,Ind] = setCartesianAndHexSampling(fhex,zeros(2*R,2*R),1,true);
    
else
    fileId    = num2str(file_ind);
    outSubFold= fileId;
    f_path = strcat(outFolder,strcat(strcat('set0',num2str(set_ind)),strcat('_',outSubFold)));
    mkdir(f_path);
    fxy    = imread(fileName);
    fxy    = fxy(:,:,1);
    N1     = size(fxy,2);
    N2     = size(fxy,1);
    if( N1 ~= N2 )
       disp('The width and height for the image file do not match. Please provide a square size image.');
       return;
    end
    R      = N2/2;
    N      = 3*R*R;
    S      = N*ones(6,1);                 % <-- scaling factors for amplitude of each frequency point.
    % Use mat2gray to scale the image between 0 and 1
    fxy    = mat2gray(fxy(:,:));
    fhex   = zeros(N,1);
    
    % hexagonally sample the grayscale image fxy and store the result in fhex,
    %   obtain its indices Ind, its HDFFT version fhexP,
    %   and the cartesian coordinates (cartX and cartY) in which to visualize the results.
    [fhex,fhexP,cartX,cartY,Ind] = setCartesianAndHexSampling(fhex,fxy,1,false);

    
end
    out_path = strcat(f_path,'/');
    %graphics_toolkit gnuplot
    % plot the hexagonally sampled image in cartesian coordinates
    fig1 = figure('Visible','off');
    scatter(cartX,cartY,16.5,fhex,'o','filled');
    colormap(gray);
    axis image off; %axis('tic','off');
    set(gcf,'units','pixel'); set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(fig1, strcat(strcat(out_path,'HEX_SAMP_IMG' ),'.png'),'png');
    
    % -------------------------------------------------------------
    % -------------------------------------------------------------
    % Obtain the forward hexagonal DFT
    FHEXP = HDFFT(fhexP,-1);
    FHEX  = changeOrdering(FHEXP);

    cartX = zeros(N,1);
    cartY = zeros(N,1);
    h           = 1.0;
    oX          = 0.5*double(R) +1;
    ind         = 1;
    for r1=0:(2*R-1)
       for r2=0:(2*R-1)
          curX = oX + h*(r1 - 0.5*r2);
          curY = 2*R - h*0.5*sqrt(3.0)*r2;
          if( r1-r2>=-R && r1-r2<R)
             cartX(ind) =  curX;
             cartY(ind) =  2*R - curY;
             ind = ind + 1;
          end
       end
    end
    %  plot the real part of the Hexagonal DFT (FHEX).
    maxVal     = max(log(abs(real(FHEX)) + 1));
    fig2_1 = figure('Visible','off');
    scatter(cartX,cartY,16.5,sign(real(FHEX)).*log(abs(real(FHEX)) + 1),'o','filled');
    colormap(redGreenBlackMap);
    caxis([-maxVal, maxVal]);
    axis image off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(fig2_1, strcat(strcat(out_path,'HDFT_REAL_PART' ),'.png'),'png');
    
    %  plot the imaginary part of the Hexagonal DFT (FHEX)
    fig2_2 = figure('Visible','off');
    scatter(cartX,cartY,16.5,sign(imag(FHEX)).*log(abs(imag(FHEX)) + 1),'o','filled');
    colormap(redGreenBlackMap);
    caxis([-maxVal, maxVal]);
    axis image off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(fig2_2, strcat(strcat(out_path,'HDFT_IMAGINARY_PART' ),'.png'),'png');
    
    
    % -------------------------------------------------------------
    % -------------------------------------------------------------
    % Band-pass filter in the hexagonal frequency domain
    % considering frequencies centered at each sampled grid frequency
    % (in both linear and exponential sampling).
    
    
    
    % declare one array for visualizing each sampling type 
    %    (exp, linear fixed and linear varied respectively).
    FhexBandPass   = zeros(size(FHEX));
    FhexBandPass(:,2) = zeros(size(FHEX));
    %    always consider DC (average of all frequencies)
    FhexBandPass(Ind(R+1,R+1),1) = FHEX(Ind(R+1,R+1));
    FhexBandPass(Ind(R+1,R+1),2) = FHEX(Ind(R+1,R+1));
    % declare one array for the inverse Hex DFT for each sampling
    %    type (exp, linear fixed and linear varied respectively).
    FhexBandPassP   = zeros(size(FHEXP));
    FhexBandPassP(:,:,2) = zeros(size(FHEXP));
    %    always consider DC (average of all frequencies)
    FhexBandPassP(1,1,1) = FHEXP(1,1);
    FhexBandPassP(1,1,2) = FHEXP(1,1);
    
    % Produce the filtering with pass-band grid cells (activating
    %   six kronecer-delta frequency amplitudes in the corners of
    %   an hexagon corresponds to summing the resulting grid
    %   patterns in the spatial domain).
    freqPoints  = zeros(N,2);% -- array for visualizing the filters
    totalNgridC = zeros(2,1);% -- array for counting the number of used grid cells
    %   --sampType indexes the sampling type (1 - Exponential, 2 - Linear),
    %   --omega is the frequency in cycles per hexagonal region of support
    %        (cycles per image).
    for sampType = 1:2
       for  f=1:nModules
          if( sampType == 1 )
             % Exponential sampling with fixed bandwidth (EF).
             omega = expbase(expBaseIndex) ^ f;
          else
             omega = f*(expbase(expBaseIndex) ^ (nModules)) / nModules;
          end
          xCurrMod = round(omega*cos(angles(f)*pi / 180.0));
          yCurrMod = round(omega*sin(angles(f)*pi / 180.0));
          
          % DC component for visualization of filters in array freqPoints
          freqPoints(Ind(R+0+1,R+0+1),sampType) = 1;
          if( totalNgridC(sampType,1) == 0 )
             totalNgridC(sampType,1) = 1;
          end
          % Perform the actual filtering with circular frequency bands sampled
          %    according to the filtering type (sampType) in the hexagonal domain.
          freqPointsModule = zeros(N);
          
          for r1=(-R):(0+R-1)
             for r2=(-R):(1*R-1)
               % conditions for delimitation of circular frequency bands according
               %    to the current frequency omega and the current bandwidth bw
               %    centered around omega.
               condC1 = ((r1-0.5*r2) - xCurrMod)^2 + ((0.5*sqrt(3)*r2) - yCurrMod)^2 <= isotropicFreqVariance;
               if( r1-r2>=-R && r1-r2<R && condC1 )
                   % linear sampling (both fixed and varying bandwidth cases)
                   if( sampType >= 2 )
                       if(totalNgridC(sampType,1) < totalNgridC(1,1) )
                          freqPoints(Ind(R+r1+1,R+r2+1),sampType) = 1;
                          freqPointsModule(Ind(R+r1+1,R+r2+1)) = 1;
                          totalNgridC(sampType,1) = totalNgridC(sampType,1) + 1;
                          FhexBandPass(Ind(R+r1+1,R+r2+1),sampType) = FHEX(Ind(R+r1+1,R+r2+1));
                          if(r2 >= 0)
                              if( r1 >= 0 )
                                  FhexBandPassP(r1+1,r2+1,sampType) = FHEXP(r1+1,r2+1);
                              else
                                  FhexBandPassP(r1+3*R+1,r2+1,sampType)  = FHEXP(r1+3*R+1,r2+1);
                              end
                          else
                              FhexBandPassP(r1+2*R+1,r2+R+1,sampType) = FHEXP(r1+2*R+1,r2+R+1);
                          end
                       end
                   % exponential sampling
                   else
                       freqPoints(Ind(R+r1+1,R+r2+1),sampType) = 1;
                       freqPointsModule(Ind(R+r1+1,R+r2+1)) = 1;
                       totalNgridC(sampType,1) = totalNgridC(sampType,1) + 1;
                       FhexBandPass(Ind(R+r1+1,R+r2+1),sampType) = FHEX(Ind(R+r1+1,R+r2+1));
                       if(r2 >= 0)
                          if( r1 >= 0 )
                             FhexBandPassP(r1+1,r2+1,sampType) = FHEXP(r1+1,r2+1);
                          else
                             FhexBandPassP(r1+3*R+1,r2+1,sampType) = FHEXP(r1+3*R+1,r2+1);
                          end
                       else
                          FhexBandPassP(r1+2*R+1,r2+R+1,sampType) = FHEXP(r1+2*R+1,r2+R+1);
                       end
                       
                   end
                   
                   % activate the rest of the five frequency points for the
                   %   current grid cell by rotating the current points by an
                   %   hexagonal rotation matrix Rot.
                   %   ------------------------------------------------------
                   theta =  pi / 3.0 ;
                   Rot = [cos(theta)+(sin(theta)/sqrt(3))  (2/sqrt(3))*sin(theta); -(2/sqrt(3))*sin(theta)   cos(theta)-(sin(theta)/sqrt(3)) ];
                   k = [ r1  r2];
                   repeated = false;
                   for j=1:5
                       k = k * Rot;
                       freqPoints(Ind(R+round(k(1))+1,R+round(k(2))+1),sampType)   = 1;
                       if( freqPointsModule(Ind(R+round(k(1))+1,R+round(k(2))+1)) == 1 )
                          repeated = true;
                       end
                       FhexBandPass(Ind(R+round(k(1))+1,R+round(k(2))+1),sampType) = FHEX(Ind(R+round(k(1))+1,R+round(k(2))+1));
                       if(round(k(2)) >= 0)
                           if( round(k(1)) >= 0 )
                                FhexBandPassP(round(k(1))+1,round(k(2))+1,sampType) = FHEXP(round(k(1))+1,round(k(2))+1);
                           else
                               FhexBandPassP(round(k(1))+3*R+1,round(k(2))+1,sampType)  = FHEXP(round(k(1))+3*R+1,round(k(2))+1);
                           end
                       else
                           FhexBandPassP(round(k(1))+2*R+1,round(k(2))+R+1,sampType) = FHEXP(round(k(1))+2*R+1,round(k(2))+R+1);
                       end
                   end
                   if( repeated == true )
                      totalNgridC(sampType,1) = totalNgridC(sampType,1) - 1;
                   end
                   %   ------------------------------------------------------
                   
               end
             end
          end % --END of actual filtering in hex domain
       end    % --END of sampType cycle
       
       disp(' Total number of grid cells (all entries should be the same): ');
       disp(totalNgridC);
       % plot the filter (frequency bands) in cartesian coordinates
       figFreqF = figure('Visible','off');
       scatter(cartX,cartY,16.5,real(freqPoints(:,sampType)),'o','filled');
       colormap(redGreenBlackMap);
       caxis([-1, 1]);
       axis image off; set(gca,'LooseInset',get(gca,'TightInset'));
       figTitle = strcat(out_path,'FILTER_REAL_PART' );
       figTitle = strcat(figTitle,'_nGridModules_');
       figTitle = strcat(figTitle,num2str(nModules));
       figTitle = strcat(figTitle,'_sampType_');figTitle = strcat(figTitle,num2str(sampType));
       figTitle = strcat(figTitle,'.png');
       saveas(figFreqF, figTitle,'png');
       
       % plot the compressed frequency domain Image in cartesian coordinates
       % -----------------------------------------------------------------
       %    --------REAL PART--------
       figFreqFT = figure('Visible','off');
       scatter(cartX,cartY,16.5,sign(real(FhexBandPass(:,sampType))).*log(abs(real(FhexBandPass(:,sampType))) + 1),'o','filled');
       colormap(redGreenBlackMap);
       caxis([-maxVal, maxVal]);
       axis image off;set(gca,'LooseInset',get(gca,'TightInset'));
       figTitle = strcat(out_path,'COMPRESSED_FREQUENCY_REAL' );
       figTitle = strcat(figTitle,'_nGridModules_');
       figTitle = strcat(figTitle,num2str(nModules));
       figTitle = strcat(figTitle,'_sampType_');figTitle = strcat(figTitle,num2str(sampType));
       figTitle = strcat(figTitle,'.png');
       saveas(figFreqFT, figTitle,'png');
       %    --------IMAGINARY PART---
       figFreqFT = figure('Visible','off');
       scatter(cartX,cartY,16.5,sign(imag(FhexBandPass(:,sampType))).*log(abs(imag(FhexBandPass(:,sampType))) + 1),'o','filled');
       colormap(redGreenBlackMap);
       caxis([-maxVal, maxVal]);
       axis image off;set(gca,'LooseInset',get(gca,'TightInset'));
       figTitle = strcat(out_path,'COMPRESSED_FREQUENCY_IMAGINARY' );
       figTitle = strcat(figTitle,'_nGridModules_');
       figTitle = strcat(figTitle,num2str(nModules));
       figTitle = strcat(figTitle,'_sampType_');figTitle = strcat(figTitle,num2str(sampType));
       figTitle = strcat(figTitle,'.png');
       saveas(figFreqFT, figTitle,'png');
       
    end
    
    %outFile = fopen(strcat(out_path,'RESULTS.txt'),'a');
    outFile = fopen(strcat(outFolder,outFileResults),'a');
    
    % Obtain inverse hexagonal DFT
    % ----------------------------------------------------------------------
    % ----------------------------------------------------------------------
    for indexHDFT=1:2
       % obtain inverse
       [IFhexP] = HDFFT(FhexBandPassP(:,:,indexHDFT),1);
       % chage to visualizing mode
       IFhex = changeOrdering(IFhexP);
       
       % plot the reconstructed image in the spatial domain in cartesian coordinates
       fig4 = figure('Visible','off');
       scatter(cartX,cartY,16.5,real(IFhex),'o','filled');
       colormap(gray);
       axis image off;set(gca,'LooseInset',get(gca,'TightInset'));
       figTitle = strcat(out_path,'RECONSTRUCTED_INVERSE' );
       figTitle = strcat(figTitle,'_sampType_');figTitle = strcat(figTitle,num2str(indexHDFT));
       figTitle = strcat(figTitle,'.png');
       saveas(fig4, figTitle,'png');
       %
       % Obtain correlation of the current type of
       %    sampling with the original image
       corr_coef = corr( real(IFhex) , real(fhex) );
       % Obtain Mean-Squared Error with the original image
       %mse_val = immse( real(IFhex) , real(fhex) );
       mse_val = (norm(IFhex(:)-fhex(:)).^2)/numel(IFhex);
       % ------------------------------------------------------
       % Print the results to a text file----------------------
       if( indexHDFT == 1 )
          fprintf(outFile,'IMAGE ID: %s\n',fileName);
          fprintf(outFile,'--EXP BASE = %f , BANDWIDTH RADIUS = %f \n', expbase(expBaseIndex),isotropicFreqVariance);
          fprintf(outFile,'--------EXPONENTIAL SAMPLING (%d grid cells)------\n',totalNgridC(1,1));
          fprintf(outFile,'--------MSE  = %f\n',mse_val);
          fprintf(outFile,'--------Corr = %f\n',corr_coef);
       else
          fprintf(outFile,'--------LINEAR SAMPLING (%d grid cells)----------\n',totalNgridC(2,1));
          fprintf(outFile,'--------MSE  = %f\n',mse_val);
          fprintf(outFile,'--------Corr = %f\n',corr_coef);
       end
    end
    

end
end
