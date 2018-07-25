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


% Function to hexagonally sample a grayscale image and obtain cartesian
% coordinates for visualization purposes.
%
% Input:  array fhex with either zeros (zeros(N,1)) or with noise values
%         where the hexagonally sampled values will be stored, array fxy
%         with the grayscale image to hexagonally sample and with the same
%         width and height (if noise = true, then fhex is assumed to have
%         noise and no hexagonal sampling is required), scalar h to set
%         the scale for the hexagonal sampling coordinates (default should
%         be set to 1) and logical noise to know  if fhex is filled with
%         noise (noise=true) or not (noise=false in which case hexagonal
%         sampling of fxy is carried out).
%
% Ouput:  array fhex with the hexagonally sampled image,
%         array fhexP with the version of fhex formated for the HDFFT,
%         cartesian coordinate values cartX and cartY for the horizonal
%         and vertical axes respectively and array of indices Ind to
%         have a 2D hexagonal index (r1,r2) for the values in fhex.
function [fhex,fhexP,cartX,cartY,Ind] = setCartesianAndHexSampling(fhex,fxy,h,noise)
   % obtain dimensions from hexagonal region and declare memory
   R = int16( sqrt( size(fhex,1) / 3 ) );
   N = 3*R*R;
   fhexP       = zeros(3*R,R);
   Ind         = zeros(2*R,2*R);
   cartX = zeros(N,1);
   cartY = zeros(N,1);
   % set cartesian coordinates (cartX and cartY), indices (Ind)
   % and fhexP array for HDFFT mode.
   if( noise == true )
      curX        = 0;  curY        = 0;
      oX          = 0.5*double(R) +1 ;
      ind         = 1;
      for r1=0:(2*R-1)
         for r2=0:(2*R-1)
            curX = oX + h*(double(r1) - 0.5*double(r2));
            curY = 2*double(R) - h*0.5*sqrt(3.0)*double(r2);
            if( curX >= 1 && curX <= 2*double(R) && curY <= 2*double(R) && curY >= 1)
               if( r1-r2>=-R && r1-r2<R)
                  if(r2 >= R)
                      if( r1 >= R )
                          fhexP(r1-R+1,r2-R+1) = fhex(ind);
                      else
                          fhexP(r1+2*R+1,r2-R+1) = fhex(ind);
                      end
                  else
                      fhexP(r1+R+1,r2+1) = fhex(ind);
                  end
                  Ind(floor(r1)+1,floor(r2)+1) = ind;
                  cartX(ind) = curX;
                  cartY(ind) = (2*double(R) - curY);
                  ind = ind + 1;
               end
            end
         end
      end
   % set cartesian coordinates (cartX and cartY), indices (Ind)
   % hexagonally sample the image fhex based on the grayscale image
   % in fxy, and obtain fhexP array for HDFFT mode.
   else
      curX        = 0;  curY        = 0;
      incX        = 0;  incY        = 0;
      interpolVal = 0;
      oX          = 0.5*double(R) +1 ;
      ind         = 1;
      for r1=0:(2*R-1)
         for r2=0:(2*R-1)
            curX = oX + h*(double(r1) - 0.5*double(r2));
            curY = 2*double(R) - h*0.5*sqrt(3.0)*double(r2);
            if( curX >=1 && curX <= 2*double(R) && curY<=2*double(R) && curY >= 1)
               % Sample hexagonal image (using bilinear interpolation)
               incX = curX - floor(curX);
               incY = ceil(curY) - curY;
               interpolVal = (1.0-incX)*(1.0-incY) * fxy( ceil(curY), floor(curX) );
               if( floor(curX)+1 <= 2*double(R) )
                  interpolVal = interpolVal + incX*(1.0-incY) * fxy( ceil(curY), floor(curX)+1 );
               end
               if( ceil(curY)-1 >= 1 )
                  interpolVal = interpolVal + (1.0-incX)*incY * fxy( ceil(curY)-1, floor(curX) );
               end
               if( floor(curX)+1 <= 2*double(R) && ceil(curY)-1 >= 1 )
                  interpolVal = interpolVal + incX*incY * fxy( ceil(curY)-1, floor(curX)+1 );
               end
               
               if( r1-r2>=-R && r1-r2<R)
                  fhex(ind) = interpolVal;
                  if(r2 >= R)
                      if( r1 >= R )
                          fhexP(r1-R+1,r2-R+1) = fhex(ind);
                      else
                          fhexP(r1+2*R+1,r2-R+1) = fhex(ind);
                      end
                  else
                      fhexP(r1+R+1,r2+1) = fhex(ind);
                  end
                  Ind(floor(r1)+1,floor(r2)+1) = ind;
                  cartX(ind) = curX;
                  cartY(ind) = (2*double(R) - curY);
                  ind = ind + 1;
               end
            end
            
         end
      end
   end
