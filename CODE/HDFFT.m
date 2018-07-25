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


% Hexagonal DFT (naive slow version, which will be updated to a fast version
%                in a future release).
%
% Input:  hexagonal array x to be Fourier-transformed
%         of size (3*R,R) where R is the radius of the
%         hexagonal region of support, and scalar d
%         to indicate if a forward transform is required (d=-1)
%         or an inverse transform is required (d=1).
%
% Ouput:  hexagonal array X Fourier-transformed of size (3*R,R).
function [X] = HDFFT(x,d)
   R = size(x,2);
   X = zeros(3*R,R);
   Z = 1;
   if( d == 1 )
       Z = 3*R*R;
   end
   for k1=0:(floor(3*R/2)-1)
    for k2=0:(floor(R/2)-1)
       F = HDFFT_R(x(1:2:(3*R),1:2:R),d,k1,k2);
       G = HDFFT_R(x(1:2:(3*R),2:2:R),d,k1,k2);
       H = HDFFT_R(x(2:2:(3*R),1:2:R),d,k1,k2);
       I = HDFFT_R(x(2:2:(3*R),2:2:R),d,k1,k2);
       WG = exp( d * i * (2*pi/(3*R)) * (2*k2-k1) ) * G;
       WH = exp( d * i * (2*pi/(3*R)) * (2*k1-k2) ) * H;
       WI = exp( d * i * (2*pi/(3*R)) *  (k1+k2)  ) * I;
       X(k1+1,k2+1) = (F + WG + WH + WI) / Z;
       X(k1+1+floor(3*R/2),k2+1) = (F - WG + WH - WI) / Z;
       X(k1+1+R,k2+1+floor(R/2)) = (F + WG - WH - WI) / Z;
       X(mod(floor(k1+(5*R/2)),3*R)+1,k2+1+floor(R/2)) = (F - WG - WH + WI) / Z;
    end
   end


function [X] = HDFFT_R(x,d,k1,k2)
   R = size(x,2);
   
   if( R > 0 )  %
      % compute HDFT by the naive definition
      X = zeros(1,1);%zeros(3*R,R);
          for r1=0:(3*R-1)
            for r2=0:(R-1)
               s1     = (2*r1-r2)*(2*k1-k2) / (3*R);
               s2     = r2*k2 / R;
               expVar = exp(d*i*pi* ( s1 + s2 ) );
               X(1,1) = X(1,1) + x(r1+1,r2+1)*expVar;
            end
          end
      return
   end
