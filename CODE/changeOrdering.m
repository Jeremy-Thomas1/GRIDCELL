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


% Function to change the ordering of the hexagonal array elements
% between origin-at-center mode (for visualization) and HDFFT mode.
%
% Input:  array x_in whose order needs to be changed.
%
% Ouput:  array x_out with the order changed. If x_in is in
%         origin-at-center mode, x_out is a 2D array compatible with
%         the HDFFT, and if x_in is in HDFFT mode, x_out is a 1D
%         array in origin-at-center mode.
function [x_out] = changeOrdering(x_in)
   % change from origin-at-center mode to HDFFT mode.
   if( size(x_in,2) == 1 )
      R = int16( sqrt( size(x_in,1) / 3 ) );
      x_out = zeros(3*R,R);
      ind         = 1;
      for r1=0:(2*R-1)
       for r2=0:(2*R-1)
             if( r1-r2>=-R && r1-r2<R)
                if(r2 >= R)
                    if( r1 >= R )
                        x_out(r1-R+1,r2-R+1) = x_in(ind);
                    else
                        x_out(r1+2*R+1,r2-R+1) = x_in(ind);
                    end
                else
                    x_out(r1+R+1,r2+1) = x_in(ind);
                end
                ind = ind + 1;
             end
       end
      end
   % change from HDFFT mode to origin-at-center mode.
   else
      R = size(x_in,2);
      x_out = zeros(3*R*R,1);
      ind         = 1;
      for r1=0:(2*R-1)
       for r2=0:(2*R-1)
             if( r1-r2>=-R && r1-r2<R)
                if(r2 >= R)
                    if( r1 >= R )
                        x_out(ind) = x_in(r1-R+1,r2-R+1);
                    else
                        x_out(ind) = x_in(r1+2*R+1,r2-R+1);
                    end
                else
                    x_out(ind) = x_in(r1+R+1,r2+1);
                end
                ind = ind + 1;
             end
       end
      end
   end
