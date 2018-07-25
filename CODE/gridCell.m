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


% Function to simulate one grid field
% on the spatial domain.
% Input:  horizontal and vertical offset n01 and n02 in
%         cartesian coordinates, orientation angle theta
%         (in radians), frequency omega for the pattern (integer)
%         specified in number of cycles across the hexagonal region
%         of support, 6x2 array fS with values to scale the value
%         of each frequency point ( the default input value
%         should be simply false if no rescaling is required ),
%         1D array S with values to scale the amplitude of each
%         frequency point and where the phase shift will be applied,
%         radius R of circle that contains the hexagonal region of
%         support and standard deviation sigma for noise to add to
%         each frequency point in hexagonal coordinates (drawn from
%         a normal distribution with mean=0 and sd=sigma).
% Ouput:  x complex signal with the grid cell pattern in the spatial
%         domain (stored in a 1D array), and complex resulting input
%         signal X in the frequency domain.
function [x,X] = gridCell(n01,n02,theta,omega,fS,S,R,sigma)
   % set offset (r0) in hexagonal coordinates for phase shift.
   n  = 6;
   r0 = zeros(2,1);
   r0(1) = n01 + (n02 / sqrt(3.0));
   r0(2) = (2 * n02 / sqrt(3.0));
   %disp('current offset in hexagonal coordinates:')
   %disp(r0)
   
   % define frequency domain points
   % -------------------------------------------------------------
   % -------------------------------------------------------------
   % define 6 points on corners of a regular
   % hexagon on the hexagonal frequency domain.
   k  = zeros(n,2);
   noise_v = zeros(n,2);
   if( sigma > 0 )
      for j=1:n
         noise_v(j,1) = normrnd(0,sigma);
         noise_v(j,2) = normrnd(0,sigma);
      end
   end
   if( fS == false )
	   k(1,1) =  omega + noise_v(1,1);
	   k(1,2) =  0             + noise_v(1,2);
	   % rescaling of room condition
	   %k(1,1) =  2*omega + 0 + noise_v(1,1);
	   %k(1,2) =  0        + noise_v(1,2);
	   
	   
	   k(2,1) =  0              + noise_v(2,1);
	   k(2,2) =  omega  + noise_v(2,2);
	   %  rescaling of room condition
	   %k(2,1) =  0.5*omega        + noise_v(2,1);
	   %k(2,2) =  2*omega - 0 + noise_v(2,2);
	   %k(2,1) =  -0.5*omega        + noise_v(2,1);
	   %k(2,2) =  1*omega - 0 + noise_v(2,2);
	   
	   k(3,1) =  omega + noise_v(3,1);
	   k(3,2) =  omega + noise_v(3,2);
	   %  rescaling of room condition
	   %k(3,1) =  1*omega + 0.5*omega - 0 + noise_v(3,1);
	   %k(3,2) =  2*omega - 0 + noise_v(3,2);
	   %k(3,1) =  1.5*omega + noise_v(3,1);
	   %k(3,2) =  1*omega - 0 + noise_v(3,2);
	   
	   %          reflection about origin of previous 3 points.
	   k(4,1) = -omega + noise_v(4,1);
	   k(4,2) =  0             + noise_v(4,2);
	   %  rescaling of room condition
	   %k(4,1) = -2*omega - 0 + noise_v(4,1);
	   %k(4,2) = 0        + noise_v(4,2);
	   
	   k(5,1) =  0             + noise_v(5,1);
	   k(5,2) = -omega + noise_v(5,2);
	   %  rescaling of room condition
	   %k(5,1) =  -0.5*omega        + noise_v(5,1);
	   %k(5,2) = -2*omega +0 + noise_v(5,2);
	   %k(5,1) =  0.5*omega        + noise_v(5,1);
	   %k(5,2) = -1*omega + noise_v(5,2);
	   
	   k(6,1) = -omega + noise_v(6,1);
	   k(6,2) = -omega + noise_v(6,2);
	   %  rescaling of room condition
	   %k(6,1) = -1*omega - 0.5*omega +0 + noise_v(6,1);
	   %k(6,2) = -2*omega +0 +noise_v(6,2);
	   %k(6,1) = -1.5*omega + noise_v(6,1);
	   %k(6,2) = -1*omega +noise_v(6,2);
   else
           for j=1:n
              k(j,1) = fS(j,1)*omega + noise_v(j,1);
              k(j,2) = fS(j,2)*omega + noise_v(j,2);
           end
   end
   
   
   % rotate the frequency points by angle theta
   % on hexagonal coordinates.
   Rot = [cos(theta)+(sin(theta)/sqrt(3))  (2/sqrt(3))*sin(theta); -(2/sqrt(3))*sin(theta)   cos(theta)-(sin(theta)/sqrt(3)) ];
   for j=1:n      
      k(j,:) = k(j,:) * Rot;
      % Also set the phase for each frequency point
      S(j,1) = exp(-i * (pi / R) * ( (4/3)*(r0(1)-0.5*r0(2))*(k(j,1)-0.5*k(j,2)) + (r0(2)*k(j,2)) ) ) * S(j,1);
   end
   
   X     = zeros(3*R*R,1);  %  <-- array for visualizing origin at center
   ind         = 1;
   for r1=(-R):(R-1)
    for r2=(-R):(R-1)
          if( r1-r2>=-R && r1-r2<R )
            for j=1:n
             if( round(k(j,1))==r1 && round(k(j,2))==r2 )
                X(ind)  = S(j,1);
             end
            end
            ind = ind + 1;
          end
    end
   end
   XP    = changeOrdering(X);    %  <-- array for computing the HDFFT
   % -------------------------------------------------------------
   % -------------------------------------------------------------

   % inverse Fourier Transform the hexagonal frequency points
   xp  =  HDFFT(XP,1);
   x   =  changeOrdering(xp);
   
