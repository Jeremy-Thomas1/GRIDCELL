----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
---- This file describes the files that are part of the Grid Cell Hexagonal   ----
---- Fourier Transform (GCHFT) software for simulating grid fields and other  ----
---- spatial-domain signals based on an Hexagonal Fourier Transform.          ----

---- The GCHFT is free software: you can redistribute it and/or modify        ----
---- it under the terms of the GNU General Public License as published by     ----
---- the Free Software Foundation, either version 3 of the License, or        ----
---- (at your option) any later version.                                      ----
----                                                                          ----
---- The GCHFT software is distributed in the hope that it will be useful,    ----
---- but WITHOUT ANY WARRANTY; without even the implied warranty of           ----
---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            ----
---- GNU General Public License for more details.                             ----
----                                                                          ----
---- You should have received a copy of the GNU General Public License        ----
---- along with this software.  If not, see <http://www.gnu.org/licenses/>.   ----
----                                                                          ----
---- Copyright 2018 Ulises Rodriguez Dominguez and Jeremy B. Caplan.          ----
----------------------------------------------------------------------------------



The following is the list of demo files for the simulations (for each demo file
an output folder needs to be specified):
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
----demoGC_Noise.m		Demo file that simulates a distorted grid field, and
				it can also be used to simulate a regular grid field.
----demoGC_Rescale.m		Demo file that simulates an change in the scale
				(or frequency) along one direction of the grid field.
----demoGridPlaceCell.m		Demo file that simulates a place field based on
				a set of grid fields.
----demoLossyCompression.m	Demo file that simulates a population of grid cells
				(set of frequency components) performing filtering
				on either white noise or natural images (lossy
                                compression of frequency information). This file
				requires an input image to be specified (see the file
				for the details). Exponential frequency sampling is
				compared against two different types of linear
				frequency sampling.


The following is the list of function files used by the demo files:
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------

----changeOrdering.m		Function file used by the demoLossyCompression.m
				and the gridCell.m files.
----gridCell.m			Function file to simulate one grid field. This
				file is used by all the demo files except for the
				demoLossyCompression.m file.
----HDFFT.m			Function file that implements Mersereau's 1979
				Hexagonal Fourier Transform algorithm. This
				function will have further optimizations for
				future implementations. This file is used by all
				the demo files (directly or through the gridCell.m
				file).
----setCartesianAndHexSampling.m
				Function file to perform hexagonal sampling of
				grayscale images and obtain proper coordinates
				for visualization purposes. This file is used
				by the demoLossyCompression.m file.


