# wavespectrum
Generation of a composed wave spectrum for a given set Hm0 and He10 values

The purpose of the program is to produce a 1D wave energy spectrum based op (expected) 
Hm0 and He10 values that can be supllied to a ship motion calculation program.

Although for the Northsea a Jonswap spectrum seems a logical choice, from a large number 
of wave spectra measured in the past it appeared that the measured combinations of Hm0 and He10
are, in many cases, not in the range of a Jonswap spectrum.

Thus as a solution the superposition of two Gamma wave spectra has been chosen. The desired Hm0 
provides the condition for the sum of both M0. The peak periods of both gamma spectra are adjusted 
to ensure that the desired He10 is covered between them, and finally the ratio of the two spectra 
is ieterated to arrive at the He10 value (within 1 cm).

Note - because this is my first Julia program many improvements will be possible...
