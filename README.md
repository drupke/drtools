## OVERVIEW

`drtools` is an IDL repository of general-purpose software.

It contains a number of routines for multicore processing in
IDL. These include four routines written by Robert da Silva which are
no longer available online, so I have bundled them in this repo in
subdirectory `rds`.

`s7nad` subdirectory contains plotting routines for [Rupke et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.4748R/abstract).

## REQUIREMENTS

IDL v8.0 or higher (tested with v8.5)

IDL libraries:
- [IDL Astronomy User's Library](http://idlastro.gsfc.nasa.gov)
- [Coyote](http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD), for graphics AND undefine.pro
  - or from the [GitHub repository](https://github.com/davidwfanning/idl-coyote/tree/master/coyote)
- SPLIT_FOR requires X_CHKFIL from [XIDL](https://github.com/profxj/xidl)

Note that the IDL Astronomy User's Library ships with some Coyote
routines. However, it's not clear how well these libraries keep track
of each other, so it may be preferable to download each package
separately and delete the redundant routines that ship within other
packages.

## USAGE

As needed.

## QUESTIONS? BUGS? WANT TO MODIFY THE CODE?

Feel free to contact David Rupke at drupke@gmail.com with questions,
bug reports, etc.

Modifications are encouraged, but subject to the license.

## LICENSE AND COPYRIGHT

Copyright (C) 2016--2021 David S. N. Rupke

These programs are free software: you can redistribute them and/or
modify them under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License or any later version.

These programs are distributed in the hope that they will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with these programs.  If not, see http://www.gnu.org/licenses/.
