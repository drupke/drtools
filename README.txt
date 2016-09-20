-------------------------------------------------------------------------
OVERVIEW
-------------------------------------------------------------------------

DRTOOLS is an IDL repository of general-purpose software. The idea is
for it to be an ever-expanding toolbox from my routine catalog, as I
get around to properly commenting and licensing old routines or make
new ones. How well this works remains to be seen, but I have seeded it
with four routines.

-------------------------------------------------------------------------
REQUIREMENTS
-------------------------------------------------------------------------

IDL v8.0 or higher (tested with v8.5)

IDL libraries:
- IDL Astronomy User's Library, for various routines
  http://idlastro.gsfc.nasa.gov
- Coyote, for graphics AND undefine.pro
   http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD
  [or from the GitHub repository:
   https://github.com/davidwfanning/idl-coyote/tree/master/coyote]

Note that the IDL Astronomy User's Library ships with some Coyote
routines. However, it's not clear how well these libraries keep track
of each other, so it may be preferable to download each package
separately and delete the redundant routines that ship within other
packages.

-------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------

As needed.

-------------------------------------------------------------------------
QUESTIONS? BUGS? WANT TO MODIFY THE CODE?
-------------------------------------------------------------------------

Feel free to contact David Rupke at drupke@gmail.com with questions,
bug reports, etc.

Modifications are encouraged, but subject to the license.

-------------------------------------------------------------------------
LICENSE AND COPYRIGHT
-------------------------------------------------------------------------

Copyright (C) 2016 David S. N. Rupke

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
