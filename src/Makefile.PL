### This file is part of SunlightLB - a constraint-based modelling
### package.

### Copyright (C) 2005-2010 Unilever UK Central Resources Ltd.

### SunlightLB is free software: you can redistribute it and/or
### modify it under the terms of the GNU General Public License as
### published by the Free Software Foundation, either version 2 of the
### License, or (at your option) any later version.

### SunlightLB is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.

### You should have received a copy of the GNU General Public License
### along with SunlightLB.  If not, see <http://www.gnu.org/licenses/>.

use ExtUtils::MakeMaker;
use Env qw(LB_INC LB_LIB);

WriteMakefile(
    NAME => 'SunlightLB',
    LIBS => [$LB_LIB],
    INC => $LB_INC,
    OBJECT => 'SunlightLB_wrap.o',
);
