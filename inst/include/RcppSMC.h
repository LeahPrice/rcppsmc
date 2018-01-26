
// This file is part of RcppSMC.
//
// RcppSMC is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppSMC is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppSMC.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>

#include <adaptMethods.h>
#include <blockpfgaussianopt.h>
#include <history.h>
#include <LinReg.h>
#include <LinReg_LA_adapt.h>
#include <LinReg_LA.h>
#include <moveset.h>
#include <nonLinPMMH.h>
#include <pflineart.h>
#include <pfnonlinbs.h>
#include <population.h>
#include <sampler.h>
#include <smc-exception.h>
#include <smctc.h>
#include <staticModelAdapt.h>
