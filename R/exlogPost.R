#######################################################################
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.

exlogPost <- function(par, type){
  dm <- c(7,2,9,5,3,3,3,5,8,10,5,10,3,3,6,14,5,7,1,2,11)
  dm <- dm[type]
  if(length(par) != dm)
    stop("par of wrong length")

  if(type == 1){
    out <- .Fortran("bwlpst", as.double(par), out = double(1))$out
  }
  if(type == 2){
    out <- .Fortran("bxlpst", as.double(par), out = double(1))$out
  }
  if(type == 3){
    out <- .Fortran("ctlpst", as.double(par), out = double(1))$out
  }
  if(type == 4){
    out <- .Fortran("exlpst", as.double(par), out = double(1))$out
  }
  if(type == 5){
    out <- .Fortran("fllpst", as.double(par), out = double(1))$out
  }
  if(type == 6){
    out <- .Fortran("hrlpst", as.double(par), out = double(1))$out
  }
  if(type == 7){
    out <- .Fortran("k3lpst", as.double(par), out = double(1))$out
  }
  if(type == 8){
    out <- .Fortran("k5lpst", as.double(par), out = double(1))$out
  }
  if(type == 9){
    out <- .Fortran("k8lpst", as.double(par), out = double(1))$out
  }
  if(type == 10){
    out <- .Fortran("lblpst", as.double(par), out = double(1))$out
  }
  if(type == 11){
    out <- .Fortran("lglpst", as.double(par), out = double(1))$out
  }
  if(type == 12){
    out <- .Fortran("lmlpst", as.double(par), out = double(1))$out
  }
  if(type == 13){
    out <- .Fortran("mtlpst", as.double(par), out = double(1))$out
  }
  if(type == 14){
    out <- .Fortran("nllpst", as.double(par), out = double(1))$out
  }
  if(type == 15){
    out <- .Fortran("nrlpst", as.double(par), out = double(1))$out
  }
  if(type == 16){
    out <- .Fortran("pglpst", as.double(par), out = double(1))$out
  }
  if(type == 17){
    out <- .Fortran("phlpst", as.double(par), out = double(1))$out
  }
  if(type == 18){
    out <- .Fortran("prlpst", as.double(par), out = double(1))$out
  }
  if(type == 19){
    out <- .Fortran("pslpst", as.double(par), out = double(1))$out
  }
  if(type == 20){
    out <- .Fortran("rdlpst", as.double(par), out = double(1))$out
  }
  if(type == 21){
    out <- .Fortran("tnlpst", as.double(par), out = double(1))$out
  }
  out
}
