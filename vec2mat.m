% ## Copyright (C) 2001 Laurent Mazet
% ##
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 2 of the License, or
% ## (at your option) any later version.
% ##
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ##
% ## You should have received a copy of the GNU General Public License
% ## along with this program; if not, write to the Free Software
% ## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% 
% ## -*- texinfo -*-
% ## @deftypefn {Function File} {@var{m} = } vec2mat (@var{v}, @var{c})
% ## @deftypefnx {Function File} {@var{m} = } vec2mat (@var{v}, @var{c}, @var{d})
% ## @deftypefnx {Function File} [{@var{m}, @var{add}] = } vec2mat (@var{..})
% ##
% ## Converts the vector @var{v} into a @var{c} column matrix with row priority
% ## arrangement and with the final column padded with the value @var{d} to the
% ## correct length. By default @var{d} is 0. The amount of padding added to
% ## the matrix is returned in @var{add}.
% ## @end deftypefn
% 
% ## 2001-02-02
% ##   initial release

function [M, d] = vec2mat (V, c, val)

  switch (nargin)
    case 1,
      M = V;
      return;
    case 2,
      val = 0;
    case 3,
      val = val;
    otherwise
      error ('sage: [M, add] = vec2mat (V, c [, d])');
  end

  V = V.';
  V = V(:);

  r = ceil (length (V) / c);

  d = r * c - length (V);
  if (d ~= 0)
    V = [ V ; val*ones(d, 1) ];
  end

  M = reshape (V, c, r)';

end