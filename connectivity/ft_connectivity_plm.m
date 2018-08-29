function [p] = ft_connectivity_plm(input, varargin)

% FT_CONNECTIVITY_PLM computes the phase linearity measurement from a cell
% array of time-domain data, where each cell is an epoch
%
% Use as
%   [p] = ft_connectivity_plm(input, ...)
%
% The input data input should be organized as a cell-array of nchan x ntime signals
%
% Additional optional input arguments come as key-value pairs:
%   bandwidth	=	scalar, half-bandwidth parameter: the frequency range
%			across which to integrate
%   fsample     =       sampling frequency, needed to convert bandwidth to number of bins
%
% The output p contains the phase lag index in the [0, 1] range.
% The output p is organized as a 3D matrix of nchan x  nchan x ntime signals
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2018, Fabio Baselice, Pierpaolo Sorrentino, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% the sequence of steps is as follows:
%  - Hilbert transformation
%  - multiply with complex conjugate
%  - fft
%  - convert bandwidth parameter to number of bins
%  - integrate over bandwidth

% for k1=1:50
% for k2=1:3
% input{k1,k2}=randn(1024, 1);
% end
% end

B = ft_getopt(varargin, 'bandwidth');
fs = ft_getopt(varargin, 'f_sample');

% NOTE BY JM: if the user inputs data with different length trials, the fft per trial is going 
% to have different frequency resolutions, which is not good. Better to throw an error in that 
% case.
nsmp = cellfun('size', input, 2);
assert(all(nsmp==nsmp(1)), 'currently there is no support for input, where the trials are of different length'); 

for k = 1:numel(input)
  input{k} = hilbert(input{k});
end

nchan=size(input,1);
ntime=size(input,2);
trial_length=length(input{1,1});
ph_min=0.1;
f=(fs/trial_length)*(0:(trial_length-1));
f_integr=(abs(f)<B) | (abs(f-fs)<B);
p=zeros(nchan, nchan, ntime);

for ktime=1:ntime
    for kchan1=1:(nchan-1)
        for kchan2=(kchan1+1):nchan
            temp=fft(input{kchan1, ktime}.*conj(input{kchan2, ktime}));    % NOTE BY FB: The inner cycle can be vectorized
            temp(1)=temp(1).*(abs(angle(temp(1)))>ph_min);  % Volume conduction suppression
            temp=(abs(temp)).^2;
            p_temp=sum(temp(f_integr))./sum(temp);
            p(kchan1, kchan2, ktime)=p_temp;
            p(kchan2, kchan1, ktime)=p_temp;
        end
    end
end
