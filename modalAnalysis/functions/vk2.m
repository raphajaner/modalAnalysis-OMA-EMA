%% Description: [x,bw,T,xr] = vk2(y,f,fs,r,filtord);
% ------------------------------------------------------------------------------------------------------------
% Function for solving: x=inv(r*r*A'*A+I)*C*y, where C = diag(exp(-j*w(1)),...,exp(-j*w(N))), 
% w = 2*pi*cumsum(f)*dt. This function appears in second generation Vold-Kalman filtering. 
% This function is for extracting a single order at a time (not simultaneous orders).
% ------------------------------------------------------------------------------------------------------------
% Input:
%   y    := N x 1, data vector
%   f    := N x 1, frequency vector [Hz] (y and f have the same length)
%   fs   := 1 x 1, sampling frequency [Hz]
%   r    := 1 x 1, weighting factor for the structural equation (typical value in 100's or 1000's)
%   filtord := 1 x 1, filter order, one less than the number of poles.  Can take values 1 or 2.
%
% Output:
%   x    := N x 1, output order vector minimizing the objective function
%         J=r*r*x'A'*A*x+(y'-x'C')*(y-Cx) (complex)
%         x is the zero-peak complex amplitude envelope
%   bw   := Nord x 1, -3 dB filter bandwidth for each filter [Hz]
%   T    := Nord x 1, 10% - 90% transition time for each filter
%   xr   := N x 1, reconstructed signals for each order (complex)
% ------------------------------------------------------------------------------------------------------------
% Function for solving: x=inv(r*r*A'*A+I)*C*y,
% where C = diag(exp(-j*w(1)),...,exp(-j*w(N))),
% w = 2*pi*cumsum(f)*dt.
% This function appears in second generation Vold-Kalman filtering.
% This function is for extracting a single order at a time (not simultaneous orders).
%
% References:
% J. Tuma, Vold-Kalman order tracking filtration, PDF presentation slides online at
% http://homel.vsb.cz/~tum52/index.php?page=download
% ------------------------------------------------------------------------------------------------------------
% Original version by J. Tuma
% Modified by Scot McNeill, April 2011.
% Copyright (c) 2011, Scot McNeill
% ------------------------------------------------------------------------------------------------------------
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% ------------------------------------------------------------------------------------------------------------
function [x,bw,T,xr] = vk2(y,f,fs,r,filtord)

if ~isvector(y)
 error('y must be a vector.');
end
y=y(:); N = length(y);
if ~isvector(f)
 error('f must be a vector.');
end
f=f(:); N2 = length(f);
if N2 ~= N
 error('length(f) must = length(y).');
end
if ~ismember(filtord,[1,2])
 error('filtord must be element of [1,2].');
end
%
dt = 1/fs;
if filtord == 1
 NR = N-2;
 e = ones(NR,1);
 A = spdiags([e -2*e e],0:2,NR,N);
 bw = fs/(2*pi)*(1.58*r.^-0.5);
 T = 2.85*r.^0.5;
elseif filtord == 2
 NR = N-3;
 e = ones(NR,1);
 A = spdiags([e -3*e 3*e -e],0:3,NR,N);
 bw = fs/(2*pi)*(1.70*r.^-(1/3));
 T =  2.80*r.^(1/3);
else
 error('filtord must be element of [1,2].');
end
AA = r*r*A'*A + speye(N);
jay=sqrt(-1);
ejth = exp(jay*2*pi*cumsum(f)*dt);
yy = conj(ejth).*y;
x = 2*(AA\yy);
xr = real(x.*ejth);
end