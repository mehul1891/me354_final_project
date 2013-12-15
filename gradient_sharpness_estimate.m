%=========================================================================%
% GRADIENT BASED SHARPNESS METRIC          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : To use a gradient based metric to deterrmine the
% level of blurriness that an image has.
% Contact information    : deman@stanford.edu & moswal@stanford.edu
% This function os based on the work of Tolga Birdal
% <http://www.tbirdal.me>
%
% Copyright (c) 2011, Tolga Birdal
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
%       the documentation and/or other materials provided with the 
%       distribution
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
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS
% u     = Image which is going to be evaluated for blurriness
%
% OUTPUT OPTIONS
% sharpness = A gradient based metric that quantifies the level of
% burriness. The higher the value the sharper the image is.
%=========================================================================%

function [sharpness]=gradient_sharpness_estimate(u)

[ux, uy]=gradient(u);

S=sqrt(ux.*ux+uy.*uy);
sharpness=sum(sum(S))./(numel(ux));

end