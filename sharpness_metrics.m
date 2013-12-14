%=========================================================================%
% SHARPNESS METRIC COMPARISON          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Function used to condense all the sharpness
% metrics
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS
% u     = Image which is going to be evaluated for blurriness
% ref   = Reference image used in the similarity index. If there is no
% reference image put "0".
%
% OUTPUT OPTIONS
% gradient = A gradient based metric that quantifies the level of
% burriness. The smaller the value the sharper the image is.
%
% simind   = A metric for the structural similarity index. The higher the
% value the sharper the image is. Comparing two of the same images give you 
% a value of 1.
%
% just_noticeable = This metric calculates the just noticeable blurr based 
% on the edge width and a treshhold of perception. The smaller the value 
% the sharper the image is.
%
%=========================================================================%

function [gradient, simind, just_noticeable] = sharpness_metrics(u, ref)

[gradient]=gradient_sharpness_estimate(u);
[just_noticeable] = JNBM_compute(u);

if size(u) == size(ref);
    [simind ssim_map] = ssim_index(u, ref);
else
    simind = 0;
    disp('No reference image found')
end





