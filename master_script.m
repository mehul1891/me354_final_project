%=========================================================================%
% MASTER SCRIPT          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Image retrieval, sharpness metric, test metric
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS

% 1. kernel_function     = 'ext' to use lpfilter.m
%                        = 'int' to use our developed Gauss convolution
%                        kernel
%
% 2. boundary_condition  = 'periodic' to use periodic boundary conditions
%                        = 'none' to not use any bundary condition
%
% 4. plot                = 'yes' or 'no'
% 5. draw_kernel         = 'yes' or 'no'


%=========================================================================%

%% Clear previous cache
clear all; close all; clc;

%% Images that this code is trained and tested with

% Test image 1: peppers.png
im1 = imread('peppers.png');

% Test image 2: lina.jpg
im2 = imread('lena.png');

% Test image 3: cameraman.tif
im3 = imread('cameraman.tif');

% Test image 4: ssphere.jpg
im4 = imread('ssphere.jpg');


