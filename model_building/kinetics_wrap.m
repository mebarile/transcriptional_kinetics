% script for the compilation of *_syms files with AMICI
% run it once per model, changinng lines 12 and 16

clear; close all; clc

%% COMPILATION

[exdir,~,~]=fileparts(which('kinetics_wrap.m'))

% compile the model
tic;
amiwrap('kinetics_leave_one','kinetics_model1_syms',exdir)
t_wrap = toc

% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/kinetics_model1']))