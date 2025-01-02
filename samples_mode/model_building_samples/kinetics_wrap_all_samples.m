% script for the compilation of *_syms files with AMICI
clear; close all; clc

%% COMPILATION

[exdir,~,~]=fileparts(which('kinetics_wrap_all_samples.m'));

% compile the model
tic;
amiwrap('kinetics_model8_all_samples','kinetics_model8_all_samples_syms',exdir)
t_wrap = toc

% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/kinetics_model8_all_samples']))
