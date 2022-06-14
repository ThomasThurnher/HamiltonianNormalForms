function install
% run this file first to install all external packages  and
% also the main software appropriately




maindirSSM = strcat(fileparts(mfilename('fullpath')),'ext/SSMTool/');

extdir = 'ext';
cocoinstall = fullfile(maindirSSM, extdir , 'coco','startup.m');
run(cocoinstall);

addpath(fullfile(maindirSSM, extdir,'combinator'));

addpath(fullfile(maindirSSM, extdir, 'tensor_toolbox'));

addpath(fullfile(maindirSSM, extdir, 'Wrappers'));

addpath(fullfile(maindirSSM, 'src'));

addpath(fullfile(maindirSSM, 'src','misc'));
end


