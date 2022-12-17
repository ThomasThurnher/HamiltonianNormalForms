function install
% run this file first to install all external packages  and
% also the main software appropriately


addpath(fullfile('src/TransformingDS'))

addpath(fullfile('src'))

addpath(fullfile('src/misc'))

maindirSSM = strcat(fileparts(mfilename('fullpath')),'/ext/SSMTool/');

extdir = 'ext';

addpath(fullfile(maindirSSM, extdir,'combinator'));

addpath(fullfile(maindirSSM, extdir, 'tensor_toolbox'));

addpath(fullfile(maindirSSM, extdir, 'Wrappers'));

addpath(fullfile(maindirSSM, 'src'));

addpath(fullfile(maindirSSM, 'src','misc'));
end


