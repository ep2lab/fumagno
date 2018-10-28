%{
This is an example user simrc, which partially overwrittes the settings
of the default simrc file when running the preprocessor

INPUT AND OUTPUT
* data: structure (optional as input; can be initially empty).
  Parameters will be added to it by running this function. 
%} 
function data = akiles2d_simrc(data)

%% General
data.akiles2d.simdir = tempname; % directory where simulation files will be saved

%% Initial guess
data.guess = load(fullfile('fixtures/akiles2d_guessfile.mat')); % path to file where initial guesses for h,phiz,?? are stored. If empty, default preprocessor values will be used. 
