% Using MATLAB Engine API for Python
%https://towardsdatascience.com/matlab-function-in-python-739c473c8176
matlabroot = 'C:\Program Files\MATLAB\R2021b';

cd (fullfile(matlabroot,'extern','engines','python'))
system('python setup.py install')