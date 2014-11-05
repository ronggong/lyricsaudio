function [ success_value ] = batchFolder_mat2arff( input_path, output_name, relation )
%   BATCHFOLDER_MAT2ARFF function is for batch processing all .mat files
%   within a folder to convert them into an .arff file for use with WEKA
%   machine learning software
%
%   (c) Jens Bandener 2011 [Jens.Bandener@ruhr-uni-bochum.de]
%                                       script version v0.1 (05-dec-2011)
%
%   .mat file conditions:
%   ------------------------------------------------------------------
%   requires following variables stored inside the .mat file.
%   [string]        AudioType
%   [double matrix] featureMatrix
%   AudioType containing the attribute for the features in this .mat file.
%   featureMatrix containing the feature data as a 2-dim or 3-dim double
%   precision matrix
%   All other variables stored in the .mat file will not be used in this
%   version.
%
%

current_path = pwd;
cd(input_path);
file_list = dir('*.mat');
num_files = size(file_list,1);
% scan files for differen genres and create genres vector
genre_string = {};
for x=1:num_files
    load(file_list(x).name);
    genre_string = [genre_string;AudioType];
end

genre_unique = unique(genre_string);

genre_string = genre_unique{1};
for ii = 2: length(genre_unique)
    genre_string = [genre_string ',' genre_unique{ii}];
end


% iterate over all feature files and call mat2arff function to create the
% arff file
cd(current_path);
for x=1:num_files
    %mat2arff([input_path '/' file_list(x).name],output_path,'featureMatrix',relation,genre_string);
    mat2arff([input_path '/' file_list(x).name], [input_path '/' output_name],relation,genre_string);
end

cd(current_path);
success_value = 1;
end

