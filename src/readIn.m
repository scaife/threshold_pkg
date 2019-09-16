function output = readIn(fName)
% Description: Imports rhessys data output. NOTE: .daily files need to have
% a .txt extension. 
%
% Input: 
%   1. fName: file name as string.
%
% Output:
%   1. output: output of all the data.
%
% Example: output = rhessysReadIn("control.txt"); 
%
% Author:   Charles Scaife
% Date:     Oct. 5, 2015

delimiter = ','; 

output = readtable(fName, 'Delimiter', delimiter);
