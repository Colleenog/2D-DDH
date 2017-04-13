function [ output ] = normalize( input )
%NORMALIZE Summary of this function goes here
%   Detailed explanation goes here

minval = min(min(min(min(min(input)))));

input = input - minval;

maxval = max(max(max(max(max(input)))));

output = input./maxval;


end

