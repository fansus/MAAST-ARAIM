function [ y ] = modnormcdf( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

y = normcdf(x);
y(y>.5)=1;


end

