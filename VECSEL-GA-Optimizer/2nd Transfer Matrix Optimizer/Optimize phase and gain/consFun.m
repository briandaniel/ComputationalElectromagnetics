function [ cons, consEq ] = consFun( x, grateSize )
%CONSFUN Summary of this function goes here
%   Detailed explanation goes here

cons = (sum(x) - grateSize).*1e7;
consEq = [];

end

