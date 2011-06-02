function [ x,u ] = flat_shelf_ice(A,rhoi,rhoo,g,x0,x1,n_points, u0, h0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = x0 : (x1-x0)/n_points : x1;

freeboard = (1- rhoi/rhoo)*h0;

u = u0 + (x-x0)*A*(1/4*rhoi*g*freeboard)^3;


end

