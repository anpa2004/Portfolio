close all; clear; clc;

const = getConst();

tspan = [0 5];
 
initial_conditions = [const.r0(1);const.v0;const.r0(2);const.v0;const.m0]