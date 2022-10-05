#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(adt_const_params)]
#![feature(split_array)]

pub const TEMP_ARY: [f64; 41] = [
    0.001, 0.0012589, 0.0015848, 0.0019952, 0.0025118, 0.0031622, 0.0039810, 0.0050118, 0.0063095,
    0.0079432, 0.01, 0.012589, 0.015848, 0.019952, 0.025118, 0.031622, 0.039810, 0.050118,
    0.063095, 0.079432, 0.1, 0.12589, 0.15848, 0.19952, 0.25118, 0.31622, 0.39810, 0.50118,
    0.63095, 0.79432, 1.0, 1.2589, 1.5848, 1.9952, 2.5118, 3.1622, 3.9810, 5.0118, 6.3095, 7.9432,
    10.0,
];

pub mod action;
pub mod cluster;
pub mod export;
pub mod field;
pub mod lattice;
pub mod metropolis;
