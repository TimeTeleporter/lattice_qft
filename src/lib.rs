#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(adt_const_params)]
#![feature(split_array)]
#![feature(is_some_with)]

pub const LONG_TEMP_ARY: [f64; 81] = [
    0.000001,
    0.0000012589,
    0.0000015848,
    0.0000019952,
    0.0000025118,
    0.0000031622,
    0.0000039810,
    0.0000050118,
    0.0000063095,
    0.0000079432,
    0.00001,
    0.000012589,
    0.000015848,
    0.000019952,
    0.000025118,
    0.000031622,
    0.000039810,
    0.000050118,
    0.000063095,
    0.000079432,
    0.0001,
    0.00012589,
    0.00015848,
    0.00019952,
    0.00025118,
    0.00031622,
    0.00039810,
    0.00050118,
    0.00063095,
    0.00079432,
    0.001,
    0.0012589,
    0.0015848,
    0.0019952,
    0.0025118,
    0.0031622,
    0.0039810,
    0.0050118,
    0.0063095,
    0.0079432,
    0.01,
    0.012589,
    0.015848,
    0.019952,
    0.025118,
    0.031622,
    0.039810,
    0.050118,
    0.063095,
    0.079432,
    0.1,
    0.12589,
    0.15848,
    0.19952,
    0.25118,
    0.31622,
    0.39810,
    0.50118,
    0.63095,
    0.79432,
    1.0,
    1.2589,
    1.5848,
    1.9952,
    2.5118,
    3.1622,
    3.9810,
    5.0118,
    6.3095,
    7.9432,
    10.0,
    12.589,
    15.848,
    19.952,
    25.118,
    31.622,
    39.810,
    50.118,
    63.095,
    79.432,
    100.0,
];

pub const TEMP_ARY: [f64; 41] = [
    0.001, 0.0012589, 0.0015848, 0.0019952, 0.0025118, 0.0031622, 0.0039810, 0.0050118, 0.0063095,
    0.0079432, 0.01, 0.012589, 0.015848, 0.019952, 0.025118, 0.031622, 0.039810, 0.050118,
    0.063095, 0.079432, 0.1, 0.12589, 0.15848, 0.19952, 0.25118, 0.31622, 0.39810, 0.50118,
    0.63095, 0.79432, 1.0, 1.2589, 1.5848, 1.9952, 2.5118, 3.1622, 3.9810, 5.0118, 6.3095, 7.9432,
    10.0,
];

pub const REL_TEMP_ARY: [f64; 31] = [
    0.01, 0.012589, 0.015848, 0.019952, 0.025118, 0.031622, 0.039810, 0.050118, 0.063095, 0.079432,
    0.1, 0.12589, 0.15848, 0.19952, 0.25118, 0.31622, 0.39810, 0.50118, 0.63095, 0.79432, 1.0,
    1.2589, 1.5848, 1.9952, 2.5118, 3.1622, 3.9810, 5.0118, 6.3095, 7.9432, 10.0,
];

pub const CLUSTER_RESULTS_PATH: &str = "./data/cluster_results.csv";
pub const CLUSTER_BINNING_PATH: &str = "./data/cluster_binning.csv";
pub const METRPLS_RESULTS_PATH: &str = "./data/metropolis_results.csv";
pub const METRPLS_BINNING_PATH: &str = "./data/metropolis_binning.csv";

pub mod algorithm;
pub mod computation;
pub mod error;
pub mod export;
pub mod field;
pub mod heightfield;
pub mod lattice;
pub mod observable;

extern crate csv;
extern crate nalgebra;
extern crate rand;
extern crate rayon;
extern crate serde;
extern crate varpro;

// A method to pause a sim until enter
use std::process::Command;

pub fn pause() {
    let _ = Command::new("cmd.exe").arg("/c").arg("pause").status();
}
