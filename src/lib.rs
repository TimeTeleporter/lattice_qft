#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(adt_const_params)]
#![feature(split_array)]
#![feature(let_chains)]
#![feature(result_option_inspect)]

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

pub const SMALL_TEMP_ARY: [f64; 21] = [
    0.02, 0.025118, 0.031622, 0.039810, 0.050118, 0.063095, 0.079432, 0.1, 0.12589, 0.15848,
    0.19952, 0.25118, 0.31622, 0.39810, 0.50118, 0.63095, 0.79432, 1.0, 1.2589, 1.5848, 2.0,
];

pub const INVESTIGATE_ARY: [f64; 10] = [0.2, 0.24, 0.28, 0.33, 0.39, 0.46, 0.55, 0.64, 0.76, 0.9];

pub const INVESTIGATE_ARY2: [f64; 10] = [0.2, 0.22, 0.25, 0.27, 0.3, 0.33, 0.37, 0.41, 0.45, 0.5];

pub const INVESTIGATE_ARY3: [f64; 30] = [
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35,
    0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
];

pub const RESULTS_PATH: &str = "./data/results.csv";
pub const RESULTS_FIT_PATH: &str = "./data/results_fit.csv";
pub const RESULTS_CORR_PATH: &str = "./data/results_corr.csv";
pub const RESULTS_COMP_PATH: &str = "./data/results_comp.csv";
pub const PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/";

pub mod algorithm;
pub mod computation;
pub mod export;
pub mod fields;
pub mod kahan;
pub mod lattice;
pub mod outputdata;

extern crate csv;
extern crate rand;
extern crate rayon;
extern crate serde;

// A method to pause a sim until enter
use std::process::Command;

pub fn pause() {
    let _ = Command::new("cmd.exe").arg("/c").arg("pause").status();
}

#[test]
fn get_log_ary() {
    const TESTING: bool = true;

    // Define the array boundries
    const LOWER: f64 = 0.2;
    const UPPER: f64 = 0.5;
    const STEPS: usize = 10;

    // Convert them to log scale
    let lower_log: f64 = LOWER.log10();
    let upper_log: f64 = UPPER.log10();

    if TESTING {
        dbg!(lower_log);
        dbg!(upper_log);
    }
    let diff: f64 = upper_log - lower_log;
    let step_size: f64 = diff / ((STEPS - 1) as f64);

    let ary: [f64; STEPS] = core::array::from_fn(|i| i)
        .map(|step| lower_log + step_size * (step as f64))
        .map(|x| f64::powf(10.0, x))
        .map(|x| (x * 100.0).round() / 100.0);

    if TESTING {
        dbg!(ary);
    };
}

#[test]
fn get_ary() {
    const TESTING: bool = true;

    // Define the array boundries
    const LOWER: f64 = 0.2;
    const UPPER: f64 = 0.3;
    const STEPS: usize = 10;

    if TESTING {
        dbg!(LOWER);
        dbg!(UPPER);
    }
    let diff: f64 = UPPER - LOWER;
    let step_size: f64 = diff / ((STEPS) as f64);

    let ary: [f64; STEPS] = core::array::from_fn(|i| i)
        .map(|step| LOWER + step_size * (step as f64))
        .map(|x| f64::powf(10.0, x))
        .map(|x| (x * 100.0).round() / 100.0);

    if TESTING {
        dbg!(ary);
    };
}
