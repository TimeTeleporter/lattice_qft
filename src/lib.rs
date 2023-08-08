#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(adt_const_params)]
#![feature(split_array)]
#![feature(let_chains)]
#![feature(result_option_inspect)]
#![feature(array_zip)]

pub const INVESTIGATE_ARY: [f64; 10] = [0.2, 0.24, 0.28, 0.33, 0.39, 0.46, 0.55, 0.64, 0.76, 0.9];
pub const INVESTIGATE_ARY2: [f64; 10] = [0.2, 0.22, 0.25, 0.27, 0.3, 0.33, 0.37, 0.41, 0.45, 0.5];
pub const INVESTIGATE_ARY3: [f64; 30] = [
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35,
    0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
];
pub const INVESTIGATE_ARY5: [f64; 17] = [
    0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285,
    0.29, 0.295, 0.3,
];

pub const INVESTIGATE_ARY6: [f64; 19] = [
    0.21, 0.22, 0.225, 0.23, 0.2325, 0.235, 0.2375, 0.24, 0.2425, 0.245, 0.2475, 0.25, 0.2525,
    0.255, 0.2575, 0.26, 0.265, 0.27, 0.28,
];

pub const INVESTIGATE_ARY7_16: [f64; 21] = [
    0.25, 0.2525, 0.255, 0.2575, 0.26, 0.2625, 0.265, 0.2675, 0.27, 0.2725, 0.275, 0.2775, 0.28,
    0.2825, 0.285, 0.2875, 0.29, 0.2925, 0.295, 0.2975, 0.3,
];

pub const INVESTIGATE_ARY7_24: [f64; 21] = [
    0.24, 0.2425, 0.245, 0.2475, 0.25, 0.2525, 0.255, 0.2575, 0.26, 0.2625, 0.265, 0.2675, 0.27,
    0.2725, 0.275, 0.2775, 0.28, 0.2825, 0.285, 0.2875, 0.29,
];

pub const INVESTIGATE_ARY7_36: [f64; 21] = [
    0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23, 0.2325, 0.235, 0.2375, 0.24,
    0.2425, 0.245, 0.2475, 0.25, 0.2525, 0.255, 0.2575, 0.26,
];

pub const INVESTIGATE_ARY7_54: [f64; 17] = [
    0.20, 0.2025, 0.205, 0.2075, 0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23,
    0.2325, 0.235, 0.2375, 0.24,
];

pub const INVESTIGATE_ARY8: [f64; 20] = [
    0.1, 0.12, 0.14, 0.16, 0.18, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52,
    0.54, 0.56, 0.58, 0.6,
];

pub const INVESTIGATE_ARY9: [f64; 10] = [0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29];

pub const INVESTIGATE_ARY10_36: [f64; 17] = [
    0.21, 0.2125, 0.215, 0.2175, 0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23,
    0.2325, 0.235, 0.2375, 0.24,
];

pub const INVESTIGATE_ARY10_54: [f64; 22] = [
    0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23, 0.2325, 0.235, 0.21, 0.2125,
    0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23, 0.2325, 0.235,
];

/// The main simulation data
pub const RESULTS_PATH: &str = "./data/results.csv";

// Correlation function and statistics
pub const CORR_FN_PATH_INCOMPLETE: &str = "./data/correlation_data/";
pub const CORR_FN_STATS_PATH_INCOMPLETE: &str = "./data/correlation_stats/";

// Compounded correlation functions
pub const COMP_RESULTS_PATH: &str = "./data/results_comp.csv";
pub const COMP_CORR_FN_STATS_PATH_INCOMPLETE: &str = "./data/comp_correlation_stats/";

// Plots
pub const DIFFERENCE_PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/difference_";

// Action and energy
pub const ACTION_DATA_PATH_INCOMPLETE: &str = "./data/action_data/action_";
pub const ENERGY_DATA_PATH_INCOMPLETE: &str = "./data/energy_data/energy_";

pub const ENERGY_RESULTS_PATH: &str = "./data/comp_energy_results.csv";
pub const COMP_ENERGY_STATS_PATH_INCOMPLETE: &str = "./data/comp_energy_stats/";
pub const STRING_TENSION_RESULTS_PATH: &str = "./data/string_tension_results.csv";

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

use num_complex::Complex64;

pub fn pause() {
    let _ = Command::new("cmd.exe").arg("/c").arg("pause").status();
}

pub fn calculate_correlation_length(corr_fn: &Vec<f64>, p1: f64, p2: f64) -> (f64, f64) {
    let g1: Complex64 = discrete_fourier_transform(corr_fn, p1);
    let g2: Complex64 = discrete_fourier_transform(corr_fn, p2);

    let cos1: f64 = p1.cos();
    let cos2: f64 = p2.cos();

    let cosh: Complex64 = (g1 * cos1 - g2 * cos2).fdiv(g1 - g2);

    assert!(!cosh.is_nan());

    let ma: Complex64 = cosh.acosh();

    (ma.re, ma.im)
}

pub fn calculate_correlation_length_errors(
    corr_fn: &Vec<f64>,
    corr_fn_err: &Vec<f64>,
    p1: f64,
    p2: f64,
) -> ((f64, f64), (f64, f64)) {
    //let corr_fn_err: &Vec<f64> = &corr_fn_err.into_iter().map(|x| x * 0.1).collect();
    let g1: Complex64 = discrete_fourier_transform(corr_fn, p1);
    let g2: Complex64 = discrete_fourier_transform(corr_fn, p2);

    let cos1: f64 = p1.cos();
    let cos2: f64 = p2.cos();

    let x: Complex64 = g1 * cos1 - g2 * cos2;
    let y: Complex64 = g1 - g2;
    let z: Complex64 = x.fdiv(y);

    let ma: Complex64 = z.acosh();

    let ma_re_err = corr_fn_err
        .iter()
        .enumerate()
        .map(|(pos, u)| {
            let dg1_du: Complex64 = Complex64::exp(Complex64::i() * (pos as f64) * p1);
            let dg2_du: Complex64 = Complex64::exp(Complex64::i() * (pos as f64) * p2);

            const COMPLEX64_ONE: Complex64 = Complex64::new(1.0, 0.0);

            let df_du: Complex64 = ((g1 * dg2_du - dg1_du * g2) * (cos1 - cos2)).fdiv(
                y * y * Complex64::sqrt(z - COMPLEX64_ONE) * Complex64::sqrt(z + COMPLEX64_ONE),
            );
            (df_du * u).norm_sqr()
        })
        .sum::<f64>()
        .sqrt();

    ((ma.re, ma.im), (ma_re_err, 0.0))
}

fn discrete_fourier_transform(corr_fn: &Vec<f64>, momentum: f64) -> Complex64 {
    corr_fn
        .iter()
        .enumerate()
        .map(|(pos, g_x)| Complex64::exp(Complex64::i() * (pos as f64) * momentum) * g_x)
        .sum()
}
