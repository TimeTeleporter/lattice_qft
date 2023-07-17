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

/// The main simulation data
pub const RESULTS_PATH: &str = "./data/results.csv";
pub const ACTION_DATA_PATH: &str = "./data/action.csv";
pub const CORR_FN_PATH_INCOMPLETE: &str = "./data/correlation_data/";

pub const RESULTS_SECOND_MOMENT_PATH: &str = "./data/results_second_moment.csv";
pub const RESULTS_COMP_PATH: &str = "./data/results_comp.csv";
pub const RESULTS_COMP_FIT_PATH: &str = "./data/results_comp_fit.csv";
pub const RESULTS_COMP_SECOND_MOMENT_PATH: &str = "./data/results_comp_second_moment.csv";
pub const COMP_CORR_FN_PATH_INCOMPLETE: &str = "./data/comp_correlation_data/";
pub const DIFFERENCE_PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/difference_";
pub const ENERGY_PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/energy_";

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
