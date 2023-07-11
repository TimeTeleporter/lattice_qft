#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(adt_const_params)]
#![feature(split_array)]
#![feature(let_chains)]
#![feature(result_option_inspect)]
#![feature(array_zip)]

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

pub const INVESTIGATE_ARY3_10: [f64; 300] = [
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35,
    0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21,
    0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37,
    0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23,
    0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
    0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25,
    0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41,
    0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27,
    0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43,
    0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
    0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45,
    0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31,
    0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47,
    0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33,
    0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35,
    0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21,
    0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37,
    0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
];

pub const INVESTIGATE_ARY3_6: [f64; 180] = [
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35,
    0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21,
    0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37,
    0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23,
    0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
    0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25,
    0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41,
    0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27,
    0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43,
    0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
    0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45,
    0.46, 0.47, 0.48, 0.49,
];

pub const INVESTIGATE_ARY3_3: [f64; 90] = [
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35,
    0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21,
    0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37,
    0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.2, 0.21, 0.22, 0.23,
    0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
    0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
];

pub const INVESTIGATE_ARY4: [f64; 18] = [
    0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38,
    0.39, 0.40,
];

pub const INVESTIGATE_ARY5: [f64; 6] = [0.22, 0.23, 0.24, 0.25, 0.26, 0.27];

pub const INVESTIGATE_ARY5_3: [f64; 18] = [
    0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.22, 0.23, 0.24, 0.25,
    0.26, 0.27,
];

pub const RESULTS_PATH: &str = "./data/results.csv";
pub const RESULTS_FIT_PATH: &str = "./data/results_fit.csv";
pub const RESULTS_SECOND_MOMENT_PATH: &str = "./data/results_second_moment.csv";
pub const RESULTS_COMP_PATH: &str = "./data/results_comp.csv";
pub const RESULTS_COMP_FIT_PATH: &str = "./data/results_comp_fit.csv";
pub const RESULTS_COMP_SECOND_MOMENT_PATH: &str = "./data/results_comp_second_moment.csv";
pub const CORR_FN_PATH_INCOMPLETE: &str = "./data/plot_data/correlation_data/";
pub const COMP_CORR_FN_PATH_INCOMPLETE: &str = "./data/plot_data/comp_correlation_data/";
pub const DIFFERENCE_PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/difference_data/";
pub const ENERGY_PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/energy_data/";
pub const OTHER_PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/other_data/";

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
    let (g1_re, g1_im): (f64, f64) = discrete_fourier_transform(corr_fn, p1);
    let (g2_re, g2_im): (f64, f64) = discrete_fourier_transform(corr_fn, p2);

    let g1: Complex64 = Complex64::new(g1_re, g1_im);
    let g2: Complex64 = Complex64::new(g2_re, g2_im);

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
    let (g1_re, g1_im): (f64, f64) = discrete_fourier_transform(corr_fn, p1);
    let (g2_re, g2_im): (f64, f64) = discrete_fourier_transform(corr_fn, p2);
    let (g1_re_err, g1_im_err): (f64, f64) = discrete_fourier_transform(corr_fn_err, p1);
    let (g2_re_err, g2_im_err): (f64, f64) = discrete_fourier_transform(corr_fn_err, p2);

    let (g1_re_err, g1_im_err): (f64, f64) = (g1_re_err.abs(), g1_im_err.abs());
    let (g2_re_err, g2_im_err): (f64, f64) = (g2_re_err.abs(), g2_im_err.abs());

    let g1: Complex64 = Complex64::new(g1_re, g1_im);
    let g2: Complex64 = Complex64::new(g2_re, g2_im);
    let g1_err: Complex64 = Complex64::new(g1_re_err, g1_im_err);
    let g2_err: Complex64 = Complex64::new(g2_re_err, g2_im_err);

    let cos1: f64 = p1.cos();
    let cos2: f64 = p2.cos();

    let z: Complex64 = (g1 * cos1 - g2 * cos2).fdiv(g1 - g2);
    assert!(!z.is_nan());

    // Error propagation for z:
    // Define x as the numerator of cosh
    let Complex64 { re: x_re, im: x_im } = g1 * cos1 - g2 * cos2;
    let Complex64 {
        re: x_re_err,
        im: x_im_err,
    } = g1_err * cos1.abs() + g2_err * cos2.abs();
    // and y as the denominator
    let Complex64 { re: y_re, im: y_im } = g1 - g2;
    let Complex64 {
        re: y_re_err,
        im: y_im_err,
    } = g1_err + g2_err;
    let y_squared: f64 = y_re * y_re + y_im * y_im;

    // Because the quotient is a holomorphic function, we can use the identities
    let dz_re_dx_re: f64 = y_re / y_squared;
    let dz_im_dx_im: f64 = dz_re_dx_re;
    let dz_re_dx_im: f64 = y_im / y_squared;
    let dz_im_dx_re: f64 = -dz_re_dx_im;

    let dz_re_dy_re: f64 =
        (x_re * (y_im * y_im - y_re * y_re) - 2.0 * x_im * y_re * y_im) / (y_squared * y_squared);
    let dz_im_dy_im: f64 = dz_re_dy_re;
    let dz_re_dy_im: f64 =
        (x_im * (y_re * y_re - y_im * y_im) - 2.0 * x_re * y_re * y_im) / (y_squared * y_squared);
    let dz_im_dy_re: f64 = -dz_re_dy_im;

    let z_re_err: f64 = dz_re_dx_re.abs() * x_re_err
        + dz_re_dx_im.abs() * x_im_err
        + dz_re_dy_re.abs() * y_re_err
        + dz_re_dy_im.abs() * y_im_err;

    let z_im_err: f64 = dz_im_dx_re.abs() * x_re_err
        + dz_im_dx_im.abs() * x_im_err
        + dz_im_dy_re.abs() * y_re_err
        + dz_im_dy_im.abs() * y_im_err;

    const COMPLEX64_ONE: Complex64 = Complex64 { re: 1.0, im: 0.0 };
    const COMPLEX64_I: Complex64 = Complex64 { re: 0.0, im: 1.0 };
    let dacosh_dz_re: Complex64 =
        COMPLEX64_ONE.fdiv((z + COMPLEX64_ONE).sqrt() * (z - COMPLEX64_ONE).sqrt());
    let dacosh_dz_im: Complex64 =
        COMPLEX64_I.fdiv((z + COMPLEX64_ONE).sqrt() * (z - COMPLEX64_ONE).sqrt());

    let ma_re_err: f64 = dacosh_dz_re.re.abs() * z_re_err + dacosh_dz_im.re.abs() * z_im_err;
    let ma_im_err: f64 = dacosh_dz_re.im.abs() * z_re_err + dacosh_dz_im.im.abs() * z_im_err;

    let ma: Complex64 = z.acosh();

    ((ma.re, ma.im), (ma_re_err, ma_im_err))
}

fn discrete_fourier_transform(corr_fn: &Vec<f64>, momentum: f64) -> (f64, f64) {
    let real: f64 = corr_fn
        .iter()
        .enumerate()
        .map(|(x, g_x)| g_x * f64::cos(momentum * (x as f64)))
        .sum();
    let imaginary: f64 = corr_fn
        .iter()
        .enumerate()
        .map(|(x, g_x)| g_x * f64::sin(momentum * (x as f64)))
        .sum();
    (real, imaginary)
}
