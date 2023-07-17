use std::error::Error;

use lattice_qft::export::{CsvData, FitResult};
use nalgebra::DVector;
use varpro::{
    prelude::*,
    solvers::levmar::{LevMarProblemBuilder, LevMarSolver},
};

const MAX_T: usize = 54;
const CORRELATION_LENGHT: f64 = 55.0;
const M: f64 = 1.0 / CORRELATION_LENGHT;
const A: f64 = 300.0;
const B: f64 = -300.0;

fn main() {
    let ary: [f64; MAX_T] = core::array::from_fn(|i| i as f64);
    let corr_fn: Vec<f64> = ary
        .into_iter()
        .map(|x| f64::cosh(M * (x - (MAX_T as f64) / 2.0)) * A + B)
        .collect();

    let path: &str = &(lattice_qft::CORR_FN_PATH_INCOMPLETE.to_owned() + &"dummy_corr_fn.csv");

    if let Err(err) = corr_fn.clone().overwrite_csv(path) {
        eprint!("{}", err);
    };

    let corr_fn: Vec<f64> = match f64::fetch_csv_data(path, false) {
        Ok(corr_fn) => corr_fn,
        Err(err) => {
            eprintln!("{}", err);
            return;
        }
    };

    const P1: f64 = 2.0 * std::f64::consts::PI / (MAX_T as f64);
    let inverse_correlation_length: f64 =
        lattice_qft::calculate_correlation_length(&corr_fn, P1, P1 * 2.0).0;
    println!(
        "Second moment: m12 = {}, corr12 = {}",
        inverse_correlation_length,
        1.0 / inverse_correlation_length
    );
    let inverse_correlation_length: f64 =
        lattice_qft::calculate_correlation_length(&corr_fn, 2.0 * P1, P1 * 3.0).0;
    println!(
        "Second moment: m23 = {}, corr23 = {}",
        inverse_correlation_length,
        1.0 / inverse_correlation_length
    );
    let inverse_correlation_length: f64 =
        lattice_qft::calculate_correlation_length(&corr_fn, P1, P1 * 3.0).0;
    println!(
        "Second moment: m13 = {}, corr13 = {}",
        inverse_correlation_length,
        1.0 / inverse_correlation_length
    );

    if let Err(err) = nonlin_regression(corr_fn).and_then(|fit| {
        println!(
            "The fit yielded: m = {}, a = {}, b = {}, corr_fit = {}",
            fit.m, fit.a, fit.b, fit.corr
        );
        Ok(())
    }) {
        eprintln!("{}", err);
    };
}

fn nonlin_regression(corr_fn: Vec<f64>) -> Result<FitResult, Box<dyn Error>> {
    match std::panic::catch_unwind(move || {
        let y_values = corr_fn;
        let n: usize = y_values.len();

        let n: f64 = n as f64;

        // Edit the data
        let Some(x_max): Option<f64> = y_values.clone().into_iter().reduce(f64::max) else {
            panic!("Unable to find the maximum");
        };

        let y_values: Vec<f64> = y_values.into_iter().map(|x| x_max - x).collect();

        let x_values: Vec<f64> = y_values.iter().enumerate().map(|(x, _)| x as f64).collect();

        let x = DVector::from_vec(x_values);
        let y = DVector::from_vec(y_values);

        // The function to fit
        fn nonlin_fn(x: &DVector<f64>, m: f64, n: f64) -> DVector<f64> {
            x.map(|x| (m * (x - n / 2.0)).cosh())
        }

        let nonlin_fn_n_given = move |x: &DVector<f64>, m: f64| nonlin_fn(x, m, n);

        // Derivative with respect to m
        fn nonlin_fn_dm(x: &DVector<f64>, m: f64, n: f64) -> DVector<f64> {
            x.map(|x| (m * (x - n / 2.0)).sinh() * (x - n / 2.0))
        }

        let nonlin_fn_dm_n_given = move |x: &DVector<f64>, m: f64| nonlin_fn_dm(x, m, n);

        /* Derivative with respect to n
        fn nonlin_fn_dn(x: &DVector<f64>, m: f64, n: f64) -> DVector<f64> {
            x.map(|x| -1.0 * (m * (x - n / 2.0)).sinh() * (m / 2.0))
        }*/

        let model = SeparableModelBuilder::<f64>::new(&["m"])
            .function(&["m"], nonlin_fn_n_given)
            .partial_deriv("m", nonlin_fn_dm_n_given)
            .invariant_function(|x| DVector::from_element(x.len(), 1.0))
            .independent_variable(x)
            .initial_parameters(vec![0.1])
            .build()
            .unwrap();

        let problem = LevMarProblemBuilder::new(model)
            .observations(y)
            .build()
            .unwrap();

        // This thing panics, we need to catch it.
        let (solved_problem, report) = LevMarSolver::new()
            .with_xtol(f64::EPSILON)
            .minimize(problem);
        if !report.termination.was_successful() {
            panic!("termination was not successful");
        }

        let alpha = solved_problem.params();
        let coeff = solved_problem
            .linear_coefficients()
            .ok_or("Unable to compute the coefficients")
            .unwrap();

        // FitResult::new(index, alpha[0], alpha[1], coeff[0], coeff[1])
        FitResult::new(0, alpha[0], n, coeff[0], coeff[1])
    }) {
        Ok(fit) => Ok(fit),
        Err(err) => {
            let err_msg = match (err.downcast_ref(), err.downcast_ref::<String>()) {
                (Some(&s), _) => s,
                (_, Some(s)) => &**s,
                _ => "<No panic message>",
            };
            Err(format!("nonlin regression paniced: {}", err_msg).into())
        }
    }
}
