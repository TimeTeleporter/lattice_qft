use lattice_qft::computation::ComputationSummary;
use lattice_qft::export::{get_correlation_fn, CsvData, FitResult};
use nalgebra::DVector;
use std::error::Error;
use varpro::prelude::*;
use varpro::solvers::levmar::{LevMarProblemBuilder, LevMarSolver};

const SYMMETRIZE: bool = true;

fn main() {
    // Read the result data
    let results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(lattice_qft::RESULTS_PATH, true) {
            Ok(res) => res,
            Err(err) => {
                eprint!("{}", err);
                return;
            }
        };

    let fitted: Vec<FitResult> = results
        .into_iter()
        .filter(|res| res.correlation_data)
        .filter_map(|res| {
            nonlin_regression(res.index, lattice_qft::CORRELATION_PLOT_PATH_INCOMPLETE)
                .map_err(|err| {
                    eprint!("{}", err);
                })
                //.or::<Box<dyn Error>>(Ok(FitResult::new(index, f64::NAN, f64::NAN, f64::NAN, f64::NAN)))
                .ok()
        })
        .collect();

    if let Err(err) = fitted.overwrite_csv(lattice_qft::RESULTS_FIT_PATH) {
        eprint!("{}", err);
    }
}

fn nonlin_regression(index: usize, incomplete_path: &str) -> Result<FitResult, Box<dyn Error>> {
    match std::panic::catch_unwind(move || {
        let mut y_values = get_correlation_fn(index, incomplete_path).unwrap();
        let n: usize = y_values.len();

        if SYMMETRIZE {
            let mut new_y_values: Vec<f64> = Vec::new();
            new_y_values.push(y_values[0]);
            for i in 1..n {
                new_y_values.push((y_values[i] + y_values[n - i]) / 2.0)
            }
            assert_eq!(n, new_y_values.len());
            y_values = new_y_values;
        }

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
        FitResult::new(index, alpha[0], n, coeff[0], coeff[1])
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
