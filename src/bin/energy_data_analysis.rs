#![feature(result_option_inspect)]
#![feature(drain_filter)]
#![feature(let_chains)]

use std::error::Error;

use lattice_qft::{
    computation::ComputationSummary,
    export::{CsvData, EnergyData, StringTension},
    kahan::WelfordsAlgorithm64,
};
use nalgebra::{DMatrix, DVector};
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
use rayon::prelude::*;

// Parameters
const BIN_SIZES: [usize; 6] = [1, 2, 4, 8, 16, 32];
const BOOTSTRAPPING_ENSEMBLES: u64 = 100;
const SEEDED: bool = true;

fn main() {
    binning_resampling(
        lattice_qft::RESULTS_PATH,
        lattice_qft::COMP_ENERGY_STATS_PATH_INCOMPLETE,
        Some(lattice_qft::ENERGY_RESULTS_PATH),
    );

    string_tension_fit(
        lattice_qft::ENERGY_RESULTS_PATH,
        lattice_qft::COMP_ENERGY_STATS_PATH_INCOMPLETE,
        lattice_qft::STRING_TENSION_RESULTS_PATH,
    );
}

fn string_tension_fit(
    energy_results_path: &str,
    comp_energy_stats_path_incomplete: &str,
    string_tension_results_path: &str,
) {
    let results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(energy_results_path, true) {
            Ok(res) => res,
            Err(err) => {
                eprintln!("Reading energy results: {}", err);
                return;
            }
        };

    let mut string_tension_results: Vec<StringTension> = Vec::new();

    let mut results: Vec<ComputationSummary> = results
        .into_iter()
        .filter(|summary| summary.comptype.ends_with("Wilson Metropolis Simulations"))
        .rev()
        .collect();

    while let Some(key_summary) = results.pop() {
        let local_results: Vec<(u64, EnergyData)> = results
            .drain_filter(|summary| summary.temp == key_summary.temp && summary.t == key_summary.t)
            .filter_map(|summary| {
                summary
                    .comptype
                    .split_ascii_whitespace()
                    .skip(1)
                    .next()
                    .and_then(|width_x_height| width_x_height.split_once("x"))
                    .ok_or(Into::<Box<dyn Error>>::into("No width found"))
                    .and_then(|(first, _)| first.parse::<u64>().map_err(|err| err.into()))
                    .and_then(|width| {
                        let path: &str = &(comp_energy_stats_path_incomplete.to_owned()
                            + &"energy_"
                            + &summary.index.to_string()
                            + &"_stats.csv");
                        EnergyData::fetch_csv_data(path, true).map(|energy| {
                            energy
                                .into_iter()
                                .map(|energy| (width, energy))
                                .collect::<Vec<(u64, EnergyData)>>()
                        })
                    })
                    .ok()
            })
            .flatten()
            .collect();

        let mut local_string_tension_results: Vec<StringTension> = BIN_SIZES
            .into_par_iter()
            .map(|bin_size| {
                let (x_values, y_values): (Vec<f64>, Vec<f64>) = local_results
                    .iter()
                    .filter(|(_, energy)| energy.bin_size == bin_size as u64)
                    .map(|(width, energy)| (*width as f64, energy.energy))
                    .unzip();

                let string_tension: f64 = match regression(x_values, y_values) {
                    Ok(val) => val,
                    Err(err) => {
                        eprintln!("String tension calc: {}", err);
                        0.0
                    }
                };

                StringTension::new(
                    key_summary.index,
                    key_summary.temp,
                    key_summary.t,
                    bin_size as u64,
                    string_tension,
                )
            })
            .collect();

        string_tension_results.append(&mut local_string_tension_results);
    }

    if let Err(err) = string_tension_results.overwrite_csv(string_tension_results_path) {
        eprintln!("Writing string tension results: {}", err);
    }
}

fn regression(x_values: Vec<f64>, y_values: Vec<f64>) -> Result<f64, Box<dyn Error>> {
    let x = DVector::from_vec(x_values);
    let y = DVector::from_vec(y_values);

    let x_1 = x.clone();
    let x_2 = x.map(|x| x.powi(-2));
    let x_3 = x.map(|_| 1.0);

    let mat = DMatrix::from_columns(&[x_1, x_2, x_3]);

    let res = (mat.transpose() * mat.clone())
        .try_inverse()
        .ok_or::<Box<dyn Error>>("Inverse failed".into())?
        * mat.transpose()
        * y;

    assert_eq!(res.len(), 3);

    Ok(res[0])
}

fn binning_resampling(
    results_path: &str,
    comp_energy_stats_path_incomplete: &str,
    comp_results_path: Option<&str>,
) {
    // Reading the results data
    let results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(results_path, true) {
            Ok(res) => res,
            Err(err) => {
                eprintln!("Reading results for binned resampling: {}", err);
                return;
            }
        };

    let _old_comp_results: Vec<ComputationSummary> = comp_results_path
        .and_then(|path| {
            ComputationSummary::fetch_csv_data(path, true)
                .inspect_err(|err| eprintln!("Reading results for binned resampling: {}", err))
                .ok()
        })
        .unwrap_or(Vec::new());

    let mut results: Vec<ComputationSummary> = results
        .into_iter()
        .filter(|summary| summary.energy_data)
        .rev()
        .collect();

    let mut comp_results: Vec<ComputationSummary> = Vec::new();

    while let Some(mut key_summary) = results.pop() {
        let mut local_results: Vec<ComputationSummary> = results
            .drain_filter(|summary| {
                comp_results_path.is_some()
                    && summary.temp == key_summary.temp
                    && summary.t == key_summary.t
                    && summary.comptype == key_summary.comptype
            })
            .collect();

        local_results.reverse();
        local_results.push(key_summary.clone());
        local_results.reverse();

        let compund_data_amount: usize = local_results.len();

        if true
        /* !old_comp_results
        .par_iter()
        .filter(|summary| summary.temp == key_summary.temp && summary.t == key_summary.t)
        .any(|summary| {
            // Get the number of simulations and compare to the now red data
            summary
                .comptype
                .split_ascii_whitespace()
                .next()
                .and_then(|first_number| first_number.parse::<usize>().ok())
                .is_some_and(|first_number| first_number == compund_data_amount)
        })*/
        {
            // Fetching the correlation functions
            let local_energy_data: Vec<Vec<f64>> = local_results
                .into_par_iter()
                .filter_map(|summary| {
                    let path: &str = &(lattice_qft::ENERGY_DATA_PATH_INCOMPLETE.to_owned()
                        + &summary.index.to_string()
                        + &".csv");
                    let energy_data: Result<Vec<f64>, Box<dyn Error>> =
                        f64::fetch_csv_data(path, false);
                    energy_data
                        .inspect_err(|err| eprintln!("Fetching energy data: {}", err))
                        .ok()
                })
                .collect();

            let bin_observations: Vec<EnergyData> = BIN_SIZES
                .into_par_iter()
                .map(|bin_size| {
                    // Binning the correlation functions
                    let binned_local_energy_data: Vec<(Vec<&[f64]>, usize)> = local_energy_data
                        .iter()
                        .map(|single_sim_energy_data| {
                            let single_sim_binned_energy_data: Vec<&[f64]> =
                                single_sim_energy_data.chunks_exact(bin_size).collect();

                            // Quick check to see if we catch all correlation functions
                            let bin_samples_amount: usize = single_sim_energy_data.len() / bin_size;
                            assert_eq!(single_sim_binned_energy_data.len(), bin_samples_amount);

                            (single_sim_binned_energy_data, bin_samples_amount)
                        })
                        .collect();

                    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

                    let bootstrapping_ensembles: Vec<f64> = (0..BOOTSTRAPPING_ENSEMBLES)
                        .into_par_iter()
                        .map(|seed| {
                            let mut rng: StdRng = if SEEDED {
                                StdRng::seed_from_u64(seed)
                            } else {
                                StdRng::from_entropy()
                            };

                            // Choosing the simulations
                            let binned_local_resampled_energy_data: Vec<&(Vec<&[f64]>, usize)> =
                                binned_local_energy_data
                                    .iter()
                                    .filter_map(|_| binned_local_energy_data.choose(&mut rng))
                                    .collect();

                            let single_ensemble_energy_data_stats: WelfordsAlgorithm64 =
                                binned_local_resampled_energy_data
                                    .iter()
                                    .map(|(binned_energy_values, bin_samples_amount)| {
                                        let single_sim_resampled_corr_fns: Vec<&f64> = (0
                                            ..*bin_samples_amount)
                                            .into_iter()
                                            .filter_map(|_| binned_energy_values.choose(&mut rng))
                                            .cloned()
                                            .flatten()
                                            .collect();
                                        single_sim_resampled_corr_fns
                                    })
                                    .flatten()
                                    .cloned()
                                    .fold(WelfordsAlgorithm64::new(), |mut acc, value| {
                                        acc.update(value);
                                        acc
                                    });

                            single_ensemble_energy_data_stats.get_mean()
                        })
                        .collect();

                    let mean_stats: WelfordsAlgorithm64 = bootstrapping_ensembles.into_iter().fold(
                        WelfordsAlgorithm64::new(),
                        |mut acc, value| {
                            acc.update(value);
                            acc
                        },
                    );

                    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                    EnergyData::new(
                        key_summary.index,
                        bin_size as u64,
                        mean_stats.get_mean(),
                        mean_stats.get_sd_unbiased(),
                    )
                })
                .collect();

            // Writing the energy data
            let path: &str = &(comp_energy_stats_path_incomplete.to_owned()
                + &"energy_"
                + &key_summary.index.to_string()
                + &"_stats.csv");
            if let Err(err) = bin_observations.overwrite_csv(path) {
                eprintln!("Writing energy stats: {}", err);
            }
        }

        key_summary.comptype = compund_data_amount.to_string() + &" " + &key_summary.comptype + "s";

        comp_results.push(key_summary);
    }

    if let Some(path) = comp_results_path {
        if let Err(err) = comp_results.overwrite_csv(path) {
            eprintln!("Writing compounded energy results: {}", err);
        }
    }
}
