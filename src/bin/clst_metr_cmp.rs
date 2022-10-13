#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use lattice_qft::{
    cluster::cluster_simulation3d,
    export::{BinData, CsvData, SimResult},
    lattice::Lattice3d,
    metropolis::metropolis_simulation3d,
};

const MAX_X: usize = 2;
const MAX_Y: usize = 2;
const MAX_T: usize = 2;

const CLUSTER_RESULTS_PATH: &str = "data/clst_metr_cmp/cluster_results.csv";
const CLUSTER_BINS_PATH: &str = "data/clst_metr_cmp/cluster_bins.csv";
const METROPOLIS_RESULTS_PATH: &str = "data/clst_metr_cmp/metropolis_results.csv";
const METROPOLIS_BINS_PATH: &str = "data/clst_metr_cmp/metropolis_bins.csv";

const BURNIN: usize = 100_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 100_000_000;
const NUM_BINS: usize = 21;

fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    let results: Vec<(
        (SimResult, [BinData; NUM_BINS]),
        (SimResult, [BinData; NUM_BINS]),
    )> = lattice_qft::LONG_TEMP_ARY
        .par_iter()
        //.filter(|&&temp| temp == 0.001)
        .map(|&temp| {
            let (clst_result, mut clst_data) =
                cluster_simulation3d(&lattice, temp, BURNIN, ITERATIONS);

            let mut binning_step: u32 = 0;
            let clst_bins: [BinData; NUM_BINS] = [(); NUM_BINS].map(|_| {
                let bin: BinData = clst_data.get_bin_data(temp, binning_step);
                clst_data = clst_data.binning();
                binning_step += 1;
                bin
            });

            let (metr_result, mut metr_data) =
                metropolis_simulation3d(&lattice, temp, BURNIN, ITERATIONS);

            let mut binning_step: u32 = 0;
            let metr_bins: [BinData; NUM_BINS] = [(); NUM_BINS].map(|_| {
                let bin: BinData = metr_data.get_bin_data(temp, binning_step);
                metr_data = metr_data.binning();
                binning_step += 1;
                bin
            });

            ((clst_result, clst_bins), (metr_result, metr_bins))
        })
        .collect();

    for ((clst_result, clst_bins), (metr_result, metr_bins)) in results {
        if let Err(err) = clst_result.read_write_csv(CLUSTER_RESULTS_PATH) {
            eprint!("{}", err);
        }
        if let Err(err) = metr_result.read_write_csv(METROPOLIS_RESULTS_PATH) {
            eprint!("{}", err);
        }
        for bin in clst_bins {
            if let Err(err) = bin.read_write_csv(CLUSTER_BINS_PATH) {
                eprint!("{}", err);
            }
        }
        for bin in metr_bins {
            if let Err(err) = bin.read_write_csv(METROPOLIS_BINS_PATH) {
                eprint!("{}", err);
            }
        }
    }
}
