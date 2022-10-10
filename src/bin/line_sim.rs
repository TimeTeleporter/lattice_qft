use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use lattice_qft::{
    export::CsvData,
    lattice::Lattice,
    metropolis::{metropolis_simulation, MetropolisSimResult},
    LONG_TEMP_ARY,
};

const DIM: usize = 1;
const SIZE_ARY: [usize; DIM] = [2];
const SIZE: usize = SIZE_ARY[0];

const BURNIN: usize = 100_000;
const ITERATIONS: usize = 100_000_000;

const PATH: &str = "data/line_sim/line_sim.csv";

fn main() {
    let lattice: Lattice<DIM, SIZE> = Lattice::new(SIZE_ARY);

    println!("{:?}", lattice);

    let results: Vec<MetropolisSimResult> = LONG_TEMP_ARY
        .par_iter()
        .map(|&temp| {
            let (result, _) = metropolis_simulation(&lattice, temp, BURNIN, ITERATIONS);
            result
        })
        .collect();

    for res in results {
        if let Err(err) = res.read_write_csv(PATH) {
            eprint!("{}", err);
        }
    }
}
