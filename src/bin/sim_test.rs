#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use lattice_qft::{
    action::Action,
    export::{clean_csv, CsvData},
    field::Field3d,
    lattice::Lattice3d,
    simulation::Simulation3d,
    TEMP_ARY,
};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;
const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points

// Parameters for the Wilson loop
const WIDTH: usize = 1;
const HEIGHT: usize = 1;

const TEST_PATH: &str = "data/test_results.csv";

const TEST_RANGE: usize = 32;
const PERMUTATIONS: usize = TEST_RANGE.pow(SIZE as u32 - 1); // 16 ^ 7 = 268’435’456
const BOUDARY: i8 = TEST_RANGE as i8 - 1;

/// Datatype to save and read simulation output.
#[derive(Debug, Default, Serialize, Deserialize)]
struct TestData {
    temp: f64,
    test_range: usize,
    test_observable: f64,
}

impl TestData {
    fn calculate_configurations(
        lattice: &Lattice3d<TEST_X, TEST_Y, TEST_T>,
        temp: f64,
    ) -> TestData {
        print!("Started to calculate the observable over the range");
        println!(" -{} to +{}.", TEST_RANGE / 2, TEST_RANGE / 2);
        println!("Those are {} field configurations.", PERMUTATIONS);

        let mut field: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::new(lattice);

        field.values[7] = (TEST_RANGE / 2) as i8;

        // The partition function is the sum over all Bolzmann weights
        let mut partfn: f64 = 0.0;

        // The observable is the sum over all weights (Bolzmann times observable),
        // devided by the partition function.
        let mut test: f64 = 0.0;

        for _ in 0..PERMUTATIONS {
            'updateconfig: for index in 0..SIZE {
                match field.values[index] {
                    x if x < BOUDARY => {
                        field.values[index] = field.values[index] + 1;
                        break 'updateconfig;
                    }
                    x if x == BOUDARY => {
                        field.values[index] = 0;
                    }
                    _ => {
                        panic!("config entry out of bounds.");
                    }
                }
            }
            let field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::from_field(field.clone());
            //field.print_values_formated();
            let bolz: f64 = (-field.action_observable(temp)).exp();
            test = test + (field.wilson_loop_observable(temp, WIDTH, HEIGHT) * bolz);
            partfn = partfn + bolz;
        }

        TestData {
            temp,
            test_range: TEST_RANGE,
            test_observable: test / partfn,
        }
    }
}

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    if let Err(err) = clean_csv(TEST_PATH) {
        eprint!("Cleaning Error: {}", err)
    };

    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    // Initialise the simulations
    let sims: Vec<Simulation3d> = Vec::new();
    let name: String =
        format!("{TEST_X}x{TEST_Y}x{TEST_T} Metropolis Simulation, {WIDTH}x{HEIGHT} Wilson loop");
    sims.push(Simulation3d::new(
        name, sim_type, observable, &lattice, TEMP, 0, 0,
    ));

    // Parallel over all temp data
    let data: Vec<TestData> = TEMP_ARY
        .par_iter()
        .map(|&temp| {
            let config = TestData::calculate_configurations(&lattice, temp);

            println!("{}: TestData calculated", temp);

            config
        })
        .collect();

    for entry in data {
        if let Err(err) = entry.read_write_csv(TEST_PATH) {
            eprint!("TestData Error: {}", err);
        };
    }
}
