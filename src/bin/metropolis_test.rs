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
};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;
const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points

const TEST_PATH: &str = "data/test_data32.csv";

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

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    const TEMP_ARY: [f64; 41] = [
        0.001, 0.0012589, 0.0015848, 0.0019952, 0.0025118, 0.0031622, 0.0039810, 0.0050118,
        0.0063095, 0.0079432, 0.01, 0.012589, 0.015848, 0.019952, 0.025118, 0.031622, 0.039810,
        0.050118, 0.063095, 0.079432, 0.1, 0.12589, 0.15848, 0.19952, 0.25118, 0.31622, 0.39810,
        0.50118, 0.63095, 0.79432, 1.0, 1.2589, 1.5848, 1.9952, 2.5118, 3.1622, 3.9810, 5.0118,
        6.3095, 7.9432, 10.0,
    ];

    if let Err(err) = clean_csv(TEST_PATH) {
        eprint!("Cleaning Error: {}", err)
    };

    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

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
            let bolz: f64 = (-field.lattice_action(temp)).exp();
            test = test + (field.lattice_action(temp) * bolz);
            partfn = partfn + bolz;
        }

        TestData {
            temp,
            test_range: TEST_RANGE,
            test_observable: test / partfn,
        }
    }
}
