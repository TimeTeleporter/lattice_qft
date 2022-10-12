#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use std::{io, io::prelude::*, time::Instant};

use lattice_qft::{action::Action, cluster::Cluster, field::Field3d, lattice::Lattice3d};

fn main() {
    // The dimensions of the lattice
    const MAX_X: usize = 10;
    const MAX_Y: usize = 10;
    const MAX_T: usize = 10;

    const TEMP: f64 = 0.01;

    let time = Instant::now();

    // Initialize a lattice with the given dimensions
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    // Set a field on the lattice points with random values form -128 to 127
    let field: Field3d<i8, MAX_X, MAX_Y, MAX_T> = Field3d::random(&lattice);

    // Transform it to a field that can have i32 values (-2_147_483_648 to 2_147_483_647)
    let mut field: Field3d<i32, MAX_X, MAX_Y, MAX_T> = Field3d::from_field(field);

    // The simulation loop
    let mut sweeps: usize = 0;
    loop {
        // Print the values of the field
        field.print_values_formated();

        // Calculate the action and the time that program has run
        println!(
            "Action: {}, Sweeps: {}, Time since startup: {}",
            field.action_observable(),
            sweeps,
            time.elapsed().as_millis(),
        );

        // Expect input to look at the data
        pause();

        // Run a nuber of sweeps
        // This gets inaccurate for below 3 sweeps, as pause() cannot be isntantly triggered again.
        for _ in 0..5 {
            field.cluster_sweep(TEMP);
            sweeps += 1;
            field.normalize_random();
        }
    }
}

/// Pauses the program until a key is pressed.
#[allow(dead_code)]
fn pause() {
    let mut stdin = io::stdin();
    let mut stdout = io::stdout();

    // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
    write!(stdout, "Press any key to continue...").unwrap();
    stdout.flush().unwrap();

    // Read a single byte and discard
    let _ = stdin.read(&mut [0u8]).unwrap();
}
