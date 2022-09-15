#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
//#![allow(dead_code)]

use std::io;
use std::io::prelude::*;
use std::time::Instant;

use field3d::Field3d;
use lattice3d::Lattice3d;
use observable::Action;

#[allow(unused_imports)]
use metropolis::Metropolis;

#[allow(unused_imports)]
use cluster::Cluster;

mod cluster;
mod field3d;
mod lattice3d;
mod metropolis;
mod observable;

fn main() {
    // The dimensions of the lattice
    const MAX_X: usize = 10;
    const MAX_Y: usize = 10;
    const MAX_T: usize = 10;

    let time = Instant::now();
    // Initialize a lattice with the given dimensions
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    // Set a field on the lattice points with random values form -128 to 127
    let field: Field3d<i8, MAX_X, MAX_Y, MAX_T> = Field3d::random(&lattice);

    // Transform it to a field that can have i32 values (-2_147_483_648 to 2_147_483_647)
    let mut field: Field3d<i32, MAX_X, MAX_Y, MAX_T> = Field3d::from_field(field);

    // The simulation loop
    loop {
        // Print the values of the field
        field.print_values_formated();

        // Calculate the action and the time that program has run
        println!(
            "Action: {}, Time since startup: {}",
            field.action_observable(),
            time.elapsed().as_millis()
        );

        // Expect input to look at the data
        pause();

        // Run a nuber of sweeps
        for _ in 0..1 {
            for _ in 0..10 {
                field.metropolis_sweep();
                field.normalize_random();
            }
        }
    }
}

/// Pauses the program until a key is pressed.
fn pause() {
    let mut stdin = io::stdin();
    let mut stdout = io::stdout();

    // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
    write!(stdout, "Press any key to continue...").unwrap();
    stdout.flush().unwrap();

    // Read a single byte and discard
    let _ = stdin.read(&mut [0u8]).unwrap();
}
