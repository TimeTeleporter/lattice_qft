#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
//#![allow(dead_code)]

use std::time::Instant;

use field3d::{ConvertField, Field3d};
use lattice3d::Lattice3d;
use metropolis::Metropolis3d;

mod field3d;
mod lattice3d;
mod metropolis;

fn main() {
    const MAX_X: usize = 10;
    const MAX_Y: usize = 10;
    const MAX_T: usize = 3;

    let time = Instant::now();

    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::default();

    let field: Field3d<i8, MAX_X, MAX_Y, MAX_T> = Field3d::random(&lattice);

    let mut field: Field3d<i32, MAX_X, MAX_Y, MAX_T> = Field3d::from_field(field);

    loop {
        field.print_values_formated();

        println!(
            "Action: {}, Time since startup: {}",
            field.calculate_action(),
            time.elapsed().as_millis()
        );

        pause();

        for _ in 0..1 {
            for _ in 0..1 {
                field.metropolis_sweep();
            }
        }
    }

    //field.print_values_formated();
}

use std::io;
use std::io::prelude::*;

fn pause() {
    let mut stdin = io::stdin();
    let mut stdout = io::stdout();

    // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
    write!(stdout, "Press any key to continue...").unwrap();
    stdout.flush().unwrap();

    // Read a single byte and discard
    let _ = stdin.read(&mut [0u8]).unwrap();
}
