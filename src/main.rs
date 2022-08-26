#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
//#![allow(dead_code)]

use std::time::Instant;

use field3d::{ConvertField, Field3d};
use lattice3d::Lattice3d;
use metropolis::Metropolis;

mod field3d;
mod lattice3d;
mod metropolis;

fn main() {
    const MAX_X: usize = 100;
    const MAX_Y: usize = 100;
    const MAX_T: usize = 100;

    let time = Instant::now();

    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::default();

    let m_field: Field3d<i8, MAX_X, MAX_Y, MAX_T> = Field3d::random(&lattice);

    let mut m_field: Field3d<f64, MAX_X, MAX_Y, MAX_T> = Field3d::from_field(m_field);

    m_field.print_values_formated();

    loop {
        for _ in 0..10 {
            m_field.metropolis_sweep();
        }

        println!("{}", m_field.calculate_action());
    }

    //m_field.print_values_formated();

    println!(
        "Program duration: {} milliseconds",
        time.elapsed().as_millis()
    );
}
