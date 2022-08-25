#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![allow(dead_code)]

use field3d::Field3d;
use lattice3d::Lattice3d;

const MAX_X: usize = 10;
const MAX_Y: usize = 10;
const MAX_T: usize = 10;

const E: f64 = 10.0;

mod field3d;
mod lattice3d;

fn main() {
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::default();

    let mut m_field: Field3d<i8, MAX_X, MAX_Y, MAX_T> = Field3d::random(&lattice);

    m_field.print_values_formated();
}
