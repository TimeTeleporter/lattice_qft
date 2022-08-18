#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

const MAX_X: usize = 10;
const MAX_Y: usize = 10;
const MAX_T: usize = 10;

mod form;
mod indexed_lattice3d;
mod lattice3d;

use indexed_lattice3d::IndexedLattice3d;
use lattice3d::LatticeValues;

fn main() {
    let lattice = IndexedLattice3d::<i32, MAX_X, MAX_Y, MAX_T>::default();

    lattice.print_values_formated();

    //println!("{:?}", lattice.indices.0);
}
