#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

const MAX_X: usize = 10;
const MAX_Y: usize = 10;
const MAX_T: usize = 10;

mod forms3d;
mod indexed_lattice3d;
mod lattice3d;

use forms3d::Codifferential;

use crate::{forms3d::ZeroForm, lattice3d::LatticeValues};

fn main() {
    let mut wilson = ZeroForm::<i32, MAX_X, MAX_Y, MAX_T>::default();

    wilson.print_values_formated();
    //println!("{:?}", lattice.indices.0);

    wilson.values.0[255] = 1;

    wilson.print_values_formated();

    let links = wilson.increase_degree();

    links.print_values_formated();

    //todo!("wilson.insert_wilson_loop()");
}
