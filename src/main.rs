#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use forms3d::{OneForm, ZeroForm};
use lattice3d::{Directions, Lattice3d};

const MAX_X: usize = 100;
const MAX_Y: usize = 100;
const MAX_T: usize = 100;

mod forms3d;
mod lattice3d;

fn main() {
    let tiny_lattice = Lattice3d::<3, 3, 3>::default();

    let tiny_zero_form = ZeroForm::<i32, 3, 3, 3>::from(&tiny_lattice);

    tiny_zero_form.print_values_formated();

    let lattice = Lattice3d::<MAX_X, MAX_Y, MAX_T>::default();

    let mut zero_form = ZeroForm::<i32, MAX_X, MAX_Y, MAX_T>::from(&lattice);

    zero_form.set_line_to_value(3, Directions::T, (0, 1, 2), (0, 1, 0));

    //zero_form.print_values_formated();

    //println!("Then the one form:");

    let _one_form = OneForm::<i32, MAX_X, MAX_Y, MAX_T>::from(&zero_form);

    //one_form.print_values_formated();
}
