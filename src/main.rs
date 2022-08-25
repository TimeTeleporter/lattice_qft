#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use forms3d::{Form, OneForm};
use lattice3d::{Directions, Lattice3d};
use metropolis::Action;

const MAX_X: usize = 10;
const MAX_Y: usize = 10;
const MAX_T: usize = 10;

const E: f64 = 10.0;

mod forms3d;
mod lattice3d;
mod metropolis;

fn main() {
    let lattice = Lattice3d::<MAX_X, MAX_Y, MAX_T>::default();

    let one_form = OneForm::<i32, MAX_X, MAX_Y, MAX_T>::from_lattice(&lattice);

    let action = Action::<MAX_X, MAX_Y, MAX_T>::from_one_form(&one_form, E);

    action.print_values_formated();

    println!("{:?}", action.get_lattice_action());
}
