#![feature(generic_const_exprs)]

// Number of spacetime dimensions of our lattice
pub const DIMENSIONS: usize = 3;

// The size of the spacetime volume
pub const LATTICE_LENGTH: usize = 2;
pub const LATTICE_WIDTH: usize = 2;
pub const LATTICE_HEIGHT: usize = 2;

mod form;
mod lattice;

use form::Form;
use lattice::Lattice;

type MyInt = i32;

type ZeroForm = Form<MyInt, DIMENSIONS, 0>;
type OneForm = Form<MyInt, DIMENSIONS, 1>;
type TwoForm = Form<MyInt, DIMENSIONS, 2>;
type ThreeForm = Form<MyInt, DIMENSIONS, 3>;

fn main() {
    let mut m_field = Lattice::<OneForm, LATTICE_LENGTH, LATTICE_WIDTH, LATTICE_HEIGHT>::new();

    println!("{:?}", m_field);
}
