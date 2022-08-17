#![feature(generic_const_exprs)]

// Number of spacetime dimensions of our lattice
pub const DIMENSIONS: usize = 3;

const CUBE_SIDE: usize = 5;

// The size of the spacetime volume
pub const MAX_X: usize = CUBE_SIDE;
pub const MAX_Y: usize = CUBE_SIDE;
pub const MAX_T: usize = CUBE_SIDE;

mod form;
mod lattice3d;

use form::Form;
use lattice3d::{Direction, Lattice3d};

// The type of integer used in the simulation. Standard would be i32, can be increased to i64 or higher if needed.
type MyInt = i32;

// The possible forms in three dimensions.
// TODO! Might want to simplify to not use forms.
type ZeroForm = Form<MyInt, DIMENSIONS, 0>;
type OneForm = Form<MyInt, DIMENSIONS, 1>;
type TwoForm = Form<MyInt, DIMENSIONS, 2>;
type ThreeForm = Form<MyInt, DIMENSIONS, 3>;

fn main() {
    let int: MyInt = MyInt::default();
}
