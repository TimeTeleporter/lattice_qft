#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

// Number of spacetime dimensions of our lattice
pub const DIMENSIONS: usize = 3;

const CUBE_SIDE: usize = 5;

// The size of the spacetime volume
pub const MAX_X: usize = CUBE_SIDE;
pub const MAX_Y: usize = CUBE_SIDE;
pub const MAX_T: usize = CUBE_SIDE;

mod form;
mod indexed_lattice3d;
mod lattice3d;

use indexed_lattice3d::IndexedLattice3d;

// The type of integer used in the simulation. Standard would be i32, can be increased to i64 or higher if needed.
type MyInt = i32;

fn main() {
    let lattice = IndexedLattice3d::<MyInt, 3, 3, 3>::default();

    println!("{:?}", lattice.indices.0);
}
