//! A three dimensional spacetime lattice saved to an array.
//!
//! The x coordinate is periodic from 0 to MAX_X, the same for y and t.
//!
//! Now consider the example of a <4,5,3> 3-dimensional lattice. We index the entries as follows:
//!
//!           y            t=0      y            t=1      y            t=2
//!         4 ^  16 17 18 19        ^  36 37 38 39        ^  56 57 58 59  
//!         3 |  12 13 14 15        |  32 33 34 35        |  52 53 54 55  
//!         2 |  8  9  10 11        |  28 29 30 31        |  48 49 50 51  
//!         1 |  4  5  6  7         |  24 25 26 27        |  44 45 46 47  
//!         0 |  0  1  2  3         |  20 21 22 23        |  40 41 42 43  
//!           .-----------> x       .-----------> x       .-----------> x
//!              0  1  2  3            0  1  2  3            0  1  2  3
//!

use crate::lattice::{Lattice, LatticeCoords};

/// The Lattice3d datatype controls lattice indices in order to aid initialize data on the lattice.
/// For each data entry it has 6 neighbours, which we save in a array.
pub struct Lattice3d<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(Lattice<3_usize>);

impl<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Lattice3d<MAX_X, MAX_Y, MAX_T> {
    pub fn get_neighbours_array(&self, index: usize) -> [usize; 3 * 2_usize] {
        self.0.values[index]
    }

    pub fn new() -> Self {
        Lattice3d(Lattice::<3>::new([MAX_X, MAX_Y, MAX_T]))
    }

    pub fn calc_index_from_coords(&self, coords: LatticeCoords<3>) -> usize {
        self.0.calc_index_from_coords(coords)
    }

    #[cfg(test)]
    #[allow(dead_code)]
    pub fn calc_coords_from_index(&self, index: usize) -> LatticeCoords<3> {
        self.0.calc_coords_from_index(index)
    }
}

#[test]
/// This tests the <4,5,3> example for the correct neighbours.
fn test_index_coordinates_conversion() {
    let lattice: Lattice3d<4, 5, 3> = Lattice3d::new();
    let coords: LatticeCoords<3> = lattice.calc_coords_from_index(29);
    let index: usize = lattice.calc_index_from_coords(coords);

    assert_eq!(coords.to_array(), [1, 2, 1]);
    assert_eq!(index, 29);
}
