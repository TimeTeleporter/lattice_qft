//! A three dimensional spacetime lattice saved to an array.
//!
//! The x coordinate is periodic from 0 to MAX_X, the same for y and t.

use std::ops::Deref;

/// D-Dimensional coordinate array to simplify the conversation to and from index.
#[derive(Debug, Clone, Copy)]
pub struct LatticeCoords<const D: usize>([usize; D]);

impl<const D: usize> LatticeCoords<D> {
    fn cleanup(mut self, size: [usize; D]) -> Self {
        for i in 0_usize..D {
            self.0[i] = self.0[i] % size[i];
        }
        self
    }

    /// Change the coordinates in one direction by one.
    fn move_one(mut self, size: [usize; D], direction: usize) -> Self {
        let mut is_forward: bool = false;
        if direction < D {
            is_forward = true;
        }
        let direction: usize = direction % D;
        self.0[direction] = if is_forward {
            self.0[direction] + 1_usize
        } else {
            self.0[direction] + (size[direction] - 1_usize)
        };
        self.cleanup(size)
    }

    #[cfg(test)]
    #[allow(dead_code)]
    pub fn into_array(self) -> [usize; D] {
        self.0
    }

    pub fn new(array: [usize; D]) -> Self {
        Self(array)
    }
}

/// The Lattice datatype controls lattice indices in order to aid initialize
/// data on the lattice. It saves for each lattice point index the indices of
/// its neighbours.
///
/// The first ```i < D``` entries of ```self.values``` are the neighbours in
/// positive coordinate direction, the entries ```D <= i < D * 2_usize``` are
/// the neighbours in negative coordinate direction.
#[derive(Debug)]
pub struct Lattice<const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    pub size: [usize; D],
    pub values: Vec<[usize; D * 2_usize]>,
}

impl<const D: usize, const SIZE: usize> Lattice<D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(size: [usize; D]) -> Self {
        let values = vec![[0; D * 2_usize]; SIZE];
        let empty_lattice: Lattice<D, SIZE> = Lattice { size, values };

        let mut values = vec![[0; D * 2_usize]; SIZE];

        for (index, neighbours) in values.iter_mut().enumerate() {
            let coords = empty_lattice.calc_coords_from_index(index);
            for direction in 0_usize..(D * 2_usize) {
                neighbours[direction] =
                    empty_lattice.calc_index_from_coords(coords.move_one(size, direction));
            }
        }

        Lattice { size, values }
    }

    /// Returns all neighbours of the given index.
    pub fn get_neighbours_array(&self, index: usize) -> [usize; D * 2_usize] {
        self.values[index]
    }

    /// Returns all neighbours of the given index in all positive coordinate directions.
    pub fn pos_neighbours_array(&self, index: usize) -> [usize; D] {
        let (&ary, _) = self.get_neighbours_array(index).split_array_ref::<D>();
        ary
    }

    pub fn calc_index_from_coords(&self, coords: LatticeCoords<D>) -> usize {
        let coords = coords.cleanup(self.size);
        let mut sum: usize = 0_usize;
        for i in 0_usize..D {
            let mut product: usize = coords.0[i];
            for j in 0_usize..i {
                product = product * self.size[j];
            }
            sum = sum + product;
        }
        sum
    }

    pub fn calc_coords_from_index(&self, index: usize) -> LatticeCoords<D> {
        let mut coords_ary: [usize; D] = [0; D];

        for i in 0_usize..D {
            let mut product: usize = 1_usize;
            for j in 0_usize..i {
                product = product * self.size[j];
            }
            coords_ary[i] = (index / product) % self.size[i];
        }

        LatticeCoords(coords_ary)
    }
}

/// The Lattice3d datatype controls lattice indices in order to aid initialize data on the lattice.
/// For each data entry it has 6 neighbours, which we save in a array.
pub struct Lattice3d<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    Lattice<3, { MAX_X * MAX_Y * MAX_T }>,
)
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Lattice3d<MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn new() -> Self {
        Lattice3d(Lattice::<3, { MAX_X * MAX_Y * MAX_T }>::new([
            MAX_X, MAX_Y, MAX_T,
        ]))
    }
}

// By implementing deref we can utilize all the functions of the tuple element
// without implementing it ourselves.
impl<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Deref
    for Lattice3d<MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type Target = Lattice<3, { MAX_X * MAX_Y * MAX_T }>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[test]
fn test_index_coordinates_conv() {
    const MEASURES: [usize; 3] = [4, 5, 3];
    const SIZE: usize = MEASURES[0] * MEASURES[1] * MEASURES[2];
    let lattice: Lattice<3, SIZE> = Lattice::new(MEASURES);
    let index: usize = 29;
    let coords = lattice.calc_coords_from_index(index);
    assert_eq!(coords.into_array(), [1, 2, 1]);
    assert_eq!(lattice.calc_index_from_coords(coords), index);
}

/// Now consider the example of a <4,5,3> 3-dimensional lattice. We index the entries as follows:
///
///           y            t=0      y            t=1      y            t=2
///         4 ^  16 17 18 19        ^  36 37 38 39        ^  56 57 58 59  
///         3 |  12 13 14 15        |  32 33 34 35        |  52 53 54 55  
///         2 |  8  9  10 11        |  28 29 30 31        |  48 49 50 51  
///         1 |  4  5  6  7         |  24 25 26 27        |  44 45 46 47  
///         0 |  0  1  2  3         |  20 21 22 23        |  40 41 42 43  
///           .-----------> x       .-----------> x       .-----------> x
///              0  1  2  3            0  1  2  3            0  1  2  3
///
#[test]
fn test_neighbour_index_array_example() {
    let lattice: Lattice<3, { 4 * 5 * 3 }> = Lattice::new([4, 5, 3]);
    let index: usize = 19;
    let test: [usize; 6] = [16, 3, 39, 18, 15, 59];

    assert_eq!(lattice.values[index], test);
}

/// This tests if for all lattice sites in all direction if the neighbours
/// neighbour is oneself.
#[test]
fn test_big_dim_neighbour_list() {
    const D: usize = 5;
    const SIZE_ARY: [usize; D] = [90, 234, 25, 3, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2] * SIZE_ARY[3] * SIZE_ARY[4];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    for (index, neighbours) in lattice.values.iter().enumerate() {
        for (direction, &neighbour) in neighbours.iter().enumerate() {
            assert_eq!(
                index,
                lattice.values[neighbour][(direction + D) % (D * 2_usize)]
            );
        }
    }
}
