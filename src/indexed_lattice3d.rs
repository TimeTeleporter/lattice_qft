use crate::lattice3d::{Direction, Lattice, Lattice3d, LatticeValues};

/// An optimized Lattice3d, where we create a list of neighbouring indices.
/// As such, we can simplify the calling of neighbouring values.
/// The neighbours are saved in the 'indices' field as follows:
///
///
///             [1] [5]      y   
///              | /         |
///         [3]-[x]-[0]      o - x
///            / |          /
///         [2] [4]        t
///
///
#[derive(Debug)]
pub struct IndexedLattice3d<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub values: Lattice3d<T, MAX_X, MAX_Y, MAX_T>,
    pub indices: Lattice3d<[usize; 6], MAX_X, MAX_Y, MAX_T>,
}

// We implement the Default trait
impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Default
    for IndexedLattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn default() -> IndexedLattice3d<T, MAX_X, MAX_Y, MAX_T> {
        let values = Lattice3d::<T, MAX_X, MAX_Y, MAX_T>::default();
        let mut indices = Lattice3d::<[usize; 6], MAX_X, MAX_Y, MAX_T>::default();

        for index in 0..MAX_X * MAX_Y * MAX_T {
            indices.0[index][0] = values.next_neighbour_index(index, Direction::X);
            indices.0[index][1] = values.next_neighbour_index(index, Direction::Y);
            indices.0[index][2] = values.next_neighbour_index(index, Direction::T);
            indices.0[index][3] = values.prev_neighbour_index(index, Direction::X);
            indices.0[index][4] = values.prev_neighbour_index(index, Direction::Y);
            indices.0[index][5] = values.prev_neighbour_index(index, Direction::T);
        }

        IndexedLattice3d { values, indices }
    }
}

// We implement the Lattice trait
impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Lattice<MAX_X, MAX_Y, MAX_T>
    for IndexedLattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn prev_neighbour_index(&self, index: usize, direction: Direction) -> usize {
        let entry: usize = match direction {
            Direction::X => 3,
            Direction::Y => 4,
            Direction::T => 5,
        };

        self.indices.0[index][entry]
    }

    fn next_neighbour_index(&self, index: usize, direction: Direction) -> usize {
        let entry: usize = match direction {
            Direction::X => 0,
            Direction::Y => 1,
            Direction::T => 2,
        };

        self.indices.0[index][entry]
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    LatticeValues<T, MAX_X, MAX_Y, MAX_T> for IndexedLattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: std::fmt::Debug,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn print_values(&self) {
        self.values.print_values()
    }

    fn get_value(&self, index: usize) -> &T {
        self.values.get_value(index)
    }
}

#[test]
fn test_indexed_neighbour() {
    let lattice = IndexedLattice3d::<i32, 4, 5, 3>::default();
    assert_eq!(lattice.indices.0[19], [16, 3, 39, 18, 15, 59]);
}

#[test]
fn test_indexed_lattice3d_default_values() {
    let lattice = IndexedLattice3d::<i32, 2, 2, 3>::default();
    assert_eq!(lattice.values.0, [0_i32; 12]);
}
