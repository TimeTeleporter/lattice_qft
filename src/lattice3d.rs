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
//! This example gets tested in the tests below.

/// Lattice operation trait concerning the conversion between indices and coordinates
trait Lattice<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> {
    fn prev_neighbour_index(&self, index: usize, direction: Direction) -> usize {
        let (mut x, mut y, mut t) = Self::get_coordinates_from_index(index);
        match direction {
            Direction::X => x = (x + MAX_X - 1) % MAX_X,
            Direction::Y => y = (y + MAX_Y - 1) % MAX_Y,
            Direction::T => t = (t + MAX_T - 1) % MAX_T,
        }

        Self::get_index_from_coordinates(x, y, t)
    }

    fn next_neighbour_index(&self, index: usize, direction: Direction) -> usize {
        let (mut x, mut y, mut t) = Self::get_coordinates_from_index(index);
        match direction {
            Direction::X => x = (x + 1) % MAX_X,
            Direction::Y => y = (y + 1) % MAX_Y,
            Direction::T => t = (t + 1) % MAX_T,
        }

        Self::get_index_from_coordinates(x, y, t)
    }

    fn get_index_from_coordinates(x: usize, y: usize, t: usize) -> usize {
        MAX_X * MAX_Y * t + MAX_X * y + x
    }

    fn get_coordinates_from_index(index: usize) -> (usize, usize, usize) {
        (
            Self::get_x_from_index(index),
            Self::get_y_from_index(index),
            Self::get_t_from_index(index),
        )
    }

    fn get_x_from_index(index: usize) -> usize {
        index % MAX_X
    }

    fn get_y_from_index(index: usize) -> usize {
        (index / MAX_X) % MAX_Y
    }

    fn get_t_from_index(index: usize) -> usize {
        index / (MAX_X * MAX_Y) % MAX_T
    }
}

/// As many directions as there are dimensions.
pub enum Direction {
    X,
    Y,
    T,
}

// This is the unoptimized lattice datatype. We optimize it further below.
#[derive(Debug)]
pub struct Lattice3d<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    pub [T; MAX_X * MAX_Y * MAX_T],
)
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Default
    for Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn default() -> Lattice3d<T, MAX_X, MAX_Y, MAX_T> {
        Lattice3d([(); MAX_X * MAX_Y * MAX_T].map(|_| T::default()))
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Lattice<MAX_X, MAX_Y, MAX_T>
    for Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn get_value_from_coordinates(&self, x: usize, y: usize, t: usize) -> &T {
        self.get_value(Self::get_index_from_coordinates(x, y, t))
    }

    pub fn get_value(&self, index: usize) -> &T {
        &self.0[index]
    }

    /// Initializes a new 3-dimensional lattice with default values for the given types.
    pub fn new() -> Lattice3d<T, MAX_X, MAX_Y, MAX_T> {
        Lattice3d::<T, MAX_X, MAX_Y, MAX_T>::default()
    }
}

#[test]
fn test_index_coordinates_conversion() {
    type MyLattice = Lattice3d<i32, 4, 5, 3>;
    let (x, y, t) = MyLattice::get_coordinates_from_index(29);
    let index = MyLattice::get_index_from_coordinates(x, y, t);

    assert_eq!(x, 1);
    assert_eq!(y, 2);
    assert_eq!(t, 1);
    assert_eq!(index, 29);
}

#[test]
fn test_neighbour_index() {
    let lattice = Lattice3d::<i32, 4, 5, 3>::default();
    let center: usize = 19;

    assert_eq!(lattice.next_neighbour_index(center, Direction::X), 16);
    assert_eq!(lattice.next_neighbour_index(center, Direction::Y), 3);
    assert_eq!(lattice.next_neighbour_index(center, Direction::T), 39);
    assert_eq!(lattice.prev_neighbour_index(center, Direction::X), 18);
    assert_eq!(lattice.prev_neighbour_index(center, Direction::Y), 15);
    assert_eq!(lattice.prev_neighbour_index(center, Direction::T), 59);
}

/// An optimized Lattice3d, where we create a list of neighbouring indices.
/// As such, we can simplify the calling of neighbouring values.
/// The neighbours are saved in the 'indices' field as follows:
///
///             [1] [5]      y   
///              | /         |
///         [3]-[x]-[0]      o - x
///            / |          /
///         [2] [4]        t
///
#[derive(Debug)]
pub struct IndexedLattice3d<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub values: Lattice3d<T, MAX_X, MAX_Y, MAX_T>,
    pub indices: Lattice3d<[usize; 6], MAX_X, MAX_Y, MAX_T>,
}

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

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Lattice<MAX_X, MAX_Y, MAX_T>
    for IndexedLattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
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
