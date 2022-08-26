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

#[derive(Debug)]
pub struct Element3d {
    pub neighbours: [usize; 6],
}

impl Default for Element3d {
    fn default() -> Self {
        Element3d {
            neighbours: [0_usize; 6],
        }
    }
}

/// As many directions as there are dimensions.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Directions {
    X,
    Y,
    T,
}

/// Implements the From trait to get the dimension from a usize index.
impl From<usize> for Directions {
    fn from(u: usize) -> Self {
        match u % 3 {
            0 => Directions::X,
            1 => Directions::Y,
            2 => Directions::T,
            _ => panic!("Tried to acces fourth dimension."),
        }
    }
}

#[allow(dead_code)]
impl Directions {
    fn into_usize(self) -> usize {
        match self {
            Directions::X => 0,
            Directions::Y => 1,
            Directions::T => 2,
        }
    }

    pub fn as_array() -> [Directions; 3] {
        [Self::X, Self::Y, Self::T]
    }
}

/// The Lattice3d datatype controls lattice indices in order to aid initialize data on the lattice.
/// For each data entry it has 6 neighbours, which we save in a array.
#[derive(Debug)]
pub struct Lattice3d<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    pub Vec<Element3d>,
)
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Lattice3d<MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    // Methods to utilize the saved data by reference
    pub fn get_next_neighbour_index(&self, index: usize, direction: Directions) -> usize {
        self.0[index].neighbours[direction.into_usize()]
    }

    pub fn get_prev_neighbour_index(&self, index: usize, direction: Directions) -> usize {
        self.0[index].neighbours[direction.into_usize() + 3]
    }

    pub fn get_index_from_coordinates(&self, (x, y, t): (usize, usize, usize)) -> usize {
        Self::calculate_index_from_coordinates(x, y, t)
    }

    // Methods to calculate the coordinates for the array index and vice versa.
    pub fn calculate_index_from_coordinates(x: usize, y: usize, t: usize) -> usize {
        MAX_X * MAX_Y * (t % MAX_T) + MAX_X * (y % MAX_Y) + (x % MAX_X)
    }

    fn calculate_coordinates_from_index(index: usize) -> (usize, usize, usize) {
        (
            Self::calculate_x_from_index(index),
            Self::calculate_y_from_index(index),
            Self::calculate_t_from_index(index),
        )
    }

    fn calculate_x_from_index(index: usize) -> usize {
        index % MAX_X
    }

    fn calculate_y_from_index(index: usize) -> usize {
        (index / MAX_X) % MAX_Y
    }

    fn calculate_t_from_index(index: usize) -> usize {
        index / (MAX_X * MAX_Y) % MAX_T
    }

    // Methods to get the neighbouring indices.
    fn calculate_next_neighbour_index(index: usize, direction: Directions) -> usize {
        let (mut x, mut y, mut t) = Self::calculate_coordinates_from_index(index);
        match direction {
            Directions::X => x = (x + 1) % MAX_X,
            Directions::Y => y = (y + 1) % MAX_Y,
            Directions::T => t = (t + 1) % MAX_T,
        }

        Self::calculate_index_from_coordinates(x, y, t)
    }

    fn calculate_prev_neighbour_index(index: usize, direction: Directions) -> usize {
        let (mut x, mut y, mut t) = Self::calculate_coordinates_from_index(index);
        match direction {
            Directions::X => x = (x + MAX_X - 1) % MAX_X,
            Directions::Y => y = (y + MAX_Y - 1) % MAX_Y,
            Directions::T => t = (t + MAX_T - 1) % MAX_T,
        }

        Self::calculate_index_from_coordinates(x, y, t)
    }
}

// Implementing the Default trait in order to have a constructor.
impl<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Default
    for Lattice3d<MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn default() -> Lattice3d<MAX_X, MAX_Y, MAX_T> {
        let elements = vec![(); MAX_X * MAX_Y * MAX_T];

        let mut elements: Vec<Element3d> = elements.iter().map(|_| Element3d::default()).collect();

        for (index, element) in elements.iter_mut().enumerate() {
            for i in 0..3 {
                element.neighbours[i] =
                    Self::calculate_next_neighbour_index(index, Directions::from(i));
                element.neighbours[i + 3] =
                    Self::calculate_prev_neighbour_index(index, Directions::from(i));
            }
        }

        Lattice3d(elements)
    }
}

#[test]
fn test_index_coordinates_conversion() {
    let (x, y, t) = Lattice3d::<4, 5, 3>::calculate_coordinates_from_index(29);
    let index = Lattice3d::<4, 5, 3>::calculate_index_from_coordinates(x, y, t);

    assert_eq!((x, y, t), (1, 2, 1));
    assert_eq!(index, 29);
}

#[test]
fn test_neighbour_index_array() {
    let lattice = Lattice3d::<4, 5, 3>::default();
    let index: usize = 19;
    let test: [usize; 6] = [16, 3, 39, 18, 15, 59];

    assert_eq!(lattice.0[index].neighbours, test);
}

#[test]
fn test_directions_usize_conversions() {
    let indices: Vec<usize> = vec![0, 1, 2];
    let directions: Vec<Directions> = vec![Directions::X, Directions::Y, Directions::T];

    let directions_2: Vec<Directions> = indices.iter().map(|&x| Directions::from(x)).collect();

    assert_eq!(directions, directions_2);

    let indices_2: Vec<usize> = directions
        .into_iter()
        .map(|direction| direction.into_usize())
        .collect();

    assert_eq!(indices, indices_2);
}
