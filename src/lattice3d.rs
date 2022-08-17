/// A three dimensional spacetime lattice.
// This might need reconsideration, as one is able to simplify it to a single array.
#[derive(Debug)]
pub struct Lattice3d<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    pub [[[T; MAX_T]; MAX_Y]; MAX_X],
);

/// As many directions as there are dimensions.
pub enum Direction {
    X,
    Y,
    T,
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Default
    for Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
{
    fn default() -> Lattice3d<T, MAX_X, MAX_Y, MAX_T> {
        let ary = [(); MAX_X].map(|_| [(); MAX_Y].map(|_| [(); MAX_T].map(|_| T::default())));
        Lattice3d(ary)
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default + Copy,
{
    pub fn prev_neighbour_value(&self, x: usize, y: usize, t: usize, direction: Direction) -> T {
        let (x, y, t) = match direction {
            Direction::X => (self.prev_index(x, direction), y, t),
            Direction::Y => (x, self.prev_index(y, direction), t),
            Direction::T => (x, y, self.prev_index(t, direction)),
        };
        self.value(x, y, t)
    }

    fn prev_index(&self, index: usize, direction: Direction) -> usize {
        if index == 0 {
            self.get_max(direction) - 1
        } else {
            index - 1
        }
    }

    pub fn next_neighbour_value(&self, x: usize, y: usize, t: usize, direction: Direction) -> T {
        let (x, y, t) = match direction {
            Direction::X => (self.next_index(x, direction), y, t),
            Direction::Y => (x, self.next_index(y, direction), t),
            Direction::T => (x, y, self.next_index(t, direction)),
        };
        self.value(x, y, t)
    }

    fn next_index(&self, index: usize, direction: Direction) -> usize {
        (index + 1) % self.get_max(direction)
    }

    fn get_max(&self, direction: Direction) -> usize {
        match direction {
            Direction::X => MAX_X,
            Direction::Y => MAX_Y,
            Direction::T => MAX_T,
        }
    }

    pub fn value(&self, x: usize, y: usize, t: usize) -> T {
        self.0[x % MAX_X][y % MAX_X][t % MAX_T]
    }

    pub fn new() -> Lattice3d<T, MAX_X, MAX_Y, MAX_T> {
        Lattice3d::<T, MAX_X, MAX_Y, MAX_T>::default()
    }
}
