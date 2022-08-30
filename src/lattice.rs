//! An implementation of a n-dimensional periodic euclidean spacetime lattice.

#[derive(Debug, Clone, Copy)]
pub struct LatticeCoords<const D: usize>([usize; D]);

impl<const D: usize> LatticeCoords<D> {
    fn cleanup(mut self, size: [usize; D]) -> Self {
        for i in 0_usize..D {
            self.0[i] = self.0[i] % size[i];
        }
        self
    }

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
    pub fn to_array(self) -> [usize; D] {
        self.0
    }

    pub fn new(array: [usize; D]) -> Self {
        Self(array)
    }
}

#[derive(Debug)]
pub struct Lattice<const D: usize>
where
    [(); D * 2_usize]:,
{
    pub size: [usize; D],
    pub values: Vec<[usize; D * 2_usize]>,
}

impl<const D: usize> Lattice<D>
where
    [(); D * 2_usize]:,
{
    pub fn new(size: [usize; D]) -> Self {
        let num_indices: usize = size.iter().product();

        let values = vec![[0; D * 2_usize]; num_indices];
        let empty_lattice: Lattice<D> = Lattice { size, values };

        let mut values = vec![[0; D * 2_usize]; num_indices];

        for (index, neighbours) in values.iter_mut().enumerate() {
            let coords = empty_lattice.calc_coords_from_index(index);
            for direction in 0_usize..(D * 2_usize) {
                neighbours[direction] =
                    empty_lattice.calc_index_from_coords(coords.move_one(size, direction));
            }
        }

        Lattice { size, values }
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

#[test]
fn test_index_coordinates_conv() {
    const SIZE: [usize; 3] = [4, 5, 3];
    let lattice: Lattice<3> = Lattice::new(SIZE);
    let index: usize = 29;
    let coords = lattice.calc_coords_from_index(index);
    assert_eq!(coords.0, [1, 2, 1]);
    assert_eq!(lattice.calc_index_from_coords(coords), index);
}

#[test]
fn test_neighbour_index_array() {
    let lattice: Lattice<3> = Lattice::new([4, 5, 3]);
    let index: usize = 19;
    let test: [usize; 6] = [16, 3, 39, 18, 15, 59];

    assert_eq!(lattice.values[index], test);
}
