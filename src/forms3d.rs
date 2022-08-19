use crate::lattice3d::Lattice3d;

pub struct ZeroForm<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    values: [T; MAX_X * MAX_Y * MAX_T],
    lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    ZeroForm<'_, T, MAX_X, MAX_Y, MAX_T>
where
    T: std::fmt::Debug + Copy,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn set_value(&mut self, value: T, x: usize, y: usize, t: usize) {
        let index: usize = self.lattice.get_index_from_coordinates(x, y, t);
        self.values[index] = value;
        //println!("Set to value {:?}", value);
    }

    pub fn print_values_formated(&self) {
        for t in 0..MAX_T {
            println!("t = {}", t);
            for y in 0..MAX_Y {
                print!("[");
                for x in 0..MAX_X {
                    if x == MAX_X - 1 {
                        println!("{:?} ]", self.get_value_from_coordinates(x, y, t));
                    } else {
                        print!("{:?}, ", self.get_value_from_coordinates(x, y, t));
                    }
                }
            }
        }
    }

    pub fn get_value_from_coordinates(&self, x: usize, y: usize, t: usize) -> T {
        self.get_value(self.lattice.get_index_from_coordinates(x, y, t))
    }

    pub fn get_value(&self, index: usize) -> T {
        self.values[index]
    }
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    From<&'a Lattice3d<MAX_X, MAX_Y, MAX_T>> for ZeroForm<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Default + Copy,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn from(lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>) -> Self {
        ZeroForm::<'a, T, MAX_X, MAX_Y, MAX_T> {
            values: [T::default(); MAX_X * MAX_Y * MAX_T],
            lattice,
        }
    }
}
