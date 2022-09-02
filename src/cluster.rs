use rand::{rngs::ThreadRng, Rng};

use crate::{field3d::Field, observable::Action};

pub trait Cluster: Action {
    fn cluster_sweep(&mut self);
}

impl<'a, const D: usize, const SIZE: usize> Cluster for Field<'a, i32, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn cluster_sweep(&mut self) {
        let mut rng = ThreadRng::default();

        // Initialize the bonds to be activated
        let mut bonds: Vec<[bool; D]> = vec![[false; D]; SIZE];

        // Set the plane to be reflected about
        let reflect: i32 = self.values[rng.gen_range(0..SIZE)];

        // Activate the bonds, that lie on the same side of the plane.
        for (index, neighbours) in self.lattice.values.iter().enumerate() {
            for link in 0..D {
                if (self.values[index] - reflect) * (self.values[neighbours[link]] - reflect) >= 0 {
                    bonds[index][link] = true;
                }
            }
        }
    }
}
