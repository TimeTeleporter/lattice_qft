use rand::{rngs::ThreadRng, Rng};

use crate::{field3d::Field, observable::Action};

pub trait Cluster: Action {
    const TEMP: f64;

    fn cluster_sweep(&mut self);
}

impl<'a, const D: usize, const SIZE: usize> Cluster for Field<'a, i32, D, SIZE>
where
    [(); D * 2_usize]:,
{
    const TEMP: f64 = 1.0;

    fn cluster_sweep(&mut self) {
        let mut rng = ThreadRng::default();

        // Initialize the bonds to be activated
        let mut bonds: Vec<[bool; D * 2_usize]> = vec![[false; D * 2_usize]; SIZE];

        // Set the plane to be reflected about
        let height: i32 = self.values[rng.gen_range(0..SIZE)];
        let modifier: i32 = match rng.gen::<bool>() {
            true => 1,
            false => -1,
        };
        println!("Reflection plane: {}", height);

        // Activate the bonds, that lie on the same side of the plane.
        for (index, neighbours) in self.lattice.values.iter().enumerate() {
            let reflected: i32 = 2 * height + modifier - self.values[index];
            for direction in 0..D {
                let neighbour: i32 = self.values[neighbours[direction]];
                // Check if they are both on the same side:
                if (self.values[index] - neighbour).abs() < (reflected - neighbour) {
                    println!(
                        "Ok: site: {}, neighbour: {}, height + mod: {}",
                        self.values[index],
                        neighbour,
                        height + modifier
                    );
                    // calculate the action
                    let action = Self::calculate_link_action(self.values[index], neighbour);
                    let reflected_action = Self::calculate_link_action(reflected, neighbour);
                    // If they are, set the link as true with a chance of ´P_{bond} = 1 - exp(-(S' - S))´
                    let draw: f64 = rng.gen_range(0.0..1.0);
                    let prob: f64 = 1.0 - (f64::from(action - reflected_action) * Self::TEMP).exp();
                    println!("{}: draw: {}, prob: {}", index, draw, prob);
                    if draw <= prob {
                        bonds[index][direction] = true;
                        bonds[neighbours[direction]][direction + D] = true;
                    }
                }
            }
        }

        // Build the clusters
        let mut check: Vec<bool> = vec![true; SIZE];
        for (index, &neighbours) in self.lattice.values.iter().enumerate() {
            if check[index] {
                let mut cluster: Vec<usize> = Vec::new();

                // I want to recursively go through the neighbours and the
                // neighbours neighbours and so on to build the cluster.
                push_cluster::<D>(&mut cluster, &mut check, &bonds, index, neighbours);

                // After having built the cluster, we invert  it with a 50% chance.
                let coin: bool = rng.gen();
                if coin {
                    for entry in cluster {
                        self.values[entry] = 2 * height + modifier - self.values[entry];
                    }
                }
            }
        }
    }
}

fn push_cluster<const D: usize>(
    cluster: &mut Vec<usize>,
    check: &mut Vec<bool>,
    bonds: &Vec<[bool; D * 2_usize]>,
    index: usize,
    neighbours: [usize; D * 2_usize],
) where
    [(); D * 2_usize]:,
{
    if check[index] {
        check[index] = false;
        cluster.push(index);
        for (direction, &neighbour) in neighbours.iter().enumerate() {
            if bonds[index][direction] {
                push_cluster::<D>(cluster, check, bonds, neighbour, neighbours);
            }
        }
    }
}
