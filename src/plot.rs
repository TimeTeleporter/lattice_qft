use serde::{Deserialize, Serialize};

use crate::{
    fields::{Action, Field, HeightVariable, WilsonField},
    lattice::Lattice,
};

#[derive(Debug, Clone)]
/// This enum serves as an input parameter for computations to determine if and what plots should be generated at the end
pub enum PlotType<'a, const SIZE: usize> {
    /// A plot of the energy is generated as the square sum over all spatial bonds minus the temporal bonds
    Energy(&'a Lattice<2, SIZE>),
    /// This plot yields arrows at given lattice points in cardinal directions of strength proportional to the difference of lattice points
    ElectricField(&'a Lattice<2, SIZE>),
}

/// A single datapoint
#[derive(Debug, Serialize, Deserialize)]
pub struct Plotpoint<T: HeightVariable<T>> {
    pos_x: usize,
    pos_y: usize,
    arr_x: T,
    arr_y: T,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Plot2d {
    size: [usize; 2],
    pub values: Vec<Plotpoint<i32>>,
}

impl Plot2d {
    fn new(size_ary: [usize; 2]) -> Plot2d {
        let vec: Vec<Plotpoint<i32>> = Vec::new();
        Plot2d {
            size: size_ary,
            values: vec,
        }
    }

    fn push(&mut self, value: Plotpoint<i32>) {
        self.values.push(value)
    }
}

// - Plotting ComputationFields -----------------------------------------------
#[derive(Debug, Clone)]
pub enum ComputationField<'a, T, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    Field(Field<'a, T, D, SIZE>),
    WilsonField(WilsonField<'a, T, SIZE>),
}

impl<'a, const SIZE: usize> ComputationField<'a, i32, 3, SIZE> {
    pub fn plot<const PLOTSIZE: usize>(
        self,
        plot_type: Option<PlotType<'a, PLOTSIZE>>,
    ) -> Option<Plot2d> {
        // Checking if the plot is required
        let Some(plot_type) = plot_type else {
			return None;
		};

        // Extracting the field from ComputationField
        let field: Field<i32, 3, SIZE> = match self {
            ComputationField::Field(field) => field,
            ComputationField::WilsonField(wilson) => wilson.field,
        };

        // Building a field a as a sum over all timeslices at a given position
        let plot_lattice: &Lattice<2, PLOTSIZE> = match plot_type {
            PlotType::Energy(plot_lattice) => plot_lattice,
            PlotType::ElectricField(plot_lattice) => plot_lattice,
        };

        let plot_field: Field<'a, i32, 2, PLOTSIZE> = Field::new(plot_lattice);
        let mut plot_field: Field<'a, i32, 2, PLOTSIZE> = Field::from_field(plot_field);

        for index in 0..PLOTSIZE {
            let mut sum: i32 = 0;
            let mut multiplier: usize = 0;
            while let Some(value) = field.values.get(index + multiplier * PLOTSIZE) {
                sum = sum + value;
                multiplier = multiplier + 1;
            }
            plot_field.values[index] = sum;
        }

        // Creating a new plot
        let mut plot: Plot2d = Plot2d::new(plot_lattice.size);

        match plot_type {
            PlotType::Energy(_) =>
            // For each point on the plot_lattice we create a plot entry for each
            // cardinal direction. If the value is negtive, we move it to the
            // inverse position.
            {
                for index in 0..PLOTSIZE {
                    for direction in 0..2 {
                        let pos: [usize; 2] =
                            plot_lattice.calc_coords_from_index(index).into_array();
                        let mut arrows: [i32; 2] = [0, 0];

                        let value: i32 = plot_field.bond_action(index, direction);

                        arrows[direction] = value;

                        let point = Plotpoint {
                            pos_x: pos[0],
                            pos_y: pos[1],
                            arr_x: arrows[0],
                            arr_y: arrows[1],
                        };
                        plot.push(point);
                    }
                }
            }
            PlotType::ElectricField(_) =>
            // For each point on the plot_lattice we create a plot entry for each
            // cardinal direction. If the value is negtive, we move it to the
            // inverse position.
            {
                for index in 0..PLOTSIZE {
                    for direction in 0..2 {
                        let mut pos: [usize; 2] =
                            plot_lattice.calc_coords_from_index(index).into_array();
                        let mut arrows: [i32; 2] = [0, 0];

                        let value: i32 = plot_field.bond_difference(index, direction);

                        if value < 0 {
                            pos[direction] = pos[direction] + 1;
                        }

                        arrows[direction] = value;

                        let point = Plotpoint {
                            pos_x: pos[0],
                            pos_y: pos[1],
                            arr_x: arrows[0],
                            arr_y: arrows[1],
                        };
                        plot.push(point);
                    }
                }
            }
        }

        Some(plot)
    }
}
