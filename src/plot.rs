use serde::Serialize;

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
#[derive(Debug, Serialize)]
pub struct Plotpoint<T: HeightVariable<T>, const D: usize>
where
    [usize; D]: Serialize,
    [T; D]: Serialize,
{
    pos: [usize; D],
    value: [T; D],
}

#[derive(Debug, Serialize)]
pub struct Plot3d {
    size: [usize; 3],
    values: Vec<Plotpoint<i32, 3>>,
}

impl Plot3d {
    fn new(size_ary: [usize; 3]) -> Plot3d {
        let vec: Vec<Plotpoint<i32, 3>> = Vec::new();
        Plot3d {
            size: size_ary,
            values: vec,
        }
    }

    fn push(&mut self, value: Plotpoint<i32, 3>) {
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
    ) -> Option<Plot3d> {
        let Some(plot_type) = plot_type else {
			return None;
		};

        let plot_lattice: &Lattice<2, PLOTSIZE> = match plot_type {
            PlotType::Energy(plot_lattice) => plot_lattice,
            PlotType::ElectricField(plot_lattice) => plot_lattice,
        };

        let plot_field: Field<'a, i32, 2, PLOTSIZE> = Field::new(plot_lattice);
        let mut plot_field: Field<'a, i32, 2, PLOTSIZE> = Field::from_field(plot_field);

        let field: Field<i32, 3, SIZE> = match self {
            ComputationField::Field(field) => field,
            ComputationField::WilsonField(wilson) => wilson.field,
        };

        for index in 0..PLOTSIZE {
            let mut sum = 0;
        }

        let plot: Plot3d = Plot3d::new(field.lattice.size);

        Some(plot)
    }
}
