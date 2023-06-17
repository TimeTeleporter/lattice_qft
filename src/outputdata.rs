//! This module implements data structures to save the field configuration
//! data of the Markov chain Monte Carlo in order to calculate observables
//! and create plots.

use rand::{rngs::ThreadRng, Rng};

use crate::{
    computation::FieldExport3d,
    fields::{Action, BondsFieldNew, Field, HeightVariable, WilsonField},
    kahan::KahanSummation,
    lattice::{Lattice, LatticeCoords},
};

// - OutputData ---------------------------------------------------------------

/// OutputData models the possible outputs that are able to be extracted from the computation.
///
/// Currently implemented is the following:
///
/// - Observables
///    - Action
/// - Plots
///    - Energy
///    - Differences
///
#[derive(Debug, Clone)]
pub enum OutputData<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    ActionObservable(ActionObservableNew),
    EnergyPlot(EnergyPlotOutputData<'a, D, SIZE>),
    DifferencePlot(DifferencePlotOutputData<'a, D, SIZE>),
    TestActionObservable(TestActionObservable),
    CorrelationData(CorrelationPlotOutputData<'a, D, SIZE>),
}

impl<'a, const D: usize, const SIZE: usize> OutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new_action_observable(temp: f64) -> Self {
        OutputData::ActionObservable(ActionObservableNew::new(temp))
    }

    pub fn new_energy_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData::EnergyPlot(EnergyPlotOutputData::new(lattice))
    }

    pub fn new_difference_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData::DifferencePlot(DifferencePlotOutputData::new(lattice))
    }

    pub fn new_test_action_observable(temp: f64) -> Self {
        OutputData::TestActionObservable(TestActionObservable::new(temp))
    }

    pub fn new_correlation_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData::CorrelationData(CorrelationPlotOutputData::new(lattice))
    }
}

// - OutputData - ActionObservableNew ------------------

#[derive(Debug, Clone)]
pub struct ActionObservableNew {
    value: KahanSummation<f64>,
    temp: f64,
}

impl ActionObservableNew {
    fn new(temp: f64) -> Self {
        ActionObservableNew {
            value: KahanSummation::new(),
            temp,
        }
    }
}

// - OutputData - TestActionObservable -----------------

/// The partition function is the sum over all Bolzmann weights. The observable
/// is the sum over all weights (Bolzmann times observable), devided by the
/// partition function.
#[derive(Debug, Clone)]
pub struct TestActionObservable {
    partfn: KahanSummation<f64>,
    obs: KahanSummation<f64>,
    temp: f64,
}

impl TestActionObservable {
    fn new(temp: f64) -> Self {
        TestActionObservable {
            partfn: KahanSummation::new(),
            obs: KahanSummation::new(),
            temp,
        }
    }

    pub fn wilson_update<'a, T: HeightVariable<T>, const SIZE: usize>(
        &mut self,
        wilson: &'a WilsonField<T, SIZE>,
    ) where
        WilsonField<'a, T, SIZE>: Action<T, 3>,
    {
        let bolz: f64 = (-wilson.action_observable(self.temp)).exp();
        self.obs
            .add(wilson.field.action_observable(self.temp) * bolz);
        self.partfn.add(bolz);
    }
}

// - OutputData - EnergyPlotOutputData -----------------------

#[derive(Debug, Clone)]
pub struct EnergyPlotOutputData<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    bonds: BondsFieldNew<'a, KahanSummation<f64>, D, SIZE>,
}

impl<'a, const D: usize, const SIZE: usize> EnergyPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        EnergyPlotOutputData {
            bonds: BondsFieldNew::new(lattice),
        }
    }
}

// - OutputData - DifferencePlotOutputData -------------------

#[derive(Debug, Clone)]
pub struct DifferencePlotOutputData<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    bonds: BondsFieldNew<'a, KahanSummation<f64>, D, SIZE>,
}

impl<'a, const D: usize, const SIZE: usize> DifferencePlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        DifferencePlotOutputData {
            bonds: BondsFieldNew::new(lattice),
        }
    }
}

// - OutputData - CorrelationPlotOutputData -------------------

/// Implements the averaging of correlation functions over configurations. When
/// calling the update function, the correlation length
///
/// $$C(t) = \sum_y (h_{x, t_0} - h_{y, t_0 + t})^2$$
///
/// gets calculated for each $t \in \{0, t_{max} - 1\}$, for a arbitrarily
/// chosen $x$ and $t_0$.
#[derive(Debug, Clone)]
pub struct CorrelationPlotOutputData<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    correlations: Vec<KahanSummation<f64>>,
}

impl<'a, const D: usize, const SIZE: usize> CorrelationPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        const TIME_DIMENSION: usize = 2;
        let max_t: usize = lattice.size[TIME_DIMENSION];
        let mut correlations: Vec<KahanSummation<f64>> = Vec::with_capacity(max_t);
        for _t in 0..max_t {
            correlations.push(KahanSummation::new())
        }

        CorrelationPlotOutputData {
            lattice,
            correlations,
        }
    }
}

// - UpdateOutputData ---------------------------------------------------------

pub trait UpdateOutputData<const D: usize, const SIZE: usize> {
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, rng: &mut ThreadRng)
    where
        [(); D * 2_usize]:;
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for OutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, rng: &mut ThreadRng)
    where
        [(); D * 2_usize]:,
    {
        match self {
            OutputData::ActionObservable(obs) => {
                <ActionObservableNew as UpdateOutputData<D, SIZE>>::update(obs, field, rng)
            }
            OutputData::TestActionObservable(obs) => {
                <TestActionObservable as UpdateOutputData<D, SIZE>>::update(obs, field, rng)
            }
            OutputData::EnergyPlot(plt) => {
                <EnergyPlotOutputData<'a, D, SIZE> as UpdateOutputData<D, SIZE>>::update(
                    plt, field, rng,
                )
            }
            OutputData::DifferencePlot(plt) => {
                <DifferencePlotOutputData<'a, D, SIZE> as UpdateOutputData<D, SIZE>>::update(
                    plt, field, rng,
                )
            }
            OutputData::CorrelationData(dat) => {
                <CorrelationPlotOutputData<'a, D, SIZE> as UpdateOutputData<D, SIZE>>::update(
                    dat, field, rng,
                )
            }
        }
    }
}

impl<const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for ActionObservableNew
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, _rng: &mut ThreadRng)
    where
        [(); D * 2_usize]:,
    {
        self.value.add(field.action_observable(self.temp))
    }
}

impl<const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for TestActionObservable
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, _rng: &mut ThreadRng)
    where
        [(); D * 2_usize]:,
    {
        let bolz: f64 = (-field.action_observable(self.temp)).exp();
        self.obs.add(field.action_observable(self.temp) * bolz);
        self.partfn.add(bolz);
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for EnergyPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, _rng: &mut ThreadRng) {
        for index in 0..SIZE {
            for direction in 0..D {
                let action = field.bond_action(index, direction);
                self.bonds
                    .get_value_mut(index, direction)
                    .add(action.into());
            }
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for DifferencePlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, _rng: &mut ThreadRng) {
        for index in 0..SIZE {
            for direction in 0..D {
                let diff = field.bond_difference(index, direction);
                self.bonds.get_value_mut(index, direction).add(diff.into());
            }
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for CorrelationPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, rng: &mut ThreadRng) {
        // Defining the first three dimensions
        const X_DIM: usize = 0;
        const Y_DIM: usize = 1;
        const T_DIM: usize = 2;

        // Reading the extents of the dimensions
        let max_x: usize = self.lattice.size[X_DIM];
        let max_y: usize = self.lattice.size[Y_DIM];
        let max_t: usize = self.lattice.size[T_DIM];

        // Building the array to get the index of our starting timeslice point
        let mut ary: [usize; D] = [0; D];
        ary[0] = rng.gen_range(0..max_x);
        ary[1] = rng.gen_range(0..max_y);

        // Index of the point on the starting timeslice (slice zero)
        let x_index: usize = self.lattice.calc_index_from_coords(LatticeCoords::new(ary));

        // Preparing the correlation function
        let mut correlation_function: Vec<KahanSummation<f64>> = Vec::with_capacity(max_t);
        for _t in 0..max_t {
            correlation_function.push(KahanSummation::new());
        }
        // For each lattice site, we calculate the difference squared and
        // associate it to a point of the correlation function.
        let x_value: T = field.get_value(x_index);
        for index in 0..SIZE {
            let t: usize = self.lattice.calc_coords_from_index(index).into_array()[T_DIM];
            let diff: T = x_value - field.get_value(index);
            let squared: T = diff * diff;
            correlation_function[t].add(squared.into());
        }
        // Converting into a Vec
        let correlation_function: Vec<f64> =
            correlation_function.into_iter().map(|x| x.sum()).collect();

        for (step, corr) in self.correlations.iter_mut().enumerate() {
            corr.add(correlation_function[step])
        }
    }
}

// - Observe ------------------------------------------------------------------

pub trait Observe {
    fn result(self) -> f64;
}

impl Observe for ActionObservableNew {
    fn result(self) -> f64 {
        self.value.mean()
    }
}

impl Observe for TestActionObservable {
    fn result(self) -> f64 {
        self.obs.mean() / self.partfn.mean()
    }
}

// - Plotting -----------------------------------------------------------------

pub trait Plotting<T> {
    fn plot(self) -> Vec<Vec<T>>;
}

impl<'a, const SIZE: usize> Plotting<FieldExport3d<f64>> for EnergyPlotOutputData<'a, 3, SIZE> {
    fn plot(self) -> Vec<Vec<FieldExport3d<f64>>> {
        let mut ary: Vec<Vec<FieldExport3d<f64>>> = Vec::new();
        for field in self.bonds.into_sub_fields().into_iter() {
            ary.push(Field::<'a, f64, 3, SIZE>::from_field(field).into_export());
        }
        ary
    }
}

impl<'a, const SIZE: usize> Plotting<FieldExport3d<f64>> for DifferencePlotOutputData<'a, 3, SIZE> {
    fn plot(self) -> Vec<Vec<FieldExport3d<f64>>> {
        let mut ary: Vec<Vec<FieldExport3d<f64>>> = Vec::new();
        for field in self.bonds.into_sub_fields().into_iter() {
            ary.push(Field::<'a, f64, 3, SIZE>::from_field(field).into_export());
        }
        ary
    }
}

impl<'a, const SIZE: usize> Plotting<f64> for CorrelationPlotOutputData<'a, 3, SIZE> {
    fn plot(self) -> Vec<Vec<f64>> {
        let mut ary: Vec<Vec<f64>> = Vec::new();
        ary.push(self.correlations.into_iter().map(|x| x.mean()).collect());
        ary
    }
}

#[test]
fn test_observable_array() {
    const D: usize = 3;
    const SIZE_ARY: [usize; D] = [4, 5, 6];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];
    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let temp: f64 = 1.0;
    let mut observables: Vec<OutputData<3, SIZE>> = Vec::new();
    observables.push(OutputData::new_action_observable(temp));
    observables.push(OutputData::new_energy_plot(&lattice));
    observables.push(OutputData::new_difference_plot(&lattice));
}
