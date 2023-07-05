//! This module implements data structures to save the field configuration
//! data of the Markov chain Monte Carlo in order to calculate observables
//! and create plots.

use rand::rngs::ThreadRng;

use crate::{
    computation::FieldExport3d,
    fields::{Action, BondsFieldNew, Field, HeightVariable, WilsonField},
    kahan::KahanSummation,
    lattice::Lattice,
};

pub trait UpdateOutputData<const D: usize, const SIZE: usize> {
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, rng: &mut ThreadRng)
    where
        [(); D * 2_usize]:;
}

pub trait Observe {
    fn result(self) -> f64;
}

/// This trait implements the export of OutputData. The first level vector
/// contains (usually three) vectors for each cardinal direction, which then
/// contain the [FieldExport3d] data, that describes values on the lattice,
/// such as energy or field strength.
pub trait Plotting<T> {
    fn plot(self) -> Vec<Vec<T>>;
}

/// OutputData models the possible outputs that are able to be extracted from the computation.
///
/// Currently implemented is the following plots and observables
///
/// - Action
///     - Simulations
///     - Tests
/// - Energy plots
/// - Field strength plots
/// - Correlation length observables
#[derive(Debug, Clone)]
pub struct OutputData<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    pub data: OutputDataType<'a, D, SIZE>,
    frequency: usize,
    repetitions: usize,
}

#[derive(Debug, Clone)]
pub enum OutputDataType<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    ActionObservable(ActionObservableNew),
    EnergyPlot(EnergyPlotOutputData<'a, D, SIZE>),
    DifferencePlot(DifferencePlotOutputData<'a, D, SIZE>),
    TestActionObservable(TestActionObservable),
    CorrelationData(CorrelationPlotOutputData<'a, D, SIZE>),
}

// Implementing the constructors
impl<'a, const D: usize, const SIZE: usize> OutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new_action_observable(temp: f64) -> Self {
        OutputData {
            data: OutputDataType::ActionObservable(ActionObservableNew::new(temp)),
            frequency: 1,
            repetitions: 1,
        }
    }

    pub fn new_energy_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData {
            data: OutputDataType::EnergyPlot(EnergyPlotOutputData::new(lattice)),
            frequency: 1,
            repetitions: 1,
        }
    }

    pub fn new_difference_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData {
            data: OutputDataType::DifferencePlot(DifferencePlotOutputData::new(lattice)),
            frequency: 1,
            repetitions: 1,
        }
    }

    pub fn new_test_action_observable(temp: f64) -> Self {
        OutputData {
            data: OutputDataType::TestActionObservable(TestActionObservable::new(temp)),
            frequency: 1,
            repetitions: 1,
        }
    }

    pub fn new_correlation_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData {
            data: OutputDataType::CorrelationData(CorrelationPlotOutputData::new(lattice)),
            frequency: 1,
            repetitions: 1,
        }
    }

    /// Sets the steps that are ignored until the next observation
    pub fn set_frequency(mut self, frequency: usize) -> Self {
        self.frequency = frequency;
        self
    }

    /// Sets the numbers of observations done when a observation occurs
    pub fn set_repetitions(mut self, repetitions: usize) -> Self {
        self.repetitions = repetitions;
        self
    }

    /// Returns true if the given step allows the next update.
    pub fn is_step_for_next_update(&self, step: usize) -> bool {
        step % self.frequency == 0
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for OutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, rng: &mut ThreadRng)
    where
        [(); D * 2_usize]:,
    {
        for _rep in 0..self.repetitions {
            self.data.update(field, rng);
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for OutputDataType<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, rng: &mut ThreadRng)
    where
        [(); D * 2_usize]:,
    {
        match self {
            OutputDataType::ActionObservable(obs) => {
                <ActionObservableNew as UpdateOutputData<D, SIZE>>::update(obs, field, rng)
            }
            OutputDataType::TestActionObservable(obs) => {
                <TestActionObservable as UpdateOutputData<D, SIZE>>::update(obs, field, rng)
            }
            OutputDataType::EnergyPlot(plt) => {
                <EnergyPlotOutputData<'a, D, SIZE> as UpdateOutputData<D, SIZE>>::update(
                    plt, field, rng,
                )
            }
            OutputDataType::DifferencePlot(plt) => {
                <DifferencePlotOutputData<'a, D, SIZE> as UpdateOutputData<D, SIZE>>::update(
                    plt, field, rng,
                )
            }
            OutputDataType::CorrelationData(dat) => {
                <CorrelationPlotOutputData<'a, D, SIZE> as UpdateOutputData<D, SIZE>>::update(
                    dat, field, rng,
                )
            }
        }
    }
}

/// This struct implements the action as a observable of the lattice simulation
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

impl Observe for ActionObservableNew {
    fn result(self) -> f64 {
        self.value.mean()
    }
}

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

impl Observe for TestActionObservable {
    fn result(self) -> f64 {
        self.obs.mean() / self.partfn.mean()
    }
}

/// This struct implements the calculation of data that then can be used to
/// plot the energy density of the lattice simulation.
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

impl<'a, const SIZE: usize> Plotting<FieldExport3d<f64>> for EnergyPlotOutputData<'a, 3, SIZE> {
    fn plot(self) -> Vec<Vec<FieldExport3d<f64>>> {
        let mut ary: Vec<Vec<FieldExport3d<f64>>> = Vec::new();
        for field in self.bonds.into_sub_fields().into_iter() {
            ary.push(Field::<'a, f64, 3, SIZE>::from_field(field).into_export());
        }
        ary
    }
}

/// This struct implements the calculation of data that then can be used to
/// plot the field strength of the lattice simulation.
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

impl<'a, const SIZE: usize> Plotting<FieldExport3d<f64>> for DifferencePlotOutputData<'a, 3, SIZE> {
    fn plot(self) -> Vec<Vec<FieldExport3d<f64>>> {
        let mut ary: Vec<Vec<FieldExport3d<f64>>> = Vec::new();
        for field in self.bonds.into_sub_fields().into_iter() {
            ary.push(Field::<'a, f64, 3, SIZE>::from_field(field).into_export());
        }
        ary
    }
}

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
    corr_fn: Vec<KahanSummation<f64>>,
}

impl<'a, const D: usize, const SIZE: usize> CorrelationPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        const T_DIM: usize = 2;
        let max_t: usize = lattice.size[T_DIM];
        let mut corr_fn: Vec<KahanSummation<f64>> = Vec::with_capacity(max_t);
        for _t in 0..max_t {
            corr_fn.push(KahanSummation::new())
        }

        CorrelationPlotOutputData { lattice, corr_fn }
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for CorrelationPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A, rng: &mut ThreadRng) {
        const T_DIM: usize = 2;
        let max_t: usize = self.lattice.size[T_DIM];

        // Index of the point on the starting timeslice (slice zero)
        let x_index: usize = self.lattice.get_random_lattice_point_index(rng);
        let x_t_coord: usize = self.lattice.calc_coords_from_index(x_index).into_array()[T_DIM];

        // Preparing the correlation function
        let mut correlation_function: Vec<KahanSummation<f64>> = Vec::with_capacity(max_t);
        for _ in 0..max_t {
            correlation_function.push(KahanSummation::new());
        }
        // For each lattice site, we calculate the difference squared and
        // associate it to a point of the correlation function.
        let x_value: T = field.get_value(x_index);
        for index in 0..SIZE {
            // Find the distance in t direction to the reference point
            let t: usize = (self.lattice.calc_coords_from_index(index).into_array()[T_DIM] + max_t
                - x_t_coord)
                % max_t;
            let diff: T = x_value - field.get_value(index);
            let squared: T = diff * diff;
            correlation_function[t].add(squared.into());
        }
        // Converting into a Vec
        let correlation_function: Vec<f64> =
            correlation_function.into_iter().map(|x| x.sum()).collect();

        for (step, corr) in self.corr_fn.iter_mut().enumerate() {
            corr.add(correlation_function[step])
        }
    }
}

impl<'a, const SIZE: usize> Plotting<f64> for CorrelationPlotOutputData<'a, 3, SIZE> {
    fn plot(self) -> Vec<Vec<f64>> {
        let mut ary: Vec<Vec<f64>> = Vec::new();
        ary.push(self.corr_fn.into_iter().map(|x| x.mean()).collect());
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
