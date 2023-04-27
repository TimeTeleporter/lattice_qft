//! This module implements data structures to save the field configuration
//! data of the Markov chain Monte Carlo in order to calculate observables
//! and create plots.

use crate::{
    computation::FieldExport3d,
    fields::{Action, BondsFieldNew, Field, HeightVariable, WilsonField},
    kahan::KahanSummation,
    lattice::Lattice,
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
    Observable(ObservableOutputData<D, SIZE>),
    Plot(PlotOutputData<'a, D, SIZE>),
}

impl<'a, const D: usize, const SIZE: usize> OutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new_action_observable(temp: f64) -> Self {
        OutputData::Observable(ObservableOutputData::new_action_observable(temp))
    }

    pub fn new_energy_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData::Plot(PlotOutputData::new_energy_plot(lattice))
    }

    pub fn new_difference_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData::Plot(PlotOutputData::new_difference_plot(lattice))
    }

    pub fn new_test_action_observable(temp: f64) -> Self {
        OutputData::Observable(ObservableOutputData::new_test_action_observable(temp))
    }

    pub fn new_correlation_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        OutputData::Plot(PlotOutputData::new_correlation_plot(lattice))
    }
}

// - OutputData - ObservableOutputData ----------------------------------------

#[derive(Debug, Clone)]
pub enum ObservableOutputData<const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    Action(ActionObservableNew),
    TestAction(TestActionObservable),
}

impl<const D: usize, const SIZE: usize> ObservableOutputData<D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn new_action_observable(temp: f64) -> Self {
        ObservableOutputData::Action(ActionObservableNew::new(temp))
    }

    fn new_test_action_observable(temp: f64) -> Self {
        ObservableOutputData::TestAction(TestActionObservable::new(temp))
    }
}

// - OutputData - ObservableOutputData - ActionObservableNew ------------------

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

// - OutputData - TestActionObservable - TestActionObservable -----------------

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

// - OutputData - PlotOutputData ----------------------------------------------

#[derive(Debug, Clone)]
pub enum PlotOutputData<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    Energy(EnergyPlotOutputData<'a, D, SIZE>),
    Difference(DifferencePlotOutputData<'a, D, SIZE>),
    Correlation(CorrelationPlotOutputData<'a, D, SIZE>),
}

impl<'a, const D: usize, const SIZE: usize> PlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn new_energy_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        PlotOutputData::Energy(EnergyPlotOutputData::new(lattice))
    }

    fn new_difference_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        PlotOutputData::Difference(DifferencePlotOutputData::new(lattice))
    }

    fn new_correlation_plot(lattice: &'a Lattice<D, SIZE>) -> Self {
        PlotOutputData::Correlation(CorrelationPlotOutputData::new(lattice))
    }
}

// - OutputData - PlotOutputData - EnergyPlotOutputData -----------------------

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

// - OutputData - PlotOutputData - DifferencePlotOutputData -------------------

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

// - OutputData - PlotOutputData - CorrelationPlotOutputData -------------------

/// Implements averaging of time slize differences over configurations. For
/// each configuration, we first build the time slices for each index in time
/// direction. Then we sum over all differences who are the same distance
/// appart. This sum then gets saved.
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
        let max_t: usize = lattice.size[2];
        let mut correlations: Vec<KahanSummation<f64>> = Vec::with_capacity(max_t);
        for t in 0..(max_t - 1) {
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
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:;
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for OutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        match self {
            OutputData::Observable(obs) => obs.update(field),
            OutputData::Plot(plt) => plt.update(field),
        }
    }
}

impl<const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for ObservableOutputData<D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        match self {
            ObservableOutputData::Action(act) => {
                <ActionObservableNew as UpdateOutputData<D, SIZE>>::update(act, field)
            }
            ObservableOutputData::TestAction(act) => {
                <TestActionObservable as UpdateOutputData<D, SIZE>>::update(act, field)
            }
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for PlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        match self {
            PlotOutputData::Energy(energy) => energy.update(field),
            PlotOutputData::Difference(diff) => diff.update(field),
            PlotOutputData::Correlation(corr) => corr.update(field),
        }
    }
}

impl<const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for ActionObservableNew
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A)
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
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A)
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
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A) {
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
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A) {
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
    fn update<T: HeightVariable<T>, A: Action<T, D>>(&mut self, field: &A) {
        const TIME_DIMENSION: usize = 2;
        let max_t: usize = self.lattice.size[TIME_DIMENSION];
        let mut time_slices: Vec<KahanSummation<f64>> = Vec::with_capacity(max_t);
        for _slice in 0..max_t {
            time_slices.push(KahanSummation::new());
        }
        for index in 0..SIZE {
            let t: usize = self.lattice.calc_coords_from_index(index).into_array()[TIME_DIMENSION];
            time_slices[t].add(field.get_value(index).into());
        }
        let time_slices: Vec<f64> = time_slices.into_iter().map(|x| x.sum()).collect();
        // Might want to consider using parallel iter mut here
        for (step, corr) in self.correlations.iter_mut().enumerate() {
            let step: usize = step + 1; // There is no zero step
            for start in 0..max_t {
                let offset: usize = (start + step) % max_t;
                let diff: f64 = time_slices[start] - time_slices[offset];
                let value: f64 = diff * diff;
                corr.add(value);
            }
        }
    }
}

// - Observe ------------------------------------------------------------------

pub trait Observe {
    fn result(self) -> f64;
}

impl<const D: usize, const SIZE: usize> Observe for ObservableOutputData<D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn result(self) -> f64 {
        match self {
            ObservableOutputData::Action(obs) => obs.result(),
            ObservableOutputData::TestAction(obs) => obs.result(),
        }
    }
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
        ary.push(self.correlations.into_iter().map(|x| x.sum()).collect());
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
