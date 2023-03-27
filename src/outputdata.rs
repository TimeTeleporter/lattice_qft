//! This module implements data structures to save the field configuration
//! data of the Markov chain Monte Carlo in order to calculate observables
//! and create plots.

use crate::{
    fields::{
        Action, DifferenceObservableField, EnergyObservableField, HeightVariable, ObservableField,
    },
    kahan::KahanSummation,
    lattice::Lattice,
};

// - DataOutput --------------------------------------------------------------

pub trait UpdateOutputData<const D: usize, const SIZE: usize> {
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:;
}

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
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for OutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        match self {
            OutputData::Observable(obs) => obs.update(field),
            OutputData::Plot(plt) => plt.update(field),
        }
    }
}

// - OutputData - Observable --------------------------------------------------

pub trait Observe {
    fn result(self) -> f64;
}

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

impl<const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for ObservableOutputData<D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
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

// - OutputData - Plot --------------------------------------------------------

#[derive(Debug, Clone)]
pub enum PlotOutputData<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    Energy(EnergyPlotOutputData<'a, D, SIZE>),
    Difference(DifferencePlotOutputData<'a, D, SIZE>),
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
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for PlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        match self {
            PlotOutputData::Energy(energy) => energy.update(field),
            PlotOutputData::Difference(diff) => diff.update(field),
        }
    }
}

// - OutputData - Observable - Action -----------------------------------------

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
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
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

// - OutputData - Observable - Test Action ------------------------------------

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
}

impl Observe for TestActionObservable {
    fn result(self) -> f64 {
        self.obs.mean() / self.partfn.mean()
    }
}

impl<const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE> for TestActionObservable
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        let bolz: f64 = (-field.action_observable(self.temp)).exp();
        self.obs.add(field.action_observable(self.temp) * bolz);
        self.partfn.add(bolz);
    }
}

// - OutputData - Plot - Energy -----------------------------------------------

#[derive(Debug, Clone)]
pub struct EnergyPlotOutputData<'a, const D: usize, const SIZE: usize>(
    EnergyObservableField<'a, D, SIZE>,
)
where
    [(); D * 2_usize]:;

impl<'a, const D: usize, const SIZE: usize> EnergyPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        EnergyPlotOutputData(EnergyObservableField::new(lattice))
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for EnergyPlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        self.0.update(field)
    }
}

// - OutputData - Plot - Difference -------------------------------------------

#[derive(Debug, Clone)]
pub struct DifferencePlotOutputData<'a, const D: usize, const SIZE: usize>(
    DifferenceObservableField<'a, D, SIZE>,
)
where
    [(); D * 2_usize]:;

impl<'a, const D: usize, const SIZE: usize> DifferencePlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        DifferencePlotOutputData(DifferenceObservableField::new(lattice))
    }
}

impl<'a, const D: usize, const SIZE: usize> UpdateOutputData<D, SIZE>
    for DifferencePlotOutputData<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn update<T: HeightVariable<T>, A: Action<T>>(&mut self, field: &A)
    where
        [(); D * 2_usize]:,
    {
        self.0.update(field)
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
