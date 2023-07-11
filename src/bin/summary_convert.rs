use lattice_qft::{
    computation::ComputationSummary,
    export::{CsvData, OldComputationSummary},
};

fn main() {
    if let Err(err) = OldComputationSummary::fetch_csv_data(lattice_qft::RESULTS_PATH, true)
        .and_then(|old_summaries| {
            let new_summaries: Vec<ComputationSummary> = old_summaries
                .into_iter()
                .map(|old_summary| convert_computation_summary(old_summary))
                .collect();
            Ok(new_summaries)
        })
        .and_then(|new_summaries| new_summaries.overwrite_csv(lattice_qft::RESULTS_PATH))
    {
        eprintln!("Summary conversion error: {}", err);
        return;
    };
}

fn convert_computation_summary(old_summary: OldComputationSummary) -> ComputationSummary {
    // We parse the computation summary in order to extract the iterations.
    let iterations = old_summary
        .comptype
        .clone()
        .map(|comptype| {
            comptype.split_whitespace().last().map(|in_brackets| {
                in_brackets
                    .trim_start_matches('(')
                    .trim_end_matches(')')
                    .parse::<u64>()
                    .ok()
            })
        })
        .flatten()
        .flatten();
    ComputationSummary {
        index: old_summary.index,
        d: old_summary.d,
        size: old_summary.size,
        x: old_summary.x,
        y: old_summary.y,
        t: old_summary.t,
        temp: old_summary.temp,
        comptype: old_summary.comptype,
        iterations,
        comptime: old_summary.comptime,
        action: old_summary.action,
        energy_data: old_summary.energy_data,
        difference_data: old_summary.difference_data,
        correlation_data: old_summary.correlation_data,
        corr12: old_summary.corr12,
        corr12_err: old_summary.corr12_err,
    }
}

#[test]
fn test_conversion() {
    let old_summary = OldComputationSummary {
        index: 300000,
        d: None,
        size: None,
        x: None,
        y: None,
        t: None,
        temp: None,
        comptype: Some("Metropolis Simulation (800000)".to_owned()),
        comptime: None,
        action: None,
        energy_data: false,
        difference_data: false,
        correlation_data: false,
        corr12: None,
        corr12_err: None,
    };
    let new_summary = convert_computation_summary(old_summary);
    dbg!(new_summary);
}
