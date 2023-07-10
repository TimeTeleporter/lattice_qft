use lattice_qft::{
    computation::ComputationSummary,
    export::{CsvData, OldComputationSummary},
};

fn main() {
    if let Err(err) = OldComputationSummary::fetch_csv_data(lattice_qft::RESULTS_PATH, true)
        .and_then(|old_summaries| {
            let new_summaries: Vec<ComputationSummary> = old_summaries
                .into_iter()
                .map(|old_summary| ComputationSummary {
                    index: old_summary.index,
                    d: old_summary.d,
                    size: old_summary.size,
                    x: old_summary.x,
                    y: old_summary.y,
                    t: old_summary.t,
                    temp: old_summary.temp,
                    comptype: old_summary.comptype,
                    comptime: old_summary.comptime.map(|x| x.round().abs() as u64),
                    action: old_summary.action,
                    energy_data: old_summary.energy_data,
                    difference_data: old_summary.difference_data,
                    correlation_data: old_summary.correlation_data,
                    corr12: old_summary.corr12,
                    corr12_err: None,
                })
                .collect();
            Ok(new_summaries)
        })
        .and_then(|new_summaries| new_summaries.overwrite_csv(lattice_qft::RESULTS_PATH))
    {
        eprintln!("Summary conversion error: {}", err);
        return;
    };
}
