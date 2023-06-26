#![feature(drain_filter)]

use lattice_qft::{computation::ComputationSummary, export::CsvData, kahan::KahanSummation};

fn main() {
    // Read the result data
    let mut results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(lattice_qft::RESULTS_PATH, true) {
            Ok(summary) => summary,
            Err(err) => {
                eprint!("{}", err);
                return;
            }
        };

    let mut comp: Vec<ComputationSummary> = Vec::new();

    while let Some(summary) = results.pop() {
        let temp: Option<f64> = summary.temp;
        let max_t: Option<usize> = summary.t;

        let mut corr: KahanSummation<f64> = KahanSummation::new();
        if let Some(corr12) = summary.corr12 {
            corr.add(corr12);
        } else {
            eprintln!("Index {} no correlation lenght availible!", summary.index)
        }

        results
            .drain_filter(|summary| summary.temp == temp && summary.t == max_t)
            .for_each(|summary| {
                if let Some(corr12) = summary.corr12 {
                    corr.add(corr12);
                } else {
                    eprintln!("Index {} no correlation lenght availible!", summary.index)
                }
            });

        comp.push(summary.set_correlation_length(corr.mean()))
    }

    if let Err(err) = comp.overwrite_csv(lattice_qft::RESULTS_COMP_PATH) {
        eprint!("{}", err);
    }
}
