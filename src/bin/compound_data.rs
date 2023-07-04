#![feature(drain_filter)]
#![feature(let_chains)]

use lattice_qft::{
    computation::ComputationSummary,
    export::{get_correlation_fn, CsvData},
    kahan::KahanSummation,
};

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

    // Initialize new array for summaries
    let mut summaries: Vec<ComputationSummary> = Vec::new();

    // For each uniue set of coupling constant, lattice size and algorithm,
    // average the correlation functions and correlation lengths
    while let Some(summary) = results.pop() {
        let temp: Option<f64> = summary.temp;
        let max_t: Option<usize> = summary.t;
        let index: usize = summary.index;
        let comptype: Option<String> = summary.comptype.clone();

        // Read the correlation functions
        let mut corr_fn: Option<Vec<KahanSummation<f64>>> =
            get_correlation_fn(index, lattice_qft::CORRELATION_PLOT_PATH_INCOMPLETE)
                .map_err(|err| eprintln!("Index {} no correlation function: {}", index, err))
                .ok()
                .map(|ary| {
                    ary.into_iter()
                        .map(|x| {
                            let mut kahan: KahanSummation<f64> = KahanSummation::new();
                            kahan.add(x);
                            kahan
                        })
                        .collect()
                });

        let mut corr: KahanSummation<f64> = KahanSummation::new();
        if let Some(corr12) = summary.corr12 {
            corr.add(corr12);
        } else {
            eprintln!("Index {} no correlation lenght availible!", index)
        }

        results
            .drain_filter(|entry| {
                entry.temp == temp && entry.t == max_t && entry.comptype == comptype
            })
            .for_each(|entry| {
                if let Some(corr12) = entry.corr12 {
                    corr.add(corr12);
                } else {
                    eprintln!("Index {} no correlation lenght availible!", entry.index)
                }
                if let Some(entry_corr_fn) = get_correlation_fn(entry.index, lattice_qft::CORRELATION_PLOT_PATH_INCOMPLETE)
                    .map_err(|err| {
                        eprintln!("Index {} no correlation function: {}", summary.index, err)
                    })
                    .ok() && let Some(corr_fn) = &mut corr_fn
                {
                    corr_fn
                        .iter_mut()
                        .zip(entry_corr_fn)
                        .for_each(|(corr1, corr2)| corr1.add(corr2))
                }
            });

        if let Some(corr_fn) = corr_fn {
            let corr_fn: Vec<f64> = corr_fn.into_iter().map(|kahan| kahan.mean()).collect();
            let path: &str = &(lattice_qft::COMPOUNDED_CORRELATION_PLOT_PATH_INCOMPLETE.to_owned()
                + &"correlation_"
                + &index.to_string()
                + &".csv");
            if let Err(err) = corr_fn.overwrite_csv(path) {
                eprintln!("Writing compounded correlation function: {}", err);
            }
        }

        summaries.push(summary.set_correlation_length(corr.mean()))
    }

    if let Err(err) = summaries.overwrite_csv(lattice_qft::RESULTS_COMP_PATH) {
        eprint!("{}", err);
    }
}
