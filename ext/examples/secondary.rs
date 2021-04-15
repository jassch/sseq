use algebra::module::homomorphism::ModuleHomomorphism;
use ext::secondary::*;
use ext::utils::construct;

fn main() -> error::Result<()> {
    let module_file_name: String = query::with_default("Module", "S_2", Ok);

    let max_s = query::with_default("Max s", "7", Ok);
    let max_f = query::with_default("Max f", "30", Ok);

    let res_save_file: Option<String> = query::optional("Resolution save file", Ok);
    #[cfg(feature = "concurrent")]
    let del_save_file: Option<String> = query::optional("Delta save file", Ok);

    #[cfg(feature = "concurrent")]
    let bucket = {
        let num_threads = query::with_default("Number of threads", "2", Ok);
        thread_token::TokenBucket::new(num_threads)
    };

    let resolution = construct(
        (&*module_file_name, algebra::AlgebraType::Milnor),
        res_save_file.as_deref(),
    )?;

    if !can_compute(&resolution) {
        eprintln!("Cannot compute d2 for the module {}", module_file_name);
        return Ok(());
    }

    #[cfg(not(feature = "concurrent"))]
    let deltas = compute_delta(&resolution, max_s, max_f);

    #[cfg(feature = "concurrent")]
    let deltas = compute_delta_concurrent(&resolution, max_s, max_f, &bucket, del_save_file);

    for f in 1..=max_f {
        for s in 1..(max_s - 1) {
            let t = s as i32 + f;
            let delta = &deltas[s as usize - 1];

            if delta.source().number_of_gens_in_degree(t + 1) == 0 {
                continue;
            }
            let d = delta.hom_k(t);

            for (i, entry) in d.into_iter().enumerate() {
                println!("d_2 x_({}, {}, {}) = {:?}", f, s, i, entry);
            }
        }
    }
    Ok(())
}
