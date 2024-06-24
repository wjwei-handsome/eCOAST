mod bedgraph;
mod cli;
mod gff;
use anyhow::Result;
use bedgraph::{BedGraph, DefaultReadFilter, OnlyDepthProcessor};
use clap::Parser;
use gff::Gene;
use indicatif::ParallelProgressIterator;
use linfa::traits::Transformer;
use linfa_hierarchical::HierarchicalCluster;
use linfa_kernel::{Kernel, KernelMethod};
use polars::lazy::dsl::*;
use polars::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rust_lapper::Lapper;
use std::io::{Error, Write};
use std::thread;
use std::{collections::HashMap, f64::consts::E, path::PathBuf};

const HARD_MIN_COV: i64 = 10;

fn format_df(df: DataFrame) -> Result<(Option<DataFrame>, i64)> {
    // add ln(value)
    let df = df
        .lazy()
        .with_column((col("value").log(E)).alias("log_value"))
        .collect()?;

    // add length = end - start
    let df = df
        .lazy()
        .with_column((col("end") - col("start")).alias("length"))
        .collect()?;

    // convert -Inf to 0 in `log_value`
    let df = df
        .lazy()
        .with_column(
            when(col("log_value").eq(lit(f64::NEG_INFINITY)))
                .then(lit(0.0))
                .otherwise(col("log_value"))
                .alias("log_value"),
        )
        .collect()?;

    // sum length
    let exon_length = df.column("length")?.i64()?.sum().unwrap_or(0);

    Ok((Some(df), exon_length))
}

fn df_from_bdg(bdgs: Vec<BedGraph>) -> Result<(Option<DataFrame>, i64)> {
    // if bdgs is empty
    if bdgs.is_empty() {
        return Ok((None, 0));
    }

    // convert bdgs to df
    let df = DataFrame::new(vec![
        Series::new(
            "chrom",
            bdgs.iter().map(|x| x.ref_seq.clone()).collect::<Vec<_>>(),
        ),
        Series::new(
            "start",
            bdgs.iter().map(|x| x.pos as i64).collect::<Vec<_>>(),
        ),
        Series::new("end", bdgs.iter().map(|x| x.end as i64).collect::<Vec<_>>()),
        Series::new(
            "value",
            bdgs.iter().map(|x| x.depth as i64).collect::<Vec<_>>(),
        ),
    ])?;

    format_df(df)
}

fn df_from_bdgfile(bdgfile: &str) -> Result<(Option<DataFrame>, i64)> {
    let csv_parse_option = CsvParseOptions::default().with_separator(b'\t');
    let mut df = CsvReadOptions::default()
        .with_parse_options(csv_parse_option)
        .with_has_header(false)
        .try_into_reader_with_file_path(Some(PathBuf::from(bdgfile)))?
        .finish()?;

    // if df is empty
    if df.height() == 0 {
        return Ok((None, 0));
    }

    // rename columns
    let colnames = vec!["chrom", "start", "end", "value"];
    df.set_column_names(&colnames)?;

    format_df(df)
}

fn cluster(df: DataFrame) -> anyhow::Result<(DataFrame, Vec<usize>)> {
    // reshape(-1,1) in log_value_array
    let log_value_array = df.column("log_value")?.f64()?;
    let array_len = log_value_array.len();
    let log_value_array = log_value_array.to_ndarray()?.into_shape((array_len, 1))?;

    // 创建高斯核函数
    let kernel = Kernel::params()
        .method(KernelMethod::Gaussian(1.0))
        .transform(log_value_array);

    // 使用 linfa 进行聚类分析，指定距离阈值
    let distance_threshold = 0.5;
    let model = HierarchicalCluster::default()
        .max_distance(distance_threshold)
        .transform(kernel)?;
    let label_array = model.targets();
    Ok((df, label_array.to_owned()))
}

fn convert_label(df: DataFrame, label_array: &[usize]) -> anyhow::Result<DataFrame> {
    // trans label_array to u32
    let mut u32_label_array = Vec::new();
    for label in label_array {
        u32_label_array.push(*label as u32);
    }
    // make label column
    let label_series = Series::new("label", u32_label_array);
    // add label column to df
    let df = df.lazy().with_column(label_series.lit()).collect()?;

    // gen rank column
    let k = df.column("label")?.u32()?.max().unwrap_or(0) + 1;
    let mut mean_values = HashMap::new();
    // count mean log_value for each cluster
    for i in 0..k {
        let mean_value = df
            .clone()
            .lazy()
            .filter(col("label").eq(lit(i)))
            .collect()?
            .column("log_value")?
            .f64()?
            .mean()
            .unwrap_or(0.0);
        mean_values.insert(i, mean_value);
    }
    // sort mean_values by value
    let mut mean_values_vec: Vec<_> = mean_values.iter().collect();
    mean_values_vec.sort_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal));
    let rank: Vec<_> = mean_values_vec.iter().map(|(k, _)| k).collect();
    // generate rank vec
    let mut rank_vec = Vec::new();
    for label in df.column("label")?.u32()?.into_no_null_iter() {
        rank_vec.push(rank.iter().position(|&x| **x == label).unwrap() as u32);
    }
    // add rank column
    let rank_series = Series::new("rank", rank_vec);
    let df = df.lazy().with_column(rank_series.lit()).collect()?;

    // stdout df
    // let stdout = std::io::stdout();
    // let mut csvwtr = CsvWriter::new(stdout.lock());
    // csvwtr.finish(&mut df)?;

    Ok(df)
}

fn compute_ratio(df: DataFrame, exon_length: i64) -> anyhow::Result<f64> {
    // compute high ratio by drop lowest rank
    // judge lowest rank's min value
    let low_rank_min_value = df
        .clone()
        .lazy()
        .filter(col("rank").eq(lit(0)))
        .collect()?
        .column("value")?
        .i64()?
        .min()
        .unwrap_or(0);
    let df = if low_rank_min_value > HARD_MIN_COV {
        // nothing
        df
    } else {
        // drop lowest rank
        df.lazy().filter(col("rank").neq(lit(0))).collect()?
    };
    // compute high ratio by summary length
    let high_length = df.column("length")?.i64()?.sum().unwrap_or(0);
    let high_ratio = high_length as f64 / exon_length as f64;
    Ok(high_ratio)
}

// aux function to init a nested vec
// fn init_vec<T>() -> Vec<Vec<T>> {
//     vec![Vec::new(), Vec::new()]
// }

// wrap process_sample_gene
fn process_gene_to_results(
    gene: &Gene,
    processer: Arc<OnlyDepthProcessor<DefaultReadFilter>>,
) -> Result<f64> {
    let mut gene_bdg_vec = Vec::new();
    let tid = gene.tid.clone();

    // merge overlapping exons
    let mut exon_lapper = Lapper::new(gene.exons.clone());
    exon_lapper.merge_overlaps();

    for iv in exon_lapper.intervals {
        let exon_start = iv.start;
        let exon_stop = iv.stop;
        let res = processer.process_region(&tid, exon_start as u32, exon_stop as u32)?;
        gene_bdg_vec.extend(res);
    }

    let (df, exon_length) = df_from_bdg(gene_bdg_vec)?;
    let df = match df {
        Some(df) => df,
        None => {
            return Ok(0.0);
        }
    };
    let (df, label_array) = cluster(df)?;
    let df = convert_label(df, &label_array)?;
    let ratio = compute_ratio(df, exon_length)?;
    Ok(ratio)
}

fn process_bdgs_to_result(bdg_file: &str) -> Result<f64> {
    let (df, exon_length) = df_from_bdgfile(bdg_file)?;
    let df = match df {
        Some(df) => df,
        None => {
            return Ok(0.0);
        }
    };
    let (df, label_array) = cluster(df)?;
    let df = convert_label(df, &label_array)?;
    let ratio = compute_ratio(df, exon_length)?;
    Ok(ratio)
}

fn process_gene_to_bdgs(gene: Gene, bam_path: &str, outdir: &str) -> anyhow::Result<()> {
    let read_filter = DefaultReadFilter::new(0, 0, 0);
    let bam_path = PathBuf::from(bam_path); // check it
    let depth_processer = OnlyDepthProcessor::new(bam_path, 0, read_filter);
    // write to file
    let outpath = format!("{}/{}.bedgraph", outdir, gene.name);
    let mut file = std::fs::File::create(outpath)?;
    let mut gene_bdg_vec = Vec::new();
    // merge overlapping exons
    let mut exon_lapper = Lapper::new(gene.exons);
    exon_lapper.merge_overlaps();
    for iv in exon_lapper.intervals {
        let exon_start = iv.start;
        let exon_stop = iv.stop;
        let res = depth_processer.process_region(&gene.tid, exon_start as u32, exon_stop as u32)?;
        gene_bdg_vec.extend(res);
    }
    for bdg in gene_bdg_vec {
        writeln!(file, "{}", bdg)?;
    }
    Ok(())
}

fn main() -> Result<()> {
    // parse cli
    let args = cli::Cli::parse();
    let threads = args.threads;
    let _verbose = args.verbose;

    // set rayon threads
    let bigger_stack_size = 8 * 1024 * 1024; //8M

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(bigger_stack_size)
        .build_global()?;

    match args.command {
        cli::Commands::MakeBedGraph { gff, bam, out } => {
            // 1. parse gff3
            let (genes, gene_count) = gff::parse_gff3(&gff)?;
            // 2. parallel process genes
            genes
                .into_par_iter()
                .progress_count(gene_count as u64)
                .for_each(|gene| {
                    let _ = process_gene_to_bdgs(gene, &bam, &out);
                });
        }
        cli::Commands::GenRatio {
            bedgraph_dir,
            sample_name,
            out,
        } => {
            // 1. get bedgraph files list
            let bedgraph_files = std::fs::read_dir(bedgraph_dir)?
                .map(|res| res.map(|e| e.path()))
                .collect::<Result<Vec<PathBuf>, Error>>()?;
            let bedgraph_lens = bedgraph_files.len();
            // 2. parallel process bedgraph files
            let ratio_res = bedgraph_files
                .into_par_iter()
                .progress_count(bedgraph_lens as u64)
                .map(|bdg_file| {
                    let bdg_file = bdg_file.to_str().unwrap().to_string();
                    let ratio = process_bdgs_to_result(&bdg_file).unwrap();
                    (bdg_file, ratio)
                })
                .collect::<HashMap<_, _>>();
            // 3. write to file
            let mut file = std::fs::File::create(out)?;
            // header
            writeln!(file, "gene\t{}", sample_name)?;
            for (bdg_file, ratio) in ratio_res {
                let gene_name = bdg_file
                    .split('/')
                    .last()
                    .unwrap()
                    .split('.')
                    .next()
                    .unwrap();
                writeln!(file, "{}\t{}", gene_name, ratio)?;
            }
        }
        cli::Commands::Run {
            gff,
            bam,
            out,
            sample_name,
        } => {
            // 1. parse gff3
            let (genes, gene_count) = gff::parse_gff3(&gff)?;
            // 2. parallel run each gene
            //  // Create the read filter
            let read_filter = DefaultReadFilter::new(0, 0, 0);
            let bam_path = PathBuf::from(bam); // check it
            let depth_processer = OnlyDepthProcessor::new(bam_path, 0, read_filter);
            let depth_processer = Arc::new(depth_processer);

            // // Process each gene
            let ratio_res = thread::spawn(move || {
                genes
                    .into_par_iter()
                    .progress_count(gene_count as u64)
                    .map(|gene| {
                        let depth_processer = Arc::clone(&depth_processer);
                        let ratio = process_gene_to_results(&gene, depth_processer).unwrap();
                        (gene.name, ratio)
                    })
                    .collect::<HashMap<_, _>>()
            });

            let ratio_res = ratio_res.join().unwrap();

            // 3. write output
            let mut file = std::fs::File::create(out)?;
            // header
            writeln!(file, "gene\t{}", sample_name)?;
            for (gene_name, ratio) in ratio_res {
                writeln!(file, "{}\t{}", gene_name, ratio)?;
            }
        }
    }
    Ok(())
}
