use clap::{command, ArgAction, Parser, Subcommand};

#[derive(Parser)]
#[command(name = "tmocliname")]
#[command(about = "descp")]
#[command(long_about = "long_about todo!!!")]
#[command(author, version)]
#[command(
help_template =
"{name} -- {about}\n\nVersion: {version}\n\nAuthors: {author}\
    \n\n{usage-heading} {usage}\n\n{all-args}"
) // change template more!
]
pub struct Cli {
    /// Threads, default 1
    #[arg(long, short, global = true, default_value = "1", help_heading = Some("GLOBAL"))]
    pub threads: usize,
    /// Logging level [-v: Info, -vv: Debug, -vvv: Trace, defalut: Warn].
    #[arg(short, long, global = true, action = ArgAction::Count, help_heading = Some("GLOBAL"))]
    pub verbose: u8,
    /// Subcommands
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Make Exon Coverage BedGraph from gff3 and BAM
    #[command(visible_alias = "mb", name = "make-bedgraph")]
    MakeBedGraph {
        /// Input GFF3 file
        #[arg(required = true, long, short)]
        gff: String,
        /// Input BAM file with index
        #[arg(required = true, long, short)]
        bam: String,
        /// Output Directory
        #[arg(required = true, long, short)]
        out: String,
    },

    /// Generate high-ratio table by AgglomerativeClustering
    #[command(visible_alias = "gh", name = "gen-ratio")]
    GenRatio {
        /// Input BedGraph Dir
        #[arg(required = true, long, short)]
        bedgraph_dir: String,
        /// Sample name
        #[arg(required = true, long, short)]
        sample_name: String,
        /// Output file
        #[arg(required = true, long, short)]
        out: String,
    },

    /// Directly compute high-ratio table from gff3 and BAM to avoid too many bedgraph file
    #[command(visible_alias = "r", name = "run")]
    Run {
        /// Input GFF3 file
        #[arg(required = true, long, short)]
        gff: String,
        /// Input BAM file with index
        #[arg(required = true, long, short)]
        bam: String,
        /// Output file
        #[arg(required = true, long, short)]
        out: String,
        /// Sample name
        #[arg(required = true, long, short)]
        sample_name: String,
    },
}
