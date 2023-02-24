use argh::FromArgs;

#[derive(FromArgs)]
/// Transcribe and translate DNA codons
pub struct Args {
    /// whether or not to transcribe before translation
    #[argh(switch, short = 'd')]
    dna: bool,
    /// whether or not to show full names of amino acids
    #[argh(switch, short = 'f')]
    full: bool,
    /// whether or not to show single letters
    #[argh(switch, short = 's')]
    single: bool,
    /// file to transcribe and translate
    #[argh(positional)]
    file: Option<String>,
}

fn main() {
    let args: Args = argh::from_env();
    let mapper = if args.dna {
        Codon::from_dna_char
    } else {
        Codon::from_rna_char
    };
    if let Some(file) = &args.file {
        let data = std::fs::read_to_string(file).expect("Failed to read file");
        let protein = Protein::new(data, mapper);
        protein.print(&args);
    };
    loop {
        let mut data = String::new();
        std::io::stdin().read_line(&mut data).unwrap();
        let protein = Protein::new(data, mapper);
        protein.print(&args);
    }
}

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct Protein {
    amino_acids: Vec<AminoAcid>,
}

impl Protein {
    pub fn new(data: String, mapper: fn(char) -> Codon) -> Self {
        let codons: Vec<Codon> = data
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(mapper)
            .collect();
        if codons.len() % 3 != 0 {
            panic!("codon list not divisible by 3!");
        }
        let amino_acids = codons
            .chunks_exact(3)
            .map(|c| AminoAcid::from_codons(c[0], c[1], c[2]))
            .collect();
        Self { amino_acids }
    }
    pub fn to_chars(&self) -> String {
        self.amino_acids
            .iter()
            .map(AminoAcid::to_codon_char)
            .collect()
    }
    pub fn to_names(&self) -> Vec<String> {
        self.amino_acids.iter().map(ToString::to_string).collect()
    }
    pub fn to_abbreviations(&self) -> Vec<String> {
        self.amino_acids
            .iter()
            .map(AminoAcid::to_abbreviation)
            .collect()
    }
    pub fn print(&self, args: &Args) {
        if args.full {
            let mut started = false;
            for acid in self.to_names() {
                if !started {
                    started = true;
                } else {
                    print!(", ");
                }
                print!("{acid}");
            }
            println!();
        } else if args.single {
            println!("{}", self.to_chars());
        } else {
            let mut started = false;
            for acid in self.to_abbreviations() {
                if !started {
                    started = true;
                } else {
                    print!(", ");
                }
                print!("{acid}");
            }
            println!();
        }
    }
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub enum Codon {
    A,
    U,
    C,
    G,
}

impl Codon {
    fn from_rna_char(value: char) -> Self {
        match value.to_ascii_uppercase() {
            'A' => Codon::A,
            'U' => Codon::U,
            'C' => Codon::C,
            'G' => Codon::G,
            _ => panic!("Value {value} is not a codon!"),
        }
    }
    fn from_dna_char(value: char) -> Self {
        match value.to_ascii_uppercase() {
            'A' => Codon::U,
            'T' => Codon::A,
            'C' => Codon::G,
            'G' => Codon::C,
            _ => panic!("Value {value} is not a codon!"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum AminoAcid {
    Alanine,
    Glycine,
    Methionine,
    Serine,
    Cysteine,
    Hisitidine,
    Asparagine,
    Threonine,
    AsparticAcid,
    Isoleucine,
    Proline,
    Valine,
    GlutamicAcid,
    Lysine,
    Glutamine,
    Tryptophan,
    Phenylalanine,
    Leucine,
    Arginine,
    Tyrosine,
    Stop,
}

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let out = match self {
            Self::Alanine => "Alanine",
            Self::Glycine => "Glycine",
            Self::Methionine => "Methionine",
            Self::Serine => "Serine",
            Self::Cysteine => "Cysteine",
            Self::Hisitidine => "Hisitidine",
            Self::Asparagine => "Asparagine",
            Self::Threonine => "Threonine",
            Self::AsparticAcid => "Aspartic acid",
            Self::Isoleucine => "Isoleucine",
            Self::Proline => "Proline",
            Self::Valine => "Valine",
            Self::GlutamicAcid => "Glutamic acid",
            Self::Lysine => "Lysine",
            Self::Glutamine => "Glutamine",
            Self::Tryptophan => "Tryptophan",
            Self::Phenylalanine => "Phenylalanine",
            Self::Leucine => "Leucine",
            Self::Arginine => "Arginine",
            Self::Tyrosine => "Tyrosine",
            Self::Stop => "■",
        };
        f.write_str(out)
    }
}

impl AminoAcid {
    pub fn to_abbreviation(&self) -> String {
        match self {
            Self::Alanine => "Ala",
            Self::Glycine => "Gly",
            Self::Methionine => "Met",
            Self::Serine => "Ser",
            Self::Cysteine => "Cys",
            Self::Hisitidine => "His",
            Self::Asparagine => "Asn",
            Self::Threonine => "Thr",
            Self::AsparticAcid => "Asp",
            Self::Isoleucine => "Ile",
            Self::Proline => "Pro",
            Self::Valine => "Val",
            Self::GlutamicAcid => "Glu",
            Self::Lysine => "Lys",
            Self::Glutamine => "Gln",
            Self::Tryptophan => "Trp",
            Self::Phenylalanine => "Phe",
            Self::Leucine => "Leu",
            Self::Arginine => "Arg",
            Self::Tyrosine => "Tyr",
            Self::Stop => "■",
        }
        .to_string()
    }
    pub fn to_codon_char(&self) -> char {
        match self {
            AminoAcid::Alanine => 'A',
            AminoAcid::Glycine => 'G',
            AminoAcid::Methionine => 'M',
            AminoAcid::Serine => 'S',
            AminoAcid::Cysteine => 'C',
            AminoAcid::Hisitidine => 'H',
            AminoAcid::Asparagine => 'N',
            AminoAcid::Threonine => 'T',
            AminoAcid::AsparticAcid => 'D',
            AminoAcid::Isoleucine => 'I',
            AminoAcid::Proline => 'P',
            AminoAcid::Valine => 'V',
            AminoAcid::GlutamicAcid => 'E',
            AminoAcid::Lysine => 'K',
            AminoAcid::Glutamine => 'Q',
            AminoAcid::Tryptophan => 'W',
            AminoAcid::Phenylalanine => 'F',
            AminoAcid::Leucine => 'L',
            AminoAcid::Arginine => 'R',
            AminoAcid::Tyrosine => 'Y',
            AminoAcid::Stop => '■',
        }
    }
    pub fn from_codons(first: Codon, second: Codon, third: Codon) -> Self {
        match first {
            Codon::A => Self::from_codons_a(second, third),
            Codon::G => Self::from_codons_g(second, third),
            Codon::C => Self::from_codons_c(second, third),
            Codon::U => Self::from_codons_u(second, third),
        }
    }
    fn from_codons_a(second: Codon, third: Codon) -> Self {
        match second {
            Codon::A => Self::from_codons_aa(third),
            Codon::G => Self::from_codons_ag(third),
            Codon::C => Self::Threonine,
            Codon::U => Self::from_codons_au(third),
        }
    }
    fn from_codons_aa(third: Codon) -> Self {
        match third {
            Codon::A | Codon::G => Self::Lysine,
            Codon::U | Codon::C => Self::Asparagine,
        }
    }
    fn from_codons_ag(third: Codon) -> Self {
        match third {
            Codon::A | Codon::G => Self::Arginine,
            Codon::C | Codon::U => Self::Serine,
        }
    }
    fn from_codons_au(third: Codon) -> Self {
        match third {
            Codon::A | Codon::C | Codon::U => Self::Isoleucine,
            Codon::G => Self::Methionine,
        }
    }
    fn from_codons_g(second: Codon, third: Codon) -> Self {
        match second {
            Codon::A => Self::from_codons_ga(third),
            Codon::G => Self::Glycine,
            Codon::C => Self::Alanine,
            Codon::U => Self::Valine,
        }
    }
    fn from_codons_ga(third: Codon) -> Self {
        match third {
            Codon::A | Codon::G => Self::GlutamicAcid,
            Codon::U | Codon::C => Self::AsparticAcid,
        }
    }
    fn from_codons_c(second: Codon, third: Codon) -> Self {
        match second {
            Codon::A => Self::from_codons_ca(third),
            Codon::G => Self::Arginine,
            Codon::C => Self::Proline,
            Codon::U => Self::Leucine,
        }
    }
    fn from_codons_ca(third: Codon) -> Self {
        match third {
            Codon::A | Codon::G => Self::Glutamine,
            Codon::U | Codon::C => Self::Hisitidine,
        }
    }
    fn from_codons_u(second: Codon, third: Codon) -> Self {
        match second {
            Codon::A => Self::from_codons_ua(third),
            Codon::G => Self::from_codons_ug(third),
            Codon::C => Self::Serine,
            Codon::U => Self::from_codons_uu(third),
        }
    }
    fn from_codons_ua(third: Codon) -> Self {
        match third {
            Codon::A | Codon::G => Self::Stop,
            Codon::U | Codon::C => Self::Tyrosine,
        }
    }
    fn from_codons_ug(third: Codon) -> Self {
        match third {
            Codon::A => Self::Stop,
            Codon::G => Self::Tryptophan,
            Codon::C | Codon::U => Self::Cysteine,
        }
    }
    fn from_codons_uu(third: Codon) -> Self {
        match third {
            Codon::A | Codon::G => Self::Leucine,
            Codon::U | Codon::C => Self::Phenylalanine,
        }
    }
}
