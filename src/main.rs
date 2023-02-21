use std::fmt::Display;

fn main() {
    loop {
        let mut sin = String::new();
        std::io::stdin().read_line(&mut sin).unwrap();
        let mapper = match std::env::args().nth(1) {
            None => Codon::from_dna_char,
            Some(v) => match v.as_str() {
                "--dna" => Codon::from_dna_char,
                "--rna" => Codon::from_rna_char,
                _ => panic!("Argument 1 must be `--dna` or `--rna`, got {v}!"),
            },
        };
        let codons: Vec<Codon> = sin
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(mapper)
            .collect();
        let protein: Vec<AminoAcid> = codons
            .chunks_exact(3)
            .map(|c| AminoAcid::from_codons(c[0], c[1], c[2]))
            .collect();
        println!("{protein:?}");
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

#[derive(Debug, PartialEq, Eq, Hash)]
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

impl Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let out = match self {
            AminoAcid::Alanine => "Ala",
            AminoAcid::Glycine => "Gly",
            AminoAcid::Methionine => "Met",
            AminoAcid::Serine => "Ser",
            AminoAcid::Cysteine => "Cys",
            AminoAcid::Hisitidine => "His",
            AminoAcid::Asparagine => "Asn",
            AminoAcid::Threonine => "Thr",
            AminoAcid::AsparticAcid => "Asp",
            AminoAcid::Isoleucine => "Ile",
            AminoAcid::Proline => "Pro",
            AminoAcid::Valine => "Val",
            AminoAcid::GlutamicAcid => "Glu",
            AminoAcid::Lysine => "Lys",
            AminoAcid::Glutamine => "Gln",
            AminoAcid::Tryptophan => "Trp",
            AminoAcid::Phenylalanine => "Phe",
            AminoAcid::Leucine => "Leu",
            AminoAcid::Arginine => "Arg",
            AminoAcid::Tyrosine => "Tyr",
            AminoAcid::Stop => "■",
        };
        f.write_str(out)
    }
}

impl AminoAcid {
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
