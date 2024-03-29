use anyhow::Context;
use needletail::Sequence;
use serde_json::json;
use std::path::Path;

pub mod utils;

#[derive(Debug)]
struct SeqCol {
    lengths: Vec<usize>,
    names: Vec<String>,
    sequences: Option<Vec<String>>,
}

#[allow(dead_code)]
#[derive(Debug)]
/// The configuration describing how a digest should
/// be computed.
enum DigestConfig {
    /// Compute the digest using only the required fields
    /// of "lengths", "names" and "sequences" (if present)
    RequiredOnly,
    /// Compute the digest additionally including the
    /// sorted list of (sequence name, length) pair digests
    WithSeqnameLenPairs,
}

/// The default configuration uses only the required fields
impl Default for DigestConfig {
    fn default() -> Self {
        Self::RequiredOnly
    }
}

#[allow(dead_code)]
impl SeqCol {
    /// Takes an iterator `it` over pairs of sequence names and lengths, and populates
    /// this [SeqCol] object from the contents of the iterator.
    pub fn from_sam_header<'a, I>(it: I) -> Self
    where
        I: IntoIterator<Item = (&'a [u8], usize)>,
    {
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        for (n, l) in it {
            names.push(std::str::from_utf8(n).unwrap().to_owned());
            lengths.push(l);
        }

        Self {
            lengths,
            names,
            sequences: None,
        }
    }

    pub fn from_seqcol(sc: &serde_json::Value) -> anyhow::Result<Self> {
        if let Some(seqcol) = sc.as_object() {
            let lengths = seqcol
                .get("lengths")
                .ok_or(anyhow::anyhow!("must contain lengths field"))?;
            let lengths = lengths
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as usize)
                .collect();

            let names = seqcol
                .get("names")
                .ok_or(anyhow::anyhow!("must contain names field"))?;
            let names = names
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_str().unwrap().to_owned())
                .collect();

            let sequences = seqcol.get("sequences").map(|seqs| {
                seqs.as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_str().unwrap().to_owned())
                    .collect()
            });

            Ok(Self {
                lengths,
                names,
                sequences,
            })
        } else {
            anyhow::bail!("seqcol object must be a valid JSON Object");
        }
    }

    pub fn from_fasta_file<P: AsRef<Path>>(fp: P) -> anyhow::Result<Self> {
        let mut reader = needletail::parse_fastx_file(&fp).with_context(|| {
            format!("cannot parse FASTA records from {}", &fp.as_ref().display())
        })?;
        let mut names = vec![];
        let mut lengths = vec![];
        let mut seqs = vec![];
        while let Some(record) = reader.next() {
            let seqrec = record?;
            let h = utils::sha512t24u_digest(seqrec.normalize(false).as_ref(), 24);
            seqs.push(format!("SQ.{h}"));
            names.push(std::str::from_utf8(seqrec.id())?.to_owned());
            lengths.push(seqrec.num_bases());
        }

        Ok(Self {
            lengths,
            names,
            sequences: Some(seqs),
        })
    }

    pub fn digest(&self, c: DigestConfig) -> anyhow::Result<String> {
        let mut sq_json = json!({
            "lengths" : self.lengths,
            "names" : self.names,
        });

        if let Some(v) = &self.sequences {
            sq_json["sequences"] = serde_json::Value::Array(
                v.iter()
                    .map(|x| serde_json::Value::String(x.to_string()))
                    .collect(),
            );
        };

        let mut snlp_digests = vec![];
        if let DigestConfig::WithSeqnameLenPairs = c {
            snlp_digests.reserve_exact(self.names.len());

            for (l, n) in self.lengths.iter().zip(self.names.iter()) {
                let cr = utils::canonical_rep(&json!({"length" : l, "name" : n}))?;
                snlp_digests.push(utils::sha512t24u_digest(cr.as_bytes(), 24));
            }
            snlp_digests.sort_unstable();

            sq_json["sorted_name_length_pairs"] = serde_json::Value::Array(
                snlp_digests
                    .into_iter()
                    .map(serde_json::Value::String)
                    .collect(),
            );
        }

        let mut digest_json = json!({});
        for (k, v) in sq_json.as_object().unwrap().iter() {
            let v2 = utils::canonical_rep(v)?;
            let h2 = utils::sha512t24u_digest(v2.as_bytes(), 24);
            digest_json[k] = serde_json::Value::String(h2);
        }
        let digest_str = utils::canonical_rep(&digest_json)?;
        Ok(utils::sha512t24u_digest(digest_str.as_bytes(), 24))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn from_sam_header_works() {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq1\tLN:13
@RG\tID:rg0
@PG\tID:pg0\tPN:noodles
@CO\tndls
";
        let header: noodles_sam::Header = s.parse().unwrap();

        let s = SeqCol::from_sam_header(
            header
                .reference_sequences()
                .iter()
                .map(|(k, v)| (k.as_slice(), v.length().into())),
        );
        let r = s.digest(DigestConfig::default()).unwrap();

        assert_eq!(r, "2HqWKZw8F4VY7q9sfYRM-JJ_RaMXv1eK");
    }

    #[test]
    fn from_seqcol_object() {
        let file = File::open("test_data/seqcol_obj.json").expect("can't open input seqcol file");
        let reader = BufReader::new(file);
        let sc = serde_json::from_reader(reader).unwrap();
        let s = SeqCol::from_seqcol(&sc).unwrap();
        let r = s.digest(DigestConfig::default()).unwrap();
        assert_eq!(r, "2HqWKZw8F4VY7q9sfYRM-JJ_RaMXv1eK");
    }

    #[test]
    fn from_fasta_file_works_with_default() {
        let s = SeqCol::from_fasta_file(Path::new("test_data/simple.fa")).unwrap();
        let r = s.digest(DigestConfig::default()).unwrap();
        assert_eq!(r, "E0cJxnAB5lrWXGP_JoWRNWKEDfdPUDUR");
    }

    #[test]
    fn from_fasta_file_works_with_seqname_length_pairs() {
        let s = SeqCol::from_fasta_file(Path::new("test_data/simple.fa")).unwrap();
        let r = s.digest(DigestConfig::WithSeqnameLenPairs).unwrap();
        assert_eq!(r, "bXpsYPctlKYGMvDGwmoHTUuS7ryH5miY");
    }
}
