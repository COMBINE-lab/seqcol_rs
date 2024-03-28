use sha2::{Sha512, Digest};
use base64::{Engine as _, prelude::BASE64_URL_SAFE};
use itertools::Itertools;
use serde_json::json;

#[derive(Debug)]
struct SeqCol {
    lengths: Vec<usize>,
    names: Vec<String>,
    sequences: Vec<String>
}

impl SeqCol {
    pub fn from_sam_header<'a, I>(it: I) -> Self
        where 
        I: IntoIterator<Item=(&'a [u8], usize)>
    {
        let mut names = Vec::new();
        let mut lengths= Vec::new();
        for (n, l) in it {
            names.push(std::str::from_utf8(n).unwrap().to_owned());
            lengths.push(l);
        }

        Self {
            lengths,
            names,
            sequences: vec![]
        }
    }

    pub fn digest(&self) -> anyhow::Result<String> {
        let sq_json = json!({
            "lengths" : self.lengths,
            "names" : self.names,
            "sequences" : self.sequences
        });

        let mut digest_json = json!({});
        for (k, v) in sq_json.as_object().unwrap().iter() {
            let v2 = canonical_rep(v)?;
            let h2 = sha512t24u_digest(v2.as_bytes(), 24);
            digest_json[k] = serde_json::Value::String(h2);
        }
        let digest_str = canonical_rep(&digest_json)?;
        Ok(sha512t24u_digest(digest_str.as_bytes(), 24))

    }
}

pub fn sha512t24u_digest(attr: &[u8], offset: usize) -> String {
    println!("adding {}", std::str::from_utf8(attr).unwrap());
    let digest = Sha512::digest(attr);
    println!("digest = {:?}", &digest[..offset]);
    let tdigest_b64us = BASE64_URL_SAFE.encode(&digest[..offset]);
    tdigest_b64us
}

pub fn canonical_rep(inp: &serde_json::Value) -> anyhow::Result<String> {
    if let Some(om) = inp.as_object() {
        let mut output_obj = vec![];//serde_json::Object::new();
        for k in om.keys().sorted() {
            output_obj.push((k, om[k].clone()))
        }
        let output_val: serde_json::Value = output_obj.into_iter().collect();
        Ok(serde_json::to_string(&output_val)?)
    } else {
        Ok(serde_json::to_string(inp)?)
    }
}

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam;
    use serde_json::json;

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
            header.reference_sequences().iter().map(|(k, v)| (k.as_slice(), v.length().into() )) 
        );

        println!("s = {s:?}");
        let lengths = json!({ "lengths" : s.lengths });
        let names = json!({ "names" : s.names });
        let seqs = json!({ "sequences" : s.sequences });

        let v = canonical_rep(&lengths["lengths"]).unwrap();
        println!("{}", v);
        let h = sha512t24u_digest(&v.as_bytes(), 24);
        println!("{h}");
        let v = canonical_rep(&names["names"]).unwrap();
        let h2 = sha512t24u_digest(&v.as_bytes(), 24);
        println!("{h2}");

        let v = canonical_rep(&seqs["sequences"]).unwrap();
        let h3 = sha512t24u_digest(&v.as_bytes(), 24);
        println!("{h3}");

        let j = json!({
            "lengths" : h,
            "names" : h2,
            "sequences" : h3
        });

        let v = canonical_rep(&j).unwrap();
        println!("{v}");
        let h4 = sha512t24u_digest(&v.as_bytes(), 24);
        println!("{h4}");

        let r = s.digest().unwrap();
        println!("digest = {r}");
    }


    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
