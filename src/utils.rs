use crate::constants;
use base64::{prelude::BASE64_URL_SAFE, Engine as _};
use itertools::Itertools;
use sha2::{Digest, Sha512};

#[inline(always)]
pub fn sha512t24u_digest_default(attr: &[u8]) -> String {
    sha512t24u_digest(attr, constants::DEFAULT_DIGEST_BYTES)
}

#[inline(always)]
pub fn sha512t24u_digest(attr: &[u8], offset: usize) -> String {
    let digest = Sha512::digest(attr);
    BASE64_URL_SAFE.encode(&digest[..offset])
}

#[inline(always)]
pub fn canonical_dict(inp: &serde_json::Value) -> serde_json::Value {
    match inp {
        serde_json::Value::Object(om) => {
            let mut output_obj = Vec::with_capacity(om.len());
            for k in om.keys().sorted() {
                output_obj.push((k, om[k].clone()))
            }
            let output_val: serde_json::Value = output_obj.into_iter().collect();
            output_val
        }
        _ => inp.clone(),
    }
}

#[inline(always)]
pub fn canonical_rep(inp: &serde_json::Value) -> anyhow::Result<String> {
    match inp {
        serde_json::Value::Object(om) => {
            let mut output_obj = Vec::with_capacity(om.len());
            for k in om.keys().sorted() {
                output_obj.push((k, om[k].clone()))
            }
            let output_val: serde_json::Value = output_obj.into_iter().collect();
            Ok(serde_json::to_string(&output_val)?)
        }
        serde_json::Value::Array(ar) => {
            let mut output_obj = Vec::with_capacity(ar.len());
            for e in ar.iter() {
                output_obj.push(canonical_dict(e));
            }
            let output_val: serde_json::Value = output_obj.into_iter().collect();
            Ok(serde_json::to_string(&output_val)?)
        }
        _ => Ok(serde_json::to_string(inp)?),
    }
}
