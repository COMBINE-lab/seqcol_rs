use crate::constants;
use base64::{Engine as _, prelude::BASE64_URL_SAFE};
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
pub fn canonical_rep(inp: &serde_json::Value) -> anyhow::Result<String> {
    let mut out = inp.clone();
    out.sort_all_objects();
    Ok(serde_json::to_string(&out)?)
}
