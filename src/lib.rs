use std::ffi::CStr;

use edlib_sys::edlibAlignmentToCigar;

pub mod edlib_sys;

#[derive(Debug)]
pub struct EdlibAlignResult {
    edlib_res: edlib_sys::EdlibAlignResult,
    cigar_str: Option<*const i8>,
}

impl EdlibAlignResult {
    pub fn new(res: edlib_sys::EdlibAlignResult) -> Self {
        EdlibAlignResult {
            edlib_res: res,
            cigar_str: None,
        }
    }

    pub fn cigar_str(&mut self, eqx: bool) -> Option<&'_ CStr> {
        if self.edlib_res.alignment.is_null() {
            return None;
        }

        if let Some(cigar_str) = self.cigar_str {
            return unsafe { Some(CStr::from_ptr(cigar_str) )};
        }

        let cigar_fmt = if eqx {
            edlib_sys::EdlibCigarFormat_EDLIB_CIGAR_EXTENDED
        } else {
            edlib_sys::EdlibCigarFormat_EDLIB_CIGAR_STANDARD
        };

        unsafe {
            let aln_str = edlibAlignmentToCigar(
                self.edlib_res.alignment,
                self.edlib_res.alignmentLength,
                cigar_fmt,
            );
            self.cigar_str = Some(aln_str);

            Some(CStr::from_ptr(aln_str))
        }
    }
}

impl From<edlib_sys::EdlibAlignResult> for EdlibAlignResult {
    fn from(result: edlib_sys::EdlibAlignResult) -> Self {
        EdlibAlignResult::new(result)
    }
}

impl Drop for EdlibAlignResult {
    fn drop(&mut self) {
        unsafe {
            edlib_sys::edlibFreeAlignResult(self.edlib_res);
            if let Some(cigar_str) = self.cigar_str {
                libc::free(cigar_str as *mut std::ffi::c_void);
            }
        }
    }
}

pub fn edlib_align(
    query: &[u8],
    target: &[u8],
    config: edlib_sys::EdlibAlignConfig,
) -> EdlibAlignResult {
    unsafe {
        edlib_sys::edlibAlign(
            query.as_ptr() as *const i8,
            query.len() as i32,
            target.as_ptr() as *const i8,
            target.len() as i32,
            config,
        )
        .into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edlib_align() {
        let query = b"ACGT";
        let target = b"ACGT";
        let config = edlib_sys::EdlibAlignConfig {
            k: 0,
            mode: edlib_sys::EdlibAlignMode_EDLIB_MODE_HW,
            task: edlib_sys::EdlibAlignTask_EDLIB_TASK_LOC,
            additionalEqualities: std::ptr::null(),
            additionalEqualitiesLength: 0,
        };
        let mut aln_res = edlib_align(query, target, config);
        println!("{:?}", aln_res);
        println!("{:?}", aln_res.cigar_str(true));
    }
}
