use std::cmp;
use memchr::memchr2;
use std::io::{self, BufRead, Read};

pub struct CompactBufReader<R> {
    inner: R,
    buf: Vec<u8>,
    start: usize, // index of first unread byte
    end: usize,   // one past last unread byte
}

#[allow(dead_code)]
impl<R: Read> CompactBufReader<R> {
    pub fn with_capacity(inner: R, capacity: usize) -> Self {
        let cap = if capacity == 0 { 8 * 1024 } else { capacity };
        let buf = vec![0u8; cap];
        Self { inner, buf, start: 0, end: 0 }
    }

    pub fn new(inner: R) -> Self {
        Self::with_capacity(inner, 8 * 1024)
    }

    pub fn capacity(&self) -> usize { self.buf.len() }

    #[inline]
    pub fn available(&self) -> &[u8] {
        &self.buf[self.start..self.end]
    }

    #[inline]
    pub fn filled_len(&self) -> usize {
        self.end - self.start
    }

    #[inline]
    pub fn find_delimiter(&self) -> Option<usize> {
        memchr2(b'\t', b'\n', self.available())
    }

    pub fn ensure_delimiter(&mut self) -> io::Result<bool> {
        loop {
            if self.find_delimiter().is_some() {
                return Ok(true);
            }

            // EOF handling
            if self.buf.is_empty() {
                return Ok(false);
            }

            self.compact();

            let before = self.end;

            while self.end < self.buf.len() {
                match self.inner.read(&mut self.buf[self.end..]) {
                    Ok(0) => {
                        // EOF
                        return Ok(self.filled_len() > 0);
                    }
                    Ok(n) => {
                        self.end += n;
                        break;
                    }
                    Err(ref e) if e.kind() == io::ErrorKind::Interrupted => continue,
                    Err(e) => return Err(e),
                }
            }

            // no progress possible
            if self.end == before {
                return Ok(self.filled_len() > 0);
            }
        }
    }

    /// Inherent consume: advance by `amt`.
    pub fn consume_front(&mut self, amt: usize) {
        let to_consume = cmp::min(amt, self.filled_len());
        self.start += to_consume;
        if self.start == self.end {
            self.start = 0;
            self.end = 0;
        }
    }

    pub fn compact(&mut self) {

        if self.start == 0 {
            return;
        }
        if self.start >= self.end {
            self.start = 0;
            self.end = 0;
            return;
        }
        let len = self.end - self.start;
        self.buf.copy_within(self.start..self.end, 0);
        self.start = 0;
        self.end = len;
    }

    /// Compact unread bytes to the front, then read until full or EOF.
    pub fn fill_buffer(&mut self) -> io::Result<&[u8]> {
        if self.buf.is_empty() {
            return Ok(&[]);
        }
        self.compact();

        while self.end < self.buf.len() {
            match self.inner.read(&mut self.buf[self.end..]) {
                Ok(0) => break, // EOF
                Ok(n) => self.end += n,
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            }
        }

        Ok(&self.buf[self.start..self.end])
    }

    pub fn get_ref(&self) -> &R { &self.inner }
    pub fn get_mut(&mut self) -> &mut R { &mut self.inner }
    pub fn into_inner(self) -> R { self.inner }
}

impl<R: Read> BufRead for CompactBufReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.fill_buffer()
    }

    fn consume(&mut self, amt: usize) {
        // Call the inherent method to avoid recursion.
        self.consume_front(amt);
    }
}

impl<R: Read> Read for CompactBufReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }

        if self.filled_len() == 0 {
            self.fill_buffer()?;
            if self.filled_len() == 0 {
                return Ok(0); // EOF
            }
        }

        let n = cmp::min(out.len(), self.filled_len());
        out[..n].copy_from_slice(&self.buf[self.start..self.start + n]);
        self.consume_front(n);
        Ok(n)
    }
}
